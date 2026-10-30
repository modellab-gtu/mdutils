"""Self-contained FE topology: preprocessing, parsing, charge checks.

``preprocess`` expands all #include/#ifdef via ``gmx grompp -pp`` so both the
pdb2gmx and the ParmEd route yield one standalone topology.  ``Topology``
parses only what the FE workflow needs ([ moleculetype ], [ atoms ],
[ molecules ]); every other line is kept verbatim so ``write`` round-trips
through grompp and later edits (extra moleculetypes, restraint blocks) can be
spliced in via ``lines``.
"""
from __future__ import annotations

import os
import subprocess
import tempfile
from dataclasses import dataclass, field, astuple
from pathlib import Path
from typing import Literal

from .rundir import RunDir, write_mdp


@dataclass
class Atom:
    nr: int
    type: str
    resnr: int
    res: str
    name: str
    cgnr: int
    q: float
    mass: float
    typeB: str
    qB: float
    massB: float


@dataclass
class MolType:
    name: str
    nrexcl: int = 3
    atoms: list[Atom] = field(default_factory=list)
    # line index of each atom in Topology.lines and a snapshot to detect edits
    _lines: list[int] = field(default_factory=list, repr=False)
    _snap: list[tuple] = field(default_factory=list, repr=False)
    _fresh: bool = field(default=True, repr=False)

    def charge(self, state: str = "A") -> float:
        return sum(a.qB if state == "B" else a.q for a in self.atoms)


def _clean(line: str) -> str:
    return line.split(";", 1)[0].strip()


def _parse_atom(toks: list[str]) -> Atom:
    nr, typ, resnr, res, name, cgnr = (int(toks[0]), toks[1], int(toks[2]),
                                       toks[3], toks[4], int(toks[5]))
    q = float(toks[6])
    mass = float(toks[7]) if len(toks) > 7 else 0.0
    typeB = toks[8] if len(toks) > 8 else typ
    qB = float(toks[9]) if len(toks) > 9 else q
    massB = float(toks[10]) if len(toks) > 10 else mass
    return Atom(nr, typ, resnr, res, name, cgnr, q, mass, typeB, qB, massB)


def _fmt_atom(a: Atom) -> str:
    return (f"{a.nr:6d} {a.type:>10s} {a.resnr:6d} {a.res:>6s} {a.name:>6s} "
            f"{a.cgnr:6d} {a.q:12.6f} {a.mass:10.5f} {a.typeB:>10s} "
            f"{a.qB:12.6f} {a.massB:10.5f}")


class Topology:
    """Parsed view over the raw topology lines.

    Atom edits (``MolType.atoms``) and ``molecules`` edits are flushed into
    ``_lines`` before every structural operation and in ``write``; structural
    operations re-index afterwards, so no stale indices exist.
    """

    def __init__(self):
        self._lines: list[str] = []
        self.moleculetypes: dict[str, MolType] = {}
        self.molecules: list[tuple[str, int]] = []
        self._mol_snap: list[tuple[str, int]] = []

    @property
    def lines(self) -> list[str]:
        self._sync()
        return list(self._lines)

    @classmethod
    def parse(cls, path) -> "Topology":
        top = cls()
        top._lines = Path(path).read_text().splitlines()
        top._index()
        return top

    # -- indexing / syncing -------------------------------------------------
    def _index(self) -> None:
        old = self.moleculetypes
        new: dict[str, MolType] = {}
        section = None
        cur = None
        expect_name = False
        mols: list[tuple[str, int]] = []
        for i, raw in enumerate(self._lines):
            s = _clean(raw)
            if not s or s.startswith("#"):
                continue
            if s.startswith("["):
                section = s.strip("[] \t").lower()
                expect_name = section == "moleculetype"
                continue
            toks = s.split()
            if section == "moleculetype" and expect_name:
                if toks[0] in new:
                    raise ValueError(f"duplicate moleculetype '{toks[0]}'")
                cur = old.get(toks[0]) or MolType(toks[0])
                cur.nrexcl = int(toks[1]) if len(toks) > 1 else 3
                cur._lines = []
                cur._fresh = toks[0] not in old
                new[cur.name] = cur
                expect_name = False
            elif section == "atoms" and cur is not None:
                cur._lines.append(i)
                if cur._fresh:
                    cur.atoms.append(_parse_atom(toks))
            elif section == "molecules":
                mols.append((toks[0], int(toks[1])))
        for mt in new.values():
            mt._snap = [astuple(a) for a in mt.atoms]
            mt._fresh = False
            if len(mt._lines) != len(mt.atoms):
                raise ValueError(f"atom count mismatch for {mt.name}")
        self.moleculetypes = new
        self.molecules = mols
        self._mol_snap = list(mols)

    def _section_range(self, name: str):
        """(header_idx, end_idx) of a section; end is exclusive."""
        hdr = None
        for i, raw in enumerate(self._lines):
            s = _clean(raw)
            if s.startswith("["):
                if hdr is not None:
                    return hdr, i
                if s.strip("[] \t").lower() == name:
                    hdr = i
        return (hdr, len(self._lines)) if hdr is not None else (None, None)

    def _mol_data_range(self):
        hdr, end = self._section_range("molecules")
        if hdr is None:
            return None
        idx = [i for i in range(hdr + 1, end)
               if _clean(self._lines[i]) and not _clean(self._lines[i]).startswith("#")]
        return (idx[0], idx[-1] + 1) if idx else (hdr + 1, hdr + 1)

    def _sync(self) -> None:
        for key, mt in self.moleculetypes.items():
            if key != mt.name or mt._fresh:
                raise ValueError(
                    f"moleculetype '{mt.name}' was not added through "
                    "add_moleculetype and cannot be written")
        for mt in self.moleculetypes.values():
            for k, (a, idx) in enumerate(zip(mt.atoms, mt._lines)):
                if astuple(a) != mt._snap[k]:
                    self._lines[idx] = _fmt_atom(a)
                    mt._snap[k] = astuple(a)
        if self.molecules != self._mol_snap:
            self._apply_molecules(self.molecules)

    def _apply_molecules(self, mols) -> None:
        rng = self._mol_data_range()
        if rng is None:
            raise ValueError("topology has no [ molecules ] section")
        new = [f"{n:<12s} {c}" for n, c in mols]
        self._lines[rng[0]:rng[1]] = new
        self.molecules = list(mols)
        self._mol_snap = list(mols)

    # -- public structural operations ---------------------------------------
    def add_moleculetype(self, name: str, text: str, after: str | None = None) -> MolType:
        """Insert a full ``[ moleculetype ]`` block (default: before [ system ])."""
        self._sync()
        if name in self.moleculetypes:
            raise ValueError(f"moleculetype '{name}' already exists")
        block = text.rstrip("\n").splitlines() + [""]
        if after is not None:
            mt = self.require_moleculetype(after)
            start = mt._lines[0]
            while not _clean(self._lines[start]).startswith("[") or \
                    _clean(self._lines[start]).strip("[] \t").lower() != "moleculetype":
                start -= 1
            pos = len(self._lines)
            for i in range(mt._lines[-1] + 1, len(self._lines)):
                s = _clean(self._lines[i])
                if s.startswith("[") and s.strip("[] \t").lower() in ("moleculetype", "system"):
                    pos = i
                    break
        else:
            pos, _ = self._section_range("system")
            if pos is None:
                pos = len(self._lines)
        self._lines[pos:pos] = block
        self._index()
        return self.require_moleculetype(name)

    def set_molecules(self, mols: list[tuple[str, int]]) -> None:
        self._sync()
        for n, _ in mols:
            self.require_moleculetype(n)
        self._apply_molecules(list(mols))
        self._index()

    def append_section(self, text: str) -> None:
        """Append raw text (e.g. [ intermolecular_interactions ]) at the end."""
        self._sync()
        self._lines.extend([""] + text.rstrip("\n").splitlines())
        self._index()

    # -- queries --------------------------------------------------------------
    def require_moleculetype(self, name: str) -> MolType:
        if name not in self.moleculetypes:
            raise ValueError(
                f"moleculetype '{name}' not found in topology; available: "
                f"{', '.join(self.moleculetypes)}")
        return self.moleculetypes[name]

    def charge(self, state: Literal["A", "B"] = "A") -> float:
        total = 0.0
        for name, n in self.molecules:
            total += n * self.require_moleculetype(name).charge(state)
        return total

    def write(self, path) -> None:
        self._sync()
        Path(path).write_text("\n".join(self._lines) + "\n")


def ligand_charge(top: Topology, ligand: str) -> int:
    q = top.require_moleculetype(ligand).charge("A")
    n = round(q)
    if abs(q - n) > 0.01:
        raise ValueError(
            f"ligand {ligand} charge {q:.4f} is not an integer (|dq|>0.01)")
    return int(n)


def preprocess(run: RunDir, outdir: Path) -> Path:
    """grompp -pp the run folder's topology into outdir/processed.top."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "processed.top"
    with tempfile.TemporaryDirectory(dir=outdir, prefix="pp_") as td:
        # minimal, define-free mdp: FE windows use no #define restraints
        em = Path(td) / "pp.mdp"
        write_mdp({"integrator": "steep", "nsteps": 0}, em)
        cmd = ["gmx", "grompp", "-f", str(em), "-c", str(run.md_gro),
               "-p", str(run.top), "-pp", str(out.resolve()),
               "-o", str(Path(td) / "pp.tpr"), "-maxwarn", "0"]
        r = subprocess.run(cmd, cwd=td,
                           capture_output=True, text=True,
                           env={**os.environ, "GMX_MAXBACKUP": "-1"})
    if r.returncode != 0 or not out.is_file():
        raise RuntimeError(f"grompp -pp failed:\n{r.stderr[-2000:]}")
    return out
