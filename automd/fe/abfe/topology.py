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
    def __init__(self):
        self.lines: list[str] = []
        self.moleculetypes: dict[str, MolType] = {}
        self.molecules: list[tuple[str, int]] = []
        self._mol_lines: list[int] = []
        self._mol_snap: list[tuple[str, int]] = []

    @classmethod
    def parse(cls, path) -> "Topology":
        top = cls()
        top.lines = Path(path).read_text().splitlines()
        section = None
        cur: MolType | None = None
        expect_name = False
        for i, raw in enumerate(top.lines):
            s = _clean(raw)
            if not s or s.startswith("#"):
                continue
            if s.startswith("["):
                section = s.strip("[] \t").lower()
                expect_name = section == "moleculetype"
                continue
            toks = s.split()
            if section == "moleculetype" and expect_name:
                cur = MolType(toks[0], int(toks[1]) if len(toks) > 1 else 3)
                top.moleculetypes[cur.name] = cur
                expect_name = False
            elif section == "atoms" and cur is not None:
                a = _parse_atom(toks)
                cur.atoms.append(a)
                cur._lines.append(i)
                cur._snap.append(astuple(a))
            elif section == "molecules":
                top.molecules.append((toks[0], int(toks[1])))
                top._mol_lines.append(i)
        top._mol_snap = list(top.molecules)
        return top

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
        out = list(self.lines)
        for mt in self.moleculetypes.values():
            for a, idx, snap in zip(mt.atoms, mt._lines, mt._snap):
                if astuple(a) != snap:
                    out[idx] = _fmt_atom(a)
        if self.molecules != self._mol_snap:
            if not self._mol_lines:
                raise ValueError("topology has no [ molecules ] section to edit")
            first, last = self._mol_lines[0], self._mol_lines[-1]
            new = [f"{n:<12s} {c}" for n, c in self.molecules]
            out[first:last + 1] = new
        Path(path).write_text("\n".join(out) + "\n")


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
    em = run.path / "MDP" / "em.mdp"
    with tempfile.TemporaryDirectory(dir=outdir, prefix="pp_") as td:
        if not em.is_file():
            em = Path(td) / "em.mdp"
            write_mdp({"integrator": "steep", "nsteps": 1}, em)
        cmd = ["gmx", "grompp", "-f", str(em), "-c", str(run.md_gro),
               "-r", str(run.md_gro), "-p", str(run.top), "-pp", str(out.resolve()),
               "-o", str(Path(td) / "pp.tpr"), "-maxwarn", "2"]
        r = subprocess.run(cmd, cwd=td,
                           capture_output=True, text=True,
                           env={**os.environ, "GMX_MAXBACKUP": "-1"})
    if r.returncode != 0 or not out.is_file():
        raise RuntimeError(f"grompp -pp failed:\n{r.stderr[-2000:]}")
    return out
