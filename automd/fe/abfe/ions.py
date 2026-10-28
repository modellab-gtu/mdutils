"""Co-alchemical counter-ions for charged ligands.

A ligand of net charge q is decoupled together with |q| ions of the opposite
sign, so the box stays neutral in the A state (all charges on) and in the
decoupled B state (ligand charge 0, ion charge 0).  The ions are moved to the
end of the coordinates and of ``[ molecules ]`` as a new moleculetype
``<ION>_ALCH`` whose B-state charge is 0 and whose LJ type is unchanged, so
their charge follows the Coulomb lambda.  A flat-bottom position restraint
keeps each ion near its start position.

Sign convention (verified by tests): for acetate (q = -1) the chosen ion is
Na+; afterwards ``top.charge("A") == 0`` and ``top.charge("B") == -1`` (only
the ion is perturbed in the topology) while ``fe_state_charge(top, lig, "B")``
(B charge with the ligand's own charges zeroed, as couple-moltype does) == 0.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

from .topology import Topology, ligand_charge

POSRE_R = 0.5      # nm, flat-bottom sphere radius
POSRE_K = 1000.0   # kJ mol-1 nm-2
PROTEIN_MIN_ATOMS = 100


def _walk(top: Topology):
    """Per-molecule list of (name, first_atom_index, natoms), in coordinate order."""
    out, i = [], 0
    for name, cnt in top.molecules:
        n = len(top.require_moleculetype(name).atoms)
        for _ in range(cnt):
            out.append((name, i, n))
            i += n
    return out, i


def _single_ion_charge(top: Topology, name: str):
    mt = top.require_moleculetype(name)
    if len(mt.atoms) != 1:
        return None
    q = mt.atoms[0].q
    return q if abs(q) > 0.5 else None


def fe_state_charge(top: Topology, ligand: str, state: str = "B") -> float:
    """Total charge of ``state`` with the ligand's charges zeroed (couple-moltype)."""
    mt = top.require_moleculetype(ligand)
    n_lig = sum(c for n, c in top.molecules if n == ligand)
    return top.charge(state) - n_lig * mt.charge(state)


def check_charges(top: Topology, ligand: str, tol: float = 1e-3) -> None:
    """State A neutral; decoupled B (ligand charges zeroed) neutral."""
    qa, qb = top.charge("A"), fe_state_charge(top, ligand, "B")
    if abs(qa) > tol or abs(qb) > tol:
        raise RuntimeError(
            f"charge bookkeeping failed: state A {qa:+.4f}, decoupled state B "
            f"{qb:+.4f} (both must be 0)")


def choose_ions(u: mda.Universe, top: Topology, ligand: str, q_lig: int,
                min_lig_dist: float = 1.5) -> list[int]:
    """0-based indices of |q_lig| ions of sign opposite to q_lig (distances in nm)."""
    if q_lig == 0:
        return []
    mols, natoms = _walk(top)
    if natoms != len(u.atoms):
        raise ValueError(f"topology has {natoms} atoms, coordinates {len(u.atoms)}")
    pos = u.atoms.positions.astype(float)
    box = u.dimensions
    lig = [k for name, i, n in mols if name == ligand for k in range(i, i + n)]
    if not lig:
        raise ValueError(f"ligand '{ligand}' not in [ molecules ]")
    prot = [k for name, i, n in mols
            if name != ligand and n >= PROTEIN_MIN_ATOMS for k in range(i, i + n)]
    want = 1.0 if q_lig < 0 else -1.0
    cand = [i for name, i, n in mols
            if n == 1 and (q := _single_ion_charge(top, name)) is not None
            and np.sign(q) == want]
    ok, d_ref = [], []
    if cand:
        cpos = pos[cand]
        d_lig = distance_array(cpos, pos[lig], box=box).min(axis=1) / 10.0
        ref = (distance_array(cpos, pos[prot], box=box).min(axis=1) / 10.0
               if prot else d_lig)
        for c, dl, dr in zip(cand, d_lig, ref):
            if dl >= min_lig_dist:
                ok.append(c)
                d_ref.append(dr)
    need = abs(int(q_lig))
    if len(ok) < need:
        raise RuntimeError(
            f"only {len(ok)} suitable {'cations' if want > 0 else 'anions'} at least "
            f"{min_lig_dist} nm from {ligand}, need {need}; re-prepare with --salt "
            "(or a larger box)")
    order = np.argsort(-np.asarray(d_ref), kind="stable")
    return [ok[k] for k in order[:need]]


def _atomtype_mass(top: Topology, typ: str) -> float:
    """Mass from [ atomtypes ] (atoms lines may omit the mass column)."""
    sec = None
    for raw in top.lines:
        s = raw.split(";", 1)[0].strip()
        if not s or s.startswith("#"):
            continue
        if s.startswith("["):
            sec = s.strip("[] \t").lower()
            continue
        if sec == "atomtypes":
            t = s.split()
            if t[0] == typ:
                # name [bond_type] [at.num] mass charge ptype sigma epsilon
                j = next((k for k in range(3, len(t)) if t[k] in ("A", "S", "V", "D")), None)
                if j is not None:
                    return float(t[j - 2])
    raise ValueError(f"cannot find mass of atomtype '{typ}'")


def _alch_text(top: Topology, src: str) -> str:
    a = top.require_moleculetype(src).atoms[0]
    if a.mass <= 0:
        a = type(a)(**{**a.__dict__, "mass": _atomtype_mass(top, a.type)})
    name = f"{src}_ALCH"
    return (
        "[ moleculetype ]\n"
        f"{name:<12s} 1\n\n"
        "[ atoms ]\n"
        "; nr type resnr residue atom cgnr charge mass typeB chargeB massB\n"
        f"{1:6d} {a.type:>10s} {1:6d} {a.res:>6s} {a.name:>6s} {1:6d} "
        f"{a.q:12.6f} {a.mass:10.5f} {a.type:>10s} {0.0:12.6f} {a.mass:10.5f}\n\n"
        "[ position_restraints ]\n"
        "; ai funct g(1=sphere) r(nm) k\n"
        f"{1:6d} {2:5d} {1:5d} {POSRE_R:8.3f} {POSRE_K:10.1f}\n")


def apply_coalchemical(top: Topology, gro_in, gro_out, ion_idx: list[int]) -> Topology:
    """Move ``ion_idx`` to the end of gro/[ molecules ] as ``<ION>_ALCH``.

    Edits ``top`` in place through Topology's structural operations and also
    returns it.  Writes ``gro_out`` and ``posre_ref.gro`` (same content) next
    to it.
    """
    lines = Path(gro_in).read_text().splitlines()
    n_at = int(lines[1])
    atom_lines, tail = lines[2:2 + n_at], lines[2 + n_at:]
    mols, natoms = _walk(top)
    if natoms != n_at:
        raise ValueError(f"topology has {natoms} atoms, gro has {n_at}")
    by_start = {i: (name, n) for name, i, n in mols}
    chosen = list(ion_idx)
    if len(set(chosen)) != len(chosen):
        raise ValueError("duplicate ion indices")
    for k in chosen:
        if k not in by_start or by_start[k][1] != 1:
            raise ValueError(f"atom {k} is not a single-atom molecule")

    alch_names = {}
    for k in chosen:
        src = by_start[k][0]
        alch = f"{src}_ALCH"
        alch_names[k] = alch
        if alch not in top.moleculetypes:
            # insert right after the source ion moleculetype
            top.add_moleculetype(alch, _alch_text(top, src), after=src)

    cs = set(chosen)
    keep_names, keep_lines = [], []
    for name, i, n in mols:
        if i in cs:
            continue
        keep_names.append(name)
        keep_lines.extend(atom_lines[i:i + n])
    new_names = keep_names + [alch_names[k] for k in chosen]
    new_lines = keep_lines + [atom_lines[k] for k in chosen]

    runs: list[list] = []
    for nm in new_names:
        if runs and runs[-1][0] == nm:
            runs[-1][1] += 1
        else:
            runs.append([nm, 1])
    top.set_molecules([(n, c) for n, c in runs])

    out = [f"{l[:15]}{(j + 1) % 100000:5d}{l[20:]}" for j, l in enumerate(new_lines)]
    gro_out = Path(gro_out)
    gro_out.write_text("\n".join([lines[0], f"{n_at:5d}"] + out + tail) + "\n")
    shutil.copyfile(gro_out, gro_out.parent / "posre_ref.gro")
    return top
