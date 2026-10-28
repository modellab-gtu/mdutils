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
from dataclasses import replace
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array, minimize_vectors

from .topology import Topology, ligand_charge

POSRE_R = 0.5      # nm, flat-bottom sphere radius
POSRE_K = 1000.0   # kJ mol-1 nm-2
# Heuristic: any non-ligand molecule with at least this many atoms counts as
# "protein" (water, ions and small cosolvents are far below it).
PROTEIN_MIN_ATOMS = 100
MAX_D_MIN = 1.5    # nm
BOX_FRACTION = 0.45
WATER_NAMES = ("SOL", "HOH", "WAT", "TIP3", "TIP3P", "SPC")


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
    return q if abs(abs(q) - 1.0) < 0.01 else None   # monovalent only


def fe_state_charge(top: Topology, ligand: str, state: str = "B") -> float:
    """Total charge of ``state`` with the ligand's charges zeroed (couple-moltype)."""
    mt = top.require_moleculetype(ligand)
    n_lig = sum(c for n, c in top.molecules if n == ligand)
    return top.charge(state) - n_lig * mt.charge(state)


def check_charges(top: Topology, ligand: str, tol: float = 1e-3) -> None:
    """State A neutral; decoupled B (ligand charges zeroed) neutral."""
    qa, qb = top.charge("A"), fe_state_charge(top, ligand, "B")
    if abs(qa) > tol or abs(qb) > tol:
        raise ValueError(
            f"charge bookkeeping failed: state A {qa:+.4f}, decoupled state B "
            f"{qb:+.4f} (both must be 0)")


def default_d_min(box) -> float:
    """min(1.5 nm, 0.45 x shortest box edge); ``box`` is MDAnalysis dimensions (A)."""
    return min(MAX_D_MIN, BOX_FRACTION * float(np.min(box[:3])) / 10.0)


def choose_ions(u: mda.Universe, top: Topology, ligand: str, q_lig: int,
                min_lig_dist: float | None = None) -> list[int]:
    """0-based indices of |q_lig| monovalent ions of sign opposite to q_lig.

    Ions at least d_min (default ``default_d_min``) from the ligand
    (minimum image) come first, farthest from the protein first (farthest from
    the ligand when there is no protein).  If too few qualify, the remaining
    ones are the opposite-sign ions farthest from the ligand;
    ``apply_coalchemical`` then swaps those with distant waters.
    """
    if q_lig == 0:
        return []
    d_min = default_d_min(u.dimensions) if min_lig_dist is None else min_lig_dist
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
    need = abs(int(q_lig))
    if len(cand) < need:
        raise RuntimeError(
            f"only {len(cand)} monovalent {'cations' if want > 0 else 'anions'} in the "
            f"box, need {need}; re-prepare with --salt")
    cpos = pos[cand]
    d_lig = distance_array(cpos, pos[lig], box=box).min(axis=1) / 10.0
    d_ref = (distance_array(cpos, pos[prot], box=box).min(axis=1) / 10.0
             if prot else d_lig)
    far = [k for k in range(len(cand)) if d_lig[k] >= d_min]
    far.sort(key=lambda k: -d_ref[k])
    rest = sorted((k for k in range(len(cand)) if d_lig[k] < d_min), key=lambda k: -d_lig[k])
    return [cand[k] for k in (far + rest)[:need]]


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
        a = replace(a, mass=_atomtype_mass(top, a.type))
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


def _swap_waters(top, mols, u, pos, ion_idx, lig, d_min):
    """Swap ions closer than d_min to the ligand with the farthest waters.

    Returns new positions (A).  Whole water goes to the ion's site (first atom
    on the ion), the ion to the water's first-atom site; each water used once.
    """
    box = u.dimensions
    near = [k for k in ion_idx
            if distance_array(pos[[k]], pos[lig], box=box).min() / 10.0 < d_min]
    if not near:
        return pos
    waters = [(i, n) for name, i, n in mols if name in WATER_NAMES]
    if len(waters) < len(near):
        raise RuntimeError("not enough water molecules to swap with the counter-ion")
    first = np.array([i for i, _ in waters])
    dw = distance_array(pos[first], pos[lig], box=box).min(axis=1)
    order = np.argsort(-dw, kind="stable")
    pos = pos.copy()
    for ion, w in zip(near, order):
        wi, wn = waters[w]
        ion_pos, w0 = pos[ion].copy(), pos[wi].copy()
        rel = minimize_vectors(pos[wi:wi + wn] - w0, box)   # whole water
        pos[wi:wi + wn] = ion_pos + rel
        pos[ion] = w0
    return pos


def apply_coalchemical(top: Topology, gro_in, gro_out, ion_idx: list[int],
                       ligand: str, min_lig_dist: float | None = None) -> Topology:
    """Move ``ion_idx`` to the end of gro/[ molecules ] as ``<ION>_ALCH``.

    Ions nearer than d_min to the ligand are first swapped with the water
    farthest from the ligand.  Edits ``top`` in place through Topology's
    structural operations and also returns it.  Writes ``gro_out`` and
    ``posre_ref.gro`` (same content) next to it.  Raises ValueError if an index
    is not a monovalent single-atom ion of sign opposite to the ligand or if the
    charge checks (A and decoupled B both 0) fail.
    """
    q_lig = ligand_charge(top, ligand)
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
    if len(chosen) != abs(q_lig):
        raise ValueError(f"ligand {ligand} has charge {q_lig:+d}: need {abs(q_lig)} "
                         f"ions, got {len(chosen)}")
    for k in chosen:
        if k not in by_start or by_start[k][1] != 1:
            raise ValueError(f"atom {k} is not a single-atom molecule")
        q = _single_ion_charge(top, by_start[k][0])
        if q is None or np.sign(q) == np.sign(q_lig):
            raise ValueError(f"atom {k} ({by_start[k][0]}) is not a monovalent ion "
                             f"of sign opposite to ligand charge {q_lig:+d}")

    # coordinate swap (only when an ion is too close to the ligand)
    u = mda.Universe(str(gro_in))
    d_min = default_d_min(u.dimensions) if min_lig_dist is None else min_lig_dist
    lig = [k for name, i, n in mols if name == ligand for k in range(i, i + n)]
    old = u.atoms.positions.astype(float)
    new = _swap_waters(top, mols, u, old, chosen, lig, d_min)
    for j in np.flatnonzero(np.any(new != old, axis=1)):
        x, y, z = new[j] / 10.0
        atom_lines[j] = f"{atom_lines[j][:20]}{x:8.3f}{y:8.3f}{z:8.3f}{atom_lines[j][44:]}"

    alch_names = {}
    for k in chosen:
        src = by_start[k][0]
        alch = f"{src}_ALCH"
        alch_names[k] = alch
        if alch not in top.moleculetypes:
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
    check_charges(top, ligand)

    out = [f"{l[:15]}{(j + 1) % 100000:5d}{l[20:]}" for j, l in enumerate(new_lines)]
    gro_out = Path(gro_out)
    gro_out.write_text("\n".join([lines[0], f"{n_at:5d}"] + out + tail) + "\n")
    shutil.copyfile(gro_out, gro_out.parent / "posre_ref.gro")
    return top
