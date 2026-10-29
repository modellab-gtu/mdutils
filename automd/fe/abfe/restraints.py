"""Boresch orientational restraint for the ABFE complex leg.

Selection uses MDRestraintsGenerator (MDRG) 0.2.1 on the equilibrium run
(``md_<r>.tpr/xtc``); the restraint is written as a GROMACS
``[ intermolecular_interactions ]`` block and its analytic release correction
(Boresch et al. 2003, eq. 32) is computed here.

Geometry and atom order
-----------------------
``Boresch.atoms`` = (P3, P2, P1, L1, L2, L3), 0-based global indices. In MDRG
terms L1..L3 = ``l_atoms[0..2]`` and P1..P3 = ``p_atoms[0..2]``
(P1 = CA, then the backbone C and N of that residue).

* r      = |P1-L1|                (MDRG ``bond``)
* thetaA = angle(P2, P1, L1)      (MDRG ``angles[1]``, reversed atom order)
* thetaB = angle(P1, L1, L2)      (MDRG ``angles[0]``, reversed)
* phiA   = dihedral(P3, P2, P1, L1) (MDRG ``dihedrals[2]``, reversed)
* phiB   = dihedral(P2, P1, L1, L2) (MDRG ``dihedrals[1]``, reversed)
* phiC   = dihedral(P1, L1, L2, L3) (MDRG ``dihedrals[0]``, reversed)

Reversing the atom order leaves angles and (IUPAC) dihedrals unchanged.

Indices refer to the atom ordering of ``md_<r>.tpr`` (the run being
analysed). The FE topology keeps that ordering for protein and ligand: the
only reordering (Task 5, co-alchemical ions) moves ions to the END of the
system, so protein/ligand atoms must precede any moved ion for the indices to
remain valid.

Units and sign
--------------
GROMACS topology units, all harmonic as V = 1/2 k (x - x0)^2:
bond funct 6: b0 nm, k kJ mol^-1 nm^-2; angle funct 1: theta0 deg,
k kJ mol^-1 rad^-2; dihedral funct 2 (harmonic improper): xi0 deg,
k kJ mol^-1 rad^-2. A state (bonded-lambda 0) has k = 0, B state
(bonded-lambda 1) the full k; equilibrium values are identical in A and B.
GROMACS switches [ intermolecular_interactions ] bonded terms with
bonded-lambdas (the FE MDPs set bonded-lambdas = restraint-lambdas).

``boresch_correction`` returns dG_restr,off (kJ/mol, negative): the free
energy of releasing the restraint on the decoupled ligand into the standard
volume V0. Binding: dG_bind = dG_solv,decouple - dG_cplx,decouple - dG_restr,off.
MDRG's ``BoreschRestraint.standard_state`` returns the same quantity, same
sign, in kcal/mol (force_constant=10 kcal/mol/rad^2 == our defaults).
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

KB = 0.0083144626      # kJ/mol/K
V0 = 1.66054           # nm^3, standard-state volume (1 M)
KCAL = 4.184           # kJ per kcal


@dataclass
class Boresch:
    atoms: tuple            # 6 x 0-based global index, order P3 P2 P1 L1 L2 L3
    names: tuple
    r0: float               # nm
    thetaA0: float          # deg, angle(P2, P1, L1)
    thetaB0: float          # deg, angle(P1, L1, L2)
    phiA0: float            # deg, dihedral(P3, P2, P1, L1)
    phiB0: float            # deg, dihedral(P2, P1, L1, L2)
    phiC0: float            # deg, dihedral(P1, L1, L2, L3)
    kr: float = 4184.0      # kJ/mol/nm^2
    ktheta: float = 41.84   # kJ/mol/rad^2
    kphi: float = 41.84     # kJ/mol/rad^2
    frame: int = -1         # trajectory frame the equilibrium values come from

    def __post_init__(self):
        self.atoms = tuple(int(a) for a in self.atoms)
        self.names = tuple(self.names)
        if len(self.atoms) != 6 or len(set(self.atoms)) != 6:
            raise ValueError(f"Boresch needs 6 distinct atoms, got {self.atoms}")


def boresch_correction(b: Boresch, T: float) -> float:
    """dG_restr,off in kJ/mol (Boresch 2003 eq. 32), negative for usual k.

    -kT ln[ 8 pi^2 V0 sqrt(Kr KthA KthB KphA KphB KphC)
            / (r0^2 sin thA sin thB (2 pi kT)^3) ]
    """
    kT = KB * T
    num = 8.0 * math.pi ** 2 * V0 * math.sqrt(b.kr * b.ktheta ** 2 * b.kphi ** 3)
    den = (b.r0 ** 2 * math.sin(math.radians(b.thetaA0))
           * math.sin(math.radians(b.thetaB0)) * (2.0 * math.pi * kT) ** 3)
    return -kT * math.log(num / den)


def check_collinear(b: Boresch, min_deg: float = 10.0) -> None:
    for name, th in (("thetaA", b.thetaA0), ("thetaB", b.thetaB0)):
        if th < min_deg or th > 180.0 - min_deg:
            raise ValueError(
                f"Boresch anchor angle {name} = {th:.1f} deg is within {min_deg:g} deg "
                f"of collinear (0/180); pick other anchors or use --force-restraint")


def wall_bond_lines(wall_bonds) -> list[str]:
    """funct-10 (restraint potential) lines: low=d_wall, up1, up2, kdr; no B state."""
    # B-state written equal to A: grompp -maxwarn 0 rejects the implicit copy for
    # perturbed atoms; A == B keeps the wall lambda-independent
    return [f"{a:>7d}{b:>7d}    10 {d:10.4f} {u1:10.3f} {u2:10.3f} {k:10.3f}"
            f" {d:10.4f} {u1:10.3f} {u2:10.3f} {k:10.3f}"
            for a, b, d, u1, u2, k in wall_bonds]


def wall_block(wall_bonds) -> str:
    """Standalone ``[ intermolecular_interactions ]`` with only the ion wall bonds."""
    return "\n".join(["; co-alchemical ion wall (automd.fe.abfe), funct 10, lambda-independent",
                      "[ intermolecular_interactions ]", "[ bonds ]",
                      ";    ai     aj  type    lowA    up1A    up2A    kdrA    lowB    up1B    up2B    kdrB"]
                     + wall_bond_lines(wall_bonds)) + "\n"


def intermolecular_block(b: Boresch, wall_bonds=()) -> str:
    """The single ``[ intermolecular_interactions ]`` text (append after [ molecules ]):
    Boresch terms plus optional ion-wall bonds (a, b, low, up1, up2, kdr; 1-based)."""
    i = [a + 1 for a in b.atoms]  # 1-based global
    P3, P2, P1, L1, L2, L3 = i
    out = ["; Boresch restraint (automd.fe.abfe); A: k=0 (off), B: full k",
           "[ intermolecular_interactions ]",
           "[ bonds ]",
           ";    ai     aj  type       bA       kA       bB         kB",
           f"{P1:>7d}{L1:>7d}     6 {b.r0:10.5f}      0.0 {b.r0:10.5f} {b.kr:10.3f}"]
    if wall_bonds:
        out += ["; co-alchemical ion wall, funct 10 (no B state: lambda-independent)",
                ";    ai     aj  type    lowA    up1A    up2A    kdrA    lowB    up1B    up2B    kdrB"] + wall_bond_lines(wall_bonds)
    out += ["[ angles ]",
           ";    ai     aj     ak  type     thA      kA     thB         kB"]
    for (a1, a2, a3), th in (((P2, P1, L1), b.thetaA0), ((P1, L1, L2), b.thetaB0)):
        out.append(f"{a1:>7d}{a2:>7d}{a3:>7d}     1 {th:10.4f}      0.0 {th:10.4f} "
                   f"{b.ktheta:10.4f}")
    out += ["[ dihedrals ]",
            ";    ai     aj     ak     al  type    phiA      kA    phiB         kB"]
    for (a1, a2, a3, a4), ph in (((P3, P2, P1, L1), b.phiA0), ((P2, P1, L1, L2), b.phiB0),
                                 ((P1, L1, L2, L3), b.phiC0)):
        out.append(f"{a1:>7d}{a2:>7d}{a3:>7d}{a4:>7d}     2 {ph:10.4f}      0.0 "
                   f"{ph:10.4f} {b.kphi:10.4f}")
    return "\n".join(out) + "\n"


def to_json(b: Boresch, T: float, path) -> dict:
    dg = boresch_correction(b, T)
    d = {
        "atoms_0based": list(b.atoms),
        "atoms_1based": [a + 1 for a in b.atoms],
        "order": ["P3", "P2", "P1", "L1", "L2", "L3"],
        "names": list(b.names),
        "r0_nm": b.r0, "thetaA0_deg": b.thetaA0, "thetaB0_deg": b.thetaB0,
        "phiA0_deg": b.phiA0, "phiB0_deg": b.phiB0, "phiC0_deg": b.phiC0,
        "kr_kJ_mol_nm2": b.kr, "ktheta_kJ_mol_rad2": b.ktheta,
        "kphi_kJ_mol_rad2": b.kphi, "frame": b.frame,
        "T": T, "V0_nm3": V0,
        "dG_restr_off_kJ": dg, "dG_restr_off_kcal": dg / KCAL,
        "sign": "dG_restr_off = release of restraint into V0 (negative); "
                "dG_bind = dG_solv,decouple - dG_cplx,decouple - dG_restr_off",
        "note": "atom indices follow md_<r>.tpr ordering",
    }
    Path(path).write_text(json.dumps(d, indent=2) + "\n")
    return d


def find_boresch(tpr, xtc, ligand_resname: str, T: float, *, max_frames: int = 200,
                 force: bool = False) -> Boresch:
    """Pick a Boresch restraint from an equilibrium trajectory with MDRG.

    Ligand anchors: MDRG ``find_ligand_atoms`` (lowest-RMSF heavy atoms after
    aligning on protein CA); host anchors: ``FindHostAtoms`` (protein CA within
    5-9 A of L1 over the trajectory, plus that residue's backbone C and N);
    ``FindBoreschRestraint`` ranks candidates by summed variance and picks
    the frame closest to the mean (MDRG ``min_frame``); the equilibrium values
    are the values in that frame (same convention as MDRG's writer and
    ``standard_state``). The trajectory is strided to at most ``max_frames``
    frames. ``T`` is accepted for interface symmetry (the correction is
    computed later with ``boresch_correction``). Raises
    ValueError on near-collinear anchors unless ``force``.
    """
    import warnings

    import MDAnalysis as mda
    from MDRestraintsGenerator.restraints import FindBoreschRestraint
    from MDRestraintsGenerator.search import FindHostAtoms, find_ligand_atoms

    u = mda.Universe(str(tpr), str(xtc))
    n = u.trajectory.n_frames
    step = max(1, math.ceil(n / max_frames))
    if step > 1:
        u.transfer_to_memory(step=step)
    l_sel = f"resname {ligand_resname} and not name H*"
    if len(u.select_atoms(l_sel)) < 3:
        raise ValueError(f"fewer than 3 heavy atoms match '{l_sel}'")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        l_atoms = find_ligand_atoms(u, l_selection=l_sel, p_align="protein and name CA")
        atom_set = []
        for la in l_atoms:
            fh = FindHostAtoms(u, la[0], p_selection="protein and name CA").run()
            for pa in fh.host_atoms:
                atom_set.append((la, pa))
        if not atom_set:
            raise RuntimeError("MDRestraintsGenerator found no host anchor atoms")
        fb = FindBoreschRestraint(u, atom_set).run()
    rest = fb.restraint
    f = int(rest.min_frame)
    d2 = rest.dihedrals[2].atomgroup.ix   # L1 P1 P2 P3
    d0 = rest.dihedrals[0].atomgroup.ix   # L3 L2 L1 P1
    atoms = (d2[3], d2[2], d2[1], d0[2], d0[1], d0[0])
    b = Boresch(
        atoms=atoms,
        names=tuple(f"{u.atoms[a].resname}{u.atoms[a].resid}:{u.atoms[a].name}"
                    for a in atoms),
        r0=float(rest.bond.values[f]) / 10.0,
        thetaA0=float(rest.angles[1].values[f]),
        thetaB0=float(rest.angles[0].values[f]),
        phiA0=float(rest.dihedrals[2].values[f]),
        phiB0=float(rest.dihedrals[1].values[f]),
        phiC0=float(rest.dihedrals[0].values[f]),
        frame=f * step,
    )
    if not (b.atoms[2] == rest.bond.atomgroup.ix[1]
            and b.atoms[3] == rest.bond.atomgroup.ix[0]):
        raise RuntimeError(
            f"MDRestraintsGenerator atom order changed: bond atoms "
            f"{list(rest.bond.atomgroup.ix)} do not match P1/L1 {b.atoms[2]}/{b.atoms[3]}")
    if not force:
        check_collinear(b)
    return b
