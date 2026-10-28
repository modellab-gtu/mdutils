import re
import subprocess

import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.lib.distances import distance_array

from automd.fe.abfe.ions import (apply_coalchemical, check_charges, choose_ions,
                                 fe_state_charge)
from automd.fe.abfe.mdp import fe_mdps
from automd.fe.abfe.rundir import load_rundir, write_mdp
from automd.fe.abfe.schedule import default_schedule
from automd.fe.abfe.topology import Topology, ligand_charge, preprocess


def _setup(rd, tmp_path, ligand):
    run = load_rundir(rd)
    top = Topology.parse(preprocess(run, tmp_path))
    return run, top, mda.Universe(str(run.md_gro)), ligand_charge(top, ligand)


def _lig_atoms(top, u, ligand):
    i = 0
    for n, c in top.molecules:
        k = len(top.moleculetypes[n].atoms)
        if n == ligand:
            return list(range(i, i + k * c))
        i += k * c


def test_choose_for_acetate(lig_rundir, tmp_path, gmx):
    run, top, u, q = _setup(lig_rundir, tmp_path, "ACE")
    assert q == -1
    # The fixture's only Na+ is a contact ion (0.22 nm) in a 3.3 nm box, so the
    # production 1.5 nm default cannot be met here (see test_too_few_ions);
    # exercise selection with a cutoff the fixture can satisfy.
    idx = choose_ions(u, top, "ACE", q, min_lig_dist=0.15)
    assert len(idx) == 1
    assert u.atoms[idx[0]].name.upper().startswith("NA")
    d = distance_array(u.atoms.positions[idx], u.atoms.positions[_lig_atoms(top, u, "ACE")],
                       box=u.dimensions).min() / 10
    assert d >= 0.15


def test_min_distance_1p5_synthetic(tmp_path):
    """1.5 nm cutoff honoured with minimum image: ion at x=8.0 in a 9 nm box is 1.0 nm
    from the ligand at x=0 through the periodic boundary and must be rejected."""
    p, gro = _synthetic(tmp_path)
    L = gro.read_text().splitlines()
    L[2 + 1] = L[2 + 1][:20] + f"{8.0:8.3f}" + L[2 + 1][28:]   # NA #1 at 8.0 (1.0 nm PBC)
    L[2 + 2] = L[2 + 2][:20] + f"{4.0:8.3f}" + L[2 + 2][28:]   # NA #2 at 4.0
    L[2 + 3] = L[2 + 3][:20] + f"{2.0:8.3f}" + L[2 + 3][28:]   # NA #3 at 2.0 (>=1.5)
    gro.write_text("\n".join(L) + "\n")
    u = mda.Universe(str(gro))
    top = Topology.parse(p)
    assert choose_ions(u, top, "LIG", -1) == [2]   # x=4.0 is farthest from the ligand
    with pytest.raises(RuntimeError, match="--salt"):
        choose_ions(u, top, "LIG", -3)


def test_neutral_returns_empty(lig_rundir, tmp_path, gmx):
    run, top, u, q = _setup(lig_rundir, tmp_path, "ACE")
    assert choose_ions(u, top, "ACE", 0) == []


def test_choose_cation_charged_ligand_and_protein_order(charged_lig_rundir, tmp_path, gmx):
    run, top, u, q = _setup(charged_lig_rundir, tmp_path, "LID")
    assert q == 1
    idx = choose_ions(u, top, "LID", q, min_lig_dist=0.5)
    assert len(idx) == 1 and u.atoms[idx[0]].name.upper().startswith("CL")


def test_protein_ordering_farthest_first(prolig_rundir, tmp_path, gmx):
    run, top, u, _ = _setup(prolig_rundir, tmp_path, "JZ4")
    # pretend the ligand is -2: two cations, farthest from protein first
    idx = choose_ions(u, top, "JZ4", -2)
    assert len(idx) == 2 and all(u.atoms[i].name.upper().startswith(("NA", "K")) for i in idx)
    prot = u.select_atoms("protein").positions
    d = [distance_array(u.atoms.positions[[i]], prot, box=u.dimensions).min() for i in idx]
    assert d[0] >= d[1]
    lig = _lig_atoms(top, u, "JZ4")
    assert distance_array(u.atoms.positions[idx], u.atoms.positions[lig],
                          box=u.dimensions).min() / 10 >= 1.5


def test_apply_charges(lig_rundir, tmp_path, gmx):
    run, top, u, q = _setup(lig_rundir, tmp_path, "ACE")
    idx = choose_ions(u, top, "ACE", q, min_lig_dist=0.15)
    apply_coalchemical(top, run.md_gro, tmp_path / "start.gro", idx)
    assert abs(top.charge("A")) < 1e-3
    assert abs(top.charge("B") - (-1)) < 1e-3
    assert abs(fe_state_charge(top, "ACE", "B")) < 1e-3
    check_charges(top, "ACE")
    assert top.molecules[-1] == ("NA_ALCH", 1)
    alch = top.moleculetypes["NA_ALCH"].atoms[0]
    src = top.moleculetypes["NA"].atoms[0]
    assert alch.qB == 0 and alch.q == src.q and alch.typeB == src.type == alch.type
    assert alch.massB == src.mass
    assert (tmp_path / "posre_ref.gro").read_text() == (tmp_path / "start.gro").read_text()
    txt = "\n".join(top.lines)
    assert re.search(r"\[ position_restraints \]\s*\n;[^\n]*\n\s*1\s+2\s+1\s+0\.5\d*\s+1000", txt)
    # last coordinate is the chosen ion, same xyz as in the input
    g = mda.Universe(str(tmp_path / "start.gro"))
    assert np.allclose(g.atoms[-1].position, u.atoms[idx[0]].position, atol=0.1)
    assert g.atoms[-1].name.upper().startswith("NA")


def test_charged_lid(charged_lig_rundir, tmp_path, gmx):
    run, top, u, q = _setup(charged_lig_rundir, tmp_path, "LID")
    idx = choose_ions(u, top, "LID", q, min_lig_dist=0.5)
    apply_coalchemical(top, run.md_gro, tmp_path / "start.gro", idx)
    assert dict(top.molecules).get("CL_ALCH") == 1
    assert abs(top.charge("A")) < 1e-3
    assert abs(top.charge("B") - 1) < 1e-3
    assert abs(fe_state_charge(top, "LID", "B")) < 1e-3


def _synthetic(tmp_path):
    """Repeated-name block NA 3, one ion in the middle of the block."""
    top = tmp_path / "s.top"
    top.write_text("""[ defaults ]
1 1 no 1 1
[ atomtypes ]
 X 6 12.0 0.0 A 0.3 0.1
 N 11 22.99 0.0 A 0.2 0.1
 C 17 35.45 0.0 A 0.4 0.1
[ moleculetype ]
LIG 3
[ atoms ]
 1 X 1 LIG X1 1 -1.0 12.0
[ moleculetype ]
NA 1
[ atoms ]
 1 N 1 NA NA 1 1.0 22.99
[ moleculetype ]
CL 1
[ atoms ]
 1 C 1 CL CL 1 -1.0 35.45
[ system ]
s
[ molecules ]
LIG 1
NA 3
CL 2
""")
    atoms = [("LIG", "X1", 0.0), ("NA", "NA", 1.0), ("NA", "NA", 2.0), ("NA", "NA", 3.0),
             ("CL", "CL", 4.0), ("CL", "CL", 5.0)]
    gro = ["synthetic", f"{len(atoms):5d}"]
    for i, (r, n, x) in enumerate(atoms, 1):
        gro.append(f"{i:5d}{r:<5s}{n:>5s}{i:5d}{x:8.3f}{0:8.3f}{0:8.3f}")
    gro.append("   9.00000   9.00000   9.00000")
    (tmp_path / "s.gro").write_text("\n".join(gro) + "\n")
    return top, tmp_path / "s.gro"


def test_split_repeated_block(tmp_path):
    p, gro = _synthetic(tmp_path)
    top = Topology.parse(p)
    apply_coalchemical(top, gro, tmp_path / "o.gro", [2])   # middle NA (x = 2.0)
    assert top.molecules == [("LIG", 1), ("NA", 2), ("CL", 2), ("NA_ALCH", 1)]
    xs = [float(l[20:28]) for l in (tmp_path / "o.gro").read_text().splitlines()[2:8]]
    assert xs == [0.0, 1.0, 3.0, 4.0, 5.0, 2.0]
    top.write(tmp_path / "o.top")
    t2 = Topology.parse(tmp_path / "o.top")
    assert t2.molecules == top.molecules
    assert abs(t2.charge("A")) < 1e-6 and abs(t2.charge("B") + 1) < 1e-6


def test_apply_rejects_non_ion(tmp_path):
    p, gro = _synthetic(tmp_path)
    with pytest.raises(ValueError):
        apply_coalchemical(Topology.parse(p), gro, tmp_path / "o.gro", [0, 0])
    p2 = Topology.parse(p)
    # atom 0 is a LIG with one atom here, so use a bad index instead
    with pytest.raises(ValueError):
        apply_coalchemical(p2, gro, tmp_path / "o.gro", [99])


def test_too_few_ions(tmp_path, lig_rundir, gmx):
    run, top, u, q = _setup(lig_rundir, tmp_path, "ACE")
    # no ion of the opposite sign satisfies the distance: impossible cutoff
    with pytest.raises(RuntimeError, match="--salt"):
        choose_ions(u, top, "ACE", q, min_lig_dist=1000.0)
    # a fixture with ions removed
    p, gro = _synthetic(tmp_path)
    t = Topology.parse(p)
    t.set_molecules([("LIG", 1), ("CL", 2)])
    g = mda.Universe(str(_gro_without(gro, tmp_path)))
    with pytest.raises(RuntimeError, match="--salt"):
        choose_ions(g, t, "LIG", -1)


def _gro_without(gro, tmp_path):
    L = gro.read_text().splitlines()
    keep = [L[2]] + L[6:8]
    out = tmp_path / "noion.gro"
    out.write_text("\n".join([L[0], f"{len(keep):5d}"] + keep + [L[-1]]) + "\n")
    return out


@pytest.mark.parametrize("fixture,ligand,mind", [("lig_rundir", "ACE", 0.15),
                                                  ("charged_lig_rundir", "LID", 0.5)])
def test_grompp_coalch(request, tmp_path, gmx, fixture, ligand, mind):
    rd = request.getfixturevalue(fixture)
    run, top, u, q = _setup(rd, tmp_path, ligand)
    idx = choose_ions(u, top, ligand, q, min_lig_dist=mind)
    gro = tmp_path / "start.gro"
    apply_coalchemical(top, run.md_gro, gro, idx)
    top.write(tmp_path / "fe.top")
    # The run folder's index.ndx is stale after the reorder (ions moved to the
    # end), so build the two tc-grps groups from the topology.
    g1, g2 = run.md_mdp["tc-grps"].split()
    lig = set(_lig_atoms(top, u, ligand))
    n = len(u.atoms)
    def block(name, ids):
        ids = [str(i + 1) for i in ids]
        return [f"[ {name} ]"] + [" ".join(ids[k:k + 15]) for k in range(0, len(ids), 15)]
    ndx = block(g1, sorted(lig)) + block(g2, [i for i in range(n) if i not in lig])
    (tmp_path / "index.ndx").write_text("\n".join(ndx) + "\n")
    sched = default_schedule("solvent")
    last = sched.n_states - 1
    for state in (0, last):
        mdps = fe_mdps(run.md_mdp, sched, ligand, state, 0.01, 1)
        write_mdp(mdps["prod"], tmp_path / f"prod{state}.mdp")
        r = subprocess.run(
            ["gmx", "grompp", "-f", f"prod{state}.mdp", "-c", "start.gro", "-r", "posre_ref.gro",
             "-p", "fe.top", "-n", "index.ndx", "-o", f"s{state}.tpr", "-maxwarn", "0"],
            cwd=tmp_path, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-3000:]
        assert "net charge" not in r.stderr.lower()


def test_parmed_repeated_block_split(parmed_prolig_rundir, tmp_path, gmx):
    """ParmEd-style repeated [ molecules ] entries: moving a mid-block ion splits counts."""
    run, top, u, _ = _setup(parmed_prolig_rundir, tmp_path, "JZ4")
    before = list(top.molecules)
    names = [n for n, _ in before]
    assert names.count("NA") > 1
    n_total_ions = sum(c for n, c in before if n == "NA")
    idx = choose_ions(u, top, "JZ4", -1)       # fake q=-1 just to drive selection
    apply_coalchemical(top, run.md_gro, tmp_path / "o.gro", idx)
    top.write(tmp_path / "o.top")
    t2 = Topology.parse(tmp_path / "o.top")
    assert t2.molecules == top.molecules and t2.molecules[-1] == ("NA_ALCH", 1)
    assert sum(c for n, c in t2.molecules if n == "NA") == n_total_ions - 1
    # same multiset of molecules otherwise, and atom count unchanged
    n_at = lambda t: sum(c * len(t.moleculetypes[n].atoms) for n, c in t.molecules)
    assert n_at(t2) == len(u.atoms) == len(mda.Universe(str(tmp_path / "o.gro")).atoms)
    g = mda.Universe(str(tmp_path / "o.gro"))
    assert np.allclose(g.atoms[-1].position, u.atoms[idx[0]].position, atol=0.1)
