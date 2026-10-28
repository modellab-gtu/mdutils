"""Task 6: Boresch restraint selection, topology block and analytic correction."""
import json
import math
import subprocess

import numpy as np
import pytest

from automd.fe.abfe.restraints import (Boresch, boresch_correction, check_collinear,
                                       find_boresch, intermolecular_block, to_json)

KCAL = 4.184


def _ref(**kw):
    d = dict(atoms=(0, 1, 2, 3, 4, 5), names=("P3", "P2", "P1", "L1", "L2", "L3"),
             r0=0.5, thetaA0=90.0, thetaB0=90.0, phiA0=0.0, phiB0=0.0, phiC0=0.0)
    d.update(kw)
    return Boresch(**d)


def test_defaults():
    b = _ref()
    assert b.kr == 4184.0 and b.ktheta == 41.84 and b.kphi == 41.84


def test_correction_reference():
    dg = boresch_correction(_ref(), 298.15)
    assert dg == pytest.approx(-28.5807, abs=1e-3)
    assert dg / KCAL == pytest.approx(-6.8310, abs=1e-3)


def test_correction_matches_mdrg():
    """MDRG BoreschRestraint._analytical_energy (via standard_state) returns
    kcal/mol with the same sign (dG_off < 0); force_constant=10 kcal/mol/rad^2
    means k_angle = k_dihedral = 10 kcal/mol/rad^2 = 41.84 kJ/mol/rad^2 and
    k_bond = 100*10 kcal/mol/nm^2 = 4184 kJ/mol/nm^2. Bond values are in A."""
    from types import SimpleNamespace
    from MDRestraintsGenerator.datatypes import BoreschRestraint
    fake = SimpleNamespace(bond=SimpleNamespace(values=np.array([5.0])),
                           angles=[SimpleNamespace(values=np.array([90.0])),
                                   SimpleNamespace(values=np.array([90.0]))],
                           min_frame=0)
    # standard_state(frame, force_constant, temperature) dispatches to this
    mdrg_kcal = BoreschRestraint._analytical_energy(fake, 0, 10.0, 298.15)
    ours = boresch_correction(_ref(), 298.15)
    assert ours == pytest.approx(mdrg_kcal * KCAL, abs=0.01)
    # non-trivial geometry too
    fake.bond.values[0] = 7.3
    fake.angles[0].values[0] = 63.0
    fake.angles[1].values[0] = 121.0
    mdrg2 = BoreschRestraint._analytical_energy(fake, 0, 10.0, 310.0)
    ours2 = boresch_correction(_ref(r0=0.73, thetaA0=121.0, thetaB0=63.0), 310.0)
    assert ours2 == pytest.approx(mdrg2 * KCAL, abs=0.01)


def test_correction_scaling():
    """dG_off = -kT ln(... sqrt(k) ...): 4x k_r lowers dG_off by kT ln 2."""
    weak = boresch_correction(_ref(), 300.0)
    strong = boresch_correction(_ref(kr=4 * 4184.0), 300.0)
    kT = 0.0083144626 * 300.0
    assert strong - weak == pytest.approx(-kT * math.log(2.0), abs=1e-9)


@pytest.mark.parametrize("kw", [dict(thetaA0=175.0), dict(thetaB0=5.0),
                                dict(thetaA0=170.5)])
def test_collinear_guard(kw):
    with pytest.raises(ValueError):
        check_collinear(_ref(**kw))


def test_collinear_ok():
    check_collinear(_ref(thetaA0=160.0, thetaB0=20.0))


def test_block_format():
    b = _ref(atoms=(10, 11, 12, 2000, 2001, 2002), r0=0.61234, thetaA0=100.5,
             thetaB0=80.25, phiA0=-170.0, phiB0=45.0, phiC0=179.9)
    txt = intermolecular_block(b)
    assert "[ intermolecular_interactions ]" in txt
    sec, rows = None, {"bonds": [], "angles": [], "dihedrals": []}
    for line in txt.splitlines():
        s = line.split(";")[0].strip()
        if not s:
            continue
        if s.startswith("["):
            sec = s.strip("[] ")
            continue
        rows[sec].append(s.split())
    assert len(rows["bonds"]) == 1 and len(rows["angles"]) == 2
    assert len(rows["dihedrals"]) == 3
    bond = rows["bonds"][0]
    assert bond[:3] == ["13", "2001", "6"]          # P1, L1 (1-based)
    assert float(bond[3]) == pytest.approx(0.61234, abs=1e-5)
    assert float(bond[4]) == 0.0 and float(bond[5]) == float(bond[3])
    assert float(bond[6]) == pytest.approx(4184.0)
    exp_ang = [(["12", "13", "2001"], 100.5), (["13", "2001", "2002"], 80.25)]
    for row, (ids, val) in zip(rows["angles"], exp_ang):
        assert row[:3] == ids and row[3] == "1"
        assert float(row[4]) == pytest.approx(val) and float(row[6]) == pytest.approx(val)
        assert float(row[5]) == 0.0 and float(row[7]) == pytest.approx(41.84)
    exp_dih = [(["11", "12", "13", "2001"], -170.0), (["12", "13", "2001", "2002"], 45.0),
               (["13", "2001", "2002", "2003"], 179.9)]
    for row, (ids, val) in zip(rows["dihedrals"], exp_dih):
        assert row[:4] == ids and row[4] == "2"
        assert float(row[5]) == pytest.approx(val) and float(row[7]) == pytest.approx(val)
        assert float(row[6]) == 0.0 and float(row[8]) == pytest.approx(41.84)


def test_to_json(tmp_path):
    b = _ref()
    p = tmp_path / "boresch.json"
    to_json(b, 298.15, p)
    d = json.loads(p.read_text())
    assert d["dG_restr_off_kJ"] == pytest.approx(-28.5807, abs=1e-3)
    assert d["dG_restr_off_kcal"] == pytest.approx(-6.8310, abs=1e-3)
    assert d["T"] == 298.15 and d["atoms_1based"] == [1, 2, 3, 4, 5, 6]
    assert d["r0_nm"] == 0.5


@pytest.fixture(scope="module")
def prolig_boresch():
    from pathlib import Path
    import os
    base = os.environ.get("MDUTILS_FE_TESTDATA")
    rd = Path(base) / "t1_run" if base else Path("/cta/users/akocak/maint/phase2_test/t1_run")
    if not rd.is_dir():
        pytest.skip(f"fixture folder missing: {rd}")
    return rd, find_boresch(rd / "md_0.tpr", rd / "md_0.xtc", "JZ4", 310.0)


def test_find_on_prolig(prolig_boresch):
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import calc_angles, calc_bonds, calc_dihedrals
    rd, b = prolig_boresch
    u = mda.Universe(str(rd / "md_0.tpr"), str(rd / "md_0.xtc"))
    ag = u.atoms[list(b.atoms)]
    assert all(a.name in ("N", "CA", "C") and a.residue.resname != "JZ4" for a in ag[:3])
    assert all(a.resname == "JZ4" and not a.name.startswith("H") for a in ag[3:])
    assert 0.3 < b.r0 < 1.5
    check_collinear(b)
    # equilibrium values reproduce the chosen frame with OUR atom order
    u.trajectory[b.frame]
    x = [a.position for a in ag]
    box = u.dimensions
    assert calc_bonds(x[2], x[3], box) / 10 == pytest.approx(b.r0, abs=1e-4)
    assert math.degrees(calc_angles(x[1], x[2], x[3], box)) == pytest.approx(b.thetaA0, abs=1e-3)
    assert math.degrees(calc_angles(x[2], x[3], x[4], box)) == pytest.approx(b.thetaB0, abs=1e-3)
    for i, phi in enumerate((b.phiA0, b.phiB0, b.phiC0)):
        assert math.degrees(calc_dihedrals(*x[i:i + 4], box=box)) == pytest.approx(phi, abs=1e-3)
    assert -60.0 < boresch_correction(b, 310.0) < 0.0


def test_grompp_with_block(prolig_boresch, tmp_path, gmx):
    """intermolecular_interactions appended after [ molecules ] of the define-free
    processed topology grompps with an FE MDP (states 0 and 5; the block is
    switched by bonded-lambdas)."""
    from automd.fe.abfe.mdp import fe_mdps
    from automd.fe.abfe.rundir import load_rundir, write_mdp
    from automd.fe.abfe.schedule import default_schedule
    from automd.fe.abfe.topology import Topology, preprocess
    rd, b = prolig_boresch
    run = load_rundir(rd)
    top = Topology.parse(preprocess(run, tmp_path))
    top.append_section(intermolecular_block(b))
    top.write(tmp_path / "fe.top")
    txt = (tmp_path / "fe.top").read_text()
    assert txt.index("[ molecules ]") < txt.index("[ intermolecular_interactions ]")
    sched = default_schedule("complex")
    for state in (0, 5):
        mdps = fe_mdps(run.md_mdp, sched, "JZ4", state, 0.01, 1)
        write_mdp(mdps["prod"], tmp_path / f"prod{state}.mdp")
        r = subprocess.run(
            [gmx, "grompp", "-f", f"prod{state}.mdp", "-c", str(run.md_gro), "-p", "fe.top",
             "-n", str(rd / "index.ndx"), "-o", f"s{state}.tpr", "-maxwarn", "0"],
            cwd=tmp_path, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-3000:]


def _restraint_energy(b, x, box):
    """1/2 k (x-x0)^2 summed over the 6 Boresch terms (kJ/mol), MDAnalysis geometry."""
    from MDAnalysis.lib.distances import calc_angles, calc_bonds, calc_dihedrals
    r = calc_bonds(x[2], x[3], box) / 10.0
    th = (calc_angles(x[1], x[2], x[3], box), calc_angles(x[2], x[3], x[4], box))
    ph = [calc_dihedrals(*x[i:i + 4], box=box) for i in range(3)]
    wrap = lambda d: (d + math.pi) % (2 * math.pi) - math.pi
    e = 0.5 * b.kr * (r - b.r0) ** 2
    e += sum(0.5 * b.ktheta * (t - math.radians(t0)) ** 2
             for t, t0 in zip(th, (b.thetaA0, b.thetaB0)))
    e += sum(0.5 * b.kphi * wrap(p - math.radians(p0)) ** 2
             for p, p0 in zip(ph, (b.phiA0, b.phiB0, b.phiC0)))
    return e


def test_restraint_follows_restraint_schedule(prolig_boresch, tmp_path, gmx):
    """GROMACS physics check: rerun md_0.xtc with the FE topology at complex
    state 5 (restraint fully on). The intermolecular bond/angle/dihedral terms
    are *bonded* interactions, so their dH/dlambda appears in the bonded
    component; it must equal our 1/2 k dx^2 at every frame (units, atom order,
    dihedral sign and the A=0/B=full switch all checked at once). The
    restraint energy must be ~0 at the frame the reference values came from."""
    import MDAnalysis as mda
    from automd.fe.abfe.mdp import fe_mdps
    from automd.fe.abfe.rundir import load_rundir, write_mdp
    from automd.fe.abfe.schedule import default_schedule
    from automd.fe.abfe.topology import Topology, preprocess
    rd, b = prolig_boresch
    run = load_rundir(rd)
    top = Topology.parse(preprocess(run, tmp_path))
    top.append_section(intermolecular_block(b))
    top.write(tmp_path / "fe.top")
    sched = default_schedule("complex")
    assert sched.restraint[5] == 1.0 and sched.coul[5] == 0.0
    m = fe_mdps(run.md_mdp, sched, "JZ4", 5, 0.01, 1)["prod"]
    m.update({"nstdhdl": "1", "calc-lambda-neighbors": "0"})
    write_mdp(m, tmp_path / "p5.mdp")
    for cmd in ([gmx, "grompp", "-f", "p5.mdp", "-c", str(run.md_gro), "-p", "fe.top",
                 "-n", str(rd / "index.ndx"), "-o", "p5.tpr", "-maxwarn", "0"],
                [gmx, "mdrun", "-s", "p5.tpr", "-rerun", str(rd / "md_0.xtc"),
                 "-deffnm", "rr", "-nt", "4"]):
        r = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-3000:]
    xvg = (tmp_path / "rr.xvg").read_text().splitlines()
    legends = [l.split('legend "', 1)[1] for l in xvg if l.startswith("@ s") and "legend" in l]
    col = 1 + next(i for i, s in enumerate(legends) if "bonded-lambda" in s)
    rows = [list(map(float, l.split())) for l in xvg if l and l[0] not in "#@"]
    u = mda.Universe(str(rd / "md_0.tpr"), str(rd / "md_0.xtc"))
    ag = u.atoms[list(b.atoms)]
    assert len(rows) == u.trajectory.n_frames
    for ts, row in zip(u.trajectory, rows):
        mine = _restraint_energy(b, ag.positions, u.dimensions)
        assert row[col] == pytest.approx(mine, rel=2e-3, abs=2e-3), ts.frame
        if ts.frame == b.frame:
            assert row[col] < 1e-2
        else:
            assert row[col] > 1.0
