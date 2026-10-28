import subprocess

import pytest

from automd.fe.abfe.rundir import load_rundir
from automd.fe.abfe.topology import Topology, ligand_charge, preprocess


@pytest.fixture
def pp(tmp_path, gmx):
    def _pp(rd):
        return preprocess(load_rundir(rd), tmp_path / rd.name)
    return _pp


def test_preprocess_both_routes(pp, lig_rundir, prolig_rundir, parmed_prolig_rundir):
    for rd in (lig_rundir, prolig_rundir, parmed_prolig_rundir):
        p = pp(rd)
        assert p.is_file()
        assert "#include" not in p.read_text()


def test_charge_neutral(pp, lig_rundir, prolig_rundir, parmed_prolig_rundir,
                        charged_lig_rundir):
    for rd in (lig_rundir, prolig_rundir, parmed_prolig_rundir, charged_lig_rundir):
        t = Topology.parse(pp(rd))
        assert abs(t.charge("A")) < 1e-3
        assert abs(t.charge("B")) < 1e-3


def test_ligand_charge(pp, lig_rundir, prolig_rundir, charged_lig_rundir):
    assert ligand_charge(Topology.parse(pp(lig_rundir)), "ACE") == -1
    assert ligand_charge(Topology.parse(pp(prolig_rundir)), "JZ4") == 0
    assert ligand_charge(Topology.parse(pp(charged_lig_rundir)), "LID") == 1


def test_parmed_route_names_and_repeated_molecules(pp, parmed_prolig_rundir):
    t = Topology.parse(pp(parmed_prolig_rundir))
    assert "system1" in t.moleculetypes and "JZ4" in t.moleculetypes
    names = [n for n, _ in t.molecules]
    assert names.count("NA") > 1 and isinstance(t.molecules, list)


def test_missing_ligand_moleculetype(pp, lig_rundir):
    t = Topology.parse(pp(lig_rundir))
    with pytest.raises(ValueError) as e:
        ligand_charge(t, "XYZ")
    assert "XYZ" in str(e.value) and "ACE" in str(e.value)


def test_b_state_defaults_to_a(pp, lig_rundir):
    a = Topology.parse(pp(lig_rundir)).moleculetypes["ACE"].atoms[0]
    assert (a.typeB, a.qB, a.massB) == (a.type, a.q, a.mass)


def test_write_verbatim_when_unedited(pp, lig_rundir, tmp_path):
    src = pp(lig_rundir)
    Topology.parse(src).write(tmp_path / "o.top")
    assert (tmp_path / "o.top").read_text() == src.read_text()


def test_edit_roundtrip(pp, lig_rundir, tmp_path):
    t = Topology.parse(pp(lig_rundir))
    t.moleculetypes["ACE"].atoms[0].qB = t.moleculetypes["ACE"].atoms[0].q + 0.5
    t.write(tmp_path / "e.top")
    t2 = Topology.parse(tmp_path / "e.top")
    assert abs(t2.charge("B") - t2.charge("A") - 0.5) < 1e-4


@pytest.mark.parametrize("which", ["lig_rundir", "parmed_prolig_rundir"])
def test_roundtrip_grompp(gmx, pp, tmp_path, request, which):
    rd = request.getfixturevalue(which)
    run = load_rundir(rd)
    t = Topology.parse(pp(rd))
    t.write(tmp_path / "rt.top")
    r = subprocess.run(
        [gmx, "grompp", "-f", str(rd / "MDP" / "em.mdp"), "-c", str(run.md_gro), "-r", str(run.md_gro),
         "-p", str(tmp_path / "rt.top"), "-o", str(tmp_path / "rt.tpr"),
         "-maxwarn", "2"], cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-1500:]
    assert (tmp_path / "rt.tpr").is_file()
