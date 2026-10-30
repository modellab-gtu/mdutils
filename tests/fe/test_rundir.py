import pytest

from automd.fe.abfe.rundir import load_rundir, read_mdp, write_mdp


def test_read_mdp_normalises(tmp_path):
    p = tmp_path / "a.mdp"
    p.write_text("; header\ntc_grps = A B ; trailing\ntc-grps2=X\nREF_T = 310 310\n\n")
    d = read_mdp(p)
    assert d["tc-grps"] == "A B"
    assert d["tc-grps2"] == "X"
    assert d["ref-t"] == "310 310"
    assert list(d) == ["tc-grps", "tc-grps2", "ref-t"]
    q = tmp_path / "b.mdp"
    write_mdp(d, q)
    assert read_mdp(q) == d


def test_load_lig_rundir(lig_rundir):
    r = load_rundir(lig_rundir)
    assert r.top.name == "lig_solv_gmx.top"
    assert r.temperature == 300.0
    assert r.replica == 0
    assert r.md_tpr.name == "md_0.tpr"


def test_load_prolig_rundir(prolig_rundir):
    r = load_rundir(prolig_rundir)
    assert r.top.name == "topol.top"
    assert r.ndx is not None and r.ndx.exists()


def test_load_parmed_prolig_rundir(parmed_prolig_rundir):
    r = load_rundir(parmed_prolig_rundir)
    assert r.top.name == "complex_gmx.top"


def test_missing_replica(tmp_path):
    (tmp_path / "MDP").mkdir()
    (tmp_path / "MDP" / "md.mdp").write_text("ref_t = 300\n")
    with pytest.raises(FileNotFoundError) as e:
        load_rundir(tmp_path)
    assert "md_0" in str(e.value)
    assert "run the MD replica first" in str(e.value)
