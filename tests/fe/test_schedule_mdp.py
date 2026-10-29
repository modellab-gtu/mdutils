import pytest

from automd.fe.abfe.mdp import fe_mdps, write_window_mdps
from automd.fe.abfe.rundir import load_rundir, read_mdp
from automd.fe.abfe.schedule import (Schedule, default_schedule, load_schedule,
                                     validate)
from pathlib import Path

SHORT = Path(__file__).resolve().parents[2] / "automd/fe/abfe/schedules/test_short.yaml"


def test_default_complex():
    s = default_schedule("complex")
    assert s.n_states == 36
    st = lambda i: (s.restraint[i], s.coul[i], s.vdw[i])
    assert st(0) == (0, 0, 0)
    assert st(5) == (1, 0, 0)
    assert st(15) == (1, 1, 0)
    assert st(35) == (1, 1, 1)
    validate(s)


def test_default_solvent():
    s = default_schedule("solvent")
    assert s.n_states == 31
    assert all(r == 0 for r in s.restraint)
    validate(s)


def test_validate_rejects_decreasing():
    with pytest.raises(ValueError):
        validate(Schedule([0, 1], [0, 1], [1, 0]))
    with pytest.raises(ValueError):
        validate(Schedule([0, 1], [0], [0, 1]))
    with pytest.raises(ValueError):
        validate(Schedule([0, 2], [0, 1], [0, 1]))


def test_short_yaml():
    assert load_schedule(SHORT, "complex").n_states == 4
    assert load_schedule(SHORT, "solvent").n_states == 3


def test_fe_mdp_keys(lig_rundir):
    base = load_rundir(lig_rundir).md_mdp
    s = default_schedule("solvent")
    m = fe_mdps(base, s, "ACE", 7, 3.0, 4242)
    assert set(m) == {"em", "nvt", "npt", "prod"}
    p = m["prod"]
    assert p["calc-lambda-neighbors"] == "-1"
    assert p["couple-moltype"] == "ACE"
    assert p["init-lambda-state"] == "7"
    assert p["rcoulomb"] == "1.2"
    assert p["nsteps"] == "1500000"
    assert p["nstxout-compressed"] == "5000"
    assert m["em"]["integrator"] == "steep"
    assert m["em"]["nsteps"] == "5000"
    assert "constraints" not in m["em"] and "pcoupl" not in m["em"]
    assert m["nvt"]["gen-seed"] == "4242" and m["nvt"]["ld-seed"] == "4242"
    assert m["nvt"]["pcoupl"] == "no" and m["nvt"]["continuation"] == "no"
    assert m["npt"]["pcoupl"].lower() == "c-rescale"
    assert m["npt"]["continuation"] == "yes"
    assert m["nvt"]["nsteps"] == "50000" and m["npt"]["nsteps"] == "100000"
    for k in m.values():
        assert k["free-energy"] == "yes"
        assert k["sc-alpha"] == "0.5" and k["couple-lambda0"] == "vdw-q"
        assert len(k["vdw-lambdas"].split()) == 31
        assert k["restraint-lambdas"].split()[0] == "0.0000"
        assert k["bonded-lambdas"] == k["restraint-lambdas"]


def test_write_window_mdps(lig_rundir, tmp_path):
    base = load_rundir(lig_rundir).md_mdp
    paths = write_window_mdps(tmp_path, base, default_schedule("solvent"),
                              "ACE", 3, 1.0, 1)
    assert [p.name for p in paths] == ["em.mdp", "nvt.mdp", "npt.mdp", "prod.mdp"]
    assert read_mdp(paths[3])["init-lambda-state"] == "3"


def test_barostat_always_c_rescale():
    """The spec mandates C-rescale in NPT/prod even if the run used another barostat."""
    from automd.fe.abfe.mdp import fe_mdps
    from automd.fe.abfe.schedule import default_schedule
    base = {"pcoupl": "Parrinello-Rahman", "tau-p": "2.0", "ref-t": "300 300",
            "tc-grps": "LIG Water_and_ions", "tau-t": "1 1"}
    m = fe_mdps(base, default_schedule("solvent"), "LIG", 0, 1.0, 1)
    assert m["npt"]["pcoupl"] == "C-rescale" and m["prod"]["pcoupl"] == "C-rescale"


def test_no_pull_keys_any_leg():
    from automd.fe.abfe.mdp import fe_mdps
    from automd.fe.abfe.schedule import default_schedule
    base = {"ref-t": "300 300", "tc-grps": "LIG Water_and_ions", "tau-t": "1 1"}
    m = fe_mdps(base, default_schedule("solvent"), "LIG", 3, 1.0, 1)
    assert not any(k.startswith("pull") for d in m.values() for k in d)
