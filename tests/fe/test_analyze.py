"""Task 9: abfe_analyze (alchemlyb MBAR/TI, combine legs, warnings, result files)."""
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from automd.fe.abfe import analyze as an

DATA = Path(__file__).parent / "data"
REPO = Path(__file__).resolve().parents[2]
KCAL = 4.184


def _setup_json(d: Path, legs, n_states, T=300.0, repeats=1, prod_ns=0.05):
    meta = {"ligand": "LIG", "legs": legs, "ligand_charge": 0,
            "n_states": n_states, "temperature": {l: T for l in legs},
            "prod_ns": {l: prod_ns for l in legs},
            "options": {"repeats": repeats}}
    (d / "setup.json").write_text(json.dumps(meta))


def _toy_abfe(tmp_path, complex_too=False):
    abfe = tmp_path / "ABFE"
    abfe.mkdir()
    shutil.copytree(DATA / "toy_solvent", abfe / "solvent")
    legs, ns = ["solvent"], {"solvent": 3}
    if complex_too:
        shutil.copytree(DATA / "toy_complex", abfe / "complex")
        legs, ns = ["complex", "solvent"], {"complex": 4, "solvent": 3}
        (abfe / "boresch.json").write_text(json.dumps(
            {"T": 300.0, "dG_restr_off_kJ": -28.0, "dG_restr_off_kcal": -28.0 / KCAL}))
    _setup_json(abfe, legs, ns)
    return abfe


def _lr(**kw):
    base = dict(dg_mbar=10.0, err_mbar=0.1, dg_ti=10.0, err_ti=0.1, min_overlap=0.1,
                fwd_bwd_diff=0.0, n_states=5)
    base.update(kw)
    return an.LegResult(**base)


@pytest.mark.parametrize("key", ["Coulomb", "VDW"])
def test_alchemtest_benzene(key):
    """analyze_files == the plain alchemlyb pipeline (same subsampling, MBAR)."""
    import alchemlyb
    from alchemlyb.estimators import MBAR
    from alchemlyb.parsing.gmx import extract_u_nk
    from alchemlyb.preprocessing import decorrelate_u_nk
    from alchemtest.gmx import load_benzene

    files = load_benzene()["data"][key]
    T = 300.0
    res = an.analyze_files(files, T, equil_ns=0.0)
    direct = MBAR().fit(alchemlyb.concat(
        [decorrelate_u_nk(extract_u_nk(f, T=T), method="dE") for f in files]))
    kT = an.KB * T
    assert res.n_states == len(files)
    assert res.dg_mbar == pytest.approx(direct.delta_f_.iloc[0, -1] * kT, abs=1e-6)
    assert res.err_mbar == pytest.approx(direct.d_delta_f_.iloc[0, -1] * kT, abs=1e-6)
    assert math.isfinite(res.dg_ti) and math.isfinite(res.fwd_bwd_diff)
    assert 0 < res.min_overlap <= 1


def test_alchemtest_abfe_complex_multicomponent():
    """Complex leg with coul/vdw/bonded lambdas: TI sums all components, close to MBAR."""
    from alchemtest.gmx import load_ABFE

    files = load_ABFE()["data"]["complex"]
    res = an.analyze_files(files, 300.0, equil_ns=0.0)
    assert res.n_states == 30
    assert math.isfinite(res.dg_mbar) and math.isfinite(res.dg_ti)
    assert abs(res.dg_mbar - res.dg_ti) < 2.0  # kJ/mol


def test_toy_runs():
    res = an.analyze_leg(DATA / "toy_solvent" / "rep0", 300.0, equil_ns=0.0)
    assert res.n_states == 3
    for v in (res.dg_mbar, res.dg_ti, res.err_ti, res.min_overlap):
        assert math.isfinite(v)
    # 3-state toy has ~no overlap: pymbar reports a negative squared uncertainty (nan
    # error); that must surface as warnings, not as a crash
    w = an.warnings({("solvent", 0): res})
    assert any("overlap" in s for s in w) and any("convergence" in s for s in w)


def test_toy_complex_does_not_crash():
    """Parser/plumbing fixture: lam_03 has dH ~1e11 kJ/mol; only shape/no crash."""
    res = an.analyze_leg(DATA / "toy_complex" / "rep0", 300.0, equil_ns=0.0)
    assert res.n_states == 4
    assert math.isfinite(res.dg_ti)


def test_missing_windows(tmp_path):
    rep = tmp_path / "rep0"
    shutil.copytree(DATA / "toy_solvent" / "rep0", rep)
    assert an.missing_windows(rep, 3) == []
    (rep / "lam_01" / "prod.xvg").unlink()
    assert an.missing_windows(rep, 3) == [1]
    (rep / "lam_02" / "prod.xvg").write_text("")
    assert an.missing_windows(rep, 3) == [1, 2]


def test_missing_window_refused(tmp_path, capsys):
    abfe = _toy_abfe(tmp_path)
    (abfe / "solvent" / "rep0" / "lam_01" / "prod.xvg").unlink()
    with pytest.raises(SystemExit) as e:
        an.main([str(abfe), "--equil-ns", "0"])
    assert e.value.code != 0
    err = capsys.readouterr().err
    assert "[1]" in err and "solvent" in err
    assert not (abfe / "results" / "abfe_result.txt").exists()


def test_missing_window_allowed(tmp_path):
    abfe = _toy_abfe(tmp_path)
    (abfe / "solvent" / "rep0" / "lam_01" / "prod.xvg").unlink()
    assert an.main([str(abfe), "--equil-ns", "0", "--allow-missing"]) == 0
    txt = (abfe / "results" / "abfe_result.txt").read_text()
    assert "missing" in txt and "[1]" in txt


def test_combine_signs():
    setup = {"legs": ["complex", "solvent"]}
    legs = {("solvent", 0): _lr(dg_mbar=4.0, dg_ti=4.0),
            ("complex", 0): _lr(dg_mbar=30.0, dg_ti=30.0)}
    res = an.combine(setup, legs, -28.0)
    assert res["quantity"] == "dG_bind"
    assert res["per_repeat"][0]["MBAR"][0] == pytest.approx((4 - 30 + 28) / KCAL)
    assert res["per_repeat"][0]["TI"][0] == pytest.approx((4 - 30 + 28) / KCAL)
    assert res["mean"]["MBAR"][0] == pytest.approx((4 - 30 + 28) / KCAL)

    hyd = an.combine({"legs": ["solvent"]}, {("solvent", 0): _lr(dg_mbar=4.0, dg_ti=4.0)},
                     None)
    assert hyd["quantity"] == "dG_hyd"
    assert hyd["per_repeat"][0]["MBAR"][0] == pytest.approx(-4 / KCAL)


def test_combine_mean_sem():
    setup = {"legs": ["solvent"]}
    legs = {("solvent", r): _lr(dg_mbar=v, dg_ti=v) for r, v in enumerate([4.0, 8.0, 6.0])}
    res = an.combine(setup, legs, None)
    vals = np.array([-4.0, -8.0, -6.0]) / KCAL
    assert res["mean"]["MBAR"][0] == pytest.approx(vals.mean())
    assert res["mean"]["MBAR"][1] == pytest.approx(vals.std(ddof=1) / math.sqrt(3))


def test_warning_thresholds():
    w = an.warnings({("solvent", 0): _lr(min_overlap=0.02)})
    assert any("overlap" in s for s in w)
    assert not an.warnings({("solvent", 0): _lr(min_overlap=0.04)})
    # fwd/bwd and MBAR-TI: 0.5 kcal/mol
    assert any("forward" in s for s in an.warnings({("solvent", 0): _lr(fwd_bwd_diff=0.6 * KCAL)}))
    assert not an.warnings({("solvent", 0): _lr(fwd_bwd_diff=0.4 * KCAL)})
    assert any("TI" in s for s in an.warnings({("solvent", 0): _lr(dg_ti=10.0 + 0.6 * KCAL)}))
    assert not an.warnings({("solvent", 0): _lr(dg_ti=10.0 + 0.4 * KCAL)})
    # repeat spread: 1.0 kcal/mol
    legs = {("solvent", 0): _lr(dg_mbar=0.0, dg_ti=0.0),
            ("solvent", 1): _lr(dg_mbar=1.1 * KCAL, dg_ti=1.1 * KCAL)}
    comb = an.combine({"legs": ["solvent"]}, legs, None)
    assert any("spread" in s for s in an.warnings(legs, comb))
    legs[("solvent", 1)] = _lr(dg_mbar=0.9 * KCAL, dg_ti=0.9 * KCAL)
    comb = an.combine({"legs": ["solvent"]}, legs, None)
    assert not an.warnings(legs, comb)


def test_main_writes_results(tmp_path):
    abfe = _toy_abfe(tmp_path, complex_too=True)
    assert an.main([str(abfe), "--equil-ns", "0"]) == 0
    res = abfe / "results"
    txt = (res / "abfe_result.txt").read_text()
    assert "dG_bind = dG_solvent,decouple - dG_complex,decouple - dG_restr,off" in txt
    assert "dG_restr,off" in txt and "-6.692" in txt  # -28 kJ/mol in kcal/mol
    assert "WARNING" in txt
    import csv
    rows = list(csv.DictReader((res / "abfe_result.csv").open()))
    assert list(rows[0].keys()) == ["quantity", "repeat", "estimator", "value_kcal",
                                    "error_kcal"]
    q = {r["quantity"] for r in rows}
    assert {"dG_complex_decouple", "dG_solvent_decouple", "dG_restr_off", "dG_bind"} <= q
    for leg in ("complex", "solvent"):
        assert (res / f"overlap_{leg}_rep0.png").stat().st_size > 0
        assert (res / f"convergence_{leg}_rep0.png").stat().st_size > 0


def test_main_hydration(tmp_path):
    abfe = _toy_abfe(tmp_path)
    assert an.main([str(abfe), "--equil-ns", "0"]) == 0
    txt = (abfe / "results" / "abfe_result.txt").read_text()
    assert "dG_hyd = -dG_solvent,decouple" in txt
    assert "dG_bind" not in txt.split("Sign convention")[0]


def test_wrapper_script(tmp_path):
    abfe = _toy_abfe(tmp_path)
    r = subprocess.run(["bash", str(REPO / "bin" / "abfe_analyze.sh"), str(abfe),
                        "--equil-ns", "0"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (abfe / "results" / "abfe_result.csv").is_file()
