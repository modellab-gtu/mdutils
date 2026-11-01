"""ligparam --opt (geometry before charges).

(a) default path: LIG.mol2/.frcmod/.itp byte-identical to the reference tool (the live
    /cta/capps/envs/ligparam/opt/ligparam before the switch-over; set LIGPARAM_REF to its kept copy after it)
(b) --opt aimnet2: charges sum to the formal charge, sqm ran AM1 with maxcyc=0, json records the model
(c) --opt ani2x on a charged ligand falls back to hf, and says so
(d) --opt none keeps the input coordinates
plus the hf job files and --geom handling.
"""
import json
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest

LIGPARAM_BIN = Path(__file__).resolve().parents[2] / "ligparam" / "bin" / "ligparam"

REF = Path(os.environ.get("LIGPARAM_REF", "/cta/capps/envs/ligparam/opt/ligparam/bin/ligparam"))
LIGANDS = ["aspirin.sdf", "lidocaineH.sdf", "acetate.sdf", "mephos.sdf"]
FORMAL = {"aspirin.sdf": 0, "lidocaineH.sdf": 1, "acetate.sdf": -1, "mephos.sdf": -2}
# the routes of /cta/share/ligparam/run_examples.sh
EXAMPLE_ROUTES = [("lidocaineH.sdf", ["--charge-model", "bcc"]), ("acetate.sdf", ["--charge-model", "nagl"]),
                  ("mephos.sdf", ["--ff", "openff"])]


def run(tool, *args):
    return subprocess.run([str(tool), *map(str, args)], capture_output=True, text=True)


def mol2_atoms(path):
    txt = Path(path).read_text()
    rows = txt.split("@<TRIPOS>ATOM\n", 1)[1].split("@<TRIPOS>", 1)[0].splitlines()
    return [r.split() for r in rows if r.strip()]


def sdf_xyz(path):
    from rdkit import Chem
    m = Chem.MolFromMolFile(str(path), removeHs=False)
    return m.GetConformer().GetPositions()


# ---------------------------------------------------------------- (a) default path unchanged
def sqm_routine(outdir):
    """sqm picks its diagonaliser by timing it at start-up (diag_routine auto): two runs of the SAME
    antechamber command can use different routines and differ in the 5th-6th charge decimal."""
    f = Path(outdir) / "sqm.out"
    if not f.exists():
        return None
    m = re.search(r"Using (\w+) routine \(diag_routine=(\d+)\)", f.read_text())
    return m.group(0) if m else None


def itp_sections(path):
    """The OpenFF route (Interchange) writes the .itp's [ pairs ] lines in a run-to-run varying order (the live
    tool does too: three runs, three md5s). [ pairs ] is compared as sorted lines, every other section as is."""
    return [sorted(b.splitlines()) if b.startswith("[ pairs ]") else b
            for b in re.split(r"(?m)^(?=\[ )", Path(path).read_text())]


def _compare(tmp_path, lig, args, files, attempts=6, unordered=()):
    """Reference and repo tool on the same input; byte-identical outputs required whenever both sqm runs
    chose the same diagonaliser (retried up to `attempts` times until they do)."""
    if not REF.exists():
        pytest.skip(f"reference ligparam missing: {REF}")
    seen = []
    for k in range(attempts):
        a, b = tmp_path / f"ref{k}", tmp_path / f"new{k}"
        r1 = run(REF, "-i", lig, "-o", a, *args)
        r2 = run(LIGPARAM_BIN, "-i", lig, "-o", b, *args)
        assert r1.returncode == 0, r1.stdout + r1.stderr
        assert r2.returncode == 0, r2.stdout + r2.stderr
        # identical external commands (the out dir is the only difference in the log)
        cmds = [[ln.replace(str(d), "OUT") for ln in (d / "ligparam.log").read_text().splitlines()
                 if ln.startswith("$ ")] for d in (a, b)]
        assert cmds[0] == cmds[1]
        ra, rb = sqm_routine(a), sqm_routine(b)
        seen.append((ra, rb))
        if ra != rb:
            continue
        for f in files:
            if f in unordered:
                assert itp_sections(a / f) == itp_sections(b / f), f"{f} differs ({a} vs {b})"
            else:
                assert (a / f).read_bytes() == (b / f).read_bytes(), f"{f} differs ({a} vs {b})"
        return
    pytest.fail(f"sqm never chose the same diagonaliser in {attempts} attempts: {seen}")


@pytest.mark.parametrize("name", LIGANDS)
def test_default_byte_identical(tmp_path, example, name):
    _compare(tmp_path, example(name), [], ["LIG.mol2", "LIG.frcmod", "LIG.itp"])


@pytest.mark.parametrize("name,args", EXAMPLE_ROUTES)
def test_example_routes_byte_identical(tmp_path, example, name, args):
    if "openff" in args:
        _compare(tmp_path, example(name), args, ["LIG.itp"], unordered=("LIG.itp",))
    else:
        _compare(tmp_path, example(name), args, ["LIG.mol2", "LIG.frcmod", "LIG.itp"])


def test_default_json_records_am1(tmp_path, example):
    r = run(LIGPARAM_BIN, "-i", example("acetate.sdf"), "-o", tmp_path / "o")
    assert r.returncode == 0, r.stdout + r.stderr
    assert json.loads((tmp_path / "o" / "ligparam.json").read_text())["opt"] == "am1"


def test_json_opt_null_without_am1_step(tmp_path, example):
    """nagl (and user / openff) charges have no AM1 step: opt is recorded as null."""
    r = run(LIGPARAM_BIN, "-i", example("acetate.sdf"), "-o", tmp_path / "o", "--charge-model", "nagl")
    assert r.returncode == 0, r.stdout + r.stderr
    assert json.loads((tmp_path / "o" / "ligparam.json").read_text())["opt"] is None


@pytest.mark.parametrize("args", [["--charge-model", "nagl"], ["--charge-model", "user"], ["--ff", "openff"]])
def test_opt_refused_without_am1_step(tmp_path, example, args):
    inp = tmp_path / "x.mol2" if "user" in args else example("aspirin.sdf")  # refused before the input is read
    r = run(LIGPARAM_BIN, "-i", inp, "-o", tmp_path / "o", "--opt", "none", *args)
    assert r.returncode != 0 and "--opt" in r.stderr and "Traceback" not in r.stderr


# ---------------------------------------------------------------- (b) NNP pre-optimisation
@pytest.mark.parametrize("name", ["aspirin.sdf", "acetate.sdf"])
def test_opt_aimnet2(tmp_path, example, name):
    out = tmp_path / "o"
    r = run(LIGPARAM_BIN, "-i", example(name), "-o", out, "--opt", "aimnet2")
    assert r.returncode == 0, r.stdout + r.stderr
    q = sum(float(f[8]) for f in mol2_atoms(out / "LIG.mol2"))
    assert abs(q - FORMAL[name]) < 1e-5
    sqm_in = (out / "sqm.in").read_text()
    assert "maxcyc=0" in sqm_in and "qm_theory='AM1'" in sqm_in
    assert re.search(r"Selected Hamiltonian \*\*\* *\n\| QMMM: AM1", (out / "sqm.out").read_text())
    rec = json.loads((out / "ligparam.json").read_text())
    assert rec["opt"] == "aimnet2"
    assert abs(rec["charge_sum"] - FORMAL[name]) < 1e-4
    # charges came from the relaxed geometry (charge_geom.sdf, read by antechamber) ...
    x0 = sdf_xyz(example(name))
    xq = sdf_xyz(out / "charge_geom.sdf")
    assert np.abs(xq - x0).max() > 1e-2
    assert "-i charge_geom.sdf" in (out / "ligparam.log").read_text()
    assert (out / "opt_aimnet2" / "final.xyz").is_file()
    # ... and every output keeps the input coordinates
    x1 = np.array([[float(v) for v in f[2:5]] for f in mol2_atoms(out / "LIG.mol2")])
    assert np.abs(x1 - x0).max() < 1e-3
    assert np.abs(sdf_xyz(out / "input.sdf") - x0).max() < 1e-3
    import parmed as pmd
    s = pmd.load_file(str(out / "LIG.prmtop"), str(out / "LIG.inpcrd"))
    assert np.abs(s.coordinates - x0).max() < 1e-3


# ---------------------------------------------------------------- (c) fallback to hf
def test_opt_ani2x_charged_falls_back_to_hf(tmp_path, example):
    out = tmp_path / "o"
    r = run(LIGPARAM_BIN, "-i", example("acetate.sdf"), "-o", out, "--opt", "ani2x")
    assert r.returncode == 0, r.stdout + r.stderr
    msg = r.stdout + r.stderr
    assert "ani2x is neutral-only" in msg and "falling back to HF" in msg
    assert (out / "opt_hf.gjf").is_file() and (out / "run_opt.sh").is_file()
    assert "--geom" in msg
    assert not (out / "LIG.mol2").exists()
    side = json.loads((out / "opt_hf.json").read_text())
    assert side["opt_requested"] == "ani2x" and "neutral-only" in side["opt_note"]


def fake_gaussian_log(path, mol_path):
    """A minimal Gaussian opt log (Input orientation of the given SDF, Optimization completed, Normal termination)."""
    from rdkit import Chem
    m = Chem.MolFromMolFile(str(mol_path), removeHs=False)
    dash = " " + "-" * 69
    rows = [f"{i + 1:7d}{a.GetAtomicNum():11d}{0:12d}    {x:12.6f}{y:12.6f}{z:12.6f}"
            for i, (a, (x, y, z)) in enumerate(zip(m.GetAtoms(), m.GetConformer().GetPositions()))]
    Path(path).write_text("\n".join([" Entering Gaussian System", "                          Input orientation:",
                                     dash, " Center     Atomic      Atomic             Coordinates (Angstroms)",
                                     " Number     Number       Type             X           Y           Z", dash]
                                    + rows + [dash, "    Optimization completed.",
                                              " Normal termination of Gaussian 16 at Sun Nov  1 12:00:00 2026."]) + "\n")


def test_fallback_provenance_survives_geom_rerun(tmp_path, example):
    """ani2x -> hf fallback on the Slurm path: the --geom rerun records the requested model and the reason."""
    out = tmp_path / "o"
    r = run(LIGPARAM_BIN, "-i", example("acetate.sdf"), "-o", out, "--opt", "ani2x")
    assert r.returncode == 0, r.stdout + r.stderr
    fake_gaussian_log(out / "opt_hf.log", example("acetate.sdf"))
    r = run(LIGPARAM_BIN, "-i", example("acetate.sdf"), "-o", out, "--opt", "hf", "--geom", out / "opt_hf.log")
    assert r.returncode == 0, r.stdout + r.stderr
    rec = json.loads((out / "ligparam.json").read_text())
    assert rec["opt"] == "hf" and rec["opt_requested"] == "ani2x"
    assert "neutral-only" in rec["opt_note"] and rec["geom"] == str(out / "opt_hf.log")
    assert "neutral-only" in r.stdout


def test_missing_covparam_clean_error(tmp_path, example):
    """covparam not importable (e.g. not readable for this user): a clean ligparam error, no traceback."""
    script = Path(__file__).resolve().parents[2] / "ligparam" / "ligparam.py"
    env = os.environ.get("LIGPARAM_ENV", "/cta/capps/envs/ligparam")
    code = ("import sys, runpy\n"
            "class Block:\n"
            "    def find_spec(self, name, path=None, target=None):\n"
            "        if name == 'covparam' or name.startswith('covparam.'): raise ImportError('blocked for the test')\n"
            "sys.meta_path.insert(0, Block())\n"
            f"sys.argv = ['ligparam', '-i', {str(example('aspirin.sdf'))!r}, '-o', {str(tmp_path / 'o')!r}, '--opt', 'aimnet2']\n"
            f"runpy.run_path({str(script)!r}, run_name='__main__')\n")
    r = subprocess.run([f"{env}/bin/python", "-c", code], capture_output=True, text=True,
                       env={**os.environ, "AMBERHOME": env, "PATH": f"{env}/bin:" + os.environ["PATH"],
                            "PYTHONNOUSERSITE": "1"})
    assert r.returncode != 0
    assert "Traceback" not in r.stderr
    assert "--opt aimnet2 needs covparam" in r.stderr


# ---------------------------------------------------------------- (d) none
def test_opt_none_keeps_coordinates(tmp_path, example):
    out = tmp_path / "o"
    r = run(LIGPARAM_BIN, "-i", example("aspirin.sdf"), "-o", out, "--opt", "none")
    assert r.returncode == 0, r.stdout + r.stderr
    x0 = sdf_xyz(example("aspirin.sdf"))
    x1 = np.array([[float(v) for v in f[2:5]] for f in mol2_atoms(out / "LIG.mol2")])
    assert np.abs(x1 - x0).max() < 1e-3
    assert "maxcyc=0" in (out / "sqm.in").read_text()
    assert json.loads((out / "ligparam.json").read_text())["opt"] == "none"


# ---------------------------------------------------------------- hf job files, --geom
def test_opt_hf_writes_job(tmp_path, example):
    out = tmp_path / "o"
    r = run(LIGPARAM_BIN, "-i", example("aspirin.sdf"), "-o", out, "--opt", "hf")
    assert r.returncode == 0, r.stdout + r.stderr
    gjf = (out / "opt_hf.gjf").read_text()
    assert "#P HF/6-31G* Opt=(MaxCycles=500,CalcFC,MaxStep=10) NoSymm" in gjf
    assert "%nprocshared=16" in gjf and "\n0 1\n" in gjf
    assert "ModRedundant" not in gjf
    sh = (out / "run_opt.sh").read_text()
    assert sh.startswith("#!/bin/bash -l\n")
    for s in ("#SBATCH -w cws02", "#SBATCH -c 16", "#SBATCH --export=NONE", "module load gaussian/g16-avx2",
              'GAUSS_SCRDIR=/tmp/$SLURM_JOB_ID', '"$g16root/g16/g16" opt_hf.gjf'):
        assert s in sh, s
    assert os.access(out / "run_opt.sh", os.X_OK)
    assert f"--geom {out}/opt_hf.log" in r.stdout
    assert not (out / "LIG.mol2").exists()


def test_geom_xyz(tmp_path, example):
    """--geom with an .xyz (same atom order): used as the charge geometry (AM1 maxcyc=0); outputs keep the input."""
    from rdkit import Chem
    src = example("aspirin.sdf")
    m = Chem.MolFromMolFile(str(src), removeHs=False)
    xyz = sdf_xyz(src) + np.array([3.0, -2.0, 1.0])  # a shifted copy: the outputs must not follow it
    lines = [str(m.GetNumAtoms()), "test"] + [f"{a.GetSymbol()} {x:.6f} {y:.6f} {z:.6f}"
                                             for a, (x, y, z) in zip(m.GetAtoms(), xyz)]
    g = tmp_path / "g.xyz"; g.write_text("\n".join(lines) + "\n")
    out = tmp_path / "o"
    r = run(LIGPARAM_BIN, "-i", src, "-o", out, "--opt", "hf", "--geom", g)
    assert r.returncode == 0, r.stdout + r.stderr
    x1 = np.array([[float(v) for v in f[2:5]] for f in mol2_atoms(out / "LIG.mol2")])
    assert np.abs(x1 - sdf_xyz(src)).max() < 1e-3
    assert np.abs(sdf_xyz(out / "charge_geom.sdf") - xyz).max() < 1e-3
    assert "maxcyc=0" in (out / "sqm.in").read_text()
    rec = json.loads((out / "ligparam.json").read_text())
    assert rec["opt"] == "hf" and rec["geom"] == str(g)


def test_geom_wrong_atoms(tmp_path, example):
    g = tmp_path / "g.xyz"; g.write_text("1\nx\nC 0 0 0\n")
    r = run(LIGPARAM_BIN, "-i", example("aspirin.sdf"), "-o", tmp_path / "o", "--opt", "hf", "--geom", g)
    assert r.returncode != 0 and "atoms" in r.stderr


@pytest.mark.parametrize("opt", ["am1", "none", "aimnet2"])
def test_geom_only_with_hf(tmp_path, example, opt):
    g = tmp_path / "g.xyz"; g.write_text("1\nx\nC 0 0 0\n")
    r = run(LIGPARAM_BIN, "-i", example("aspirin.sdf"), "-o", tmp_path / "o", "--geom", g, "--opt", opt)
    assert r.returncode != 0 and "--geom" in r.stderr and "Traceback" not in r.stderr


def test_opt_openff_refused(tmp_path, example):
    r = run(LIGPARAM_BIN, "-i", example("aspirin.sdf"), "-o", tmp_path / "o", "--ff", "openff", "--opt", "aimnet2")
    assert r.returncode != 0 and "--opt" in r.stderr


HF_LOG = Path(os.environ.get("LIGPARAM_HF_LOG", "/cta/users/akocak/maint/ligparam_opt_work/hfjob/asp/opt_hf.log"))


def test_geom_gaussian_log(tmp_path, example):
    """--geom with the real cws02 HF/6-31G* optimisation of aspirin (skipped if that log is not there)."""
    if not HF_LOG.is_file():
        pytest.skip(f"no Gaussian log at {HF_LOG}")
    out = tmp_path / "o"
    r = run(LIGPARAM_BIN, "-i", example("aspirin.sdf"), "-o", out, "--opt", "hf", "--geom", HF_LOG)
    assert r.returncode == 0, r.stdout + r.stderr
    assert abs(sum(float(f[8]) for f in mol2_atoms(out / "LIG.mol2"))) < 1e-5
    x1 = np.array([[float(v) for v in f[2:5]] for f in mol2_atoms(out / "LIG.mol2")])
    assert np.abs(x1 - sdf_xyz(example("aspirin.sdf"))).max() < 1e-3
    assert json.loads((out / "ligparam.json").read_text())["opt"] == "hf"
