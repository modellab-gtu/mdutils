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
    tool does too: three runs, three md5s). Compared per [ section ] as sorted lines."""
    import re as _re
    return [sorted(b.splitlines()) for b in _re.split(r"(?m)^(?=\[ )", Path(path).read_text())]


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
    # the geometry really moved, and it was put back onto the input pose (aligned)
    x0 = sdf_xyz(example(name))
    x1 = np.array([[float(v) for v in f[2:5]] for f in mol2_atoms(out / "LIG.mol2")])
    rmsd = np.sqrt(((x1 - x0) ** 2).sum(1).mean())
    assert 1e-3 < rmsd < 1.0
    assert (out / "opt_aimnet2" / "final.xyz").is_file()


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
    """--geom with an .xyz (same atom order): used as the charge geometry, AM1 not re-optimised."""
    from rdkit import Chem
    src = example("aspirin.sdf")
    m = Chem.MolFromMolFile(str(src), removeHs=False)
    xyz = sdf_xyz(src) + np.array([3.0, -2.0, 1.0])  # a shifted copy: the alignment must undo the shift
    lines = [str(m.GetNumAtoms()), "test"] + [f"{a.GetSymbol()} {x:.6f} {y:.6f} {z:.6f}"
                                             for a, (x, y, z) in zip(m.GetAtoms(), xyz)]
    g = tmp_path / "g.xyz"; g.write_text("\n".join(lines) + "\n")
    out = tmp_path / "o"
    r = run(LIGPARAM_BIN, "-i", src, "-o", out, "--opt", "hf", "--geom", g)
    assert r.returncode == 0, r.stdout + r.stderr
    x1 = np.array([[float(v) for v in f[2:5]] for f in mol2_atoms(out / "LIG.mol2")])
    assert np.abs(x1 - sdf_xyz(src)).max() < 1e-3
    assert "maxcyc=0" in (out / "sqm.in").read_text()
    rec = json.loads((out / "ligparam.json").read_text())
    assert rec["opt"] == "hf" and rec["geom"] == str(g)


def test_geom_wrong_atoms(tmp_path, example):
    g = tmp_path / "g.xyz"; g.write_text("1\nx\nC 0 0 0\n")
    r = run(LIGPARAM_BIN, "-i", example("aspirin.sdf"), "-o", tmp_path / "o", "--opt", "hf", "--geom", g)
    assert r.returncode != 0 and "atoms" in r.stderr


def test_geom_needs_opt_not_am1(tmp_path, example):
    g = tmp_path / "g.xyz"; g.write_text("1\nx\nC 0 0 0\n")
    r = run(LIGPARAM_BIN, "-i", example("aspirin.sdf"), "-o", tmp_path / "o", "--geom", g)
    assert r.returncode != 0 and "--geom" in r.stderr


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
    assert json.loads((out / "ligparam.json").read_text())["opt"] == "hf"
