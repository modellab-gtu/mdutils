"""Task 7: abfe_setup CLI, per-leg Slurm array files, bin launcher."""
import json
import os
import re
import subprocess
from pathlib import Path

import pytest

from conftest import REPO_ROOT

LAUNCHER = REPO_ROOT / "bin" / "abfe_setup.sh"


def _setup(argv):
    from automd.fe.abfe.setup import main
    return main([str(a) for a in argv])


def _windows(leg_dir, rep):
    return sorted((Path(leg_dir) / f"rep{rep}").glob("lam_*"))


def _grompp_all(gmx, leg_dir, stages=("em", "nvt", "npt", "prod")):
    """grompp every window MDP of every repeat against the leg's FE system."""
    leg_dir = Path(leg_dir)
    n = 0
    for rep in sorted(leg_dir.glob("rep*")):
        for w in sorted(rep.glob("lam_*")):
            for st in stages:
                r = subprocess.run(
                    [gmx, "grompp", "-f", f"{st}.mdp", "-c", "../../start.gro",
                     "-r", "../../posre_ref.gro", "-p", "../../topol_fe.top",
                     "-n", "../../index.ndx", "-o", "_chk.tpr", "-po", "_chk.mdp",
                     "-maxwarn", "0"],
                    cwd=w, capture_output=True, text=True,
                    env={**os.environ, "GMX_MAXBACKUP": "-1"})
                assert r.returncode == 0, f"{w}/{st}.mdp\n{r.stderr[-3000:]}"
                for f in ("_chk.tpr", "_chk.mdp"):
                    (w / f).unlink()
                n += 1
    return n


def _listing(d):
    """Every entry (dotfiles included) of a run folder, for no-write checks."""
    return sorted(os.listdir(d))


def _ndx_groups(path):
    groups, cur = {}, None
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if s.startswith("["):
            cur = s.strip("[] ")
            groups[cur] = []
        elif s and cur is not None:
            groups[cur].extend(int(x) for x in s.split())
    return groups


@pytest.fixture(scope="module")
def solvent_only(tmp_path_factory):
    """Hydration setup on acetate (-1): one co-alchemical Na+."""
    from conftest import _rundir
    import shutil
    if shutil.which("gmx") is None:
        pytest.skip("gmx not on PATH")
    lig = _rundir("t4_ace_spce_run")
    before = _listing(lig)
    out = tmp_path_factory.mktemp("solv") / "ABFE"
    assert _setup(["--solvent", lig, "--ligand", "ACE", "--schedule", "test_short",
                   "--repeats", "2", "--outdir", out]) == 0
    assert _listing(lig) == before, "setup wrote into the run folder"
    return lig, out


def test_setup_solvent_only(solvent_only, gmx):
    lig, out = solvent_only
    solv = out / "solvent"
    assert [w.name for w in _windows(solv, 0)] == ["lam_00", "lam_01", "lam_02"]
    assert [w.name for w in _windows(solv, 1)] == ["lam_00", "lam_01", "lam_02"]
    meta = json.loads((out / "setup.json").read_text())
    assert meta["ligand_charge"] == -1
    assert len(meta["alch_ions"]["solvent"]) == 1
    assert meta["n_states"] == {"solvent": 3}
    assert meta["temperature"]["solvent"] == pytest.approx(300.0)
    assert meta["temperature_override"] is None
    assert "same rule in both legs" in meta["seed_rule"]
    assert not (out / "complex").exists()
    assert not (out / "boresch.json").exists()
    for f in ("topol_fe.top", "start.gro", "posre_ref.gro", "index.ndx"):
        assert (solv / f).is_file(), f
    assert "NA_ALCH" in (solv / "topol_fe.top").read_text()
    # seeds differ per repeat and state
    seeds = {(r, s) for r in (0, 1) for s in range(3)}
    got = {(int(k.split("/")[0][3:]), int(k.split("/")[1][4:])): v
           for k, v in meta["seeds"]["solvent"].items()}
    assert set(got) == seeds and len(set(got.values())) == 6
    nvt = (solv / "rep1" / "lam_02" / "nvt.mdp").read_text()
    assert re.search(r"gen-seed\s*=\s*%d\b" % got[(1, 2)], nvt)
    # tc-grps groups exist; Water_and_ions contains the alch ion (last atom)
    g = _ndx_groups(solv / "index.ndx")
    assert "ACE" in g and "Water_and_ions" in g
    natoms = int((solv / "start.gro").read_text().splitlines()[1])
    assert natoms in g["Water_and_ions"]
    assert set(g["ACE"]).isdisjoint(g["Water_and_ions"])
    assert len(g["ACE"]) + len(g["Water_and_ions"]) == natoms
    # sbatch: solvent on cws02, CPU module
    sb = (solv / "rep0" / "run_abfe.sbatch").read_text()
    assert "#SBATCH -w cws02" in sb and "gromacs/2026.1-tmpi\n" in sb
    assert "--array=0-2" in sb and "gpu" not in sb
    # every window MDP grompps with the FE system
    assert _grompp_all(gmx, solv) == 2 * 3 * 4


def test_setup_full(prolig_rundir, jz4_lig_rundir, tmp_path, gmx, capsys):
    before = {d: _listing(d) for d in (prolig_rundir, jz4_lig_rundir)}
    out = tmp_path / "ABFE"
    args = ["--complex", prolig_rundir, "--solvent", jz4_lig_rundir, "--ligand", "JZ4",
            "--schedule", "test_short", "--repeats", "2", "--outdir", out]
    # fixtures: complex md.mdp 310 K, solvent 300 K -> refused without --temperature
    with pytest.raises(SystemExit) as e:
        _setup(args)
    assert e.value.code != 0
    assert "--temperature" in capsys.readouterr().err
    assert not out.exists()
    assert _setup(args + ["--temperature", "300"]) == 0
    for d, lst in before.items():
        assert _listing(d) == lst, f"setup wrote into {d}"
    top = (out / "complex" / "topol_fe.top").read_text()
    assert "intermolecular_interactions" in top
    assert "intermolecular_interactions" not in (out / "solvent" / "topol_fe.top").read_text()
    b = json.loads((out / "boresch.json").read_text())
    assert len(b["atoms_0based"]) == 6 and b["dG_restr_off_kJ"] < 0
    meta = json.loads((out / "setup.json").read_text())
    assert meta["n_states"] == {"complex": 4, "solvent": 3}
    assert meta["ligand_charge"] == 0
    assert meta["alch_ions"] == {"complex": [], "solvent": []}
    assert meta["temperature"] == {"complex": 300.0, "solvent": 300.0}
    assert meta["md_temperature"] == {"complex": 310.0, "solvent": 300.0}
    assert meta["temperature_override"] == 300.0
    assert b["T"] == 300.0
    for mdp in out.glob("*/rep*/lam_*/*.mdp"):
        t = mdp.read_text()
        assert re.search(r"^ref-t\s*=\s*300 300$", t, re.M), mdp
        if mdp.name == "nvt.mdp":
            assert re.search(r"^gen-temp\s*=\s*300$", t, re.M), mdp
    assert "git_describe" in meta and "versions" in meta
    assert len(_windows(out / "complex", 1)) == 4
    sub = (out / "submit_abfe.sh").read_text()
    listed = re.findall(r"\S+/run_abfe\.sbatch", sub)
    assert sorted(listed) == sorted(f"{leg}/rep{r}/run_abfe.sbatch"
                                    for leg in ("complex", "solvent") for r in (0, 1))
    for p in listed:
        assert (out / p).is_file()
    sb = (out / "complex" / "rep0" / "run_abfe.sbatch").read_text()
    assert "#SBATCH -w cws01" in sb and "--gres=gpu:1" in sb and "-c 8" in sb
    assert "gromacs/2026.1-tmpi-cuda" in sb and "CUDA_DEVICE_ORDER=PCI_BUS_ID" in sb
    g = _ndx_groups(out / "complex" / "index.ndx")
    assert "Protein_JZ4" in g and "Water_and_ions" in g
    assert _grompp_all(gmx, out / "complex") == 2 * 4 * 4
    assert _grompp_all(gmx, out / "solvent") == 2 * 3 * 4


def test_setup_charged_cation(charged_lig_rundir, tmp_path, gmx):
    """LID (+1): one co-alchemical Cl-, every window grompps."""
    out = tmp_path / "ABFE"
    assert _setup(["--solvent", charged_lig_rundir, "--ligand", "LID",
                   "--schedule", "test_short", "--repeats", "1", "--outdir", out]) == 0
    meta = json.loads((out / "setup.json").read_text())
    assert meta["ligand_charge"] == 1
    assert len(meta["alch_ions"]["solvent"]) == 1
    assert "CL_ALCH" in (out / "solvent" / "topol_fe.top").read_text()
    assert _grompp_all(gmx, out / "solvent") == 1 * 3 * 4


def test_sbatch_restart(tmp_path):
    from automd.fe.abfe.slurm import write_leg_sbatch
    leg = tmp_path / "complex"
    (leg / "rep0").mkdir(parents=True)
    p = write_leg_sbatch(leg, "complex", 4, 0, "gtu", {"node": "cws01"})
    assert p == leg / "rep0" / "run_abfe.sbatch"
    t = p.read_text()
    assert "--array=0-3" in t
    assert "-cpi prod.cpt" in t
    assert "[[ -s prod.gro ]]" in t
    for st in ("em", "nvt", "npt"):
        assert f"[[ -s {st}.gro ]]" in t
    assert t.count("-r ../../posre_ref.gro") == 4 and t.count("-n ../../index.ndx") == 4
    assert 'lam_$(printf %02d $SLURM_ARRAY_TASK_ID)' in t
    r = subprocess.run(["bash", "-n", str(p)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_truba_sbatch(tmp_path):
    from automd.fe.abfe.slurm import write_leg_sbatch
    leg = tmp_path / "solvent"
    (leg / "rep1").mkdir(parents=True)
    opts = {"account": "akocak", "partition": "orfoz", "ntasks": 1, "omp": 16,
            "gmx_module": "apps/gromacs/2024.1-oneapi2024", "time": "3-00:00:00"}
    t = write_leg_sbatch(leg, "solvent", 3, 1, "truba", opts).read_text()
    assert "gmx_mpi" in t and "#SBATCH -C weka" in t and "mpirun" not in t
    assert "#SBATCH -A akocak" in t and "#SBATCH -n 1" in t and "#SBATCH -c 16" in t
    assert "module load apps/gromacs/2024.1-oneapi2024" in t
    assert "-w cws" not in t and "--array=0-2" in t
    opts["ntasks"] = 4
    t4 = write_leg_sbatch(leg, "solvent", 3, 1, "truba", opts).read_text()
    assert "mpirun -np $SLURM_NTASKS gmx_mpi mdrun" in t4


def test_truba_cli(solvent_only, tmp_path):
    lig, _ = solvent_only
    out = tmp_path / "T"
    assert _setup(["--solvent", lig, "--ligand", "ACE", "--schedule", "test_short",
                   "--repeats", "1", "--cluster", "truba", "--truba-account", "akocak",
                   "--outdir", out]) == 0
    t = (out / "solvent" / "rep0" / "run_abfe.sbatch").read_text()
    assert "gmx_mpi" in t and "-C weka" in t and "mpirun" not in t
    assert "LD_LIBRARY_PATH" not in t   # no python in the window jobs


def test_refuses_existing_outdir(solvent_only):
    lig, out = solvent_only
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    r = subprocess.run([str(LAUNCHER), "--solvent", str(lig), "--ligand", "ACE",
                        "--schedule", "test_short", "--outdir", str(out)],
                       capture_output=True, text=True, env=env)
    assert r.returncode != 0
    assert "--overwrite" in r.stderr
    assert (out / "setup.json").is_file()   # untouched


def test_refuses_outdir_in_run_folder(lig_rundir):
    with pytest.raises(SystemExit) as e:
        _setup(["--solvent", lig_rundir, "--ligand", "ACE", "--outdir",
                Path(lig_rundir) / "ABFE"])
    assert e.value.code != 0
    assert not (Path(lig_rundir) / "ABFE").exists()


def test_complex_requires_solvent(prolig_rundir, tmp_path):
    with pytest.raises(SystemExit) as e:
        _setup(["--complex", prolig_rundir, "--ligand", "JZ4", "--outdir", tmp_path / "X"])
    assert e.value.code != 0
    assert not (tmp_path / "X").exists()


def test_launcher_help():
    r = subprocess.run([str(LAUNCHER), "-h"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    for opt in ("--complex", "--solvent", "--ligand", "--replica", "--repeats",
                "--complex-ns", "--solvent-ns", "--schedule", "--cluster",
                "--complex-node", "--solvent-node", "--truba-account", "--truba-partition",
                "--truba-ntasks", "--truba-omp", "--truba-gmx-module", "--truba-time",
                "--outdir", "--submit", "--overwrite", "--force-restraint"):
        assert opt in r.stdout, opt


def test_index_unknown_tc_group(lig_rundir, tmp_path, gmx):
    from automd.fe.abfe.setup import SetupError, build_index
    with pytest.raises(SetupError, match="Membrane"):
        build_index(Path(lig_rundir) / "md_0.gro", tmp_path / "i.ndx",
                    ["ACE", "Membrane"], "ACE")


def test_boresch_name_check(tmp_path):
    from automd.fe.abfe.restraints import Boresch
    from automd.fe.abfe.setup import SetupError, check_boresch_names
    gro = tmp_path / "x.gro"
    rows = [(1, "LEU", "N"), (1, "LEU", "CA"), (1, "LEU", "C"),
            (164, "JZ4", "C1"), (164, "JZ4", "C2"), (164, "JZ4", "C3")]
    gro.write_text("t\n6\n" + "".join(f"{r:5d}{n:<5s}{a:>5s}{i + 1:5d}   0.000   0.000   0.000\n"
                                       for i, (r, n, a) in enumerate(rows)) + "   3.0 3.0 3.0\n")
    names = ("LEU1:C", "LEU1:CA", "LEU1:N", "JZ4164:C1", "JZ4164:C2", "JZ4164:C3")
    b = Boresch(atoms=(2, 1, 0, 3, 4, 5), names=names, r0=0.5, thetaA0=90, thetaB0=90,
                phiA0=0, phiB0=0, phiC0=0)
    check_boresch_names(b, gro)
    b2 = Boresch(atoms=(1, 2, 0, 3, 4, 5), names=names, r0=0.5, thetaA0=90, thetaB0=90,
                 phiA0=0, phiB0=0, phiC0=0)
    with pytest.raises(SetupError, match="boresch.json says LEU1:C"):
        check_boresch_names(b2, gro)
    # same names, other residue number -> refused
    b3 = Boresch(atoms=(2, 1, 0, 3, 4, 5), names=names[:3] + ("JZ4165:C1",) + names[4:],
                 r0=0.5, thetaA0=90, thetaB0=90, phiA0=0, phiB0=0, phiC0=0)
    with pytest.raises(SetupError, match="JZ4165:C1"):
        check_boresch_names(b3, gro)


def test_ligand_mismatch(prolig_rundir, jz4_lig_rundir, tmp_path, capsys):
    """A solvent ligand with other partial charges is refused (no override)."""
    import shutil
    solv = tmp_path / "jz4_solv_copy"
    shutil.copytree(jz4_lig_rundir, solv, symlinks=True)
    itp = solv / "JZ4_GMX.itp"
    txt = itp.read_text()
    # move 0.01 e from C1 to C2: integer total charge kept, atoms differ
    txt = txt.replace("JZ4     C1      1 -0.09510000", "JZ4     C1      1 -0.08510000", 1)
    txt = txt.replace("JZ4     C2      2 -0.16500000", "JZ4     C2      2 -0.17500000", 1)
    assert txt != itp.read_text()
    itp.write_text(txt)
    out = tmp_path / "ABFE"
    with pytest.raises(SystemExit) as e:
        _setup(["--complex", prolig_rundir, "--solvent", solv, "--ligand", "JZ4",
                "--schedule", "test_short", "--repeats", "1", "--temperature", "300",
                "--outdir", out])
    assert e.value.code != 0
    err = capsys.readouterr().err
    assert "#1 C1" in err and "#2 C2" in err and "charge" in err
    assert not out.exists()


def test_fixture_folders_untouched():
    """Runs last in this module: no fixture run folder has MDAnalysis offsets
    files or other leftovers from the FE tests."""
    from conftest import _DEFAULTS
    for d in _DEFAULTS.values():
        if os.path.isdir(d):
            stray = [f for f in os.listdir(d) if "_offsets" in f]
            assert not stray, f"{d}: {stray}"
