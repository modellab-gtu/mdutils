import os
import shutil
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

_DEFAULTS = {
    "t4_ace_spce_run": "/cta/users/akocak/maint/phase2_test/t4_ace_spce_run",
    "t1_run": "/cta/users/akocak/maint/phase2_test/t1_run",
    "t2_run": "/cta/users/akocak/maint/phase2_test/t2_run",
    "t3_lid_truba_run": "/cta/users/akocak/maint/phase2_test/t3_lid_truba_run",
    "jz4_lig_run": "/cta/users/akocak/maint/abfe_test/fixtures/jz4_lig_run",
}


def _rundir(name):
    base = os.environ.get("MDUTILS_FE_TESTDATA")
    p = Path(base) / name if base else Path(_DEFAULTS[name])
    if not p.is_dir():
        pytest.skip(f"fixture folder missing: {p}")
    return p


@pytest.fixture
def lig_rundir():
    """Acetate ACE (-1) in SPC/E, ParmEd route, lig_solv_gmx.top."""
    return _rundir("t4_ace_spce_run")


@pytest.fixture
def prolig_rundir():
    """T4L + JZ4, pdb2gmx route, topol.top."""
    return _rundir("t1_run")


@pytest.fixture
def parmed_prolig_rundir():
    """T4L + JZ4, ParmEd route, complex_gmx.top (protein moleculetype system1)."""
    return _rundir("t2_run")


@pytest.fixture
def charged_lig_rundir():
    """Lidocaine LID (+1), topol.top, 1 CL."""
    return _rundir("t3_lid_truba_run")


@pytest.fixture
def jz4_lig_rundir():
    """JZ4 in water."""
    return _rundir("jz4_lig_run")


@pytest.fixture
def gmx():
    """Path to gmx; skip unless on PATH (module load gromacs/2026.1-tmpi)."""
    exe = shutil.which("gmx")
    if exe is None:
        pytest.skip("gmx not on PATH (module load gromacs/2026.1-tmpi)")
    return exe


def ro_traj(rd, dest, r=0):
    """Symlink md_<r>.tpr/xtc of a fixture run folder into ``dest`` and return them.

    MDAnalysis writes .<xtc>_offsets.npz/.lock next to the trajectory it opens;
    opening the symlinks keeps the fixture folders unmodified."""
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    out = []
    for ext in ("tpr", "xtc"):
        link = dest / f"md_{r}.{ext}"
        if not link.exists():
            link.symlink_to((Path(rd) / f"md_{r}.{ext}").resolve())
        out.append(link)
    return tuple(out)
