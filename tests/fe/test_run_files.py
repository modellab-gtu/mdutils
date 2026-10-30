import os
import subprocess

import pytest

from conftest import REPO_ROOT

FIXTURE = os.environ.get("MDUTILS_FE_JZ4_LIG_UNPATCHED",
                         "/cta/users/akocak/maint/abfe_test/fixtures/jz4_lig")


def test_water_and_ions_group_created(gmx):
    if not os.path.isdir(FIXTURE):
        pytest.skip(f"fixture folder missing: {FIXTURE}")
    r = subprocess.run(["bash", "-e", str(REPO_ROOT / "tests/run_files/test_water_and_ions.sh"), FIXTURE],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "PASS" in r.stdout
