import os
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
LIGPARAM = REPO_ROOT / "ligparam" / "bin" / "ligparam"
EXAMPLES = Path(os.environ.get("LIGPARAM_EXAMPLES", "/cta/share/ligparam"))


@pytest.fixture
def example():
    def get(name):
        p = EXAMPLES / name
        if not p.is_file():
            pytest.skip(f"example ligand missing: {p}")
        return p
    return get


@pytest.fixture
def ligparam_bin():
    return LIGPARAM
