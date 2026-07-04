#!/usr/bin/env bash
# Creates and populates the mdutils conda environment.
# Run with:  bash setup_env.sh
# Re-running is safe; conda skips already-installed packages.

set -euo pipefail

ENV_NAME="mdutils"
CONDA_CMD="${CONDA_EXE:-conda}"

# ── helper ────────────────────────────────────────────────────────────────────
info()  { echo "[setup] $*"; }
die()   { echo "[ERROR] $*" >&2; exit 1; }

command -v "$CONDA_CMD" &>/dev/null || die "conda not found. Install Miniconda/Anaconda first."

# ── 1. Create environment ─────────────────────────────────────────────────────
info "Step 1: creating conda environment '${ENV_NAME}' with Python 3.12"
"$CONDA_CMD" create -y -n "$ENV_NAME" python=3.12 \
    --channel conda-forge --override-channels

# ── 2. OpenBabel ──────────────────────────────────────────────────────────────
info "Step 2: installing OpenBabel"
"$CONDA_CMD" install -y -n "$ENV_NAME" openbabel \
    --channel conda-forge --override-channels

# ── 3. PDBFixer (+ OpenMM) ────────────────────────────────────────────────────
info "Step 3: installing PDBFixer (pulls in OpenMM)"
"$CONDA_CMD" install -y -n "$ENV_NAME" pdbfixer \
    --channel conda-forge --override-channels

# ── 4. GROMACS ────────────────────────────────────────────────────────────────
info "Step 4: installing GROMACS"
"$CONDA_CMD" install -y -n "$ENV_NAME" gromacs \
    --channel conda-forge --override-channels

# ── 5. AmberTools ─────────────────────────────────────────────────────────────
# conda-forge currently ships AmberTools 24.8 (the latest stable build).
# AmberTools 25 will appear on conda-forge when it is packaged upstream.
info "Step 5: installing AmberTools (latest available on conda-forge)"
"$CONDA_CMD" install -y -n "$ENV_NAME" ambertools \
    --channel conda-forge --override-channels

# ── Verify ────────────────────────────────────────────────────────────────────
info "Verifying installation"

CONDA_PREFIX="$("$CONDA_CMD" info --json | python3 -c "
import sys, json
envs = json.load(sys.stdin)['envs']
for e in envs:
    if e.endswith('${ENV_NAME}'):
        print(e); break
")"

[[ -n "$CONDA_PREFIX" ]] || die "Could not locate environment prefix for '${ENV_NAME}'"

check() {
    local cmd="$1"
    if "$CONDA_PREFIX/bin/$cmd" --version &>/dev/null 2>&1; then
        echo "  [OK] $cmd $("$CONDA_PREFIX/bin/$cmd" --version 2>&1 | head -1)"
    else
        echo "  [WARN] $cmd not found or failed"
    fi
}

check obabel
"$CONDA_PREFIX/bin/python" -c "import pdbfixer; print('  [OK] pdbfixer', pdbfixer.__version__)" 2>/dev/null \
    || echo "  [WARN] pdbfixer import failed"
check gmx
check antechamber
check tleap
check parmchk2
"$CONDA_PREFIX/bin/python" -c "import parmed; print('  [OK] parmed (amb2gro_top_gro.py dependency)', parmed.__version__)" 2>/dev/null \
    || echo "  [WARN] parmed import failed"

info "Done. Activate with:  conda activate ${ENV_NAME}"
