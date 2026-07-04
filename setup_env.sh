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

CONDA_PREFIX="$("$CONDA_CMD" env list | awk -v env="$ENV_NAME" '$1 == env {print $NF}')"
[[ -n "$CONDA_PREFIX" ]] || die "Could not locate environment prefix for '${ENV_NAME}'"

BIN="$CONDA_PREFIX/bin"
PY="$BIN/python"

ok()   { echo "  [OK]   $*"; }
warn() { echo "  [WARN] $*"; }

# obabel: version flag is -V (capital V)
if [[ -x "$BIN/obabel" ]]; then
    ok "obabel $("$BIN/obabel" -V 2>&1 | head -1)"
else
    warn "obabel not found"
fi

# pdbfixer: no __version__; check openmm version instead
if "$PY" -c "import pdbfixer, openmm; print('pdbfixer OK, openmm', openmm.__version__)" 2>/dev/null | grep -q OK; then
    ok "$("$PY" -c "import pdbfixer, openmm; print('pdbfixer + openmm', openmm.__version__)" 2>/dev/null)"
else
    warn "pdbfixer or openmm import failed"
fi

# gmx
if [[ -x "$BIN/gmx" ]]; then
    ok "$("$BIN/gmx" --version 2>&1 | grep 'GROMACS version' | xargs)"
else
    warn "gmx not found"
fi

# AmberTools binaries exit non-zero on -h, so just check they exist
for cmd in antechamber parmchk2 tleap amb2gro_top_gro.py; do
    if [[ -x "$BIN/$cmd" ]]; then
        ok "$cmd"
    else
        warn "$cmd not found"
    fi
done

# parmed Python library (used by amb2gro_top_gro.py)
if ver="$("$PY" -c "import parmed; print(parmed.__version__)" 2>/dev/null)"; then
    ok "parmed $ver"
else
    warn "parmed import failed"
fi

# numpy sanity check (catches broken installations like missing __init__.py)
if ver="$("$PY" -c "import numpy; print(numpy.__version__)" 2>/dev/null)"; then
    ok "numpy $ver"
else
    warn "numpy import failed — environment may be corrupted"
fi

info "Done. Activate with:  conda activate ${ENV_NAME}"
