# Installation Guide

## Requirements

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- Linux x86_64 (the pipeline has not been tested on macOS or Windows)
- Gaussian 16 (required only for `--charge-method resp`; all other charge methods work without it)

## Quick setup

```bash
bash setup_env.sh
conda activate mdutils
```

The script installs everything in order into a fresh `mdutils` environment. Re-running is safe.

---

## Step-by-step

### 1. Create the environment

```bash
conda create -y -n mdutils python=3.12 --channel conda-forge --override-channels
```

### 2. OpenBabel

Used for ligand file format conversion (SDF, MOL2, PDB, MOL, etc.).

```bash
conda install -y -n mdutils openbabel --channel conda-forge --override-channels
```

### 3. PDBFixer

Repairs missing residues, missing heavy atoms, and adds hydrogens at a chosen pH.
Pulls in OpenMM as a dependency.

```bash
conda install -y -n mdutils pdbfixer --channel conda-forge --override-channels
```

### 4. GROMACS

Molecular dynamics engine. The conda-forge build is CPU-only (no GPU).
For GPU-accelerated production runs, install GROMACS from source with CUDA support.

```bash
conda install -y -n mdutils gromacs --channel conda-forge --override-channels
```

Verify:

```bash
conda run -n mdutils gmx --version
```

### 5. AmberTools

Provides `antechamber`, `parmchk2`, `tleap`, and `parmed` (via `amb2gro_top_gro.py`).

> **Note:** conda-forge currently packages AmberTools 24.8. AmberTools 25 will appear
> once it is packaged upstream. The scripts in this repository are compatible with
> both 24.x and 25.x.

```bash
conda install -y -n mdutils ambertools --channel conda-forge --override-channels
```

---

## Gaussian 16 (RESP charges only)

Gaussian 16 is commercial software and cannot be installed via conda.
It must be installed separately and its executables (`g16`, `formchk`) must be on `PATH`.

The script uses G16 only when `--charge-method resp` is given.
All other charge methods (`bcc`, `abcg2`, `mul`, `cm2`, …) work without G16.

---

## Verification

```bash
conda activate mdutils
obabel --version
python -c "import pdbfixer; print('pdbfixer', pdbfixer.__version__)"
gmx --version | grep "GROMACS version"
antechamber -h 2>&1 | head -2
python -c "import parmed; print('parmed', parmed.__version__)"
```

---

## GPU-accelerated GROMACS (optional)

The conda-forge `gromacs` package runs on CPU only. For GPU production runs:

1. Load your cluster's CUDA module (e.g. `module load cuda/12.x`)
2. Build GROMACS from source following the [official instructions](https://manual.gromacs.org/documentation/current/install-guide/index.html)
3. Point `PATH` and `GMXLIB` to the custom build before running the prep script

The preparation steps (pdb2gmx, editconf, solvate, genion, grompp, energy minimization)
are CPU-bound and fast; the conda-forge build is sufficient for those.
