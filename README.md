# mdutils

Utilities for setting up and analysing molecular dynamics simulations,
developed at the Molecular Modeling Lab, Gebze Technical University.

---

## automd/proprep — Protein–Ligand GROMACS Preparation

Automated pipeline for preparing protein–ligand complexes for GROMACS MD simulations.

**Script:** `automd/proprep/prepare_complex_gmx_labready_resp.sh`

### Features

- Downloads protein from RCSB or accepts local PDB files
- Extracts and parameterises ligand with GAFF2 (AmberTools)
- Supports multiple charge methods: `bcc` (AM1-BCC), `resp` (Gaussian16 RESP), `abcg2`, `mul`, `cm2`, …
- RESP: adds `nosymm` to prevent G16 geometry reorientation and grafts original coordinates back
- Optional PDBFixer preprocessing for missing residues/atoms and protonation at chosen pH
- Runs pdb2gmx, editconf, solvate, genion, energy minimisation
- Generates production-ready `index.ndx` with `Protein_<LIGID>` temperature-coupling group

### Input modes

```bash
# 1. RCSB protein + RCSB ligand
prepare_complex_gmx_labready_resp.sh --pdb 3HTB --ligand JZ4 --charge-method bcc --charge 0

# 2. RCSB protein + local ligand
prepare_complex_gmx_labready_resp.sh --pdb 3HTB --local-ligand lig.sdf --ligand JZ4 \
    --charge-method bcc --charge 0

# 3. Local protein + local ligand (RESP with Gaussian16 optimisation)
prepare_complex_gmx_labready_resp.sh --local-protein prot_noH.pdb --local-ligand lig.mol2 \
    --ligand LIG --charge-method resp --charge 0 --resp-opt yes
```

### MDP templates

Ready-to-use MDP files are in `automd/proprep/MDP/`:

| File | Stage | Notes |
|------|-------|-------|
| `ions.mdp` | Ion placement | minimal grompp input |
| `em.mdp` | Energy minimisation | steep descent, POSRES + POSRES\_LIG + FLEXIBLE |
| `nvt.mdp` | NVT equilibration | sd integrator, 1 ns |
| `npt.mdp` | NPT equilibration | sd + C-rescale barostat, 1 ns |
| `md.mdp` | Production MD | sd + C-rescale, 100 ns, no restraints |

> **Note:** Before running NVT/NPT/MD, replace `Protein_` in `tc-grps` with your
> actual ligand name (e.g. `Protein_JZ4`). The prep script prints the exact line
> at the end of its run.

---

## Installation

### Requirements

- Linux x86_64
- Miniconda or Anaconda
- Gaussian 16 *(optional — only needed for `--charge-method resp`)*

### Quick setup

```bash
bash setup_env.sh
conda activate mdutils
```

This creates a `mdutils` conda environment with all required tools.
See [INSTALL.md](INSTALL.md) for step-by-step instructions and GPU GROMACS notes.

### Tested versions (2026-07-04)

| Package | Version |
|---------|---------|
| Python | 3.12 |
| OpenBabel | 3.1.0 |
| OpenMM | 8.5.2 |
| PDBFixer | 1.12 |
| GROMACS | 2025.4 |
| AmberTools | 24.8 |
| parmed | 4.3.1 |
| numpy | 2.5.0 |

All installed from [conda-forge](https://conda-forge.org/) via `setup_env.sh`.
