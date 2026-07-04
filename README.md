# mdutils

Utilities for setting up and analysing molecular dynamics simulations,
developed at the Molecular Modeling Lab, Gebze Technical University.

---

## automd/prolig_normalMD — Protein–Ligand Classical MD Preparation

Automated pipeline for preparing protein–ligand complexes for classical MD simulations in GROMACS or AMBER.

**Script:** `automd/prolig_normalMD/prepare_prolig_md.sh`

### Features

- Downloads protein from RCSB or accepts local PDB files
- Extracts and parameterises ligand with GAFF2 (AmberTools)
- Supports multiple charge methods: `bcc` (AM1-BCC), `resp` (Gaussian16 RESP), `abcg2`, `mul`, `cm2`, …
- RESP: adds `nosymm` to prevent G16 geometry reorientation and grafts original coordinates back
- Optional PDBFixer preprocessing for missing residues/atoms and protonation at chosen pH
- Two topology pathways selectable via `--ff`:
  - **pdb2gmx + amb2gro** (default): uses GROMACS force fields (`amber99sb-ildn`, `amber99sb`, …)
  - **tleap + parmed** (`--use-parmed`): builds the full system in tleap and converts to GROMACS; used automatically when the chosen FF has no GROMACS equivalent (e.g. `ff14SB`, `ff19SB`)
- Generates production-ready `index.ndx` with `Protein_<LIGID>` temperature-coupling group

### Force field selection

A single `--ff` flag controls both pathways:

| `--ff` value | GROMACS route | tleap+parmed route |
|---|---|---|
| `amber99sb-ildn` *(default)* | `amber99sb-ildn.ff` | `oldff/leaprc.ff99SBildn` |
| `amber99sb` | `amber99sb.ff` | `oldff/leaprc.ff99SB` |
| `amber03` | `amber03.ff` | `leaprc.protein.ff03.r1` |
| `ff14SB` | *(auto-switches to parmed)* | `leaprc.protein.ff14SB` |
| `ff19SB` | *(auto-switches to parmed)* | `leaprc.protein.ff19SB` |

If `--ff` names a force field that has no GROMACS `.ff` directory, the script automatically switches to the `--use-parmed` pathway and notifies you.

### Input modes

```bash
# 1. RCSB protein + RCSB ligand (default FF, pdb2gmx route)
prepare_prolig_md.sh --pdb 3HTB --ligand JZ4 --charge-method bcc --charge 0

# 2. RCSB protein + local ligand
prepare_prolig_md.sh --pdb 3HTB --local-ligand lig.sdf --ligand JZ4 \
    --charge-method bcc --charge 0

# 3. Local protein + local ligand (RESP with Gaussian16 optimisation)
prepare_prolig_md.sh --local-protein prot_noH.pdb --local-ligand lig.mol2 \
    --ligand LIG --charge-method resp --charge 0 --resp-opt yes

# 4. tleap+parmed pathway with ff14SB (auto-detected from --ff)
prepare_prolig_md.sh --pdb 3HTB --ligand JZ4 --charge-method bcc --charge 0 \
    --ff ff14SB

# 5. Explicitly request parmed pathway with default FF
prepare_prolig_md.sh --pdb 3HTB --ligand JZ4 --charge-method bcc --charge 0 \
    --use-parmed
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
