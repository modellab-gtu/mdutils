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
- Keeps only the first alternate location (blank or `A`) of protein atoms
- `--water` is honoured on both routes (tleap+parmed: `tip3p`, `spce`, `opc3`, `fb3`, `tip4pew`, `opc`)
- Writes a ready-to-run MD folder (see [Running the MD](#running-the-md))

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

Ready-to-use MDP files are in `automd/prolig_normalMD/MDP/`:

| File | Stage | Notes |
|------|-------|-------|
| `ions.mdp` | Ion placement | minimal grompp input |
| `em.mdp` | Energy minimisation | steep descent, POSRES + POSRES\_LIG + FLEXIBLE |
| `nvt.mdp` | NVT equilibration | sd integrator, 1 ns |
| `npt.mdp` | NPT equilibration | sd + C-rescale barostat, 1 ns |
| `md.mdp` | Production MD | sd + C-rescale, 100 ns, no restraints |

The prepare script copies these into the output folder's `MDP/` with `tc-grps` already set to
`Protein_<LIGID> Water_and_ions`. Put your own versions in `./MDP/` (where you run the script) to use them instead.

---

## automd/ligprep — Ligand + Solvent Classical MD Preparation

Pipeline for preparing a solvated ligand-only system for classical MD.

**Script:** `automd/ligprep/prepare_lig_md.sh`

### Features

- Local ligand file or RCSB PDB extraction
- Ligand parameterisation with GAFF2 (AmberTools antechamber/parmchk2)
- Same charge methods as prolig\_normalMD: `bcc`, `resp`, `abcg2`, …
- Two topology pathways via `--ff` (same logic as prolig\_normalMD):
  - **tleap + amb2gro** (default): GROMACS force fields for water/ion topology
  - **tleap + parmed** (`--use-parmed`): AMBER-native; auto-selected when `--ff` has no GROMACS equivalent

### Input modes

```bash
# 1. Local ligand file (SDF, MOL2, PDB, …)
prepare_lig_md.sh --local-ligand lig.sdf --ligand LIG --charge-method bcc --charge 0

# 2. Extract ligand from RCSB PDB
prepare_lig_md.sh --pdb 3HTB --ligand JZ4 --charge-method bcc --charge 0

# 3. RESP charges from local file
prepare_lig_md.sh --local-ligand lig.mol2 --ligand LIG \
    --charge-method resp --charge 0 --resp-opt yes

# 4. parmed pathway (AMBER ff14SB water model via tleap)
prepare_lig_md.sh --local-ligand lig.sdf --ligand LIG \
    --charge-method bcc --charge 0 --ff ff14SB
```

### MDP templates

Ready-to-use MDP files are in `automd/ligprep/MDP/`:

| File | Stage | Notes |
|------|-------|-------|
| `ions.mdp` | Ion placement | minimal grompp input |
| `em.mdp` | Energy minimisation | steep descent, POSRES\_LIG + FLEXIBLE |
| `nvt.mdp` | NVT equilibration | sd integrator, 1 ns |
| `npt.mdp` | NPT equilibration | sd + C-rescale barostat, 1 ns |
| `md.mdp` | Production MD | sd + C-rescale, 10 ns at 300 K, no restraints |

The prepare script copies these into the output folder's `MDP/` with `tc-grps` already set to
`<LIGID> Water_and_ions`. For small ligands the box padding is widened so periodic images stay at least
`--min-image` (3.2 nm) apart.

---

## Ligand input (both scripts)

The ligand's protonation state is your choice: give the ligand **as you want it simulated**.

| Source | How bond orders, charges and hydrogens are obtained |
|---|---|
| `--local-ligand lig.sdf` / `lig.mol2` *(preferred)* | used as is; must carry explicit hydrogens and formal charges |
| RCSB ligand (`--pdb ID --ligand CCD`) | bond orders and formal charges from the PDB Chemical Component Dictionary ideal SDF (RDKit), hydrogens added on the crystal coordinates; the CCD form is usually neutral |
| `--local-ligand lig.pdb` | Open Babel guesses bond orders and adds hydrogens (least reliable; a warning is printed) |

- `--charge` must equal the sum of the formal charges; the scripts stop otherwise. For an ionised ligand (carboxylate,
  protonated amine, …) give a protonated SDF with the matching `--charge`.
- `--lig-ph 7.4` re-protonates the ligand with Open Babel first (off by default; check the resulting charge).
- `parmchk2` runs with `-s 2` (GAFF2), matching the GAFF2 atom types.
- Implementation: `automd/common/ligand_input.py` (RDKit + Open Babel).

## Running the MD

Both scripts write, next to `topol.top`:

| File | Purpose |
|---|---|
| `MDP/` | em/ions/nvt/npt/md.mdp with `tc-grps` filled in |
| `run_em.sbatch` | energy minimisation (skipped if `em.gro` exists) |
| `run_md.sbatch` | NVT → NPT → production, one Slurm array task per replica |
| `submit_md.sh` | submits EM, then the replica array after it: `bash submit_md.sh` |

| Option | Default | Meaning |
|---|---|---|
| `--cluster gtu\|truba` | `gtu` | GTU ModelLab cluster or TRUBA/ARF |
| `--md-node cws01\|cws02` | `cws01` | GTU: cws01 = GPU, cws02 = CPU |
| `--replicas N` | 5 (prolig), 1 (lig) | independent replicas |
| `--em inline\|job` | `inline` (`job` with the GTU module, and always on TRUBA) | minimise now, or in `run_em.sbatch` |
| `--submit` | off | run `submit_md.sh` at the end (GTU) |
| `--truba-account/-partition/-ntasks/-omp/-gmx-module/-time` | `$USER`, `orfoz`, 110 ranks × 1 thread (prolig) or 1 × 16 (lig), `apps/gromacs/2024.1-oneapi2024`, `3-00:00:00` | TRUBA settings |

- **TRUBA:** prepare here, copy the folder to TRUBA, run `bash submit_md.sh` there. `em.tpr` is rebuilt on TRUBA,
  because its GROMACS 2024.1 cannot read tpr files from 2026.1.
- A ligand-in-water box (about 3 nm) cannot be split into many MPI domains, so `prepare_lig_md.sh` uses one MPI rank
  with OpenMP threads on TRUBA.
- `bin/mdprep-prolig` and `bin/mdprep-lig` are short wrappers for the two scripts.

## Free energy (ABFE / hydration)

Absolute binding free energies (double decoupling with a Boresch restraint) and hydration free energies, GROMACS
only (Amber is planned). Code: `automd/fe/abfe/`; launchers `bin/abfe_setup.sh` and `bin/abfe_analyze.sh`
(`-h` lists every option).

### Workflow

```bash
mdprep-prolig ... && mdprep-lig ...     # 1. prepare both MD folders (same ligand parameters)
bash submit_md.sh                       # 2. run MD replicas in each folder (the restraint comes from this MD)
abfe_setup.sh --complex PROLIG_DIR --solvent LIG_DIR --ligand JZ4 --submit   # 3. build + submit all windows
abfe_analyze.sh ABFE                    # 4. MBAR/TI, Boresch correction, warnings
# hydration free energy only (no protein): abfe_setup.sh --solvent LIG_DIR --ligand JZ4
```

- Both legs must share the same ligand parameters (atom names, types, masses, charges; checked) and temperature;
  if the two `md.mdp` files differ, give `--temperature T` (sets `ref-t`/`gen-temp` of every FE stage of both legs).
- `--replica R` picks the MD replica (default 0) whose `md_R.gro/tpr/xtc` are used.

| Option | Default | Meaning |
|---|---|---|
| `--repeats N` | 3 | independent repeats per leg |
| `--complex-ns` / `--solvent-ns` | 5 / 3 | production ns per window |
| `--schedule NAME\|FILE.yaml` | built-in | lambda schedule; `test_short` for smoke tests, or a custom YAML |
| `--cluster gtu\|truba` | `gtu` | `--complex-node` (cws01, GPU), `--solvent-node` (cws02, CPU, 8 threads); `--truba-*` as above |
| `--outdir`, `--overwrite` | `ABFE` | output folder (must be new unless `--overwrite`) |
| `--force-restraint` | off | accept a Boresch anchor angle within 10 deg of 0/180 |
| `--submit` | off | run `submit_abfe.sh` at the end (GTU) |

### Output layout

```
ABFE/
  setup.json  boresch.json  submit_abfe.sh
  complex|solvent/
    topol_fe.top  start.gro  posre_ref.gro  index.ndx
    rep<r>/run_abfe.sbatch  rep<r>/lam_NN/{em,nvt,npt,prod}.*
  results/                       (after abfe_analyze.sh)
```

Each window runs steep EM, NVT 100 ps, NPT 200 ps (C-rescale), then production; the tpr is rebuilt inside the job.
Restart a failed window with `sbatch --array=N rep<r>/run_abfe.sbatch`: finished steps are skipped and production
resumes from its checkpoint.

### Lambda schedules

- Complex, 36 states: restraint `0 0.05 0.1 0.2 0.5 1`, then Coulomb `0..1` step 0.1, then van der Waals `0..1`
  step 0.05. Solvent, 31 states: Coulomb, then van der Waals.
- Soft-core: `sc-alpha 0.5`, `sc-power 1`, `sc-r-power 6`; `couple-intramol = no`.

### Restraint

The Boresch restraint (1 distance, 2 angles, 3 dihedrals) is derived from the MD replica with
MDRestraintsGenerator, applied as a bonded-lambda (`bonded-lambdas` = `restraint-lambdas`) and switched on first in
the complex leg. Its analytic correction (release of the decoupled ligand into V0 = 1.66054 nm3) is stored in
`boresch.json`.

### Charged ligands

The net charge is compensated by co-alchemical monovalent counter-ion(s) that are decoupled with the ligand. They are
chosen at least min(1.5 nm, 0.45 x box) from the ligand (swapped with a far water molecule if needed) and held in
place by a 0.5 nm flat-bottom restraint.

### Analysis

`abfe_analyze.sh ABFE [--equil-ns 0.5] [--allow-missing] [--outdir DIR]` reads `setup.json` (temperature, number of
states), discards the equilibration, decorrelates each window and runs MBAR (main estimate) and TI per leg and
repeat, then combines the legs:

```
dG_bind = dG_solvent,decouple - dG_complex,decouple - dG_restr,off
dG_hyd  = -dG_solvent,decouple
```

(`dG_restr,off` is negative.) Output in `ABFE/results/`: `abfe_result.txt`, `abfe_result.csv`,
`overlap_<leg>_rep<r>.png`, `convergence_<leg>_rep<r>.png`. Values are given in kcal/mol and kJ/mol, as the mean and
standard error over the repeats. Warnings are printed for: minimum neighbour overlap < 0.03; |first-half - second-half|
> 0.5 kcal/mol; |MBAR - TI| > 0.5 kcal/mol; repeat spread (max - min) > 1.0 kcal/mol. Missing windows stop the
analysis unless `--allow-missing`.

### Resources

- Complex leg on cws01 (GPU): about 150 ns/day for T4L-size systems in FE windows. For example 36 windows x 5 ns x
  3 repeats = 540 ns, i.e. a few days of GPU time. Solvent leg on cws02 CPU, 8 threads per window.
- TRUBA: `--cluster truba` (same options as `mdprep-*`; one rank with OpenMP threads per window).
- Needs alchemlyb, pymbar and MDRestraintsGenerator (see [INSTALL.md](INSTALL.md)).

## legacy/

Old scripts kept for reference (not used): the former ligand/conformer utilities and the old TRUBA run script.

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
| AmberTools | 25 |
| parmed | 4.3.1 |
| numpy | 2.5.0 |

All installed from [conda-forge](https://conda-forge.org/) via `setup_env.sh`.
