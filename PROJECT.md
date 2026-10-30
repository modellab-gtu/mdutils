# mdutils — Project Aims

Automated pipelines for classical and alchemical MD simulations of protein–ligand systems, developed at the Molecular Modeling Lab, Gebze Technical University.

---

## Aim 1 — Protein + Ligand + Solvent Classical MD (proprep)

**Directory:** `automd/prolig_normalMD/`

**Status:** Nearly complete.

Full pipeline for preparing a solvated protein–ligand complex and running classical MD in GROMACS or AMBER.

- Downloads protein from RCSB or accepts local PDB files
- Ligand parameterisation with GAFF2 (antechamber/parmchk2)
- Charge methods: AM1-BCC, RESP (Gaussian16), abcg2, mul, cm2
- Optional PDBFixer preprocessing (missing residues/atoms, protonation)
- Two topology pathways selectable via `--ff`:
  - **pdb2gmx + amb2gro** — GROMACS force fields (amber99sb-ildn, …)
  - **tleap + parmed** — AMBER force fields (ff14SB, ff19SB, …); auto-selected when the chosen FF has no GROMACS equivalent
- Solvation, ionisation, energy minimisation
- GROMACS output: topology, coordinates, index, MDP templates (NVT/NPT/MD)
- AMBER output: prmtop, inpcrd, and input file templates (min/nvt/npt/md)

---

## Aim 2 — Ligand + Solvent Classical MD (ligprep)

**Directory:** `automd/ligprep/`

**Script:** `automd/ligprep/prepare_lig_md.sh`

**Status:** Complete.

Pipeline for preparing a solvated ligand-only system and running classical MD — used for computing ligand properties in solution (e.g. hydration free energy reference leg, conformational sampling).

- Local ligand file or RCSB extraction input
- Ligand parameterisation with GAFF2 (antechamber/parmchk2)
- Charge methods: AM1-BCC, RESP (Gaussian16), abcg2, mul, cm2
- Two topology pathways selectable via `--ff`:
  - **tleap + amb2gro** — GROMACS force fields (for water/ion topology)
  - **tleap + parmed** — AMBER-native; auto-selected when FF has no GROMACS equivalent
- Solvation, ionisation, energy minimisation
- GROMACS output: topology, coordinates, index, MDP templates (NVT/NPT/MD)
- AMBER output (parmed route): prmtop, inpcrd, input file templates (min/nvt/npt/md)

---

## Aim 3 — Protein + Ligand + Solvent Alchemical Free Energy (complex leg)

**Directory:** `automd/fe/abfe/` (`bin/abfe_setup.sh`, `bin/abfe_analyze.sh`)

**Status:** Implemented for GROMACS (Boresch restraint, MBAR/TI analysis); Amber planned in Phase 4a.

Alchemical decoupling of the ligand from the solvated protein–ligand complex for absolute binding free energy (ABFE) calculations.

- Ligand decoupled via staged λ schedules: Coulombic then van der Waals
- GROMACS: softcore potentials, `free-energy` MDP templates per λ window
- AMBER: `icfe`/`clambda` input templates per λ window
- Supports thermodynamic integration (TI) and multistate Bennett acceptance ratio (MBAR) analysis

---

## Aim 4 — Ligand + Solvent Alchemical Free Energy (solvent leg)

**Directory:** `automd/fe/abfe/` (solvent leg; hydration free energy with `abfe_setup.sh --solvent`)

**Status:** Implemented for GROMACS; Amber planned in Phase 4a.

Alchemical decoupling of the ligand from bulk solvent — the reference leg of the ABFE thermodynamic cycle.

- Mirrors Aim 3 λ schedule and analysis workflow
- Combined with Aim 3 results: ΔG_bind = ΔG_solvent − ΔG_complex − ΔG_restr (decoupling legs; see README)
- Supports the GROMACS alchemical framework (Amber: Phase 4a, planned)
