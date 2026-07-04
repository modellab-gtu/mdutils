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

**Status:** In progress.

Pipeline for preparing a solvated ligand-only system and running classical MD — used for computing ligand properties in solution (e.g. hydration free energy reference leg, conformational sampling).

- Ligand parameterisation with GAFF2
- Solvation in explicit water box
- GROMACS and/or AMBER output

---

## Aim 3 — Protein + Ligand + Solvent Alchemical Free Energy (complex leg)

**Directory:** `automd/prolig_normalMD/` *(to be extended)*

**Status:** Planned.

Alchemical decoupling of the ligand from the solvated protein–ligand complex for absolute binding free energy (ABFE) calculations.

- Ligand decoupled via staged λ schedules: Coulombic then van der Waals
- GROMACS: softcore potentials, `free-energy` MDP templates per λ window
- AMBER: `icfe`/`clambda` input templates per λ window
- Supports thermodynamic integration (TI) and multistate Bennett acceptance ratio (MBAR) analysis

---

## Aim 4 — Ligand + Solvent Alchemical Free Energy (solvent leg)

**Directory:** `automd/ligprep/` *(to be extended)*

**Status:** Planned.

Alchemical decoupling of the ligand from bulk solvent — the reference leg of the ABFE thermodynamic cycle.

- Mirrors Aim 3 λ schedule and analysis workflow
- Combined with Aim 3 results: ΔG_bind = ΔG_complex − ΔG_solvent
- Supports both GROMACS and AMBER alchemical frameworks
