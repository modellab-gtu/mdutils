# mdutils ABFE (GROMACS absolute binding free energy) — design

Date: 2026-10-28. Status: approved in conversation, awaiting written-spec review.
Roadmap: advanced workflows Phase 4b (fulfils mdutils PROJECT.md Aims 3 and 4).

## Goal

Add absolute binding free energy (ABFE) and hydration free energy calculations to mdutils.
Users go prepare → MD → ABFE inside one tool. Each result is a ΔG with an error bar and quality warnings.

What the user decided:
- The code lives in mdutils, under `automd/fe/`.
- The Boresch restraint is chosen from the mdutils production MD.
- Independent λ windows run as Slurm arrays (approach A). λ replica exchange may come later as an option.
- Charged ligands are handled with a co-alchemical counter-ion.
- Validation is full.

## Non-goals (v1)

- λ replica exchange (`-multidir -replex`).
- Relative FE (pmx is Phase 4c; Amber is 4a).
- Analytic finite-size corrections.
- Multiple ligands or cofactors decoupled at once.
- Amber engine.

## Inputs

These are the run folders written by `prepare_prolig_md.sh` (complex) and `prepare_lig_md.sh` (solvent), after at least one MD replica has finished.

| Item | Source | Use |
|---|---|---|
| `topol.top` (+ includes), `index.ndx` | run folder (`FINAL_TOP` of `run_files.sh`) | topology and groups |
| `md_<r>.tpr`, `md_<r>.xtc`, `md_<r>.gro` | finished replica `r` (default 0) | restraint search; start coordinates |
| `MDP/md.mdp` | run folder | thermostat, cutoff and PME settings copied into the FE MDPs |
| ligand residue name | `--ligand` (must match the `[ moleculetype ]` name) | `couple-moltype` |

`--solvent` is optional. With only `--solvent` (no `--complex`), the tool runs a hydration free energy, which is the solvent leg alone.

## Commands (installed in `bin/`, sourced helpers in `automd/fe/`)

```
abfe_setup.sh  --complex DIR --solvent DIR --ligand LIG
               [--replica 0] [--repeats 3] [--complex-ns 5] [--solvent-ns 3]
               [--schedule FILE] [--cluster gtu|truba] [--complex-node cws01] [--solvent-node cws02]
               [--truba-* as in prepare scripts] [--outdir ABFE] [--submit]
abfe_analyze.sh ABFE_DIR [--equil-ns 0.5] [--temperature from mdp]
```

Hydration only: `abfe_setup.sh --solvent DIR --ligand LIG` (no complex leg and no restraint).

## Output layout

```
ABFE/
  setup.json                 # all options, versions, ligand charge, schedule, source dirs
  boresch.json               # anchors (atom indices + names), r0/θ0/φ0, force constants, ΔG_restr (kJ and kcal)
  complex/  topol_fe.top  start.gro  index.ndx  MDP/lam_NN.mdp  rep{0..R-1}/lam_NN/  run_abfe.sbatch
  solvent/  (same, without restraints)
  submit_abfe.sh             # submits one array per leg × repeat
  results/  abfe_result.txt  abfe_result.csv  overlap_<leg>_rep<r>.png  convergence_<leg>_rep<r>.png
```

## λ schedule (default; `--schedule` YAML overrides)

Each leg uses one MDP set with full lambda vectors (`restraint-lambdas` and `bonded-lambdas` (identical), `coul-lambdas`, `vdw-lambdas`). The states are ordered as restraint on → Coulomb off → vdW off.

- **Complex (36 states):**
  - restraint: 0, 0.05, 0.1, 0.2, 0.5, 1.0;
  - Coulomb: 0 → 1 in steps of 0.1;
  - vdW: 0 → 1 in steps of 0.05.
- **Solvent (31 states):** the same Coulomb and vdW ramps, no restraint.

## Per-window protocol (one Slurm array task = one window)

1. **EM:** steepest descent.
2. **NVT:** 100 ps.
3. **NPT:** 200 ps, C-rescale barostat.
4. **Production:**
   - complex 5 ns, solvent 3 ns;
   - `nstdhdl = 100` (0.2 ps);
   - `calc-lambda-neighbors = -1`;
   - `dhdl-print-energy = total`.

FE MDP settings:
- sd integrator, 2 fs, `constraints = h-bonds`;
- `couple-moltype = LIG`, `couple-lambda0 = vdw-q`, `couple-lambda1 = none`, `couple-intramol = no`;
- `sc-alpha = 0.5`, `sc-power = 1`, `sc-r-power = 6`;
- `separate-dhdl-file = yes`;
- `refcoord-scaling = com` in the barostat stages (needed with the co-alchemical ion's position restraint).

Cutoffs, PME and tc-grps come from the run folder's `md.mdp`.

**Starting coordinates:**
- the last frame of `md_<r>`;
- the repeat seed changes `ld-seed` and gen-vel.

**Resources:**
- **gtu complex:** cws01, `--gres=gpu:1 -c 8`, so 2 windows run at a time.
- **gtu solvent:** cws02, `-c 8`.
- **TRUBA:** same conventions as `run_files.sh` (1 rank × OMP for small systems, `gmx_mpi`).

A failed window is rerun with `sbatch --array=N run_abfe.sbatch`. Steps already finished (non-empty `.gro`) are skipped.

## Boresch restraint (complex leg)

- MDRestraintsGenerator 0.2.1 `FindBoreschRestraint` runs on `md_<r>.tpr/xtc`. Ligand atoms are the `--ligand` residue; host atoms are protein backbone.
- The script writes the restraint as `[ intermolecular_interactions ]` in `topol_fe.top`. It uses a bond (type 6) plus angles (type 1) and dihedrals (type 2) with A-state force constant 0 and B-state full force constant, controlled by `bonded-lambdas`. GROMACS switches `[ intermolecular_interactions ]` bond, angle and dihedral terms with `bonded-lambdas`, not `restraint-lambdas` (verified by rerun in Task 6), so both are set to the same schedule.
- The analytic release correction (Boresch et al. 2003, eq. 32) at the run temperature is

  `ΔG_restr,off = −kT ln[ 8π² V° √(K_r K_θA K_θB K_φA K_φB K_φC) / (r0² sinθA0 sinθB0 (2πkT)³) ]`

  with V° = 1.66054 nm³. It is stored in `boresch.json`.
- **Guard:** an anchor angle within 10° of 0° or 180° is an error, unless `--force-restraint` is given.

## Charged ligands: co-alchemical counter-ion

- Ligand net charge q ≠ 0 means |q| ions of the opposite sign are perturbed alongside it:
  - q < 0: Na⁺ (or the K⁺ already in the box);
  - q > 0: Cl⁻.
- In both legs the ions are picked from those already in the box. Each must be at least d_min = min(1.5 nm, 0.45 × shortest box edge) from the ligand in `start.gro`; among those, the ions farthest from the protein are taken.
  - If too few ions qualify (for example, a single neutralising ion sitting next to the ligand in a small ligand box), the ion is swapped with the water molecule farthest from the ligand. The whole water is moved to the ion's site and the ion to the water's site, and EM removes any clash.
  - The protein boxes have only about 1.2 nm of padding, so a protein-distance cutoff is not used.
  - Opposite-sign ions always exist: the neutralising ions in the ligand box, and the 0.15 M salt in the complex box. If there are too few anyway, setup stops with an error telling the user to re-prepare with salt.
- **Topology:** the chosen ions are moved to the end of the molecule list as a new moleculetype `<ION>_ALCH`. The B-state has charge 0 and unchanged LJ, so its charge follows `coul-lambdas`. Each ion has a flat-bottom position restraint (sphere of radius 0.5 nm around its start position, k = 1000 kJ mol⁻¹ nm⁻²), so it cannot drift into the site.
- **Ligand–ion wall (validation finding):** without it the free ligand forms a contact ion pair with its co-alchemical ion in the solvent leg (min distance 0.23 nm, Coulomb-ramp overlap ~0). Each ion gets a λ-independent one-sided distance wall to the ligand.
  - It is a `[ intermolecular_interactions ]` `[ bonds ]` entry of function type 10 (restraint potential) between the ligand heavy atom nearest the ligand centre of mass in `start.gro` (the anchor, recorded with its name in `setup.json` under `alch_wall`) and the ion: low = d_wall, up1 = 50 nm, up2 = 51 nm, kdr = 1000 kJ mol⁻¹ nm⁻². The B-state parameters equal the A state, so it is λ-independent (function type 10 follows `restraint-lambdas`; with B = A the schedule does not matter).
  - d_wall = round(min(1.0, 0.8 × d_min), 3) nm. It is zero above d_wall and harmonic below.
  - In the complex leg the wall bonds share the single `[ intermolecular_interactions ]` section with the Boresch terms.
  - The pull code is not used: it aborts when the group distance exceeds 0.49 × the box edge, which free ions do.
  - Intermolecular bonded interactions need single-rank mdrun (`-ntmpi 1`, no domain decomposition), as in the current protocol; setup refuses `--truba-ntasks > 1` for a charged ligand. mdrun still builds a 1×1×1 domain decomposition on one rank and sizes the cell from the longest bonded interaction at start (the wall bond, often > 2 nm), so the leg's mdrun lines get `-rdd 1.0`; with one rank all atoms are local and the wall energy is unchanged (verified by rerun).
  - If the wall excludes more than 10 % of the box, (4/3)π d_wall³ / V, setup warns of finite-size effects. The start structure must already satisfy the wall, otherwise setup stops.
  - No free-energy correction is applied: the wall enforces the non-interacting ligand/ion reference state that the co-alchemical method assumes. `abfe_analyze` prints a note line.
- **Checks:** setup asserts that the total charge is 0 in state A and in state B. It also asserts that the ions and the ligand are both at full charge when Coulomb λ = 0 and both at zero charge when Coulomb λ = 1.

## Analysis (`abfe_analyze.sh`)

- **Estimators:** alchemlyb 2.5 with pymbar 4.0.3.
  - `extract_u_nk` / `extract_dHdl` from GROMACS `dhdl.xvg`;
  - the first `--equil-ns` is discarded, then `statistical_inefficiency` subsampling;
  - MBAR is the main estimator, TI (trapezoid) the cross-check.
- **Per leg and repeat:** ΔG, overlap matrix, and forward/backward convergence.
- **Result:** ΔG_bind = ΔG_solvent,decouple − ΔG_complex,decouple − ΔG_restr,off (signs documented in the output), in kcal/mol and kJ/mol.
  - It is reported per repeat and as mean ± SEM over repeats.
  - Hydration ΔG_hyd = −ΔG_solvent,decouple.
- **Warnings** (printed in the result file; they do not fail the run):
  - any neighbour overlap below 0.03;
  - forward and backward final values differ by more than 0.5 kcal/mol;
  - MBAR and TI differ by more than 0.5 kcal/mol;
  - repeat spread above 1 kcal/mol;
  - missing or short windows. Analysis refuses a leg that has missing windows, unless `--allow-missing` is given.

## Environment and installation

- alchemlyb 2.5.0, pymbar 4.0.3 and MDRestraintsGenerator 0.2.1 are installed into `/cta/capps/envs/ligparam`.
  - A pip/conda dry-run comes first, so no existing package changes version (MDAnalysis 2.11, openmm 8.6.1, openff, ambertools-dac).
  - If the solve changes pinned packages, fall back to a separate `/cta/capps/envs/fe-analysis`, called by the scripts through `MDUTILS_FE_PYTHON`.
- **Code:** branch `abfe` → PR to modellab-gtu/mdutils.
- **Rollout:**
  - The installed clone is updated only after merge.
  - The module change (if any, e.g. `MDUTILS_FE_PYTHON`) is staged, then promoted by the user.
  - The `/cta/share/fe/gromacs` example (T4L-benzene + hydration) is staged in maint and promoted by the user.

## Testing and validation (FULL)

1. **Unit tests** (pytest, `tests/fe/`):
   - schedule and MDP generation (vectors equal length, ordering, settings copied);
   - topology edits (`intermolecular_interactions` block, `_ALCH` ion moleculetype, net charge 0 in A and B);
   - Boresch correction against a hand-computed reference value;
   - analysis on alchemtest's GROMACS benzene dataset, matching the alchemlyb documented value.
2. **Hydration** (solvent leg only):
   - 15–20 neutral FreeSolv molecules, both `--charge-method bcc` and `abcg2`.
   - Target: MUE vs experiment about 1.0 kcal/mol for bcc, and abcg2 not worse (literature trend, MUE about 0.4–0.7). Report MUE, RMSE and R per charge model.
3. **T4L-L99A/benzene** (PDB 181L) full ABFE, 3 repeats. Target: −5.2 kcal/mol (exp.) within ±1.5 kcal/mol, and no overlap warnings.
4. **Charged smoke:** one anionic ligand (e.g. a carboxylate from the FreeSolv set) solvent leg plus a short complex test. Check net charge 0 at every state, and that ΔG is stable across 2 repeats.

Validations run one after another on the 2 GPUs. **Estimated GPU time:** T4L about 36 windows × 5 ns × 3 repeats ≈ 540 ns. That is a few days on cws01; it is the only long job.

## Risks

- **MDRestraintsGenerator 0.2.1 vs MDAnalysis 2.11** may not be compatible. Mitigation: the dry-run, and a pinned separate env if needed.
- **GROMACS 2026 `dhdl.xvg` format vs the alchemlyb parser.** Mitigation: an early smoke run on a 2-window toy before the full schedule.
- **Wall time of the complex leg.** Mitigation: `--complex-ns` and the schedule are configurable; a short test schedule ships for smoke runs.
