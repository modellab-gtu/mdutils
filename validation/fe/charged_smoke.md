# Charged-ligand smoke test (co-alchemical counter-ion)

Date: 2026-10-28 to 2026-10-30. GROMACS 2026.1 (cws02 CPU, cws01 GPU), mdutils branch `abfe`.
Work folder: `/cta/users/akocak/maint/abfe_test/charged` (README there has the exact commands).

Spec criterion (design, "Testing and validation" item 4): net charge 0 at every state, and ΔG stable across 2 repeats.
Result: **pass, after a fix** (the ligand–ion wall, below).

## 1. Acetate hydration (ACE, q = −1; one Na⁺ co-decoupled)

- **Fixture:** `phase2_test/t4_ace_spce_run` (SPC/E). The box is small: shortest edge 2.72 nm, V ≈ 28 nm³.
- **Protocol:** 31 states (Coulomb ramp 0 → 1 in 0.1 steps, then vdW in 0.05 steps), 3 ns production per window, 2 repeats, `--equil-ns 0.5`.
- **Reported ΔG_hyd:** −ΔG_decouple of the ligand plus counter-ion pair. It is not comparable to an experimental single-ion value; the test checks consistency, not accuracy.

| run | rep | ΔG_hyd MBAR (kcal/mol) | TI | min neighbour overlap | MBAR − TI |
|---|---|---|---|---|---|
| no wall (`ace_hyd`) | 0 | −181.81 ± 0.67 | −184.40 | 0.0005 | 2.59 |
| no wall | 1 | −183.87 ± 0.86 | −184.64 | 0.0002 | 0.77 |
| **with wall** (`ace_hyd_wall`) | 0 | **−182.61 ± 0.09** | −182.79 | 0.018 | 0.19 |
| **with wall** | 1 | **−182.56 ± 0.09** | −182.64 | 0.016 | 0.08 |

Repeat spread: 2.06 kcal/mol without the wall, **0.05 kcal/mol** with it.

### What was wrong without the wall

The overlap was near zero only at the two ends of the Coulomb ramp (states 0→1 and 9→10).

The Hamiltonian was consistent. Re-running the same frames with the λ0 and λ1 tprs gave identical dH/dλ (folder `diag_rerun`).

The cause was configurational. In the solvent leg the ligand is unrestrained, so it diffuses to the restrained co-alchemical Na⁺ and forms a contact ion pair: the minimum distance is 0.23 nm. How often the pair is present depends on λ:

- 13 % of frames at λ_coul = 0;
- 100 % of frames for λ_coul = 0.1–0.9;
- 0 % at λ_coul = 1.

⟨dH/dλ_coul⟩ jumped from 606 to 409 kT and from 209 to 46 kT at the two ends of the ramp.

### The fix (spec, "Ligand–ion wall")

- Each co-alchemical ion gets a λ-independent flat-bottom bond (funct 10, in `[ intermolecular_interactions ]`) to the ligand heavy atom nearest the ligand COM.
  - It is zero above d_wall = min(1.0, 0.8 d_min) nm and harmonic (k = 1000 kJ mol⁻¹ nm⁻²) below.
  - No free-energy correction is applied.
- Verified by rerun on single-rank mdrun:
  - the energy uses the minimum-image distance (ion shifted by box vectors: 19.976 kJ/mol at 0.800 nm, as expected);
  - the energy is identical across λ states;
  - it adds nothing to dH/dλ.
- Setup warns when the wall excludes more than 10 % of the box. This acetate box gives 14.8 %.

The pull code was tried first and rejected: GROMACS aborts when a pull distance exceeds 0.49 × the box.

### Remaining minor issue

Even with the wall, the last Coulomb step (λ_coul 0.9 → 1.0, states 9→10) has overlap 0.016–0.018, below the 0.03 warning threshold. MBAR and TI still agree within 0.2 kcal/mol, and the repeats within 0.05 kcal/mol, so the estimate is sound. A denser Coulomb ramp near full decoupling (e.g. extra states at 0.85 and 0.95) for charged ligands would remove the warning. This is left as a follow-up.

## 2. JZ4 phenolate (q = −1) in T4L-L99A, short complex + solvent test

- **Inputs:** 3HTB, phenolate built with RDKit, bcc charges; 1 ns MD per leg.
- **Protocol:** `test_short` schedule (4 complex states, 3 solvent states), 0.2 ns per window, 1 repeat, `--temperature 300`.
- **Scope:** mechanics only. With 4 and 3 states over 0.2 ns, overlap is ~0 by design and the ΔG values are meaningless.

| Check | Result |
|---|---|
| Ion choice | Complex Na⁺ md#33321 → start#33382 (end of topology). Solvent Na⁺ md#2241 → start#2241. |
| Net charge 0 in state A and state B | Setup's charge check passed. |
| Ligand and ion at full charge at λ_coul = 0, both at zero charge at λ_coul = 1 | Setup's charge check passed. |
| grompp `-maxwarn 0`, all windows | Every window ran, so every grompp passed. |
| One `[ intermolecular_interactions ]` section in the complex leg, holding the Boresch terms and the wall bond | Yes. |
| Ion flat-bottom posres | Energy 0 in all complex windows, so the ion stayed within its 0.5 nm sphere. Solvent λ0: mean 0.5, max 6.3 kJ/mol, i.e. at most about 0.11 nm outside the sphere. |
| Wall | Anchor C7 in both legs, d_wall 1.0 nm. The solvent box gets the >10 % warning (18.4 %). |

### Bug found and fixed

The first solvent λ0 window aborted at NVT start: "no domain decomposition for 1 ranks … minimum cell size 2.29 nm".

- **Cause:** GROMACS builds a 1×1×1 domain decomposition even on one rank, and sizes the cell from the longest bonded interaction at start. That is the wall bond, about 2.1 nm.
- **Fix:** mdrun lines of legs with a wall get `-rdd 1.0` (commit 6cb6005). With one rank every atom is local; rerun energies of the PBC test frames are identical with and without `-rdd`.
- Setup refuses multi-rank TRUBA runs for charged ligands.
- The window was rerun (job 4478) and completed.

## Verdict

The co-alchemical machinery works:
- charge is neutral at every state;
- the ion stays restrained;
- the wall prevents ion pairing;
- the acetate repeats agree within 0.05 kcal/mol, and MBAR agrees with TI.

Follow-up: denser Coulomb spacing near λ_coul = 1 for charged ligands.
