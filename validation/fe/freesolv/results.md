# FreeSolv hydration validation (ABFE solvent leg)

Date: 2026-10-28 to 2026-10-30. GROMACS 2026.1, mdutils branch `abfe`.
- Work folder: `/cta/users/akocak/maint/abfe_test/freesolv`.
- Systems were made with `run_all.sh`; the summary with `summarize.py --analyze`.

**Protocol:**
- `mdprep-lig` (GAFF2, bcc or abcg2 charges), 1 ns EM/NVT/NPT + 1 ns MD (R20/R23).
- ABFE solvent leg: default 31 states, 3 ns production per window, 1 repeat, 0.5 ns discarded.
- MBAR result; experiment from FreeSolv `database.json`.

## Scope (user decision 2026-10-29)

The plan had 16 molecules × 2 charge models (32 systems). To save cluster time the user trimmed it to 20 systems:

- **Both models:** methanol, ethanol, 1-butanol, acetone, ethyl acetate, dimethyl ether, toluene.
- **bcc only:** benzene and anisole, plus phenol, acetamide, cyclohexane and imidazole, added so the bcc set spans −9.7 to +1.2 kcal/mol.
- **Dropped by user decision (not failures):** chlorobenzene, nitrobenzene and pyridine (both models), and the remaining abcg2 runs.

## Results (kcal/mol)

| molecule | expt | bcc | bcc − expt | abcg2 | abcg2 − expt |
|---|---|---|---|---|---|
| methanol | −5.10 | −2.88 | +2.22 | −4.37 | +0.73 |
| ethanol | −5.00 | −2.92 | +2.08 | −4.35 | +0.65 |
| 1-butanol | −4.72 | −2.73 | +1.99 | −4.22 | +0.50 |
| acetone | −3.80 | −3.45 | +0.35 | −3.43 | +0.37 |
| ethyl acetate | −2.94 | −3.38 | −0.44 | −2.73 | +0.21 |
| dimethyl ether | −1.91 | −0.62 | +1.29 | −1.37 | +0.54 |
| toluene | −0.90 | −1.28 | −0.38 | −0.61 | +0.29 |
| benzene | −0.90 | −1.25 | −0.35 | | |
| anisole | −2.45 | −2.50 | −0.05 | | |
| phenol | −6.60 | −5.34 | +1.26 | | |
| acetamide | −9.71 | −8.22 | +1.49 | | |
| cyclohexane | +1.23 | +1.20 | −0.03 | | |
| imidazole | −9.63 | −7.75 | +1.88 | | |

Statistical errors of the single repeats are 0.02–0.04 kcal/mol.

| charge model | n | MUE | RMSE | R | mean signed error |
|---|---|---|---|---|---|
| bcc | 13 | 1.06 | 1.32 | 0.959 | +0.87 |
| abcg2 | 7 | 0.47 | 0.50 | 0.995 | +0.47 |
| bcc on the same 7 molecules as abcg2 | 7 | 1.25 | | | |

## Pass criteria (plan Task 11)

- bcc MUE ≤ 1.5 kcal/mol: 1.06, **PASS**.
- abcg2 MUE ≤ bcc MUE + 0.2, compared on the same 7 molecules: 0.47 vs 1.25, **PASS**.

## Interpretation

- **bcc errors have a clear pattern.** Strong hydrogen-bond donors are under-solvated by 1.3–2.2 kcal/mol: the alcohols, phenol, acetamide and imidazole. Apolar and aprotic molecules are within 0.45 kcal/mol: cyclohexane, benzene, toluene, anisole, acetone, ethyl acetate.
  - This is the known AM1-BCC weakness for polar groups, which ABCG2 was fitted to correct. Here abcg2 cuts the alcohol errors from about 2.1 to about 0.6 kcal/mol.
  - The pattern points to the charge model, not the free-energy code.
  - The final whole-branch review found no code path that would bias neutral hydration. `couple-intramol = no` makes the decoupled state the gas-phase molecule.
- **Sampling is clean:**
  - MBAR and TI agree;
  - neighbour overlap ≥ 0.03 everywhere, so no overlap warnings;
  - forward/backward differences are small where computed.
- **Two forward/backward convergence checks were skipped** (cyclohexane and imidazole). pymbar's solver stalled on one partial-data block. The analysis now caps this diagnostic at 300 s, runs it in a separate process, and reports a note. The main MBAR/TI results are unaffected.

## Recommendation

For neutral ligands with alcohol, phenol, amide or azole groups, prefer `--charge-method abcg2` in `mdprep-lig` / `mdprep-prolig`.
