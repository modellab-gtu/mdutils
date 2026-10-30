# FreeSolv hydration free energies: bcc vs abcg2

dG_hyd = -dG_solvent,decouple (MBAR, 1 repeat, error = statistical error of the single repeat); kcal/mol; experiment from FreeSolv database.json.

| res | molecule | FreeSolv id | expt | bcc calc | bcc err | abcg2 calc | abcg2 err |
|---|---|---|---|---|---|---|---|
| F01 | methanol | mobley_1636752 | -5.10 +- 0.60 | -2.88 +- 0.02 | +2.22 | -4.37 +- 0.02 | +0.73 |
| F02 | ethanol | mobley_2310185 | -5.00 +- 0.60 | -2.92 +- 0.02 | +2.08 | -4.35 +- 0.03 | +0.65 |
| F03 | 1-butanol | mobley_1019269 | -4.72 +- 0.60 | -2.73 +- 0.03 | +1.99 | -4.22 +- 0.03 | +0.50 |
| F04 | acetone | mobley_3867265 | -3.80 +- 0.60 | -3.45 +- 0.03 | +0.35 | -3.43 +- 0.03 | +0.37 |
| F05 | ethyl acetate | mobley_6973347 | -2.94 +- 0.60 | -3.38 +- 0.03 | -0.44 | -2.73 +- 0.03 | +0.21 |
| F06 | dimethyl ether | mobley_7015518 | -1.91 +- 0.60 | -0.62 +- 0.02 | +1.29 | -1.37 +- 0.02 | +0.54 |
| F07 | benzene | mobley_3053621 | -0.90 +- 0.20 | -1.25 +- 0.03 | -0.35 | n/a | n/a |
| F08 | toluene | mobley_1873346 | -0.90 +- 0.20 | -1.28 +- 0.03 | -0.38 | -0.61 +- 0.03 | +0.29 |
| F09 | phenol | mobley_20524 | -6.60 +- 0.20 | -5.34 +- 0.03 | +1.26 | n/a | n/a |
| F10 | anisole | mobley_7295828 | -2.45 +- 0.60 | -2.50 +- 0.04 | -0.05 | n/a | n/a |
| F11 | chlorobenzene | mobley_7608462 | -1.12 +- 0.60 | n/a | n/a | n/a | n/a |
| F12 | nitrobenzene | mobley_4193752 | -4.12 +- 0.60 | n/a | n/a | n/a | n/a |
| F13 | pyridine | mobley_296847 | -4.69 +- 0.60 | n/a | n/a | n/a | n/a |
| F14 | acetamide | mobley_8048190 | -9.71 +- 0.60 | -8.22 +- 0.03 | +1.49 | n/a | n/a |
| F15 | cyclohexane | mobley_2689721 | 1.23 +- 0.60 | 1.20 +- 0.03 | -0.03 | n/a | n/a |
| F16 | imidazole | mobley_7735340 | -9.63 +- 0.60 | -7.75 +- 0.03 | +1.88 | n/a | n/a |

## Statistics vs experiment

| charge model | n | MUE | RMSE | R | mean signed error |
|---|---|---|---|---|---|
| bcc | 13 | 1.06 | 1.32 | 0.959 | +0.87 |
| abcg2 | 7 | 0.47 | 0.50 | 0.995 | +0.47 |

## Pass criteria

- bcc MUE <= 1.5 kcal/mol: 1.06 -> PASS
- abcg2 MUE <= bcc MUE + 0.2: 0.47 vs 1.26 -> PASS
- (provisional: only 13 bcc and 7 abcg2 of 16 molecules done)

## Warnings per system

- F01_bcc (methanol): |calc - expt| = 2.22 > 2 kcal/mol
- F02_bcc (ethanol): |calc - expt| = 2.08 > 2 kcal/mol
- F15_bcc (cyclohexane): solvent rep0 forward/backward convergence failed (TimeoutError: forward/backward MBAR fits exceeded 300 s)
- F16_bcc (imidazole): solvent rep0 forward/backward convergence failed (TimeoutError: forward/backward MBAR fits exceeded 300 s)

## Missing results

F07_abcg2, F09_abcg2, F10_abcg2, F11_bcc, F11_abcg2, F12_bcc, F12_abcg2, F13_bcc, F13_abcg2, F14_abcg2, F15_abcg2, F16_abcg2
