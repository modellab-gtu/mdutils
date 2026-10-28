# ABFE test fixtures

Fixtures are mdutils run folders (each with md_0.tpr/xtc/gro, MDP/md.mdp, run_md.sbatch, index.ndx),
located by `tests/fe/conftest.py` and overridable via `MDUTILS_FE_TESTDATA` (a directory containing
the folder names). Tests skip when a folder is missing.

- t4_ace_spce_run: acetate ACE (-1), SPC/E, ParmEd route, lig_solv_gmx.top (made in phase2_test)
- t1_run: T4L + JZ4, pdb2gmx route, topol.top
- t2_run: T4L + JZ4, ParmEd route, complex_gmx.top
- t3_lid_truba_run: lidocaine +1, topol.top
- jz4_lig_run: JZ4 in water, made with
  `mdprep-lig --local-ligand jz4_h.sdf --ligand JZ4 --charge 0 --md-node cws02 --replicas 1 --em job --outdir jz4_lig`
  then copied to jz4_lig_run with nvt/npt nsteps=5000, md nsteps=10000, array 0-0, `bash submit_md.sh`.
  Note: this run's index.ndx was hand-patched with a `[ Water_and_ions ]` group before the run_files.sh fix
  (the ion-free box had none, so grompp failed). With the fix (commit "run_files: create Water_and_ions tc
  group when the box has no ions") a fresh `mdprep-lig` run needs no patch. The unpatched original is kept as
  fixtures/jz4_lig and is the input of tests/fe/test_run_files.py.

## dhdl fixtures for the analysis (tests/fe/test_analyze.py)

- toy_solvent/rep0/lam_00..02/prod.xvg: ACE (-1) solvent leg (t4_ace_spce_run), schedule
  `test_short` (3 states; (coul, vdw, bonded, restraint) = (0,0,0,0) (1,0,0,0) (1,1,0,0)),
  50 ps production, T = 300 K. Slurm job 3585 on cws02.
- toy_complex/rep0/lam_00..03/prod.xvg: T4L + JZ4 complex leg (t1_run), schedule `test_short`
  (4 states; (0,0,0,0) (0,0,1,1) (1,0,1,1) (1,1,1,1)), 50 ps production, T = 300 K set with
  `abfe_setup.sh --temperature 300` although t1_run's MDP/md.mdp is at 310 K. Slurm job 3588 on cws01 (GPU).
  The Boresch restraint dH/dl is in the bonded-lambda column; the restraint-lambda column is 0.
- Both were made by the Task 8 toy runs (abfe_setup.sh + the generated run_abfe.sbatch,
  GROMACS 2026.1; run dirs /cta/users/akocak/maint/abfe_test/{toy,toy_complex}); the prod.xvg files
  are copied whole (251 frames, nstdhdl = 100, no decimation).
- These are **parser/plumbing fixtures, not ΔG/convergence fixtures**: 3-4 states and 50 ps give
  essentially no neighbour overlap (complex lam_03 has ΔH to the other states of 1e7-1e11 kJ/mol),
  so MBAR/TI numbers on them are meaningless. Tests on them check only shape, finiteness and
  warnings, and must use `--equil-ns 0`. Estimator correctness is tested on alchemtest's GROMACS
  benzene and ABFE sets.
