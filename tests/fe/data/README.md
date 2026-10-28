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
