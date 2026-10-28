# legacy/

Kept for reference; not used by the prepare scripts and not installed on the cluster.

- `ligprep_utils/`: the old ligand-preparation and conformer scripts (`ligPrep.py`, `runConfGen.py`, …). They hard-code
  `$HOME/Desktop` paths. For conformers use Auto3D, CREST or DeepConf instead (on the GTU cluster:
  `/cta/share/conformers`).
- `runMD_orfoz.slurm`: the old TRUBA run script. The prepare scripts now write `run_em.sbatch`, `run_md.sbatch` and
  `submit_md.sh` for GTU (`--md-node cws01|cws02`) or TRUBA (`--cluster truba`).
