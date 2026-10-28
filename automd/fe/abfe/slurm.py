"""Slurm array files for the ABFE windows (one array per leg x repeat).

Conventions mirror ``automd/common/run_files.sh``:

* gtu: partition ``main``; cws01 = GPU node (``--gres=gpu:1``,
  ``gromacs/2026.1-tmpi-cuda``, ``CUDA_DEVICE_ORDER=PCI_BUS_ID``), any other
  node = CPU (``gromacs/2026.1-tmpi``); one thread-MPI rank x ``-c`` threads.
* truba: ``-C weka``, ``gmx_mpi``; with one rank ``gmx_mpi mdrun`` is started
  directly (mpirun would pin the rank and all its threads to one core).

The tpr files are always rebuilt inside the job (TRUBA's GROMACS cannot read
2026.1 tpr files). Every path in the job is relative to the repeat folder, so
an ABFE folder can be copied to another cluster as is.
"""
from __future__ import annotations

from pathlib import Path

GTU_DEFAULTS = {"cpus": 8, "time": "3-00:00:00"}
TRUBA_DEFAULTS = {"account": "", "partition": "orfoz", "ntasks": 1, "omp": 16,
                  "gmx_module": "apps/gromacs/2024.1-oneapi2024", "time": "3-00:00:00"}
GPU_NODES = ("cws01",)


def _header(leg: str, n_states: int, repeat: int, cluster: str, o: dict):
    """(#SBATCH lines, environment block, mdrun for EM, mdrun for MD)."""
    common = [f"#SBATCH -J abfe_{leg}_r{repeat}"]
    if cluster == "truba":
        o = {**TRUBA_DEFAULTS, **{k: v for k, v in o.items() if v is not None}}
        if not o["account"]:
            raise ValueError("TRUBA needs an account (--truba-account)")
        ntasks, omp = int(o["ntasks"]), int(o["omp"])
        head = common + [f"#SBATCH -p {o['partition']}", f"#SBATCH -A {o['account']}",
                         "#SBATCH -N 1", f"#SBATCH -n {ntasks}", f"#SBATCH -c {omp}",
                         "#SBATCH -C weka", f"#SBATCH -t {o['time']}"]
        load = (f"module load {o['gmx_module']}\n"
                f"export GMX_MAXBACKUP=-1 OMP_NUM_THREADS={omp}\n"
                "GMX=gmx_mpi")
        if ntasks == 1:
            # one rank: start gmx_mpi directly (mpirun would bind the rank, and all
            # its threads, to one core)
            mdrun = f"gmx_mpi mdrun -ntomp {omp}"
        else:
            mdrun = f"mpirun -np $SLURM_NTASKS gmx_mpi mdrun -ntomp {omp}"
        return head, load, mdrun, mdrun
    if cluster != "gtu":
        raise ValueError(f"cluster must be gtu or truba, got {cluster!r}")
    o = {**GTU_DEFAULTS, **{k: v for k, v in o.items() if v is not None}}
    node = o.get("node")
    if not node:
        raise ValueError("gtu needs a node (cws01/cws02)")
    gpu = node in GPU_NODES
    head = common + ["#SBATCH -p main", f"#SBATCH -w {node}", "#SBATCH --ntasks=1"]
    if gpu:
        head.append("#SBATCH --gres=gpu:1")
    head += [f"#SBATCH -c {int(o['cpus'])}", f"#SBATCH -t {o['time']}"]
    load = (f"module load {'gromacs/2026.1-tmpi-cuda' if gpu else 'gromacs/2026.1-tmpi'}\n"
            "export GMX_MAXBACKUP=-1 CUDA_DEVICE_ORDER=PCI_BUS_ID\n"
            "GMX=gmx")
    mdrun_em = "gmx mdrun -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK"
    # perturbed non-bonded pairs run on the CPU, the rest (and PME) on the GPU;
    # bonded/update stay automatic (sd cannot update on the GPU)
    mdrun = mdrun_em + (" -nb gpu -pme gpu" if gpu else "")
    return head, load, mdrun_em, mdrun


def write_leg_sbatch(leg_dir, leg: str, n_states: int, repeat: int, cluster: str,
                     node_opts: dict) -> Path:
    """Write ``<leg_dir>/rep<repeat>/run_abfe.sbatch``: one array task per window."""
    if n_states < 1:
        raise ValueError("n_states must be >= 1")
    rep_dir = Path(leg_dir) / f"rep{repeat}"
    rep_dir.mkdir(parents=True, exist_ok=True)
    head, load, mdrun_em, mdrun = _header(leg, n_states, repeat, cluster, dict(node_opts))
    head_txt = "\n".join(head)
    text = f"""#!/bin/bash
{head_txt}
#SBATCH -o abfe_%A_%a.out
#SBATCH --array=0-{n_states - 1}
# ABFE {leg} leg, repeat {repeat}: one array task per lambda window,
# EM -> NVT (new velocities) -> NPT -> production. A step whose .gro exists is
# skipped and production resumes from prod.cpt, so a failed window is rerun with
#   sbatch --array=N run_abfe.sbatch        (from this folder)
set -e
{load}
cd "$SLURM_SUBMIT_DIR"
[[ -f run_abfe.sbatch && -d lam_00 ]] || {{ echo "submit from the rep folder (bash submit_abfe.sh does)"; exit 1; }}
w=lam_$(printf %02d $SLURM_ARRAY_TASK_ID)
cd "$w"
echo "window $w on $(hostname), $(date)"

if [[ -s em.gro ]]; then echo "em.gro exists: skipping EM"; else
    $GMX grompp -f em.mdp -c ../../start.gro -r ../../posre_ref.gro -p ../../topol_fe.top -n ../../index.ndx -o em.tpr -maxwarn 0
    {mdrun_em} -deffnm em
fi
if [[ -s nvt.gro ]]; then echo "nvt.gro exists: skipping NVT"; else
    $GMX grompp -f nvt.mdp -c em.gro -r ../../posre_ref.gro -p ../../topol_fe.top -n ../../index.ndx -o nvt.tpr -maxwarn 0
    {mdrun} -deffnm nvt
fi
if [[ -s npt.gro ]]; then echo "npt.gro exists: skipping NPT"; else
    $GMX grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r ../../posre_ref.gro -p ../../topol_fe.top -n ../../index.ndx -o npt.tpr -maxwarn 0
    {mdrun} -deffnm npt
fi
if [[ -s prod.gro ]]; then echo "prod.gro exists: window $w is done"; else
    # keep the tpr of an interrupted run so -cpi continues it
    [[ -s prod.tpr && -s prod.cpt ]] || \\
        $GMX grompp -f prod.mdp -c npt.gro -t npt.cpt -r ../../posre_ref.gro -p ../../topol_fe.top -n ../../index.ndx -o prod.tpr -maxwarn 0
    {mdrun} -deffnm prod -cpi prod.cpt
fi
echo "window $w finished, $(date)"
"""
    path = rep_dir / "run_abfe.sbatch"
    path.write_text(text)
    return path


def write_submit(outdir, sbatches: list[str]) -> Path:
    """``submit_abfe.sh``: submit every repeat array (paths relative to outdir)."""
    lst = "\n".join(f"    {s}" for s in sbatches)
    text = f"""#!/bin/bash
# Submit one Slurm array (one task per lambda window) per leg x repeat.
set -e
cd "$(dirname "$(readlink -f "$0")")"
SBATCHES=(
{lst}
)
for s in "${{SBATCHES[@]}}"; do
    d=$(dirname "$s")
    jid=$(cd "$d" && sbatch --parsable "$(basename "$s")")
    echo "$d: job $jid"
done
"""
    p = Path(outdir) / "submit_abfe.sh"
    p.write_text(text)
    p.chmod(0o755)
    return p
