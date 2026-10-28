# shellcheck shell=bash
# Shared by prepare_prolig_md.sh and prepare_lig_md.sh (sourced): writes the MD run folder content.
#   MDP/          em/ions/nvt/npt/md.mdp with tc-grps filled in for this system
#   run_em.sbatch energy minimisation (skipped if em.gro already exists)
#   run_md.sbatch NVT -> NPT -> production per replica (Slurm array, one task per replica)
#   submit_md.sh  submits EM, then the replica array after it
# Settings (set by the calling script): CLUSTER gtu|truba, MD_NODE cws01|cws02, REPLICAS, TRUBA_* and the
# source MDP files (MDP source dir in RUN_MDP_SRC).

write_run_files() {
    local tcgrps="$1"   # e.g. "Protein_JZ4 Water_and_ions"
    local head load mdrun mdrun_em ndx=""
    [[ -s index.ndx ]] && ndx="-n index.ndx"

    # every tc group must exist in index.ndx; a ligand whose residue name GROMACS files under another group
    # (e.g. ACE/NME are protein caps in residuetypes.dat) gets its own group from the residue name
    if [[ -s index.ndx ]]; then
        local g
        for g in $tcgrps; do
            grep -q "^\[ ${g} \]" index.ndx && continue
            if [[ "$g" == "$LIGID" ]]; then
                local n; n=$(grep -c '^\[' index.ndx)
                printf 'r %s\nname %s %s\nq\n' "$LIGID" "$n" "$LIGID" | gmx make_ndx -f "$EM_COORD" -n index.ndx -o index.ndx >/dev/null 2>&1
                grep -q "^\[ ${g} \]" index.ndx && log "Added index group [ ${g} ] (residue name ${LIGID})" && continue
            fi
            log "Warning: tc group [ ${g} ] is not in index.ndx; edit tc-grps in MDP/*.mdp"
        done
    fi

    mkdir -p MDP
    local f
    for f in em ions nvt npt md; do
        [[ -f "${RUN_MDP_SRC}/${f}.mdp" ]] || continue
        sed -E "s/^(tc[-_]grps[[:space:]]*=).*/\1 ${tcgrps}/" "${RUN_MDP_SRC}/${f}.mdp" > "MDP/${f}.mdp"
    done
    log "MDP files with tc-grps = ${tcgrps} written to MDP/"

    if [[ "$CLUSTER" == "truba" ]]; then
        head="#SBATCH -p ${TRUBA_PARTITION}
#SBATCH -A ${TRUBA_ACCOUNT}
#SBATCH -N 1
#SBATCH -n ${TRUBA_NTASKS}
#SBATCH -c ${TRUBA_OMP:-1}
#SBATCH -C weka
#SBATCH -t ${TRUBA_TIME}"
        load="module load ${TRUBA_GMX_MODULE}
export GMX_MAXBACKUP=-1 OMP_NUM_THREADS=${TRUBA_OMP:-1}   # set explicitly: the default 4 threads per rank made 110-rank runs ~3x slower
GMX=gmx_mpi"
        if (( TRUBA_NTASKS == 1 )); then
            # one rank: start gmx_mpi directly (mpirun would bind the rank, and all its threads, to one core)
            mdrun="gmx_mpi mdrun -ntomp ${TRUBA_OMP:-1}"
        else
            mdrun="mpirun -np \$SLURM_NTASKS gmx_mpi mdrun -ntomp ${TRUBA_OMP:-1}"
        fi
        mdrun_em="$mdrun"
    else
        local gpu=0; [[ "$MD_NODE" == "cws01" ]] && gpu=1
        head="#SBATCH -p main
#SBATCH -w ${MD_NODE}
#SBATCH --ntasks=1
$( (( gpu )) && printf '#SBATCH --gres=gpu:1\n#SBATCH -c 8' || printf '#SBATCH -c 32')
#SBATCH -t 7-00:00:00"
        load="module load $( (( gpu )) && echo gromacs/2026.1-tmpi-cuda || echo gromacs/2026.1-tmpi)
export GMX_MAXBACKUP=-1 CUDA_DEVICE_ORDER=PCI_BUS_ID
GMX=gmx"
        mdrun='gmx mdrun -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK'
        mdrun_em="$mdrun"
        (( gpu )) && mdrun="$mdrun -nb gpu -pme gpu -bonded gpu"   # update stays automatic: the sd integrator cannot update on the GPU
    fi

    cat > run_em.sbatch <<EOF
#!/bin/bash
#SBATCH -J em_${LIGID}
${head}
#SBATCH -o em_%j.out
set -e
${load}
cd "\$SLURM_SUBMIT_DIR"
[[ -s em.gro ]] && { echo "em.gro exists: nothing to do"; exit 0; }
# em.tpr is rebuilt here so the tpr matches this machine's GROMACS (TRUBA cannot read 2026.1 tpr files)
\$GMX grompp -f MDP/em.mdp -c ${EM_COORD} -r ${EM_COORD} -p ${FINAL_TOP} -o em.tpr -maxwarn 1
${mdrun_em} -deffnm em
EOF

    cat > run_md.sbatch <<EOF
#!/bin/bash
#SBATCH -J md_${LIGID}
${head}
#SBATCH -o md_%A_%a.out
#SBATCH --array=0-$((REPLICAS - 1))
# One array task per replica: NVT (new velocities) -> NPT -> production.
set -e
${load}
cd "\$SLURM_SUBMIT_DIR"
r=\${SLURM_ARRAY_TASK_ID:-0}
[[ -s em.gro ]] || { echo "em.gro missing: run run_em.sbatch first (or bash submit_md.sh)"; exit 1; }
\$GMX grompp -f MDP/nvt.mdp -c em.gro -r em.gro -p ${FINAL_TOP} ${ndx} -o nvt_\$r.tpr -maxwarn 1
${mdrun} -deffnm nvt_\$r
\$GMX grompp -f MDP/npt.mdp -c nvt_\$r.gro -r nvt_\$r.gro -t nvt_\$r.cpt -p ${FINAL_TOP} ${ndx} -o npt_\$r.tpr -maxwarn 1
${mdrun} -deffnm npt_\$r
\$GMX grompp -f MDP/md.mdp -c npt_\$r.gro -t npt_\$r.cpt -p ${FINAL_TOP} ${ndx} -o md_\$r.tpr -maxwarn 1
${mdrun} -deffnm md_\$r
EOF

    cat > submit_md.sh <<'EOF'
#!/bin/bash
# Energy minimisation (if needed), then the replica array after it.
cd "$(dirname "$0")"
if [[ -s em.gro ]]; then
    sbatch run_md.sbatch
else
    # if EM fails, Slurm cancels the MD array instead of leaving it pending forever
    jid=$(sbatch --parsable run_em.sbatch) && sbatch --dependency=afterok:${jid} --kill-on-invalid-dep=yes run_md.sbatch
fi
EOF
    chmod +x submit_md.sh
    log "Run files: run_em.sbatch, run_md.sbatch (${REPLICAS} replica(s), $( [[ "$CLUSTER" == truba ]] && echo "TRUBA ${TRUBA_PARTITION}" || echo "$MD_NODE")), submit_md.sh"
    if [[ "$CLUSTER" == "truba" ]]; then
        log "TRUBA: copy this folder to TRUBA and run 'bash submit_md.sh' there"
    elif (( SUBMIT )); then
        log "Submitting: $(bash submit_md.sh | tr '\n' ' ')"
    fi
}
