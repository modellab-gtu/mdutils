#!/bin/bash
# FreeSolv hydration validation: 16 molecules x {bcc, abcg2}, prepare and submit everything.
# Idempotent: finished stages are skipped, jobs already in the queue are reused as dependencies.
#
#   bash run_all.sh                 prepare + submit all 32 systems (full protocol)
#   bash run_all.sh --smoke F01     one system, test_short schedule, 0.02 ns solvent windows, in smoke/
#   bash run_all.sh --only F05      only this residue (both methods), full protocol
#   options: --methods "bcc abcg2"  --work DIR
#
# Per system: mdprep-lig (--em job, 1 replica) -> sbatch run_em -> (afterok) run_md (MD shortened to 1 ns)
#   -> (afterok) setup job: abfe_setup.sh --solvent ... --repeats 1 --submit (runs on cws02).
# Node split: odd molecules solvent leg on cws02 (CPU), even on cws01 (GPU); MD runs on the same node.
# Job ids: $WORK/jobs.tsv (key, stage, jobid, time).
set -u
HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
WORK=/cta/users/akocak/maint/abfe_test/freesolv
PY=/cta/capps/envs/ligparam/bin/python
METHODS="bcc abcg2"; ONLY=""; SMOKE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --smoke) SMOKE="$2"; shift 2 ;;
        --only) ONLY="$2"; shift 2 ;;
        --methods) METHODS="$2"; shift 2 ;;
        --work) WORK="$2"; shift 2 ;;
        *) echo "unknown option $1"; exit 2 ;;
    esac
done
JOBS="$WORK/jobs.tsv"
MD_NS_STEPS=500000          # 1 ns at dt = 0.002
[[ -n "$SMOKE" ]] && { ONLY="$SMOKE"; SCHED=(--schedule test_short --solvent-ns 0.02); BASE="$WORK/smoke"; } \
    || { SCHED=(--solvent-ns 3); BASE="$WORK"; }
mkdir -p "$BASE"; touch "$JOBS"

# ligparam python + gromacs (as the other mdutils tools on chn)
module load gromacs/2026.1-tmpi >/dev/null 2>&1
export PATH=/cta/capps/envs/ligparam/bin:$PATH PYTHONNOUSERSITE=1 AMBERHOME=/cta/capps/envs/ligparam
unset PYTHONPATH; export LD_LIBRARY_PATH=/cta/capps/envs/ligparam/lib:${LD_LIBRARY_PATH:-}
export LC_ALL=C

alive() { [[ -n "${1:-}" ]] && [[ -n "$(squeue -h -j "$1" 2>/dev/null)" ]]; }
last_job() { awk -F'\t' -v k="$1" -v s="$2" '$1==k && $2==s {j=$3} END {print j}' "$JOBS"; }
record() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$(date +%FT%T)" >> "$JOBS"; }

process() {   # res name smiles idx method
    local res="$1" name="$2" smi="$3" idx="$4" m="$5"
    local node=cws02; (( idx % 2 == 0 )) && node=cws01
    local d="$BASE/${res}_${m}" key="${SMOKE:+smoke/}${res}_${m}"
    mkdir -p "$d"
    # 1. prepare (runs on chn, minutes)
    if [[ ! -f "$d/.prep_done" ]]; then
        rm -rf "$d/run" "$d/lig.sdf"
        "$PY" -I "$HERE/build_sdf.py" "$smi" "$d/lig.sdf" || { echo "$key: SDF failed"; return 1; }
        ( cd "$d" && bash "$REPO/bin/mdprep-lig" --local-ligand lig.sdf --ligand "$res" --charge 0 \
            --charge-method "$m" --replicas 1 --em job --md-node "$node" --outdir run > prep.log 2>&1 ) \
            || { echo "$key: mdprep-lig failed (see $d/prep.log)"; return 1; }
        [[ -s "$d/run/run_md.sbatch" && -s "$d/run/solv_ions.gro" ]] || { echo "$key: incomplete run folder"; return 1; }
        grep -Eq '^dt[[:space:]]*=[[:space:]]*0\.002' "$d/run/MDP/md.mdp" || { echo "$key: unexpected dt in md.mdp"; return 1; }
        sed -i -E "s/^nsteps[[:space:]]*=.*/nsteps                   = ${MD_NS_STEPS}   ; 1 ns (FreeSolv validation)/" "$d/run/MDP/md.mdp"
        touch "$d/.prep_done"
    fi
    # setup job script (abfe_setup on a compute node: ligparam python + gromacs only)
    cat > "$d/abfe_setup_job.sbatch" <<EOS
#!/bin/bash
#SBATCH -J fsset_${res}${m:0:1}
#SBATCH -p main
#SBATCH -w cws02
#SBATCH --ntasks=1
#SBATCH -c 2
#SBATCH -t 04:00:00
#SBATCH -o setup_%j.out
set -e
module load gromacs/2026.1-tmpi
export PATH=/cta/capps/envs/ligparam/bin:\$PATH PYTHONNOUSERSITE=1 AMBERHOME=/cta/capps/envs/ligparam LC_ALL=C
unset PYTHONPATH; export LD_LIBRARY_PATH=/cta/capps/envs/ligparam/lib:\${LD_LIBRARY_PATH:-}
cd "$d"
[[ -s $d/run/md_0.gro ]] || { echo "MD did not finish"; exit 1; }
bash $REPO/bin/abfe_setup.sh --solvent $d/run --ligand $res --replica 0 --repeats 1 ${SCHED[*]} \\
    --solvent-node $node --outdir $d/ABFE --submit
touch setup_done
EOS
    local jem jmd jset dep
    # 2. EM
    if [[ -s "$d/run/em.gro" ]]; then jem=""
    else jem="$(last_job "$key" em)"
        alive "$jem" || { jem=$(cd "$d/run" && sbatch --parsable run_em.sbatch) || return 1; record "$key" em "$jem"; echo "$key: EM job $jem"; }
    fi
    # 3. MD
    if [[ -s "$d/run/md_0.gro" ]]; then jmd=""
    else jmd="$(last_job "$key" md)"
        if ! alive "$jmd"; then
            dep=(); [[ -n "$jem" ]] && dep=(--dependency="afterok:$jem")
            jmd=$(cd "$d/run" && sbatch --parsable "${dep[@]}" run_md.sbatch) || return 1
            record "$key" md "$jmd"; echo "$key: MD job $jmd (after ${jem:-nothing})"
        fi
    fi
    # 4. ABFE setup + submit
    if [[ -f "$d/setup_done" ]]; then echo "$key: setup done"; return 0; fi
    jset="$(last_job "$key" setup)"
    if ! alive "$jset"; then
        dep=(); [[ -n "$jmd" ]] && dep=(--dependency="afterok:$jmd")
        # a half-written ABFE folder from a failed setup would block abfe_setup (needs a new folder)
        [[ -f "$d/ABFE/submit_abfe.sh" ]] || rm -rf "$d/ABFE"
        jset=$(cd "$d" && sbatch --parsable "${dep[@]}" abfe_setup_job.sbatch) || return 1
        record "$key" setup "$jset"; echo "$key: setup job $jset (after ${jmd:-nothing})"
    fi
}

idx=0
while IFS=, read -r res name fid iupac smiles expt dexpt; do
    [[ "$res" == res ]] && continue
    idx=$((idx + 1))
    [[ -n "$ONLY" && "$res" != "$ONLY" ]] && continue
    for m in $METHODS; do
        process "$res" "$name" "$smiles" "$idx" "$m" </dev/null || echo "FAILED: ${res}_${m}"
    done
done < "$HERE/molecules.csv"
