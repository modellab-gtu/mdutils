#!/bin/bash
# Test: write_run_files creates the Water_and_ions tc group when the box has no ions.
# Usage: test_water_and_ions.sh <fixture dir: unpatched jz4_lig (index.ndx without Water_and_ions)>
# Needs gmx on PATH.
set -euo pipefail
FIX="${1:?fixture dir}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE=/cta/users/akocak/maint/abfe_test/pytest-tmp
mkdir -p "$BASE"
W=$(mktemp -d "$BASE/wai.XXXXXX")
trap 'rm -rf "$W"' EXIT
cp -r "$FIX"/. "$W"/
cd "$W"
rm -rf MDP run_em.sbatch run_md.sbatch submit_md.sh
! grep -q '^\[ Water_and_ions \]' index.ndx || { echo "FAIL: fixture index already has Water_and_ions"; exit 1; }

log() { echo "LOG: $*"; }
export GMX_MAXBACKUP=-1
CLUSTER=gtu MD_NODE=cws02 REPLICAS=1 SUBMIT=0
LIGID=JZ4 EM_COORD=solv_ions.gro FINAL_TOP=topol.top
RUN_MDP_SRC="$FIX/MDP"
# shellcheck source=/dev/null
source "$REPO/automd/common/run_files.sh"
write_run_files "JZ4 Water_and_ions" | tee out.log

grep -q '^\[ Water_and_ions \]' index.ndx || { echo "FAIL: index.ndx lacks [ Water_and_ions ]"; exit 1; }
grep -q 'Added index group \[ Water_and_ions \]' out.log || { echo "FAIL: no log line"; exit 1; }
! grep -q 'Warning: tc group' out.log || { echo "FAIL: warning emitted"; exit 1; }
gmx grompp -f MDP/nvt.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -n index.ndx -o t.tpr -maxwarn 5 >grompp.log 2>&1 \
    || { tail -20 grompp.log; echo "FAIL: grompp"; exit 1; }
echo "PASS"
