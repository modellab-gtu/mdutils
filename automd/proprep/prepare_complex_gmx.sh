#!/usr/bin/env bash
set -Eeuo pipefail

###############################################################################
# prepare_complex_gmx.sh
#
# Unified protein-ligand GROMACS preparation workflow.
# Supports single or dual ligand, RCSB download or local pre-prepared files.
#
# Modes (auto-detected from flags):
#   RCSB + 1 ligand  : -p PDB_ID -l RESNAME [--lig-chain C --lig-resi N]
#   RCSB + 2 ligands : -p PDB_ID --lig1 RES --lig2 RES [--lig1-chain/resi ...]
#   Local + 1 ligand : --local-protein FILE --local-ligand FILE --ligand ID
#   Local + 2 ligands: --local-protein FILE --local-lig1 FILE --lig1-id ID \
#                      --local-lig2 FILE --lig2-id ID
#
# Features:
#   - Auto-selects ligand copy by chain/residue number (RCSB mode)
#   - Optional PDBFixer preprocessing (RCSB mode)
#   - GAFF2/AmberTools parameterization per ligand
#   - GAFF2 atomtype deduplication for dual-ligand systems
#   - Builds topol.top from scratch (protein + ligand(s) + water + ions)
#   - Solvates, ionizes, minimizes
###############################################################################

# ── Defaults ──────────────────────────────────────────────────────────────────
PROID=""
PROID_GIVEN=0
LIG1ID=""
LIG2ID=""
LOCAL_PROTEIN_PDB=""
LOCAL_LIG1_PDB=""
LOCAL_LIG2_PDB=""
INPUT_MODE="rcsb"
TWO_LIG=0
LIG1CHAIN=""
LIG1RESI=""
LIG2CHAIN=""
LIG2RESI=""
FF="amber99sb-ildn"
WATER="tip3p"
BOXTYPE="dodecahedron"
BOXDIST="1.0"
IONS_MDP=""
EM_MDP=""
OUTDIR=""
KEEP=1
USE_PDBFIXER=0
PH=7.0
SCRIPT_START_DIR="$(pwd)"
MDPDIR="${SCRIPT_START_DIR}/MDP"

usage() {
cat <<EOF
Usage:
  $(basename "$0") [options]

Modes (auto-detected):
  RCSB + 1 ligand  : -p PDB_ID -l RESNAME
  RCSB + 2 ligands : -p PDB_ID --lig1 RES --lig2 RES
  Local + 1 ligand : --local-protein FILE --local-ligand FILE --ligand ID
  Local + 2 ligands: --local-protein FILE --local-lig1 FILE --lig1-id ID \\
                     --local-lig2 FILE --lig2-id ID

Common options:
  -p, --pdb ID              PDB ID from RCSB (required in RCSB mode)
  -l, --ligand, --lig1 RES  Ligand 1 residue name / molecule ID (required in RCSB mode)
      --lig-chain, --lig1-chain CHAIN
      --lig-resi,  --lig1-resi  RESI   Select specific ligand instance (RCSB)
      --lig2 RESNAME               Ligand 2 residue name — enables 2-ligand mode
      --lig2-chain CHAIN
      --lig2-resi  RESI
  -f, --ff NAME             Force field (default: ${FF})
  -w, --water NAME          Water model (default: ${WATER})
  -b, --boxtype TYPE        Box type (default: ${BOXTYPE})
  -d, --distance NM         Solvent padding nm (default: ${BOXDIST})
      --ions-mdp FILE       ions.mdp (default: MDP/ions.mdp)
      --em-mdp FILE         em.mdp   (default: MDP/em.mdp)
      --outdir DIR          Output directory
      --cleanup             Remove intermediates after run
      --pdbfixer            Run PDBFixer (RCSB mode only)
      --ph VALUE            pH for PDBFixer (default: ${PH})
  -h, --help                Show this help

Local-mode options:
      --local-protein FILE  Pre-cleaned protein PDB (ATOM records only)
      --local-ligand,
      --local-lig1 FILE     Ligand 1 PDB (HETATM records, one molecule)
      --lig1-id ID          Molecule ID for ligand 1
      --local-lig2 FILE     Ligand 2 PDB — enables 2-ligand mode
      --lig2-id ID          Molecule ID for ligand 2

Examples:
  # RCSB, 1 ligand
  $(basename "$0") -p 3HTB -l JZ4 --pdbfixer

  # RCSB, 2 ligands
  $(basename "$0") -p 5D41 --lig1 ANP --lig2 57N --lig1-chain A --lig1-resi 1102

  # Local, 1 ligand
  $(basename "$0") --local-protein egfr_wt.pdb --local-ligand atp.pdb --ligand ATP

  # Local, 2 ligands (project use)
  $(basename "$0") \\
    --local-protein structures/prepared/egfr_wt_eai001_anp.pdb \\
    --local-lig1 ligands/atp/atp_placed.pdb --lig1-id ATP \\
    --local-lig2 ligands/eai001/eai001_placed.pdb --lig2-id EAI001 \\
    --outdir md/S09_wt_eai001
EOF
}

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
die() { echo "Error: $*" >&2; exit 1; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"; }
optional_cmd() { command -v "$1" >/dev/null 2>&1; }
cleanup_on_error() { echo "Workflow failed." >&2; }
trap cleanup_on_error ERR

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--pdb)                    PROID="$2"; PROID_GIVEN=1; shift 2 ;;
        -l|--ligand|--lig1)          LIG1ID="$2";            shift 2 ;;
        --lig1-id)                   LIG1ID="$2";            shift 2 ;;
        --lig-chain|--lig1-chain)    LIG1CHAIN="$2";         shift 2 ;;
        --lig-resi|--lig1-resi)      LIG1RESI="$2";          shift 2 ;;
        --lig2)                      LIG2ID="$2"; TWO_LIG=1; shift 2 ;;
        --lig2-id)                   LIG2ID="$2";            shift 2 ;;
        --lig2-chain)                LIG2CHAIN="$2";         shift 2 ;;
        --lig2-resi)                 LIG2RESI="$2";          shift 2 ;;
        --local-protein)             LOCAL_PROTEIN_PDB="$2"; INPUT_MODE="local"; shift 2 ;;
        --local-ligand|--local-lig1) LOCAL_LIG1_PDB="$2";   INPUT_MODE="local"; shift 2 ;;
        --local-lig2)                LOCAL_LIG2_PDB="$2";   TWO_LIG=1; INPUT_MODE="local"; shift 2 ;;
        -f|--ff)                     FF="$2";                shift 2 ;;
        -w|--water)                  WATER="$2";             shift 2 ;;
        -b|--boxtype)                BOXTYPE="$2";           shift 2 ;;
        -d|--distance)               BOXDIST="$2";           shift 2 ;;
        --ions-mdp)                  IONS_MDP="$2";          shift 2 ;;
        --em-mdp)                    EM_MDP="$2";            shift 2 ;;
        --outdir)                    OUTDIR="$2";            shift 2 ;;
        --cleanup)                   KEEP=0;                 shift ;;
        --pdbfixer)                  USE_PDBFIXER=1;         shift ;;
        --ph)                        PH="$2";                shift 2 ;;
        -h|--help)                   usage; exit 0 ;;
        *)                           die "Unknown option: $1" ;;
    esac
done

# ── Validate arguments ────────────────────────────────────────────────────────
# In local mode, PROID is not meaningful unless user explicitly passed -p.
# Use "LOCAL" as placeholder so outdir/logfile names are not misleading.
[[ "$INPUT_MODE" == "local" && "$PROID_GIVEN" == "0" ]] && PROID="LOCAL"

if [[ "$INPUT_MODE" == "local" ]]; then
    [[ -n "$LOCAL_PROTEIN_PDB" ]] || die "Local mode requires --local-protein"
    [[ -n "$LOCAL_LIG1_PDB"    ]] || die "Local mode requires --local-ligand / --local-lig1"
    [[ -f "$LOCAL_PROTEIN_PDB" ]] || die "Protein PDB not found: $LOCAL_PROTEIN_PDB"
    [[ -f "$LOCAL_LIG1_PDB"    ]] || die "Ligand 1 PDB not found: $LOCAL_LIG1_PDB"
    if (( TWO_LIG )); then
        [[ -n "$LOCAL_LIG2_PDB" ]] || die "2-ligand local mode requires --local-lig2"
        [[ -n "$LIG2ID"          ]] || die "2-ligand local mode requires --lig2-id"
        [[ -f "$LOCAL_LIG2_PDB"  ]] || die "Ligand 2 PDB not found: $LOCAL_LIG2_PDB"
    fi
else
    [[ -n "$PROID"  ]] || die "RCSB mode requires -p / --pdb PDB_ID"
    [[ -n "$LIG1ID" ]] || die "RCSB mode requires -l / --lig1 RESNAME"
    if (( TWO_LIG )); then
        [[ -n "$LIG2ID" ]] || die "2-ligand RCSB mode requires --lig2 RESNAME"
    fi
fi

(( TWO_LIG == 0 )) || [[ "$LIG1ID" != "$LIG2ID" ]] || \
    die "Ligand IDs must be distinct (got '${LIG1ID}' for both)"

# ── Derived paths ─────────────────────────────────────────────────────────────
if [[ -z "$OUTDIR" ]]; then
    (( TWO_LIG )) && OUTDIR="prep_${PROID}_${LIG1ID}_${LIG2ID}" \
                  || OUTDIR="prep_${PROID}_${LIG1ID}"
fi
[[ -n "$IONS_MDP" ]] || IONS_MDP="${MDPDIR}/ions.mdp"
[[ -n "$EM_MDP"   ]] || EM_MDP="${MDPDIR}/em.mdp"

# ── Dependency check ──────────────────────────────────────────────────────────
for cmd in awk grep sed gmx obabel antechamber parmchk2 tleap tr sort uniq head tail wc; do
    require_cmd "$cmd"
done
[[ "$INPUT_MODE" == "rcsb" ]] && require_cmd wget
require_cmd amb2gro_top_gro.py

IONS_MDP="$(realpath "$IONS_MDP")"
EM_MDP="$(realpath "$EM_MDP")"
[[ -f "$IONS_MDP" ]] || die "ions.mdp not found: $IONS_MDP"
[[ -f "$EM_MDP"   ]] || die "em.mdp not found: $EM_MDP"

# Resolve local input paths to absolute before cd into output dir
if [[ "$INPUT_MODE" == "local" ]]; then
    LOCAL_PROTEIN_PDB="$(realpath "$LOCAL_PROTEIN_PDB")"
    LOCAL_LIG1_PDB="$(realpath "$LOCAL_LIG1_PDB")"
    (( TWO_LIG )) && [[ -n "$LOCAL_LIG2_PDB" ]] && \
        LOCAL_LIG2_PDB="$(realpath "$LOCAL_LIG2_PDB")"
fi

mkdir -p "$OUTDIR"
cd "$OUTDIR"

if (( TWO_LIG )); then
    LOGFILE="prepare_${PROID}_${LIG1ID}_${LIG2ID}.log"
else
    LOGFILE="prepare_${PROID}_${LIG1ID}.log"
fi
exec > >(tee -a "$LOGFILE") 2>&1

# ── File names ────────────────────────────────────────────────────────────────
RAW_PDB="${PROID}.pdb"
WORK_PDB="input_working.pdb"
PROTEIN_PDB="protein.pdb"
LIG1_PDB="${LIG1ID}.pdb"
LIG2_PDB="${LIG2ID}.pdb"
COMPLEX_PDB="complex.pdb"

PROTEIN_GRO="protein.gro"
PROTEIN_TOP="topol_protein.top"
PROTEIN_ITP="topol_protein.itp"

MERGED_ATOMTYPES_ITP="merged_GAFF2_atomtypes.itp"
ACTIVE_ATOMTYPES_ITP=""   # set after ligand prep (lig1 itp or merged)
FINAL_TOP="topol.top"
COMPLEX_GRO="complex.gro"
INDEX_NDX="index.ndx"

declare -a PROTEIN_MOLNAMES=()
LIG1_MOLNAME=""
LIG2_MOLNAME=""

# ── Functions ─────────────────────────────────────────────────────────────────

setup_inputs() {
    if [[ "$INPUT_MODE" == "local" ]]; then
        log "Using local PDB files"
        cp "$LOCAL_PROTEIN_PDB" "$PROTEIN_PDB"
        cp "$LOCAL_LIG1_PDB"    "$LIG1_PDB"
        if (( TWO_LIG )); then
            cp "$LOCAL_LIG2_PDB" "$LIG2_PDB"
            cat "$PROTEIN_PDB" "$LIG1_PDB" "$LIG2_PDB" > "$COMPLEX_PDB"
        else
            cat "$PROTEIN_PDB" "$LIG1_PDB" > "$COMPLEX_PDB"
        fi
    else
        log "Downloading ${PROID} from RCSB"
        wget -O "$RAW_PDB" "https://files.rcsb.org/download/${PROID}.pdb"
        [[ -s "$RAW_PDB" ]] || die "Downloaded PDB is empty"
        cp "$RAW_PDB" "$WORK_PDB"
    fi
}

maybe_run_pdbfixer() {
    (( USE_PDBFIXER == 1 )) || return 0
    if [[ "$INPUT_MODE" == "local" ]]; then
        log "PDBFixer skipped in local mode (files are pre-prepared)"
        return 0
    fi
    log "PDBFixer requested"

    if optional_cmd pdbfixer; then
        log "Running pdbfixer CLI"
        pdbfixer "$WORK_PDB" --output="$WORK_PDB.fixed" --add-atoms=heavy --keep-heterogens=all
        [[ -s "$WORK_PDB.fixed" ]] || die "pdbfixer CLI did not produce output"
        mv "$WORK_PDB.fixed" "$WORK_PDB"
        return 0
    fi

    if python3 - <<'PY' >/dev/null 2>&1
import importlib.util, sys
sys.exit(0 if importlib.util.find_spec("pdbfixer") else 1)
PY
    then
        log "Running pdbfixer via Python"
        python3 - "$WORK_PDB" "$PH" <<'PY'
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile
inp = sys.argv[1]; ph = float(sys.argv[2])
fixer = PDBFixer(filename=inp)
fixer.findMissingResidues(); fixer.findMissingAtoms()
fixer.addMissingAtoms(); fixer.addMissingHydrogens(ph)
with open(inp + ".fixed", "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
PY
        [[ -s "$WORK_PDB.fixed" ]] || die "Python pdbfixer did not produce output"
        mv "$WORK_PDB.fixed" "$WORK_PDB"
        return 0
    fi

    die "PDBFixer requested but neither CLI nor Python package is available"
}

# Extract one ligand from WORK_PDB by residue name, auto-detecting chain/resi.
# Usage: extract_one_ligand LIG_ID CHAIN_VAR RESI_VAR OUT_PDB
# CHAIN_VAR / RESI_VAR are names of global variables (nameref); may be empty
# on entry and will be set to the auto-detected values.
extract_one_ligand() {
    local lig_id="$1"
    local -n _chain="$2"
    local -n _resi="$3"
    local out_pdb="$4"

    awk -v lig="$lig_id" '
        /^HETATM/ {
            resn  = substr($0,18,3); alt = substr($0,17,1)
            chain = substr($0,22,1); resi = substr($0,23,4)
            gsub(/ /,"",resn); gsub(/ /,"",alt)
            gsub(/ /,"",chain); gsub(/ /,"",resi)
            if (resn == lig && (alt == "" || alt == "A")) print chain "\t" resi
        }' "$WORK_PDB" | sort -u > "${lig_id}_sites.tsv"

    [[ -s "${lig_id}_sites.tsv" ]] || die "No HETATM records for ${lig_id} in ${WORK_PDB}"
    log "Ligand ${lig_id} site(s):"; cat "${lig_id}_sites.tsv"

    if [[ -n "$_chain" || -n "$_resi" ]]; then
        [[ -n "$_chain" && -n "$_resi" ]] || die "--lig?-chain and --lig?-resi must be used together"
        grep -qP "^${_chain}\t${_resi}$" "${lig_id}_sites.tsv" \
            || die "Site chain=${_chain} resi=${_resi} not found for ${lig_id}"
    else
        _chain="$(awk 'NR==1{print $1}' "${lig_id}_sites.tsv")"
        _resi="$(awk  'NR==1{print $2}' "${lig_id}_sites.tsv")"
        [[ $(wc -l < "${lig_id}_sites.tsv") -gt 1 ]] && \
            log "Multiple ${lig_id} copies; using chain=${_chain} resi=${_resi}"
    fi
    log "Selected ${lig_id}: chain=${_chain} resi=${_resi}"

    awk -v lig="$lig_id" -v sc="$_chain" -v sr="$_resi" '
        /^HETATM/ {
            resn  = substr($0,18,3); alt = substr($0,17,1)
            chain = substr($0,22,1); resi = substr($0,23,4)
            gsub(/ /,"",resn); gsub(/ /,"",alt)
            gsub(/ /,"",chain); gsub(/ /,"",resi)
            if (resn == lig && chain == sc && resi == sr && (alt == "" || alt == "A")) print
        }' "$WORK_PDB" > "$out_pdb"
    [[ -s "$out_pdb" ]] || die "Extraction produced empty file for ${lig_id}"
}

detect_and_extract_rcsb() {
    [[ "$INPUT_MODE" == "rcsb" ]] || return 0
    log "Extracting protein ATOM records"
    grep '^ATOM  ' "$WORK_PDB" > "$PROTEIN_PDB"
    [[ -s "$PROTEIN_PDB" ]] || die "protein.pdb is empty"

    extract_one_ligand "$LIG1ID" LIG1CHAIN LIG1RESI "$LIG1_PDB"
    if (( TWO_LIG )); then
        extract_one_ligand "$LIG2ID" LIG2CHAIN LIG2RESI "$LIG2_PDB"
        cat "$PROTEIN_PDB" "$LIG1_PDB" "$LIG2_PDB" > "$COMPLEX_PDB"
    else
        cat "$PROTEIN_PDB" "$LIG1_PDB" > "$COMPLEX_PDB"
    fi
}

prepare_protein() {
    log "Running pdb2gmx"
    gmx pdb2gmx -f "$PROTEIN_PDB" -o "$PROTEIN_GRO" -p "$PROTEIN_TOP" \
        -i posre_protein.itp -ff "$FF" -water "$WATER" -ignh -merge all
    [[ -s "$PROTEIN_TOP" ]] || die "Protein topology not created"
    [[ -s "$PROTEIN_GRO" ]] || die "Protein gro not created"

    log "Extracting protein include block"
    sed -n '/moleculetype/,/endif/{/endif/!p;/endif/p;}' "$PROTEIN_TOP" > "$PROTEIN_ITP"

    log "Detecting protein molecule names"
    mapfile -t PROTEIN_MOLNAMES < <(
        awk '
            BEGIN { in_mol = 0 }
            /^\[ *molecules *\]/ { in_mol = 1; next }
            !in_mol { next }
            /^[[:space:]]*$/ || /^[[:space:]]*;/ || /^[[:space:]]*\[/ { next }
            { print $1 }
        ' "$PROTEIN_TOP"
    )
    (( ${#PROTEIN_MOLNAMES[@]} > 0 )) || die "Could not detect protein molecule name(s)"
    log "Protein molecule names: ${PROTEIN_MOLNAMES[*]}"
}

# Parameterize one ligand with GAFF2/AmberTools.
# Usage: prepare_ligand LIG_ID LIG_PDB MOLNAME_VAR
# Writes all output files named <LIG_ID>_* in the current directory.
# MOLNAME_VAR receives the moleculetype name detected from the resulting ITP.
prepare_ligand() {
    local lig_id="$1"
    local lig_pdb="$2"
    local -n _molname_out="$3"

    local h_mol2="${lig_id}_h.mol2"
    local ac_mol2="AC_${lig_id}.mol2"
    local ac_frcmod="AC_${lig_id}.frcmod"
    local ac_prmtop="AC_${lig_id}.prmtop"
    local ac_inpcrd="AC_${lig_id}.inpcrd"
    local gmx_top="${lig_id}_GMX.top"
    local gmx_gro="${lig_id}_GMX.gro"
    local gmx_pdb="${lig_id}_GMX.pdb"
    local atomtypes_itp="${lig_id}_GMX_GAFF.itp"
    local mol_itp="${lig_id}_GMX.itp"
    local posre_itp="posre_${lig_id}.itp"

    log "[${lig_id}] Generating mol2 with hydrogens"
    obabel "$lig_pdb" -O "$h_mol2" -h
    [[ -s "$h_mol2" ]] || die "[${lig_id}] Open Babel produced no output"

    log "[${lig_id}] Running antechamber (AM1-BCC charges)"
    antechamber -i "$h_mol2" -fi mol2 -o "$ac_mol2" -fo mol2 -c bcc -s 2

    log "[${lig_id}] Running parmchk2"
    parmchk2 -i "$ac_mol2" -f mol2 -o "$ac_frcmod"

    log "[${lig_id}] Running tleap"
    tleap -f - <<EOF
source leaprc.gaff2
lig = loadmol2 ${ac_mol2}
loadamberparams ${ac_frcmod}
saveamberparm lig ${ac_prmtop} ${ac_inpcrd}
quit
EOF
    [[ -s "$ac_prmtop" ]] || die "[${lig_id}] prmtop missing after tleap"
    [[ -s "$ac_inpcrd"  ]] || die "[${lig_id}] inpcrd missing after tleap"

    log "[${lig_id}] Converting Amber topology to GROMACS"
    amb2gro_top_gro.py -p "$ac_prmtop" -c "$ac_inpcrd" -t "$gmx_top" -g "$gmx_gro" -b "$gmx_pdb"
    [[ -s "$gmx_top" ]] || die "[${lig_id}] GROMACS top missing"
    [[ -s "$gmx_gro" ]] || die "[${lig_id}] GROMACS gro missing"

    log "[${lig_id}] Splitting topology"
    sed -n '/atomtypes/,/moleculetype/{/moleculetype/b;p}' "$gmx_top" > "$atomtypes_itp"
    sed -n '/moleculetype/,/system/{/system/b;p}'          "$gmx_top" > "$mol_itp"
    [[ -s "$atomtypes_itp" ]] || die "[${lig_id}] atomtypes ITP is empty"
    [[ -s "$mol_itp"        ]] || die "[${lig_id}] mol ITP is empty"

    cat >> "$mol_itp" <<EOF

; Include Position restraint file
#ifdef POSRES
#include "${posre_itp}"
#endif
EOF

    log "[${lig_id}] Generating position restraints"
    echo 0 | gmx genrestr -f "$gmx_gro" -o "$posre_itp" -fc 1000 1000 1000

    _molname_out="$(
        awk '/^\[ *moleculetype *\]/ {
            if (getline <= 0) exit
            if (getline <= 0) exit
            print $1; exit
        }' "$mol_itp"
    )"
    if [[ -z "$_molname_out" || "$_molname_out" == ";" ]]; then
        log "Warning: [${lig_id}] moleculetype name not detected, using ${lig_id}"
        _molname_out="$lig_id"
    fi
    log "[${lig_id}] Molecule name: ${_molname_out}"
}

# Merge GAFF2 atomtypes from both ligands, deduplicating by type name.
# Required because GAFF2 ligands share many type names (c3, ca, ha, …) and
# GROMACS will error on duplicate [ atomtypes ] entries.
merge_atomtypes() {
    log "Merging GAFF2 atomtypes (deduplicating)"
    local itp1="${LIG1ID}_GMX_GAFF.itp"
    local itp2="${LIG2ID}_GMX_GAFF.itp"
    {
        grep -E '^\[|^;' "$itp1" | head -2
        { grep -Ev '^\[|^;|^[[:space:]]*$' "$itp1"; \
          grep -Ev '^\[|^;|^[[:space:]]*$' "$itp2"; } | awk '!seen[$1]++'
    } > "$MERGED_ATOMTYPES_ITP"
    [[ -s "$MERGED_ATOMTYPES_ITP" ]] || die "Merged atomtypes file is empty"
    local n_types
    n_types=$(grep -Evc '^\[|^;|^[[:space:]]*$' "$MERGED_ATOMTYPES_ITP" || true)
    log "Merged atomtypes: ${n_types} unique GAFF2 types"
}

validate_topology_names() {
    log "Validating topology naming"
    grep -q '^\[ *moleculetype *\]' "$PROTEIN_ITP"      || die "Protein ITP lacks [ moleculetype ]"
    grep -q '^\[ *moleculetype *\]' "${LIG1ID}_GMX.itp" || die "Lig1 ITP lacks [ moleculetype ]"
    [[ -s "$ACTIVE_ATOMTYPES_ITP" ]]                     || die "Atomtypes ITP is empty"

    [[ "$LIG1_MOLNAME" == "$LIG1ID" ]] || \
        log "Warning: lig1 label (${LIG1ID}) != moleculetype (${LIG1_MOLNAME}); topology uses ${LIG1_MOLNAME}"

    if (( TWO_LIG )); then
        grep -q '^\[ *moleculetype *\]' "${LIG2ID}_GMX.itp" || die "Lig2 ITP lacks [ moleculetype ]"
        [[ "$LIG2_MOLNAME" == "$LIG2ID" ]] || \
            log "Warning: lig2 label (${LIG2ID}) != moleculetype (${LIG2_MOLNAME}); topology uses ${LIG2_MOLNAME}"
        [[ "$LIG1_MOLNAME" != "$LIG2_MOLNAME" ]] || \
            die "Both ligands resolved to moleculetype '${LIG1_MOLNAME}' — GROMACS requires unique names"
    fi
}

build_final_topology() {
    log "Building final topol.top"
    {
        echo "; Include forcefield parameters"
        echo "#include \"${FF}.ff/forcefield.itp\""
        echo "#include \"./${ACTIVE_ATOMTYPES_ITP}\""
        echo
        echo "#include \"./${PROTEIN_ITP}\""
        echo "#include \"${LIG1ID}_GMX.itp\""
        (( TWO_LIG )) && echo "#include \"${LIG2ID}_GMX.itp\""
        echo
        echo "; Include water topology"
        case "$WATER" in
            tip3p)   echo "#include \"${FF}.ff/tip3p.itp\""   ;;
            tip4p)   echo "#include \"${FF}.ff/tip4p.itp\""   ;;
            tip4pew) echo "#include \"${FF}.ff/tip4pew.itp\"" ;;
            spc)     echo "#include \"${FF}.ff/spc.itp\""     ;;
            spce)    echo "#include \"${FF}.ff/spce.itp\""    ;;
            *)       log "Warning: unknown water '${WATER}', defaulting to tip3p"
                     echo "#include \"${FF}.ff/tip3p.itp\""   ;;
        esac
        echo
        cat <<'EOF'
#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif
EOF
        echo
        echo "; Include topology for ions"
        echo "#include \"${FF}.ff/ions.itp\""
        echo
        echo "[ system ]"
        echo "; Name"
        (( TWO_LIG )) && echo "Protein-2lig complex in water" \
                      || echo "Protein-ligand in water"
        echo
        echo "[ molecules ]"
        echo "; Compound        #mols"
        for mol in "${PROTEIN_MOLNAMES[@]}"; do
            printf "%-18s %5d\n" "$mol" 1
        done
        printf "%-18s %5d\n" "$LIG1_MOLNAME" 1
        (( TWO_LIG )) && printf "%-18s %5d\n" "$LIG2_MOLNAME" 1
    } > "$FINAL_TOP"
    [[ -s "$FINAL_TOP" ]] || die "Failed to build ${FINAL_TOP}"
}

# Merge protein + ligand GRO files. Handles 2-way (1 lig) and 3-way (2 ligs).
# GRO format: line 1 = title, line 2 = atom count, lines 3..N+2 = atoms, last = box.
merge_gro_files() {
    local pro_atm lig1_atm total
    pro_atm=$(sed  -n '2p' "$PROTEIN_GRO"      | tr -d '[:space:]')
    lig1_atm=$(sed -n '2p' "${LIG1ID}_GMX.gro" | tr -d '[:space:]')
    [[ "$pro_atm"  =~ ^[0-9]+$ ]] || die "Invalid protein atom count in GRO"
    [[ "$lig1_atm" =~ ^[0-9]+$ ]] || die "Invalid ligand 1 atom count in GRO"

    if (( TWO_LIG )); then
        log "Merging protein + ${LIG1ID} + ${LIG2ID} GRO files"
        local lig2_atm
        lig2_atm=$(sed -n '2p' "${LIG2ID}_GMX.gro" | tr -d '[:space:]')
        [[ "$lig2_atm" =~ ^[0-9]+$ ]] || die "Invalid ligand 2 atom count in GRO"
        total=$((pro_atm + lig1_atm + lig2_atm))
        log "Atom counts: protein=${pro_atm}, ${LIG1ID}=${lig1_atm}, ${LIG2ID}=${lig2_atm}, total=${total}"
        {
            sed -n '1p' "$PROTEIN_GRO"
            echo "$total"
            sed -n '3,$p' "$PROTEIN_GRO"      | head -n "$pro_atm"
            sed -n '3,$p' "${LIG1ID}_GMX.gro" | head -n "$lig1_atm"
            sed -n '3,$p' "${LIG2ID}_GMX.gro" | head -n "$lig2_atm"
            tail -n 1 "$PROTEIN_GRO"
        } > "$COMPLEX_GRO"
    else
        log "Merging protein + ${LIG1ID} GRO files"
        total=$((pro_atm + lig1_atm))
        log "Atom counts: protein=${pro_atm}, ${LIG1ID}=${lig1_atm}, total=${total}"
        {
            sed -n '1p' "$PROTEIN_GRO"
            echo "$total"
            sed -n '3,$p' "$PROTEIN_GRO"      | head -n "$pro_atm"
            sed -n '3,$p' "${LIG1ID}_GMX.gro" | head -n "$lig1_atm"
            tail -n 1 "$PROTEIN_GRO"
        } > "$COMPLEX_GRO"
    fi
    [[ -s "$COMPLEX_GRO" ]] || die "complex.gro not created"
}

generate_index() {
    log "Generating index file"
    gmx make_ndx -f "$COMPLEX_GRO" -o "$INDEX_NDX" <<EOF
q
EOF
    [[ -s "$INDEX_NDX" ]] || die "index.ndx was not created"
}

solvate_and_minimize() {
    log "Creating simulation box"
    gmx editconf -f "$COMPLEX_GRO" -o newbox.gro -bt "$BOXTYPE" -d "$BOXDIST"

    log "Solvating"
    gmx solvate -cp newbox.gro -cs spc216.gro -p "$FINAL_TOP" -o solv.gro

    log "Preparing ions"
    gmx grompp -f "$IONS_MDP" -c solv.gro -p "$FINAL_TOP" -o ions.tpr

    log "Adding ions"
    echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p "$FINAL_TOP" \
        -pname NA -nname CL -neutral

    local maxwarn=1
    (( TWO_LIG )) && maxwarn=2

    log "Preparing energy minimization"
    gmx grompp -f "$EM_MDP" -c solv_ions.gro -p "$FINAL_TOP" -o em.tpr -maxwarn "$maxwarn"

    log "Running energy minimization"
    gmx mdrun -v -deffnm em
}

cleanup_files() {
    [[ "$KEEP" -eq 1 ]] && return 0
    log "Cleaning intermediate files"
    local files=(
        "${LIG1ID}_h.mol2"  "AC_${LIG1ID}.mol2"  "AC_${LIG1ID}.frcmod"
        "AC_${LIG1ID}.prmtop" "AC_${LIG1ID}.inpcrd"
        "${LIG1ID}_GMX.top" "${LIG1ID}_GMX.pdb"
        "${LIG1ID}_sites.tsv" "$WORK_PDB" ligand_sites.tsv
    )
    if (( TWO_LIG )); then
        files+=(
            "${LIG2ID}_h.mol2"  "AC_${LIG2ID}.mol2"  "AC_${LIG2ID}.frcmod"
            "AC_${LIG2ID}.prmtop" "AC_${LIG2ID}.inpcrd"
            "${LIG2ID}_GMX.top" "${LIG2ID}_GMX.pdb"
            "${LIG2ID}_sites.tsv"
        )
    fi
    rm -f "${files[@]}"
}

# ── Main ──────────────────────────────────────────────────────────────────────
log "Mode: INPUT=${INPUT_MODE}, LIGANDS=$(( TWO_LIG + 1 ))"

setup_inputs
maybe_run_pdbfixer
detect_and_extract_rcsb
prepare_protein
prepare_ligand "$LIG1ID" "$LIG1_PDB" LIG1_MOLNAME

if (( TWO_LIG )); then
    prepare_ligand "$LIG2ID" "$LIG2_PDB" LIG2_MOLNAME
    merge_atomtypes
    ACTIVE_ATOMTYPES_ITP="$MERGED_ATOMTYPES_ITP"
else
    ACTIVE_ATOMTYPES_ITP="${LIG1ID}_GMX_GAFF.itp"
fi

validate_topology_names
build_final_topology
merge_gro_files
generate_index
solvate_and_minimize
cleanup_files

log "All steps completed successfully"
log "  Protein molecule(s)      : ${PROTEIN_MOLNAMES[*]}"
log "  Ligand 1 (${LIG1ID}) moltype : ${LIG1_MOLNAME}"
(( TWO_LIG )) && log "  Ligand 2 (${LIG2ID}) moltype : ${LIG2_MOLNAME}"
log "  Output directory         : ${OUTDIR}"
log "  Minimized structure      : em.gro"
log "  Topology                 : topol.top"
log "  Index                    : index.ndx"
