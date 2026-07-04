#!/usr/bin/env bash
set -Eeuo pipefail

###############################################################################
# prepare_lig_md.sh
#
# Ligand + solvent classical MD preparation pipeline.
# Parameterises a ligand with GAFF2 and prepares a solvated system for MD.
#
# Two topology pathways:
#   - amb2gro  (default): tleap ligand-only → amb2gro_top_gro.py → GROMACS
#   - parmed   (--use-parmed): tleap full solvated system → parmed → GROMACS
#     (also writes AMBER input files for sander/pmemd)
#
# Input modes:
#   1) Local ligand file:   --local-ligand lig.sdf --ligand LIG
#   2) Extract from RCSB:  --pdb 3HTB --ligand JZ4
###############################################################################

LIGID="LIG"
LIGID_USER_SET=0
LOCAL_LIGAND_PDB=""
LOCAL_LIGAND_GIVEN=0
INPUT_MODE="local"
PROID=""
LIGCHAIN=""
LIGRESI=""
FF="amber99sb-ildn"
WATER="tip3p"
BOXTYPE="dodecahedron"
BOXDIST="1.0"
IONS_MDP=""
EM_MDP=""
OUTDIR=""
VERBOSE=0
RUN_SOLVATE_MINIMIZE=1

USE_PARMED=0

CHARGE_METHOD="bcc"
LIG_CHARGE=0
LIG_MULT=1
RESP_OPT="yes"
NPROC=32
MEM="20GB"
BASIS_OPT="B3LYP/6-31G*"
BASIS_ESP="HF/6-31G*"

SCRIPT_START_DIR="$(pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MDPDIR="${SCRIPT_START_DIR}/MDP"
MDPDIR_BUNDLED="${SCRIPT_DIR}/MDP"

usage() {
cat <<EOF_USAGE
Usage:
  $(basename "$0") [options]

Input modes:
  1) Local ligand file (default):
    $(basename "$0") --local-ligand lig.sdf --ligand LIG --charge-method bcc --charge 0

  2) Extract ligand from RCSB PDB:
    $(basename "$0") --pdb 3HTB --ligand JZ4 --charge-method bcc --charge 0

General options:
  -l, --ligand RES             Ligand residue/molecule name (required), default: ${LIGID}
      --local-ligand FILE      Local ligand file (PDB, SDF, MOL2, …)
  -p, --pdb ID                 RCSB PDB ID for ligand extraction
      --lig-chain CHAIN        Ligand chain ID (RCSB mode, selects specific copy)
      --lig-resi RESI          Ligand residue number (RCSB mode)
  -f, --ff NAME                Force field for GROMACS water/ion topology, default: ${FF}
                               (also controls auto-routing: AMBER-only FFs → parmed pathway)
  -w, --water NAME             Water model, default: ${WATER}
  -b, --boxtype TYPE           Box type, default: ${BOXTYPE}
  -d, --distance NM            Solvent box padding in nm, default: ${BOXDIST}
      --ions-mdp FILE          ions.mdp path (default: MDP/ions.mdp if present, else bundled)
      --em-mdp FILE            em.mdp path  (default: MDP/em.mdp  if present, else bundled)
      --outdir DIR             Output directory, default: prep_lig_<LIG>_<method>
      --verbose                Keep all intermediate files (default: remove intermediates)
      --skip-solvate-minimize  Stop after topology, coordinates, and index
      --use-parmed             Use tleap+parmed pathway (auto-selected for AMBER-only FFs)

Ligand charge options:
      --charge-method METHOD   resp (Gaussian16 RESP) or any antechamber -c value
                               (bcc, abcg2, mul, cm2, …), default: ${CHARGE_METHOD}
  -q, --charge INT             Net ligand charge, default: ${LIG_CHARGE}
  -m, --mult INT               Spin multiplicity for RESP/Gaussian, default: ${LIG_MULT}
      --resp-opt yes|no        For RESP: optimize with G16 before ESP, default: ${RESP_OPT}
      --nproc INT              Gaussian nprocshared, default: ${NPROC}
      --mem MEM                Gaussian memory, default: ${MEM}
      --basis-opt LEVEL        Gaussian optimization level, default: ${BASIS_OPT}
      --basis-esp LEVEL        Gaussian ESP level, default: ${BASIS_ESP}

Examples:
  # AM1-BCC charges, local SDF file
  $(basename "$0") \\
    --local-ligand lig.sdf \\
    --ligand LIG \\
    --charge-method bcc \\
    --charge 0 \\
    --outdir prep_lig_bcc

  # RESP charges, extract ligand JZ4 from RCSB 3HTB
  $(basename "$0") \\
    --pdb 3HTB --ligand JZ4 \\
    --charge-method resp \\
    --charge 0 --mult 1 --resp-opt yes \\
    --nproc 32 --mem 40GB \\
    --outdir prep_JZ4_resp

  # parmed pathway (AMBER ff14SB water model via tleap)
  $(basename "$0") \\
    --local-ligand lig.mol2 --ligand LIG \\
    --charge-method bcc --charge 0 \\
    --ff ff14SB \\
    --outdir prep_lig_parmed
EOF_USAGE
}

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
die() { echo "Error: $*" >&2; exit 1; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"; }
optional_cmd() { command -v "$1" >/dev/null 2>&1; }
cleanup_on_error() { echo "Workflow failed." >&2; }
trap cleanup_on_error ERR

is_yes() {
    case "$(echo "${1:-}" | tr '[:upper:]' '[:lower:]')" in
        yes|y|true|1) return 0 ;;
        *) return 1 ;;
    esac
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--pdb)            PROID="$2"; INPUT_MODE="rcsb"; shift 2 ;;
        -l|--ligand)         LIGID="$2"; LIGID_USER_SET=1; shift 2 ;;
        --local-ligand)      LOCAL_LIGAND_PDB="$2"; LOCAL_LIGAND_GIVEN=1; shift 2 ;;
        --lig-chain)         LIGCHAIN="$2"; shift 2 ;;
        --lig-resi)          LIGRESI="$2"; shift 2 ;;
        -f|--ff)             FF="$2"; shift 2 ;;
        -w|--water)          WATER="$2"; shift 2 ;;
        -b|--boxtype)        BOXTYPE="$2"; shift 2 ;;
        -d|--distance)       BOXDIST="$2"; shift 2 ;;
        --ions-mdp)          IONS_MDP="$2"; shift 2 ;;
        --em-mdp)            EM_MDP="$2"; shift 2 ;;
        --outdir)            OUTDIR="$2"; shift 2 ;;
        --verbose)           VERBOSE=1; shift ;;
        --skip-solvate-minimize) RUN_SOLVATE_MINIMIZE=0; shift ;;
        --use-parmed)        USE_PARMED=1; shift ;;
        --charge-method)     CHARGE_METHOD="$(echo "$2" | tr '[:upper:]' '[:lower:]')"; shift 2 ;;
        -q|--charge)         LIG_CHARGE="$2"; shift 2 ;;
        -m|--mult)           LIG_MULT="$2"; shift 2 ;;
        --resp-opt)          RESP_OPT="$(echo "$2" | tr '[:upper:]' '[:lower:]')"; shift 2 ;;
        --nproc)             NPROC="$2"; shift 2 ;;
        --mem)               MEM="$2"; shift 2 ;;
        --basis-opt)         BASIS_OPT="$2"; shift 2 ;;
        --basis-esp)         BASIS_ESP="$2"; shift 2 ;;
        -h|--help)           usage; exit 0 ;;
        *)                   die "Unknown option: $1" ;;
    esac
done

# ── Input validation ──────────────────────────────────────────────────────────
if [[ "$INPUT_MODE" == "local" ]] && (( LOCAL_LIGAND_GIVEN == 0 )); then
    die "Provide --local-ligand FILE or --pdb PDBID to specify the ligand source"
fi

if (( LOCAL_LIGAND_GIVEN )); then
    [[ -n "$LOCAL_LIGAND_PDB" ]] || die "--local-ligand requires a file path"
    [[ -f "$LOCAL_LIGAND_PDB" ]]  || die "Local ligand file not found: $LOCAL_LIGAND_PDB"
    LOCAL_LIGAND_PDB="$(realpath "$LOCAL_LIGAND_PDB")"
fi

[[ "$LIGID_USER_SET" -eq 1 ]] || die "--ligand RESNAME is required"

case "$RESP_OPT" in
    yes|y|no|n) ;;
    *) die "--resp-opt must be yes or no" ;;
esac

[[ "$LIG_CHARGE" =~ ^-?[0-9]+$ ]] || die "--charge must be an integer, got: $LIG_CHARGE"
[[ "$LIG_MULT"   =~ ^[0-9]+$  ]]  || die "--mult must be a positive integer, got: $LIG_MULT"
[[ "$NPROC"      =~ ^[0-9]+$  ]]  || die "--nproc must be a positive integer, got: $NPROC"

if [[ -n "$PROID" ]]; then
    [[ -z "$OUTDIR" ]] && OUTDIR="prep_${PROID}_${LIGID}_${CHARGE_METHOD}"
else
    [[ -z "$OUTDIR" ]] && OUTDIR="prep_lig_${LIGID}_${CHARGE_METHOD}"
fi

if [[ -z "$IONS_MDP" ]]; then
    if [[ -f "${MDPDIR}/ions.mdp" ]]; then IONS_MDP="${MDPDIR}/ions.mdp"
    else IONS_MDP="${MDPDIR_BUNDLED}/ions.mdp"; fi
fi
if [[ -z "$EM_MDP" ]]; then
    if [[ -f "${MDPDIR}/em.mdp" ]]; then EM_MDP="${MDPDIR}/em.mdp"
    else EM_MDP="${MDPDIR_BUNDLED}/em.mdp"; fi
fi

for cmd in awk grep sed gmx obabel antechamber parmchk2 tleap tr sort uniq head tail wc realpath; do
    require_cmd "$cmd"
done
[[ "$INPUT_MODE" == "rcsb" ]] && require_cmd wget

# If the chosen FF has no GROMACS force field directory, auto-switch to parmed.
_gmx_prefix=$(gmx --version 2>/dev/null | awk '/Data prefix:/{print $NF}')
if [[ ! -d "${_gmx_prefix}/share/gromacs/top/${FF}.ff" ]]; then
    if [[ "$USE_PARMED" -eq 1 ]]; then
        :
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Note: '${FF}' is not a GROMACS force field — automatically switching to tleap+parmed pathway"
        USE_PARMED=1
    fi
fi
unset _gmx_prefix

[[ "$USE_PARMED" -eq 0 ]] && require_cmd amb2gro_top_gro.py

if [[ "$CHARGE_METHOD" == "resp" ]]; then
    require_cmd g16
fi

IONS_MDP="$(realpath "$IONS_MDP")"
EM_MDP="$(realpath "$EM_MDP")"

if [[ "$RUN_SOLVATE_MINIMIZE" -eq 1 ]]; then
    [[ -f "$IONS_MDP" ]] || die "ions.mdp not found: $IONS_MDP"
    [[ -f "$EM_MDP"   ]] || die "em.mdp not found: $EM_MDP"
fi

mkdir -p "$OUTDIR"
cd "$OUTDIR"

LOGFILE="prepare_lig_${LIGID}_${CHARGE_METHOD}.log"
exec > >(tee -a "$LOGFILE") 2>&1

if [[ "$INPUT_MODE" == "rcsb" ]]; then
    log "Input mode:              RCSB PDB extraction (${PROID})"
else
    log "Input mode:              local ligand file"
fi
log "Ligand name:             $LIGID"
log "Charge method:           $CHARGE_METHOD"
log "Ligand net charge:       $LIG_CHARGE"
log "Ligand multiplicity:     $LIG_MULT"
if [[ "$CHARGE_METHOD" == "resp" ]]; then
    log "RESP Gaussian optimize:  $RESP_OPT"
    log "RESP opt level:          $BASIS_OPT"
    log "RESP ESP level:          $BASIS_ESP"
    log "Gaussian resources:      nproc=$NPROC mem=$MEM"
fi
log "Topology pathway:        $( [[ "$USE_PARMED" -eq 1 ]] && echo "tleap+parmed" || echo "tleap+amb2gro" )"
log "Force field:             $FF"
log "Water model:             $WATER"
log "Box type:                $BOXTYPE"
log "Box padding (nm):        $BOXDIST"
log "Output directory:        $(pwd)"
log "Solvate/minimize:        $RUN_SOLVATE_MINIMIZE"

# ── File name variables ───────────────────────────────────────────────────────
RAW_PDB="${PROID:-ligand_rcsb}.pdb"
LIGAND_PDB="${LIGID}.pdb"

LIG_H_MOL2="${LIGID}_h.mol2"
LIG_XYZ="${LIGID}.xyz"
LIG_COORD="${LIGID}.coord"
AC_MOL2="AC_${LIGID}.mol2"
AC_FRCMOD="AC_${LIGID}.frcmod"
AC_PRMTOP="AC_${LIGID}.prmtop"
AC_INPCRD="AC_${LIGID}.inpcrd"

LIG_GMX_TOP="${LIGID}_GMX.top"
LIG_GMX_GRO="${LIGID}_GMX.gro"
LIG_GMX_PDB="${LIGID}_GMX.pdb"

LIG_ATOMTYPES_ITP="${LIGID}_GMX_GAFF.itp"
LIG_MOL_ITP="${LIGID}_GMX.itp"
LIG_POSRE_ITP="posre_${LIGID}.itp"

FINAL_TOP="topol.top"
INDEX_NDX="index.ndx"
COMPLEX_GRO="$LIG_GMX_GRO"

SEL_CHAIN=""
SEL_RESI=""
LIG_MOLNAME=""

# ── Utility functions ─────────────────────────────────────────────────────────
check_g16_normal() {
    local logfile="$1"
    [[ -s "$logfile" ]] || die "Gaussian log missing/empty: $logfile"
    if ! grep -q "Normal termination" "$logfile"; then
        die "Gaussian did not terminate normally. Check: $logfile"
    fi
}

check_charge_sum_mol2() {
    local mol2="$1"
    [[ -s "$mol2" ]] || die "Charge-check mol2 missing: $mol2"
    awk '
        /^@<TRIPOS>ATOM/ {in_atom=1; next}
        /^@<TRIPOS>/ && in_atom {in_atom=0}
        in_atom && NF >= 9 {s += $9; n++}
        END {
            if (n == 0) { printf "  Charge check: no MOL2 atoms parsed\n"; exit 0 }
            printf "  MOL2 atom count = %d\n", n
            printf "  Sum charge     = %.6f\n", s
        }
    ' "$mol2"
}

copy_local_ligand() {
    local src="$LOCAL_LIGAND_PDB"
    local ext="${src##*.}"; ext="${ext,,}"
    if [[ "$ext" == "pdb" ]]; then
        cp "$src" "$LIGAND_PDB"
    else
        log "Local ligand is ${ext} format; converting to PDB via obabel"
        obabel "$src" -O "$LIGAND_PDB" 2>/dev/null \
            || die "obabel failed to convert local ligand (${ext}) to PDB"
        [[ -s "$LIGAND_PDB" ]] || die "obabel produced empty PDB from local ligand"
    fi
    log "Local ligand copied to ${LIGAND_PDB}"
}

detect_ligand_sites() {
    if (( LOCAL_LIGAND_GIVEN )); then
        log "Local ligand provided: skipping RCSB site detection"
        return 0
    fi

    awk -v lig="$LIGID" '
        /^HETATM/ {
            resn  = substr($0,18,3)
            alt   = substr($0,17,1)
            chain = substr($0,22,1)
            resi  = substr($0,23,4)
            gsub(/ /,"",resn); gsub(/ /,"",alt); gsub(/ /,"",chain); gsub(/ /,"",resi)
            if (resn == lig && (alt == "" || alt == "A")) print chain " " resi
        }
    ' "$RAW_PDB" | sort -u > ligand_sites.tsv

    [[ -s ligand_sites.tsv ]] || die "No ligand ${LIGID} found in ${RAW_PDB}"
    log "Detected ligand site(s):"
    cat ligand_sites.tsv

    if [[ -n "$LIGCHAIN" || -n "$LIGRESI" ]]; then
        [[ -n "$LIGCHAIN" && -n "$LIGRESI" ]] || die "Use both --lig-chain and --lig-resi together"
        awk -v c="$LIGCHAIN" -v r="$LIGRESI" '($1==c && $2==r){found=1} END{exit !found}' ligand_sites.tsv \
            || die "Requested ligand chain/resi not found"
        SEL_CHAIN="$LIGCHAIN"
        SEL_RESI="$LIGRESI"
    else
        SEL_CHAIN="$(awk 'NR==1{print $1}' ligand_sites.tsv)"
        SEL_RESI="$(awk 'NR==1{print $2}' ligand_sites.tsv)"
        if [[ $(wc -l < ligand_sites.tsv) -gt 1 ]]; then
            log "Multiple ligand copies found; defaulting to first: chain=${SEL_CHAIN} resi=${SEL_RESI}"
        fi
    fi
    log "Selected ligand: ${LIGID}, chain=${SEL_CHAIN}, resi=${SEL_RESI}"
}

load_ligand() {
    if [[ "$INPUT_MODE" == "rcsb" ]]; then
        log "Downloading ${PROID} from RCSB"
        wget -O "$RAW_PDB" "https://files.rcsb.org/download/${PROID}.pdb"
        [[ -s "$RAW_PDB" ]] || die "Downloaded PDB is empty"
        detect_ligand_sites
        if (( LOCAL_LIGAND_GIVEN )); then
            copy_local_ligand
        else
            log "Extracting ligand ${LIGID} from ${RAW_PDB}"
            awk -v lig="$LIGID" -v schain="$SEL_CHAIN" -v sresi="$SEL_RESI" '
                /^HETATM/ {
                    resn  = substr($0,18,3); alt   = substr($0,17,1)
                    chain = substr($0,22,1); resi  = substr($0,23,4)
                    gsub(/ /,"",resn); gsub(/ /,"",alt); gsub(/ /,"",chain); gsub(/ /,"",resi)
                    if (resn == lig && chain == schain && resi == sresi && (alt == "" || alt == "A")) print
                }
            ' "$RAW_PDB" > "$LIGAND_PDB"
            [[ -s "$LIGAND_PDB" ]] || die "Ligand extraction from ${RAW_PDB} failed"
            log "Ligand extracted to ${LIGAND_PDB}"
        fi
    else
        copy_local_ligand
    fi
}

graft_coords_to_mol2() {
    local ref_mol2="$1"
    local target_mol2="$2"
    log "Grafting original coordinates from ${ref_mol2} into ${target_mol2}"
    python3 - "$ref_mol2" "$target_mol2" <<'PYEOF'
import sys

ref_path, tgt_path = sys.argv[1], sys.argv[2]

def mol2_atom_coords(path):
    coords = []
    in_atom = False
    with open(path) as f:
        for line in f:
            if line.startswith('@<TRIPOS>ATOM'):
                in_atom = True; continue
            if line.startswith('@<TRIPOS>') and in_atom:
                break
            if in_atom and line.strip():
                p = line.split()
                coords.append((float(p[2]), float(p[3]), float(p[4])))
    return coords

orig = mol2_atom_coords(ref_path)

lines = open(tgt_path).readlines()
out = []
in_atom = False
idx = 0
for line in lines:
    if line.startswith('@<TRIPOS>ATOM'):
        in_atom = True; out.append(line); continue
    if line.startswith('@<TRIPOS>') and in_atom:
        in_atom = False
    if in_atom and line.strip():
        p = line.split()
        if idx < len(orig):
            p[2], p[3], p[4] = f'{orig[idx][0]:.4f}', f'{orig[idx][1]:.4f}', f'{orig[idx][2]:.4f}'
        idx += 1
        out.append('{:>7s} {:<8s} {:>10.4f} {:>10.4f} {:>10.4f} {:<9s} {:>4s} {:<8s} {:>10.6f}\n'.format(
            p[0], p[1], float(p[2]), float(p[3]), float(p[4]),
            p[5], p[6], p[7], float(p[8])))
        continue
    out.append(line)
if idx != len(orig):
    print(f'ERROR: atom count mismatch: ref={len(orig)}, target={idx}', file=sys.stderr)
    sys.exit(1)
with open(tgt_path, 'w') as f:
    f.writelines(out)
print(f'  Coordinate graft complete: {idx} atoms restored from {ref_path}')
PYEOF
    [[ $? -eq 0 ]] || die "Coordinate grafting failed (atom count mismatch between ${1} and ${2})"
}

prepare_ligand_antechamber() {
    log "Ligand charge route: antechamber (method=${CHARGE_METHOD})"
    log "Generating ligand mol2 with hydrogens"
    obabel "$LIGAND_PDB" -O "$LIG_H_MOL2" -h
    [[ -s "$LIG_H_MOL2" ]] || die "Open Babel failed to create $LIG_H_MOL2"

    log "Running antechamber (charge method: ${CHARGE_METHOD})"
    antechamber \
        -i "$LIG_H_MOL2" -fi mol2 \
        -o "$AC_MOL2" -fo mol2 \
        -c "$CHARGE_METHOD" -s 2 \
        -rn "$LIGID" -at gaff2 \
        -nc "$LIG_CHARGE" -pf yes
    [[ -s "$AC_MOL2" ]] || die "Antechamber did not create $AC_MOL2"
    check_charge_sum_mol2 "$AC_MOL2"
}

prepare_ligand_resp() {
    log "Ligand charge route: RESP through Gaussian16 + antechamber"

    log "Saving original ligand geometry as mol2 (coordinate reference)"
    obabel "$LIGAND_PDB" -O "$LIG_H_MOL2" -h
    [[ -s "$LIG_H_MOL2" ]] || die "Open Babel failed to create $LIG_H_MOL2"

    log "Generating ligand XYZ with hydrogens for Gaussian"
    obabel "$LIGAND_PDB" -O "$LIG_XYZ" -h
    [[ -s "$LIG_XYZ" ]] || die "Open Babel failed to create $LIG_XYZ"
    tail -n +3 "$LIG_XYZ" > "$LIG_COORD"
    [[ -s "$LIG_COORD" ]] || die "Failed to create coordinate file $LIG_COORD"

    if is_yes "$RESP_OPT"; then
        log "Writing Gaussian optimization input: ${LIGID}_opt.com"
        cat > "${LIGID}_opt.com" <<EOF_OPT
%mem=${MEM}
%nprocshared=${NPROC}
%nosave
%chk=${LIGID}_opt.chk
#p opt nosymm ${BASIS_OPT}

${LIGID} optimization

${LIG_CHARGE} ${LIG_MULT}
$(cat "$LIG_COORD")


EOF_OPT

        log "Writing Gaussian ESP input from optimized checkpoint: ${LIGID}_ESP.com"
        cat > "${LIGID}_ESP.com" <<EOF_ESP
%mem=${MEM}
%nprocshared=${NPROC}
%nosave
%chk=${LIGID}_ESP.chk
#p nosymm ${BASIS_ESP} guess=read geom=check Pop=(MK) IOp(6/50=1)

${LIGID} ESP from optimized geometry

${LIG_CHARGE} ${LIG_MULT}

${LIGID}.esp

EOF_ESP

        log "Running Gaussian optimization"
        g16 < "${LIGID}_opt.com" > "${LIGID}_opt.log"
        check_g16_normal "${LIGID}_opt.log"

        log "Running Gaussian ESP from optimized checkpoint"
        [[ -s "${LIGID}_opt.chk" ]] || die "Missing Gaussian checkpoint: ${LIGID}_opt.chk"
        cp "${LIGID}_opt.chk" "${LIGID}_ESP.chk"
        g16 < "${LIGID}_ESP.com" > "${LIGID}_ESP.log"
        check_g16_normal "${LIGID}_ESP.log"
    else
        log "Writing Gaussian ESP input without optimization: ${LIGID}_ESP_noopt.com"
        cat > "${LIGID}_ESP_noopt.com" <<EOF_NOOPT
%mem=${MEM}
%nprocshared=${NPROC}
%nosave
%chk=${LIGID}_ESP_noopt.chk
#p nosymm ${BASIS_ESP} Pop=(MK) IOp(6/50=1)

${LIGID} ESP no optimization

${LIG_CHARGE} ${LIG_MULT}
$(cat "$LIG_COORD")

${LIGID}.esp

EOF_NOOPT

        log "Running Gaussian ESP without optimization"
        g16 < "${LIGID}_ESP_noopt.com" > "${LIGID}_ESP_noopt.log"
        check_g16_normal "${LIGID}_ESP_noopt.log"
    fi

    [[ -s "${LIGID}.esp" ]] || die "Gaussian did not produce ${LIGID}.esp"

    log "Running antechamber RESP using Gaussian ESP"
    antechamber \
        -i "${LIGID}.esp" -fi gesp \
        -o "$AC_MOL2" -fo mol2 \
        -c resp -s 2 \
        -rn "$LIGID" -at gaff2 \
        -nc "$LIG_CHARGE" -pf yes
    [[ -s "$AC_MOL2" ]] || die "Antechamber RESP did not create $AC_MOL2"

    graft_coords_to_mol2 "$LIG_H_MOL2" "$AC_MOL2"
    check_charge_sum_mol2 "$AC_MOL2"
}

prepare_ligand() {
    case "$CHARGE_METHOD" in
        resp) prepare_ligand_resp ;;
        *)    prepare_ligand_antechamber ;;
    esac

    log "Running parmchk2"
    parmchk2 -i "$AC_MOL2" -f mol2 -o "$AC_FRCMOD"
    [[ -s "$AC_FRCMOD" ]] || die "parmchk2 did not create $AC_FRCMOD"

    if [[ "$USE_PARMED" -eq 1 ]]; then
        LIG_MOLNAME="$LIGID"
        log "Parmed pathway: ligand charges ready; solvated system will be built by tleap+parmed"
        return 0
    fi

    log "Running tleap (ligand-only parameterization)"
    tleap -f - <<EOF_TLEAP
source leaprc.gaff2
lig = loadmol2 ${AC_MOL2}
loadamberparams ${AC_FRCMOD}
saveamberparm lig ${AC_PRMTOP} ${AC_INPCRD}
quit
EOF_TLEAP

    [[ -s "$AC_PRMTOP" ]] || die "Ligand prmtop missing"
    [[ -s "$AC_INPCRD" ]] || die "Ligand inpcrd missing"

    log "Converting ligand AMBER files to GROMACS"
    amb2gro_top_gro.py -p "$AC_PRMTOP" -c "$AC_INPCRD" -t "$LIG_GMX_TOP" -g "$LIG_GMX_GRO" -b "$LIG_GMX_PDB"

    [[ -s "$LIG_GMX_TOP" ]] || die "Ligand GROMACS top missing"
    [[ -s "$LIG_GMX_GRO" ]] || die "Ligand GROMACS gro missing"

    log "Splitting ligand topology"
    sed -n '/^\[ *atomtypes/,/^\[ *moleculetype/{/^\[ *moleculetype/b;p;}' "$LIG_GMX_TOP" > "$LIG_ATOMTYPES_ITP"
    sed -n '/^\[ *moleculetype/,/^\[ *system/{/^\[ *system/b;p;}' "$LIG_GMX_TOP" > "$LIG_MOL_ITP"
    [[ -s "$LIG_ATOMTYPES_ITP" ]] || die "Ligand atomtypes include is empty"
    [[ -s "$LIG_MOL_ITP" ]] || die "Ligand molecule include is empty"

    cat >> "$LIG_MOL_ITP" <<EOF_POSRE

; Include Position restraint file
#ifdef POSRES_LIG
#include "./${LIG_POSRE_ITP}"
#endif
EOF_POSRE

    log "Generating ligand position restraints (local indices)"
    echo "0" | gmx genrestr -f "$LIG_GMX_GRO" -o "$LIG_POSRE_ITP" -fc 1000 1000 1000
    [[ -s "$LIG_POSRE_ITP" ]] || die "Ligand position restraint file was not created"

    log "Detecting ligand molecule name from [ moleculetype ]"
    LIG_MOLNAME="$(
        awk '
            /^\[ *moleculetype *\]/ {
                if (getline <= 0) exit
                while ($0 ~ /^[[:space:]]*;/ || $0 ~ /^[[:space:]]*$/) {
                    if (getline <= 0) exit
                }
                print $1
                exit
            }
        ' "$LIG_MOL_ITP"
    )"
    if [[ -z "$LIG_MOLNAME" || "$LIG_MOLNAME" == ";" ]]; then
        log "Warning: ligand moleculetype parser returned '${LIG_MOLNAME:-EMPTY}', falling back to LIGID=${LIGID}"
        LIG_MOLNAME="$LIGID"
    fi
    log "Ligand molecule name detected: ${LIG_MOLNAME}"
}

build_final_topology() {
    log "Building final topol.top (ligand-only)"
    {
        echo "; Include forcefield parameters"
        echo "#include \"${FF}.ff/forcefield.itp\""
        echo "#include \"./${LIG_ATOMTYPES_ITP}\""
        echo
        echo "#include \"./${LIG_MOL_ITP}\""
        echo
        echo "; Include water topology"
        case "$WATER" in
            tip3p)   echo "#include \"${FF}.ff/tip3p.itp\"" ;;
            tip4p)   echo "#include \"${FF}.ff/tip4p.itp\"" ;;
            tip4pew) echo "#include \"${FF}.ff/tip4pew.itp\"" ;;
            spc)     echo "#include \"${FF}.ff/spc.itp\"" ;;
            spce)    echo "#include \"${FF}.ff/spce.itp\"" ;;
            *)       log "Warning: unknown water model '${WATER}', defaulting to tip3p.itp"
                     echo "#include \"${FF}.ff/tip3p.itp\"" ;;
        esac
        echo
        cat <<'EOF_WATER_POSRE'
#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif
EOF_WATER_POSRE
        echo
        echo "; Include topology for ions"
        echo "#include \"${FF}.ff/ions.itp\""
        echo
        echo "[ system ]"
        echo "; Name"
        echo "Ligand in water"
        echo
        echo "[ molecules ]"
        echo "; Compound        #mols"
        printf "%-18s %5d\n" "$LIG_MOLNAME" 1
    } > "$FINAL_TOP"
    [[ -s "$FINAL_TOP" ]] || die "Failed to build ${FINAL_TOP}"
}

generate_tc_index() {
    local coord="$1"
    [[ -s "$coord" ]] || die "Cannot generate index; coordinate file missing: $coord"

    log "Generating index from ${coord}"
    gmx make_ndx -f "$coord" -o "$INDEX_NDX" <<EOF_NDX
q
EOF_NDX
    [[ -s "$INDEX_NDX" ]] || die "index.ndx was not created"

    local found_lig found_wat
    found_lig="$(awk -v lig="$LIGID" -v mol="$LIG_MOLNAME" '
        /^\[/ { name=$0; gsub(/^\[ */,"",name); gsub(/ *\]$/,"",name)
                if (name==lig || name==mol) print "yes" }
    ' "$INDEX_NDX" | head -1)"
    found_wat="$(awk '
        /^\[/ { name=$0; gsub(/^\[ */,"",name); gsub(/ *\]$/,"",name)
                if (name=="Water_and_ions") print "yes" }
    ' "$INDEX_NDX" | head -1)"

    [[ "$found_lig" == "yes" ]] || log "Warning: group [${LIGID}] not found in index; check tc-grps"
    [[ "$found_wat" == "yes" ]] || log "Warning: group [Water_and_ions] not found in index; check tc-grps"

    log "Index ready: ${INDEX_NDX}"
    log "Use in MDP: tc-grps = ${LIG_MOLNAME} Water_and_ions"
}

prepare_system_parmed() {
    local lig_prmtop="lig_solv.prmtop"
    local lig_inpcrd="lig_solv.inpcrd"
    local parmed_top="lig_solv_gmx.top"
    local parmed_gro="lig_solv_gmx.gro"
    local box_cmd box_dist_a

    case "$BOXTYPE" in
        dodecahedron|octahedron) box_cmd="solvateOct" ;;
        *) box_cmd="solvateBox" ;;
    esac
    box_dist_a=$(python3 -c "print(${BOXDIST} * 10)")

    log "Building solvated ligand system with tleap (GAFF2 + TIP3P)"
    tleap -f - <<EOF_TLEAP
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams ${AC_FRCMOD}
LIG = loadMol2 ${AC_MOL2}
${box_cmd} LIG TIP3PBOX ${box_dist_a}
addions LIG Na+ 0
addions LIG Cl- 0
saveamberparm LIG ${lig_prmtop} ${lig_inpcrd}
quit
EOF_TLEAP
    [[ -s "$lig_prmtop" ]] || die "tleap did not produce ${lig_prmtop}"
    [[ -s "$lig_inpcrd" ]] || die "tleap did not produce ${lig_inpcrd}"

    log "Converting AMBER → GROMACS with parmed"
    python3 - <<EOF_PARMED
from parmed import load_file

amber = load_file('${lig_prmtop}', '${lig_inpcrd}')

for res in amber.residues:
    rn = res.name.strip()
    if rn in ('WAT', 'HOH'):
        res.name = 'SOL'
    elif rn in ('Na+', 'NA+', 'Na'):
        res.name = 'NA'
    elif rn in ('Cl-', 'CL-', 'Cl'):
        res.name = 'CL'

amber.save('${parmed_top}', format='gromacs', overwrite=True)
amber.save('${parmed_gro}', overwrite=True)
print('parmed conversion complete')
EOF_PARMED
    [[ -s "$parmed_top" ]] || die "parmed did not produce ${parmed_top}"
    [[ -s "$parmed_gro" ]] || die "parmed did not produce ${parmed_gro}"

    # Ligand position restraints: use mol2 as single-molecule file for local indices
    log "Generating ligand position restraints"
    obabel "$AC_MOL2" -O lig_posre_tmp.pdb 2>/dev/null
    [[ -s "lig_posre_tmp.pdb" ]] || die "obabel could not convert ligand mol2 to PDB for genrestr"
    echo "0" | gmx genrestr -f lig_posre_tmp.pdb -o "${LIG_POSRE_ITP}" -fc 1000 1000 1000
    rm -f lig_posre_tmp.pdb
    [[ -s "${LIG_POSRE_ITP}" ]] || die "${LIG_POSRE_ITP} was not created"

    # Inject POSRES_LIG and POSRES_WATER into topology
    log "Injecting POSRES directives into ${parmed_top}"
    python3 - "$parmed_top" "$LIGID" "$LIG_POSRE_ITP" <<'EOF_INJECT'
import re, sys

top_file, lig_name, posre_lig = sys.argv[1], sys.argv[2], sys.argv[3]

with open(top_file) as f:
    lines = f.readlines()

result = []
i = 0
current_mol = None

def posres_block(kind):
    if kind == 'lig':
        return ['\n#ifdef POSRES_LIG\n', f'#include "./{posre_lig}"\n', '#endif\n\n']
    elif kind == 'water':
        return ['\n#ifdef POSRES_WATER\n',
                '; Position restraint for each water oxygen\n',
                '[ position_restraints ]\n',
                ';  i funct       fcx        fcy        fcz\n',
                '   1    1       1000       1000       1000\n',
                '#endif\n\n']
    return []

while i < len(lines):
    line = lines[i]
    new_section = re.match(r'^\[ *(moleculetype|system)', line)

    if new_section and current_mol is not None:
        if current_mol == lig_name:
            result.extend(posres_block('lig'))
        elif current_mol == 'SOL':
            result.extend(posres_block('water'))
        current_mol = None

    result.append(line)

    if re.match(r'^\[ *moleculetype', line):
        i += 1
        while i < len(lines):
            result.append(lines[i])
            s = lines[i].strip()
            if s and not s.startswith(';'):
                current_mol = s.split()[0]
                i += 1
                break
            i += 1
        continue

    i += 1

with open(top_file, 'w') as f:
    f.writelines(result)
print(f'POSRES injected ({lig_name} → POSRES_LIG, SOL → POSRES_WATER)')
EOF_INJECT
    [[ $? -eq 0 ]] || die "POSRES injection into parmed topology failed"

    FINAL_TOP="$parmed_top"
    COMPLEX_GRO="$parmed_gro"

    generate_tc_index "$COMPLEX_GRO"
}

solvate_and_minimize() {
    if [[ "$RUN_SOLVATE_MINIMIZE" -eq 0 ]]; then
        log "Skipping solvation, ionization, and minimization by request"
        return 0
    fi

    local em_coord

    if [[ "$USE_PARMED" -eq 0 ]]; then
        log "Creating box"
        gmx editconf -f "$LIG_GMX_GRO" -o newbox.gro -bt "$BOXTYPE" -d "$BOXDIST"

        log "Solvating"
        gmx solvate -cp newbox.gro -cs spc216.gro -p "$FINAL_TOP" -o solv.gro

        log "Preparing ions"
        gmx grompp -f "$IONS_MDP" -c solv.gro -r solv.gro -p "$FINAL_TOP" -o ions.tpr -maxwarn 1

        log "Adding ions"
        echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p "$FINAL_TOP" -pname NA -nname CL -neutral

        generate_tc_index "solv_ions.gro"
        em_coord="solv_ions.gro"
    else
        # parmed pathway: tleap already solvated and ionised; index already built
        em_coord="$COMPLEX_GRO"
    fi

    log "Preparing minimization"
    gmx grompp -f "$EM_MDP" -c "$em_coord" -r "$em_coord" -p "$FINAL_TOP" -o em.tpr

    log "Running minimization"
    gmx mdrun -v -deffnm em
}

write_amber_inputs() {
    log "Writing AMBER MD input files"

    cat > amber_min.in <<'EOF_MIN'
Minimization - steepest descent then CG, ligand restrained
 &cntrl
  imin         = 1,
  ntmin        = 1,
  maxcyc       = 10000,
  ncyc         = 5000,
  ntpr         = 100,
  ntwx         = 0,
  cut          = 12.0,
  ntb          = 1,
  ntr          = 1,
  restraint_wt = 10.0,
  restraintmask = '!(:WAT,Na+,Cl-)',
 /
EOF_MIN

    cat > amber_nvt.in <<'EOF_NVT'
NVT equilibration - Langevin dynamics, ligand restrained (1 ns)
 &cntrl
  imin         = 0,
  irest        = 0,
  ntx          = 1,
  nstlim       = 500000,
  dt           = 0.002,
  ntf          = 2,
  ntc          = 2,
  cut          = 12.0,
  ntb          = 1,
  ntt          = 3,
  gamma_ln     = 1.0,
  temp0        = 310.0,
  tempi        = 310.0,
  ntpr         = 1000,
  ntwx         = 5000,
  ntwr         = 5000,
  iwrap        = 1,
  ntr          = 1,
  restraint_wt = 10.0,
  restraintmask = '!(:WAT,Na+,Cl-)',
 /
EOF_NVT

    cat > amber_npt.in <<'EOF_NPT'
NPT equilibration - Langevin dynamics + MC barostat, ligand restrained (1 ns)
 &cntrl
  imin         = 0,
  irest        = 1,
  ntx          = 5,
  nstlim       = 500000,
  dt           = 0.002,
  ntf          = 2,
  ntc          = 2,
  cut          = 12.0,
  ntb          = 2,
  ntt          = 3,
  gamma_ln     = 1.0,
  temp0        = 310.0,
  ntp          = 1,
  barostat     = 2,
  pres0        = 1.01325,
  mcbarint     = 100,
  ntpr         = 1000,
  ntwx         = 5000,
  ntwr         = 5000,
  iwrap        = 1,
  ntr          = 1,
  restraint_wt = 10.0,
  restraintmask = '!(:WAT,Na+,Cl-)',
 /
EOF_NPT

    cat > amber_md.in <<'EOF_MD'
Production MD - Langevin dynamics + MC barostat, no restraints (100 ns)
 &cntrl
  imin         = 0,
  irest        = 1,
  ntx          = 5,
  nstlim       = 50000000,
  dt           = 0.002,
  ntf          = 2,
  ntc          = 2,
  cut          = 12.0,
  ntb          = 2,
  ntt          = 3,
  gamma_ln     = 1.0,
  temp0        = 310.0,
  ntp          = 1,
  barostat     = 2,
  pres0        = 1.01325,
  mcbarint     = 100,
  ntpr         = 5000,
  ntwx         = 5000,
  ntwr         = 25000,
  iwrap        = 1,
 /
EOF_MD
}

cleanup_files() {
    if [[ "$VERBOSE" -eq 1 ]]; then return 0; fi
    log "Removing intermediates (use --verbose to keep all files)"

    rm -f "$LIGAND_PDB" "${RAW_PDB}" ligand_sites.tsv
    rm -f "$LIG_H_MOL2" "$LIG_XYZ" "$LIG_COORD" "$AC_MOL2" "$AC_FRCMOD"
    rm -f ANTECHAMBER_*.AC ANTECHAMBER_*.AC0 ATOMTYPE.INF \
          sqm.in sqm.out sqm.pdb leap.log mdinfo mdout.mdp

    if [[ "$USE_PARMED" -eq 0 ]]; then
        rm -f "$AC_PRMTOP" "$AC_INPCRD"
        rm -f "$LIG_GMX_TOP" "$LIG_GMX_PDB"
        rm -f newbox.gro solv.gro ions.tpr
        if [[ "$RUN_SOLVATE_MINIMIZE" -eq 1 ]]; then
            rm -f solv_ions.gro "$LIG_GMX_GRO" em.tpr em.edr em.log em.trr
        fi
    else
        if [[ "$RUN_SOLVATE_MINIMIZE" -eq 1 ]]; then
            rm -f "$COMPLEX_GRO" em.tpr em.edr em.log em.trr
        fi
    fi

    rm -f \#*\#
}

write_summary() {
    log "Writing preparation summary"
    {
        cat <<EOF_SUMMARY
Ligand name:           ${LIGID}
Input mode:            $( [[ "$INPUT_MODE" == "rcsb" ]] && echo "RCSB extraction (${PROID})" || echo "local file" )
Topology pathway:      $( [[ "$USE_PARMED" -eq 1 ]] && echo "tleap+parmed" || echo "tleap+amb2gro" )
Force field:           ${FF}
Water model:           ${WATER}
Charge method:         ${CHARGE_METHOD}
Ligand charge:         ${LIG_CHARGE}
Ligand multiplicity:   ${LIG_MULT}
RESP optimization:     ${RESP_OPT}

--- GROMACS output ---
Final topology:        ${FINAL_TOP}
Coordinates:           $( [[ "$RUN_SOLVATE_MINIMIZE" -eq 1 ]] && echo "em.gro" || echo "$COMPLEX_GRO" )
Index file:            ${INDEX_NDX}
Production tc-grps:    ${LIG_MOLNAME} Water_and_ions
EOF_SUMMARY
        if [[ "$USE_PARMED" -eq 1 ]]; then
            cat <<EOF_AMBER

--- AMBER output ---
Topology (prmtop):     lig_solv.prmtop
Coordinates (inpcrd):  lig_solv.inpcrd
Minimization input:    amber_min.in
NVT input:             amber_nvt.in
NPT input:             amber_npt.in
Production MD input:   amber_md.in

Example run (pmemd.cuda):
  pmemd.cuda -O -i amber_min.in -p lig_solv.prmtop -c lig_solv.inpcrd -ref lig_solv.inpcrd -o min.out -r min.rst -x min.nc
  pmemd.cuda -O -i amber_nvt.in -p lig_solv.prmtop -c min.rst          -ref lig_solv.inpcrd -o nvt.out -r nvt.rst -x nvt.nc
  pmemd.cuda -O -i amber_npt.in -p lig_solv.prmtop -c nvt.rst          -ref lig_solv.inpcrd -o npt.out -r npt.rst -x npt.nc
  pmemd.cuda -O -i amber_md.in  -p lig_solv.prmtop -c npt.rst                               -o md.out  -r md.rst  -x md.nc
EOF_AMBER
        fi
        echo
        echo "Log file:              ${LOGFILE}"
    } > preparation_summary.txt
}

# ── Main flow ─────────────────────────────────────────────────────────────────
load_ligand

if [[ "$USE_PARMED" -eq 0 ]]; then
    prepare_ligand
    build_final_topology
    generate_tc_index "$LIG_GMX_GRO"
else
    prepare_ligand
    prepare_system_parmed
    write_amber_inputs
fi

solvate_and_minimize
cleanup_files
write_summary

log "All steps completed successfully"
