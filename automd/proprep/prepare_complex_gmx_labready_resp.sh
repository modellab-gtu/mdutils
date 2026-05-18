#!/usr/bin/env bash
set -Eeuo pipefail

###############################################################################
# prepare_complex_gmx_labready_resp.sh
#
# Robust protein-ligand GROMACS preparation workflow with selectable ligand
# charge route:
#   - AM1-BCC charges through AmberTools/antechamber
#   - RESP charges through Gaussian16 ESP + AmberTools/antechamber
#
# Main features:
#   - RCSB PDB or local protein/ligand input
#   - Ligand instance selection by residue name, chain, residue number
#   - Optional PDBFixer preprocessing
#   - Protein preparation through gmx pdb2gmx
#   - Ligand parameterization with GAFF2
#   - AM1-BCC or RESP ligand charges
#   - Optional Gaussian16 optimization before RESP ESP calculation
#   - Amber ligand topology conversion to GROMACS through amb2gro_top_gro.py
#   - Final topol.top construction, solvation, ionization, minimization
#  - Post-ionization index with Protein_<LIGID> group for tc-grps
###############################################################################

PROID="3HTB"
LIGID="JZ4"
LIGID_USER_SET=0
LOCAL_PROTEIN_PDB=""
LOCAL_LIGAND_PDB=""
INPUT_MODE="rcsb"
LIGCHAIN=""
LIGRESI=""
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
RUN_SOLVATE_MINIMIZE=1

# Ligand charge/RESP controls
CHARGE_METHOD="bcc"        # bcc or resp
LIG_CHARGE=0
LIG_MULT=1
RESP_OPT="yes"             # yes or no
NPROC=32
MEM="20GB"
BASIS_OPT="B3LYP/6-31G*"
BASIS_ESP="HF/6-31G*"

SCRIPT_START_DIR="$(pwd)"
MDPDIR="${SCRIPT_START_DIR}/MDP"

usage() {
cat <<EOF_USAGE
Usage:
  $(basename "$0") [options]

Input modes:
  RCSB mode:
    $(basename "$0") --pdb 3HTB --ligand JZ4 --charge-method bcc --charge 0

  Local mode:
    $(basename "$0") --local-protein prot_noH.pdb --local-ligand lig.pdb \\
      --ligand LIG --charge-method resp --charge 0 --resp-opt yes

General options:
  -p, --pdb ID                 PDB ID from RCSB, default: ${PROID}
  -l, --ligand RES             Ligand residue/molecule name, default: ${LIGID}
      --local-protein FILE     Local protein PDB file
      --local-ligand FILE      Local ligand PDB file
      --lig-chain CHAIN        Ligand chain ID for RCSB mode
      --lig-resi RESI          Ligand residue number for RCSB mode
  -f, --ff NAME                GROMACS protein force field, default: ${FF}
  -w, --water NAME             Water model, default: ${WATER}
  -b, --boxtype TYPE           Box type, default: ${BOXTYPE}
  -d, --distance NM            Solvent distance in nm, default: ${BOXDIST}
      --ions-mdp FILE          ions.mdp path, default: \$PWD/MDP/ions.mdp
      --em-mdp FILE            em.mdp path, default: \$PWD/MDP/em.mdp
      --outdir DIR             Output directory, default: prep_<PDB>_<LIG>
      --cleanup                Remove selected intermediates
      --pdbfixer               Run PDBFixer preprocessing if available/requested
      --ph VALUE               pH for PDBFixer hydrogen addition, default: ${PH}
      --skip-solvate-minimize  Stop after complex.gro, topol.top, and index.ndx
  -h, --help                   Show this help

Ligand charge options:
      --charge-method METHOD   bcc or resp, default: ${CHARGE_METHOD}
  -q, --charge INT             Net ligand charge, default: ${LIG_CHARGE}
  -m, --mult INT               Spin multiplicity for RESP/Gaussian, default: ${LIG_MULT}
      --resp-opt yes|no        For RESP: optimize with G16 before ESP, default: ${RESP_OPT}
      --nproc INT              Gaussian nprocshared, default: ${NPROC}
      --mem MEM                Gaussian memory, default: ${MEM}
      --basis-opt LEVEL        Gaussian optimization level, default: ${BASIS_OPT}
      --basis-esp LEVEL        Gaussian ESP level, default: ${BASIS_ESP}

Examples:
  # AM1-BCC charges, local files
  $(basename "$0") \\
    --local-protein prot_noH.pdb \\
    --local-ligand lig.pdb \\
    --ligand LIG \\
    --charge-method bcc \\
    --charge 0 \\
    --outdir prep_bcc

  # RESP charges with Gaussian optimization
  $(basename "$0") \\
    --local-protein prot_noH.pdb \\
    --local-ligand lig.pdb \\
    --ligand LIG \\
    --charge-method resp \\
    --charge 0 \\
    --mult 1 \\
    --resp-opt yes \\
    --nproc 32 \\
    --mem 40GB \\
    --outdir prep_resp_opt

  # RESP charges without Gaussian optimization, ESP only on input geometry
  $(basename "$0") \\
    --local-protein prot_noH.pdb \\
    --local-ligand lig.pdb \\
    --ligand LIG \\
    --charge-method resp \\
    --charge 0 \\
    --resp-opt no \\
    --outdir prep_resp_noopt
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
        -p|--pdb) PROID="$2"; shift 2 ;;
        -l|--ligand) LIGID="$2"; LIGID_USER_SET=1; shift 2 ;;
        --local-protein) LOCAL_PROTEIN_PDB="$2"; INPUT_MODE="local"; shift 2 ;;
        --local-ligand) LOCAL_LIGAND_PDB="$2"; INPUT_MODE="local"; shift 2 ;;
        --lig-chain) LIGCHAIN="$2"; shift 2 ;;
        --lig-resi) LIGRESI="$2"; shift 2 ;;
        -f|--ff) FF="$2"; shift 2 ;;
        -w|--water) WATER="$2"; shift 2 ;;
        -b|--boxtype) BOXTYPE="$2"; shift 2 ;;
        -d|--distance) BOXDIST="$2"; shift 2 ;;
        --ions-mdp) IONS_MDP="$2"; shift 2 ;;
        --em-mdp) EM_MDP="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --cleanup) KEEP=0; shift ;;
        --pdbfixer) USE_PDBFIXER=1; shift ;;
        --ph) PH="$2"; shift 2 ;;
        --skip-solvate-minimize) RUN_SOLVATE_MINIMIZE=0; shift ;;
        --charge-method) CHARGE_METHOD="$(echo "$2" | tr '[:upper:]' '[:lower:]')"; shift 2 ;;
        -q|--charge) LIG_CHARGE="$2"; shift 2 ;;
        -m|--mult) LIG_MULT="$2"; shift 2 ;;
        --resp-opt) RESP_OPT="$(echo "$2" | tr '[:upper:]' '[:lower:]')"; shift 2 ;;
        --nproc) NPROC="$2"; shift 2 ;;
        --mem) MEM="$2"; shift 2 ;;
        --basis-opt) BASIS_OPT="$2"; shift 2 ;;
        --basis-esp) BASIS_ESP="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown option: $1" ;;
    esac
done

if [[ "$INPUT_MODE" == "local" ]]; then
    [[ -n "$LOCAL_PROTEIN_PDB" && -n "$LOCAL_LIGAND_PDB" ]] || die "Local mode requires both --local-protein and --local-ligand"
    [[ -f "$LOCAL_PROTEIN_PDB" ]] || die "Local protein PDB not found: $LOCAL_PROTEIN_PDB"
    [[ -f "$LOCAL_LIGAND_PDB"  ]] || die "Local ligand PDB not found: $LOCAL_LIGAND_PDB"

    # Important: absolute paths are required because the script later cd's into OUTDIR.
    LOCAL_PROTEIN_PDB="$(realpath "$LOCAL_PROTEIN_PDB")"
    LOCAL_LIGAND_PDB="$(realpath "$LOCAL_LIGAND_PDB")"

    PROID="PROT"
    if [[ "$LIGID_USER_SET" -eq 0 ]]; then
        LIGID="LIG"
    fi
fi

case "$CHARGE_METHOD" in
    bcc|resp) ;;
    *) die "--charge-method must be bcc or resp" ;;
esac

case "$RESP_OPT" in
    yes|y|no|n) ;;
    *) die "--resp-opt must be yes or no" ;;
esac

[[ "$LIG_CHARGE" =~ ^-?[0-9]+$ ]] || die "--charge must be an integer, got: $LIG_CHARGE"
[[ "$LIG_MULT" =~ ^[0-9]+$ ]] || die "--mult must be a positive integer, got: $LIG_MULT"
[[ "$NPROC" =~ ^[0-9]+$ ]] || die "--nproc must be a positive integer, got: $NPROC"

[[ -n "$OUTDIR" ]] || OUTDIR="prep_${PROID}_${LIGID}_${CHARGE_METHOD}"
[[ -n "$IONS_MDP" ]] || IONS_MDP="${MDPDIR}/ions.mdp"
[[ -n "$EM_MDP"   ]] || EM_MDP="${MDPDIR}/em.mdp"

for cmd in wget awk grep sed gmx obabel antechamber parmchk2 tleap tr sort uniq head tail wc realpath; do
    require_cmd "$cmd"
done
require_cmd amb2gro_top_gro.py

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

LOGFILE="prepare_${PROID}_${LIGID}_${CHARGE_METHOD}.log"
exec > >(tee -a "$LOGFILE") 2>&1

log "Input mode:              $INPUT_MODE"
log "Protein/PDB ID:          $PROID"
log "Ligand name:             $LIGID"
log "Ligand charge method:    $CHARGE_METHOD"
log "Ligand net charge:       $LIG_CHARGE"
log "Ligand multiplicity:     $LIG_MULT"
if [[ "$CHARGE_METHOD" == "resp" ]]; then
    log "RESP Gaussian optimize:  $RESP_OPT"
    log "RESP opt level:          $BASIS_OPT"
    log "RESP ESP level:          $BASIS_ESP"
    log "Gaussian resources:      nproc=$NPROC mem=$MEM"
fi
log "Output directory:        $(pwd)"
log "Solvate/minimize:        $RUN_SOLVATE_MINIMIZE"

RAW_PDB="${PROID}.pdb"
WORK_PDB="input_working.pdb"
PROTEIN_PDB="protein.pdb"
LIGAND_PDB="${LIGID}.pdb"
COMPLEX_PDB="complex.pdb"

PROTEIN_GRO="protein.gro"
PROTEIN_TOP="topol_protein.top"
PROTEIN_ITP="topol_protein.itp"

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
COMPLEX_GRO="complex.gro"
INDEX_NDX="index.ndx"

SEL_CHAIN=""
SEL_RESI=""
LIG_MOLNAME=""
declare -a PROTEIN_MOLNAMES=()

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
            if (n == 0) {
                printf "  Charge check: no MOL2 atoms parsed\n"
                exit 0
            }
            printf "  MOL2 atom count = %d\n", n
            printf "  Sum charge     = %.6f\n", s
        }
    ' "$mol2"
}

download_pdb() {
    if [[ "$INPUT_MODE" == "local" ]]; then
        log "Using local protein and ligand PDB files"
        cp "$LOCAL_PROTEIN_PDB" "$PROTEIN_PDB"
        cp "$LOCAL_LIGAND_PDB" "$LIGAND_PDB"
        cat "$PROTEIN_PDB" "$LIGAND_PDB" > "$COMPLEX_PDB"
        cp "$PROTEIN_PDB" "$WORK_PDB"
    else
        log "Downloading ${PROID} from RCSB"
        wget -O "$RAW_PDB" "https://files.rcsb.org/download/${PROID}.pdb"
        [[ -s "$RAW_PDB" ]] || die "Downloaded PDB is empty"
        cp "$RAW_PDB" "$WORK_PDB"
    fi
}

maybe_run_pdbfixer() {
    (( USE_PDBFIXER == 1 )) || return 0
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

inp = sys.argv[1]
ph = float(sys.argv[2])

fixer = PDBFixer(filename=inp)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(ph)
with open(inp + ".fixed", "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
PY
        [[ -s "$WORK_PDB.fixed" ]] || die "Python pdbfixer did not produce output"
        mv "$WORK_PDB.fixed" "$WORK_PDB"
        return 0
    fi

    die "PDBFixer requested but neither CLI nor Python package is available"
}

detect_ligand_sites() {
    if [[ "$INPUT_MODE" == "local" ]]; then
        log "Local mode: ligand selection by chain/residue is skipped"
        SEL_CHAIN=""
        SEL_RESI=""
        return 0
    fi

    awk -v lig="$LIGID" '
        /^HETATM/ {
            resn  = substr($0,18,3)
            alt   = substr($0,17,1)
            chain = substr($0,22,1)
            resi  = substr($0,23,4)
            gsub(/ /,"",resn)
            gsub(/ /,"",alt)
            gsub(/ /,"",chain)
            gsub(/ /,"",resi)
            if (resn == lig && (alt == "" || alt == "A")) print chain " " resi
        }
    ' "$WORK_PDB" | sort -u > ligand_sites.tsv

    [[ -s ligand_sites.tsv ]] || die "No ligand ${LIGID} found in ${WORK_PDB}"
    log "Detected ligand site(s):"
    cat ligand_sites.tsv

    if [[ -n "$LIGCHAIN" || -n "$LIGRESI" ]]; then
        [[ -n "$LIGCHAIN" && -n "$LIGRESI" ]] || die "Use both --lig-chain and --lig-resi"
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

extract_coordinates() {
    if [[ "$INPUT_MODE" == "local" ]]; then
        log "Local mode: using provided protein and ligand files directly"
        [[ -s "$PROTEIN_PDB" ]] || die "Local protein file copy failed"
        [[ -s "$LIGAND_PDB"  ]] || die "Local ligand file copy failed"
        cat "$PROTEIN_PDB" "$LIGAND_PDB" > "$COMPLEX_PDB"
        return 0
    fi

    log "Extracting protein ATOM records"
    grep '^ATOM  ' "$WORK_PDB" > "$PROTEIN_PDB"
    [[ -s "$PROTEIN_PDB" ]] || die "protein.pdb is empty"

    log "Extracting selected ligand"
    awk -v lig="$LIGID" -v schain="$SEL_CHAIN" -v sresi="$SEL_RESI" '
        /^HETATM/ {
            resn  = substr($0,18,3)
            alt   = substr($0,17,1)
            chain = substr($0,22,1)
            resi  = substr($0,23,4)
            gsub(/ /,"",resn)
            gsub(/ /,"",alt)
            gsub(/ /,"",chain)
            gsub(/ /,"",resi)
            if (resn == lig && chain == schain && resi == sresi && (alt == "" || alt == "A")) print
        }
    ' "$WORK_PDB" > "$LIGAND_PDB"
    [[ -s "$LIGAND_PDB" ]] || die "Ligand extraction failed"

    cat "$PROTEIN_PDB" "$LIGAND_PDB" > "$COMPLEX_PDB"
}

prepare_protein() {
    log "Running pdb2gmx"
    gmx pdb2gmx -f "$PROTEIN_PDB" -o "$PROTEIN_GRO" -p "$PROTEIN_TOP" -i posre_protein.itp -ff "$FF" -water "$WATER" -ignh -merge all
    [[ -s "$PROTEIN_TOP" ]] || die "Protein topology not created"
    [[ -s "$PROTEIN_GRO" ]] || die "Protein gro not created"

    log "Extracting protein include block"
    sed -n '/moleculetype/,/endif/{/endif/!p;/endif/p;}' "$PROTEIN_TOP" > "$PROTEIN_ITP"
    [[ -s "$PROTEIN_ITP" ]] || die "Protein include file not created"

    log "Detecting protein molecule names from [ molecules ] in ${PROTEIN_TOP}"
    mapfile -t PROTEIN_MOLNAMES < <(
        awk '
            BEGIN { in_mol = 0 }
            /^\[ *molecules *\]/ { in_mol = 1; next }
            !in_mol { next }
            /^[[:space:]]*$/ { next }
            /^[[:space:]]*;/ { next }
            /^[[:space:]]*\[/ { next }
            { print $1 }
        ' "$PROTEIN_TOP"
    )
    (( ${#PROTEIN_MOLNAMES[@]} > 0 )) || die "Could not detect protein molecule name(s)"
    log "Protein molecule names detected: ${PROTEIN_MOLNAMES[*]}"
}

prepare_ligand_bcc() {
    log "Ligand charge route: AM1-BCC"
    log "Generating ligand mol2 with hydrogens"
    obabel "$LIGAND_PDB" -O "$LIG_H_MOL2" -h
    [[ -s "$LIG_H_MOL2" ]] || die "Open Babel failed to create $LIG_H_MOL2"

    log "Running antechamber AM1-BCC"
    antechamber \
        -i "$LIG_H_MOL2" -fi mol2 \
        -o "$AC_MOL2" -fo mol2 \
        -c bcc -s 2 \
        -rn "$LIGID" -at gaff2 \
        -nc "$LIG_CHARGE" -pf yes
    [[ -s "$AC_MOL2" ]] || die "Antechamber did not create $AC_MOL2"
    check_charge_sum_mol2 "$AC_MOL2"
}

prepare_ligand_resp() {
    log "Ligand charge route: RESP through Gaussian16 + antechamber"

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
#p opt ${BASIS_OPT}

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
#p ${BASIS_ESP} guess=read geom=check Pop=(MK) IOp(6/50=1)

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
#p ${BASIS_ESP} Pop=(MK) IOp(6/50=1)

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
    check_charge_sum_mol2 "$AC_MOL2"
}

prepare_ligand() {
    case "$CHARGE_METHOD" in
        bcc)  prepare_ligand_bcc ;;
        resp) prepare_ligand_resp ;;
        *) die "Unsupported charge method: $CHARGE_METHOD" ;;
    esac

    log "Running parmchk2"
    parmchk2 -i "$AC_MOL2" -f mol2 -o "$AC_FRCMOD"
    [[ -s "$AC_FRCMOD" ]] || die "parmchk2 did not create $AC_FRCMOD"

    log "Running tleap"
    tleap -f - <<EOF_TLEAP
source leaprc.gaff2
lig = loadmol2 ${AC_MOL2}
loadamberparams ${AC_FRCMOD}
saveamberparm lig ${AC_PRMTOP} ${AC_INPCRD}
quit
EOF_TLEAP

    [[ -s "$AC_PRMTOP" ]] || die "Ligand prmtop missing"
    [[ -s "$AC_INPCRD" ]] || die "Ligand inpcrd missing"

    log "Converting ligand Amber files to GROMACS"
    amb2gro_top_gro.py -p "$AC_PRMTOP" -c "$AC_INPCRD" -t "$LIG_GMX_TOP" -g "$LIG_GMX_GRO" -b "$LIG_GMX_PDB"

    [[ -s "$LIG_GMX_TOP" ]] || die "Ligand GROMACS top missing"
    [[ -s "$LIG_GMX_GRO" ]] || die "Ligand GROMACS gro missing"

    log "Splitting ligand topology"
    sed -n '/atomtypes/,/moleculetype/{/moleculetype/b;p}' "$LIG_GMX_TOP" > "$LIG_ATOMTYPES_ITP"
    sed -n '/moleculetype/,/system/{/system/b;p}' "$LIG_GMX_TOP" > "$LIG_MOL_ITP"
    [[ -s "$LIG_ATOMTYPES_ITP" ]] || die "Ligand atomtypes include is empty"
    [[ -s "$LIG_MOL_ITP" ]] || die "Ligand molecule include is empty"

    cat >> "$LIG_MOL_ITP" <<EOF_POSRE

; Include Position restraint file
#ifdef POSRES
#include "posre_${LIGID}.itp"
#endif
EOF_POSRE

    log "Generating ligand position restraints"
    echo 0 | gmx genrestr -f "$LIG_GMX_GRO" -o "$LIG_POSRE_ITP" -fc 1000 1000 1000
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

validate_topology_names() {
    log "Validating topology naming consistency"
    [[ "$LIG_MOLNAME" == "$LIGID" ]] || log "Warning: ligand label (${LIGID}) differs from moleculetype (${LIG_MOLNAME}); [ molecules ] will use ${LIG_MOLNAME}"
    grep -q '^\[ *moleculetype *\]' "$PROTEIN_ITP" || die "Protein include file lacks [ moleculetype ]"
    grep -q '^\[ *moleculetype *\]' "$LIG_MOL_ITP" || die "Ligand include file lacks [ moleculetype ]"
    [[ -s "$LIG_ATOMTYPES_ITP" ]] || die "Ligand atomtypes include is empty"
}

build_final_topology() {
    log "Building final topol.top from scratch"
    {
        echo "; Include forcefield parameters"
        echo "#include \"${FF}.ff/forcefield.itp\""
        echo "#include \"./${LIG_ATOMTYPES_ITP}\""
        echo
        echo "#include \"./${PROTEIN_ITP}\""
        echo "#include \"./${LIG_MOL_ITP}\""
        echo
        echo "; Include water topology"
        case "$WATER" in
            tip3p) echo "#include \"${FF}.ff/tip3p.itp\"" ;;
            tip4p) echo "#include \"${FF}.ff/tip4p.itp\"" ;;
            tip4pew) echo "#include \"${FF}.ff/tip4pew.itp\"" ;;
            spc) echo "#include \"${FF}.ff/spc.itp\"" ;;
            spce) echo "#include \"${FF}.ff/spce.itp\"" ;;
            *) log "Warning: unknown water model '${WATER}', defaulting to tip3p.itp"; echo "#include \"${FF}.ff/tip3p.itp\"" ;;
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
        echo "Protein-ligand in water"
        echo
        echo "[ molecules ]"
        echo "; Compound        #mols"
        for mol in "${PROTEIN_MOLNAMES[@]}"; do
            printf "%-18s %5d\n" "$mol" 1
        done
        printf "%-18s %5d\n" "$LIG_MOLNAME" 1
    } > "$FINAL_TOP"
    [[ -s "$FINAL_TOP" ]] || die "Failed to build ${FINAL_TOP}"
}

merge_gro_files() {
    log "Merging protein and ligand GRO files"
    local pro_atm lig_atm total
    pro_atm=$(sed -n '2p' "$PROTEIN_GRO" | tr -d '[:space:]')
    lig_atm=$(sed -n '2p' "$LIG_GMX_GRO" | tr -d '[:space:]')
    [[ "$pro_atm" =~ ^[0-9]+$ ]] || die "Invalid protein atom count"
    [[ "$lig_atm" =~ ^[0-9]+$ ]] || die "Invalid ligand atom count"
    total=$((pro_atm + lig_atm))

    awk -v natoms="$total" '
        FNR==1 { if (NR==1) title=$0; next }
        FNR==2 { next }
        { lines[FILENAME, ++n[FILENAME]] = $0 }
        END {
            print title
            print natoms
            for (i=1; i<n[ARGV[1]]; i++) print lines[ARGV[1], i]
            for (i=1; i<n[ARGV[2]]; i++) print lines[ARGV[2], i]
            print lines[ARGV[1], n[ARGV[1]]]
        }
    ' "$PROTEIN_GRO" "$LIG_GMX_GRO" > "$COMPLEX_GRO"
    [[ -s "$COMPLEX_GRO" ]] || die "complex.gro not created"
}

generate_index() {
    log "Generating initial index file from unsolvated complex"
    gmx make_ndx -f "$COMPLEX_GRO" -o "$INDEX_NDX" <<EOF_NDX
q
EOF_NDX
    [[ -s "$INDEX_NDX" ]] || die "index.ndx was not created"
}

generate_tc_index_after_ions() {
    local coord="${1:-solv_ions.gro}"
    local tmp_ndx="default_after_ions.tmp.ndx"
    local tc_group="Protein_${LIGID}"
    local group_count protein_group ligand_group new_group

    [[ -s "$coord" ]] || die "Cannot generate production index; coordinate file missing: $coord"

    log "Generating default post-ionization index from ${coord}"
    gmx make_ndx -f "$coord" -o "$tmp_ndx" <<EOF_NDX_DEFAULT
q
EOF_NDX_DEFAULT
    [[ -s "$tmp_ndx" ]] || die "Temporary post-ionization index was not created"

    group_count="$(
        awk '/^\[/{n++} END{print n+0}' "$tmp_ndx"
    )"
    [[ "$group_count" =~ ^[0-9]+$ ]] || die "Could not count groups in $tmp_ndx"

    protein_group="$(
        awk '
            BEGIN { i=-1 }
            /^\[/ {
                i++
                name=$0
                gsub(/^\[ */, "", name)
                gsub(/ *\]$/, "", name)
                if (name == "Protein") {
                    print i
                    exit
                }
            }
        ' "$tmp_ndx"
    )"
    [[ -n "$protein_group" ]] || die "Could not find default [ Protein ] group in $tmp_ndx"

    ligand_group="$(
        awk -v lig="$LIGID" -v mol="$LIG_MOLNAME" '
            BEGIN { i=-1 }
            /^\[/ {
                i++
                name=$0
                gsub(/^\[ */, "", name)
                gsub(/ *\]$/, "", name)
                if (name == lig || name == mol) {
                    print i
                    exit
                }
            }
        ' "$tmp_ndx"
    )"

    if [[ -z "$ligand_group" ]]; then
        log "Warning: could not find ligand group named ${LIGID} or ${LIG_MOLNAME} in default index"
        log "Warning: falling back to make_ndx group 13, matching the manual 1|13 workflow"
        ligand_group=13
    fi

    [[ "$protein_group" =~ ^[0-9]+$ ]] || die "Invalid Protein group index: $protein_group"
    [[ "$ligand_group" =~ ^[0-9]+$ ]] || die "Invalid ligand group index: $ligand_group"

    if (( protein_group >= group_count )); then
        die "Protein group index ${protein_group} is outside default index group count ${group_count}"
    fi
    if (( ligand_group >= group_count )); then
        die "Ligand group index ${ligand_group} is outside default index group count ${group_count}"
    fi

    # make_ndx group numbering is zero-based. After one OR operation, the newly
    # created group number is the previous number of groups.
    new_group="$group_count"

    log "Creating production tc-grps index group: ${tc_group}"
    log "make_ndx operation: ${protein_group}|${ligand_group}; rename group ${new_group} to ${tc_group}"
    gmx make_ndx -f "$coord" -o "$INDEX_NDX" <<EOF_NDX_TC
${protein_group}|${ligand_group}
name ${new_group} ${tc_group}
q
EOF_NDX_TC
    [[ -s "$INDEX_NDX" ]] || die "Production index.ndx was not created"

    awk -v grp="$tc_group" '
        BEGIN { found=0 }
        /^\[/ {
            name=$0
            gsub(/^\[ */, "", name)
            gsub(/ *\]$/, "", name)
            if (name == grp) found=1
        }
        END { exit(found ? 0 : 1) }
    ' "$INDEX_NDX" || die "Production index lacks expected group [ ${tc_group} ]"

    if ! awk '
        BEGIN { found=0 }
        /^\[/ {
            name=$0
            gsub(/^\[ */, "", name)
            gsub(/ *\]$/, "", name)
            if (name == "Water_and_ions") found=1
        }
        END { exit(found ? 0 : 1) }
    ' "$INDEX_NDX"; then
        log "Warning: [ Water_and_ions ] was not found in index.ndx; check tc-grps in your NVT/NPT/MD MDP files"
    fi

    rm -f "$tmp_ndx"
    log "Production index ready: ${INDEX_NDX}"
    log "Use in MDP: tc-grps = ${tc_group} Water_and_ions"
}

solvate_and_minimize() {
    if [[ "$RUN_SOLVATE_MINIMIZE" -eq 0 ]]; then
        log "Skipping solvation, ionization, and minimization by request"
        return 0
    fi

    log "Creating box"
    gmx editconf -f "$COMPLEX_GRO" -o newbox.gro -bt "$BOXTYPE" -d "$BOXDIST"

    log "Solvating"
    gmx solvate -cp newbox.gro -cs spc216.gro -p "$FINAL_TOP" -o solv.gro

    log "Preparing ions"
    gmx grompp -f "$IONS_MDP" -c solv.gro -r solv.gro -p "$FINAL_TOP" -o ions.tpr -maxwarn 1

    log "Adding ions"
    echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p "$FINAL_TOP" -pname NA -nname CL -neutral

    generate_tc_index_after_ions "solv_ions.gro"

    log "Preparing minimization"
    gmx grompp -f "$EM_MDP" -c solv_ions.gro -r solv_ions.gro -p "$FINAL_TOP" -o em.tpr

    log "Running minimization"
    gmx mdrun -v -deffnm em
}

cleanup_files() {
    [[ "$KEEP" -eq 1 ]] && return 0
    log "Cleaning intermediates"
    rm -f \
        "$LIG_H_MOL2" "$LIG_XYZ" "$LIG_COORD" \
        "$AC_MOL2" "$AC_FRCMOD" "$AC_PRMTOP" "$AC_INPCRD" \
        "$LIG_GMX_TOP" "$LIG_GMX_PDB" \
        ligand_sites.tsv "$WORK_PDB"
}

write_summary() {
    log "Writing preparation summary"
    cat > preparation_summary.txt <<EOF_SUMMARY
Protein/PDB ID:        ${PROID}
Ligand name:           ${LIGID}
Input mode:            ${INPUT_MODE}
Force field:           ${FF}
Water model:           ${WATER}
Charge method:         ${CHARGE_METHOD}
Ligand charge:         ${LIG_CHARGE}
Ligand multiplicity:   ${LIG_MULT}
RESP optimization:     ${RESP_OPT}
Final topology:        ${FINAL_TOP}
Complex coordinates:   ${COMPLEX_GRO}
Index file:            ${INDEX_NDX}
Production tc-grps:    Protein_${LIGID} Water_and_ions
Log file:              ${LOGFILE}
EOF_SUMMARY
}

download_pdb
maybe_run_pdbfixer
detect_ligand_sites
extract_coordinates
prepare_protein
prepare_ligand
validate_topology_names
build_final_topology
merge_gro_files
generate_index
solvate_and_minimize
cleanup_files
write_summary

log "All steps completed successfully"
