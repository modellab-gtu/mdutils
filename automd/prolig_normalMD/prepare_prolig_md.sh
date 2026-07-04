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
LOCAL_LIGAND_GIVEN=0   # set when --local-ligand is given (independent of protein source)
INPUT_MODE="rcsb"      # rcsb = download protein; local = use --local-protein
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
USE_PDBFIXER=0
PH=7.0
RUN_SOLVATE_MINIMIZE=1

# Topology pathway
USE_PARMED=0               # 0 = pdb2gmx + amb2gro_top_gro.py  |  1 = full tleap + parmed

# Ligand charge/RESP controls
CHARGE_METHOD="bcc"        # resp (Gaussian16) or any antechamber -c value: bcc, abcg2, mul, cm2, …
LIG_CHARGE=0
LIG_MULT=1
RESP_OPT="yes"             # yes or no
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
  1) RCSB protein + RCSB ligand (default):
    $(basename "$0") --pdb 3HTB --ligand JZ4 --charge-method bcc --charge 0

  2) RCSB protein + local ligand:
    $(basename "$0") --pdb 3HTB --local-ligand lig.pdb --ligand JZ4 \\
      --charge-method resp --charge 0

  3) Local protein + local ligand:
    $(basename "$0") --local-protein prot_noH.pdb --local-ligand lig.pdb \\
      --ligand LIG --charge-method resp --charge 0 --resp-opt yes

General options:
  -p, --pdb ID                 PDB ID from RCSB, default: ${PROID}
  -l, --ligand RES             Ligand residue/molecule name, default: ${LIGID}
      --local-protein FILE     Local protein PDB file
      --local-ligand FILE      Local ligand PDB file
      --lig-chain CHAIN        Ligand chain ID for RCSB mode
      --lig-resi RESI          Ligand residue number for RCSB mode
  -f, --ff NAME                Protein force field (GROMACS name used by both pathways), default: ${FF}
  -w, --water NAME             Water model, default: ${WATER}
  -b, --boxtype TYPE           Box type, default: ${BOXTYPE}
  -d, --distance NM            Solvent distance in nm, default: ${BOXDIST}
      --ions-mdp FILE          ions.mdp path (default: MDP/ions.mdp if present, else bundled template)
      --em-mdp FILE            em.mdp path  (default: MDP/em.mdp  if present, else bundled template)
      --outdir DIR             Output directory, default: prep_<PDB>_<LIG>
      --verbose                Keep all intermediate files (default: remove intermediates)
      --pdbfixer               Run PDBFixer preprocessing if available/requested
      --ph VALUE               pH for PDBFixer hydrogen addition, default: ${PH}
      --skip-solvate-minimize  Stop after complex.gro, topol.top, and index.ndx
      --use-parmed             Use tleap+parmed pathway instead of pdb2gmx+amb2gro
  -h, --help                   Show this help

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
        --local-ligand)  LOCAL_LIGAND_PDB="$2"; LOCAL_LIGAND_GIVEN=1; shift 2 ;;
        --lig-chain) LIGCHAIN="$2"; shift 2 ;;
        --lig-resi) LIGRESI="$2"; shift 2 ;;
        -f|--ff) FF="$2"; shift 2 ;;
        -w|--water) WATER="$2"; shift 2 ;;
        -b|--boxtype) BOXTYPE="$2"; shift 2 ;;
        -d|--distance) BOXDIST="$2"; shift 2 ;;
        --ions-mdp) IONS_MDP="$2"; shift 2 ;;
        --em-mdp) EM_MDP="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --verbose) VERBOSE=1; shift ;;
        --pdbfixer) USE_PDBFIXER=1; shift ;;
        --ph) PH="$2"; shift 2 ;;
        --skip-solvate-minimize) RUN_SOLVATE_MINIMIZE=0; shift ;;
        --use-parmed) USE_PARMED=1; shift ;;
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

# ── Validate and resolve input sources ───────────────────────────────────────
# Three supported scenarios:
#   1. RCSB protein + RCSB ligand  : (default, --pdb --ligand)
#   2. RCSB protein + local ligand : --pdb + --local-ligand --ligand
#   3. Local protein + local ligand: --local-protein + --local-ligand --ligand

if [[ "$INPUT_MODE" == "local" ]]; then
    [[ -n "$LOCAL_PROTEIN_PDB" ]] || die "Local protein mode requires --local-protein"
    [[ -f "$LOCAL_PROTEIN_PDB" ]] || die "Local protein file not found: $LOCAL_PROTEIN_PDB"
    LOCAL_PROTEIN_PDB="$(realpath "$LOCAL_PROTEIN_PDB")"
    [[ "$LIGID_USER_SET" -eq 0 ]] && PROID="PROT"
fi

if (( LOCAL_LIGAND_GIVEN )); then
    [[ -n "$LOCAL_LIGAND_PDB" ]] || die "--local-ligand requires a file path"
    [[ -f "$LOCAL_LIGAND_PDB"  ]] || die "Local ligand file not found: $LOCAL_LIGAND_PDB"
    LOCAL_LIGAND_PDB="$(realpath "$LOCAL_LIGAND_PDB")"
    [[ "$LIGID_USER_SET" -eq 1 ]] || die "--local-ligand requires --ligand RESNAME to set the molecule ID"
fi

if [[ "$INPUT_MODE" == "local" ]] && ! (( LOCAL_LIGAND_GIVEN )); then
    die "Local protein mode also requires --local-ligand (and --ligand RESNAME)"
fi

case "$RESP_OPT" in
    yes|y|no|n) ;;
    *) die "--resp-opt must be yes or no" ;;
esac

[[ "$LIG_CHARGE" =~ ^-?[0-9]+$ ]] || die "--charge must be an integer, got: $LIG_CHARGE"
[[ "$LIG_MULT" =~ ^[0-9]+$ ]] || die "--mult must be a positive integer, got: $LIG_MULT"
[[ "$NPROC" =~ ^[0-9]+$ ]] || die "--nproc must be a positive integer, got: $NPROC"

[[ -n "$OUTDIR" ]] || OUTDIR="prep_${PROID}_${LIGID}_${CHARGE_METHOD}"
if [[ -z "$IONS_MDP" ]]; then
    if [[ -f "${MDPDIR}/ions.mdp" ]]; then
        IONS_MDP="${MDPDIR}/ions.mdp"
    else
        IONS_MDP="${MDPDIR_BUNDLED}/ions.mdp"
    fi
fi
if [[ -z "$EM_MDP" ]]; then
    if [[ -f "${MDPDIR}/em.mdp" ]]; then
        EM_MDP="${MDPDIR}/em.mdp"
    else
        EM_MDP="${MDPDIR_BUNDLED}/em.mdp"
    fi
fi

for cmd in awk grep sed gmx obabel antechamber parmchk2 tleap tr sort uniq head tail wc realpath; do
    require_cmd "$cmd"
done
[[ "$INPUT_MODE" == "rcsb" ]] && require_cmd wget

# If the chosen FF has no GROMACS force field directory, automatically use the parmed pathway.
_gmx_prefix=$(gmx --version 2>/dev/null | awk '/Data prefix:/{print $NF}')
if [[ ! -d "${_gmx_prefix}/share/gromacs/top/${FF}.ff" ]]; then
    if [[ "$USE_PARMED" -eq 1 ]]; then
        : # user already requested parmed, FF is fine for tleap
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

LOGFILE="prepare_${PROID}_${LIGID}_${CHARGE_METHOD}.log"
exec > >(tee -a "$LOGFILE") 2>&1

if [[ "$INPUT_MODE" == "local" ]] && (( LOCAL_LIGAND_GIVEN )); then
    log "Input mode:              local protein + local ligand"
elif [[ "$INPUT_MODE" == "rcsb" ]] && (( LOCAL_LIGAND_GIVEN )); then
    log "Input mode:              RCSB protein + local ligand"
else
    log "Input mode:              RCSB protein + RCSB ligand"
fi
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
log "Topology pathway:        $( [[ "$USE_PARMED" -eq 1 ]] && echo "tleap+parmed" || echo "pdb2gmx+amb2gro" )"
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

copy_local_ligand() {
    # Copy local ligand to LIGAND_PDB, converting to PDB format via obabel if needed.
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

download_pdb() {
    # Protein source
    if [[ "$INPUT_MODE" == "local" ]]; then
        log "Using local protein: ${LOCAL_PROTEIN_PDB}"
        cp "$LOCAL_PROTEIN_PDB" "$PROTEIN_PDB"
        cp "$LOCAL_PROTEIN_PDB" "$WORK_PDB"
    else
        log "Downloading ${PROID} from RCSB"
        wget -O "$RAW_PDB" "https://files.rcsb.org/download/${PROID}.pdb"
        [[ -s "$RAW_PDB" ]] || die "Downloaded PDB is empty"
        cp "$RAW_PDB" "$WORK_PDB"
    fi

    # Ligand source (independent of protein source)
    if (( LOCAL_LIGAND_GIVEN )); then
        copy_local_ligand
    fi
    # RCSB ligand is extracted later in extract_coordinates after pdbfixer runs
}

maybe_run_pdbfixer() {
    (( USE_PDBFIXER == 1 )) || return 0
    log "PDBFixer requested"

    if optional_cmd pdbfixer; then
        log "Running pdbfixer CLI (pH=${PH})"
        pdbfixer "$WORK_PDB" --output="$WORK_PDB.fixed" \
            --add-atoms=heavy --keep-heterogens=all --ph="$PH"
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
    if (( LOCAL_LIGAND_GIVEN )); then
        log "Local ligand provided: skipping RCSB site detection"
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
    # ── Protein ───────────────────────────────────────────────────────────────
    if [[ "$INPUT_MODE" == "local" ]]; then
        if (( USE_PDBFIXER )); then
            log "Propagating pdbfixer output to protein.pdb"
            cp "$WORK_PDB" "$PROTEIN_PDB"
        fi
        [[ -s "$PROTEIN_PDB" ]] || die "protein.pdb missing after local copy"
    else
        log "Extracting protein ATOM records from ${WORK_PDB}"
        grep '^ATOM  ' "$WORK_PDB" > "$PROTEIN_PDB"
        [[ -s "$PROTEIN_PDB" ]] || die "protein.pdb is empty after extraction"
    fi

    # ── Ligand ────────────────────────────────────────────────────────────────
    if ! (( LOCAL_LIGAND_GIVEN )); then
        log "Extracting selected ligand (${LIGID}) from ${WORK_PDB}"
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
    else
        [[ -s "$LIGAND_PDB" ]] || die "Local ligand file missing: ${LIGAND_PDB}"
    fi

    cat "$PROTEIN_PDB" "$LIGAND_PDB" > "$COMPLEX_PDB"
    [[ -s "$COMPLEX_PDB" ]] || die "complex.pdb is empty"
}

prepare_protein() {
    log "Running pdb2gmx"
    gmx pdb2gmx -f "$PROTEIN_PDB" -o "$PROTEIN_GRO" -p "$PROTEIN_TOP" -i posre_protein.itp -ff "$FF" -water "$WATER" -ignh -merge all
    [[ -s "$PROTEIN_TOP" ]] || die "Protein topology not created"
    [[ -s "$PROTEIN_GRO" ]] || die "Protein gro not created"

    log "Extracting protein include block"
    sed -n '/^\[ *moleculetype/,/^\[ *system/{/^\[ *system/b;p;}' "$PROTEIN_TOP" > "$PROTEIN_ITP"
    [[ -s "$PROTEIN_ITP" ]] || die "Protein include file not created"
    grep -q '^\[ *moleculetype' "$PROTEIN_ITP" || die "Protein ITP extraction failed: missing [ moleculetype ]"

    # GROMACS 2024+ places water/ion #includes before [ system ] so they land in
    # the extracted ITP. Strip them out — topol.top owns those includes.
    python3 - "$PROTEIN_ITP" <<'EOF_STRIP'
import re, sys
with open(sys.argv[1]) as f:
    txt = f.read()
txt = re.sub(r';[^\n]*[Ii]nclude water[^\n]*\n#include[^\n]*\n', '', txt)
txt = re.sub(r'#ifdef POSRES_WATER.*?#endif\n?', '', txt, flags=re.DOTALL)
txt = re.sub(r';[^\n]*[Ii]nclude[^\n]*ions[^\n]*\n#include[^\n]*ions\.itp[^\n]*\n', '', txt)
with open(sys.argv[1], 'w') as f:
    f.write(txt)
EOF_STRIP

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

graft_coords_to_mol2() {
    # Replace coordinates in target mol2 with those from ref mol2, atom-for-atom.
    # Used to restore original ligand geometry after G16 reorientation.
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
        # Rebuild line; mol2 ATOM has 9 fields
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
        log "Parmed pathway: ligand charges ready; full system will be built by tleap+parmed"
        return 0
    fi

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
    [[ -s "$LIG_POSRE_ITP" ]]     || die "Ligand position restraint file missing: ${LIG_POSRE_ITP}"
    [[ -f "posre_protein.itp" ]]  || die "Protein position restraint file missing: posre_protein.itp"
    grep -q '#ifdef POSRES_LIG' "$LIG_MOL_ITP" || die "Ligand mol ITP missing #ifdef POSRES_LIG block"
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

prepare_system_parmed() {
    local amber_ff_leaprc box_cmd box_dist_a
    local complex_prmtop="complex_solv.prmtop"
    local complex_inpcrd="complex_solv.inpcrd"
    local parmed_top="complex_gmx.top"
    local parmed_gro="complex_gmx.gro"
    local tmp_ndx="parmed_tmp.ndx"

    # Map GROMACS FF names to AmberTools leaprc equivalents.
    # Legacy FFs (ff99SBildn, ff99SB, ff96, ff99, ff03) live under oldff/.
    case "${FF%.ff}" in
        amber99sb-ildn) amber_ff_leaprc="oldff/leaprc.ff99SBildn" ;;
        amber99sb)      amber_ff_leaprc="oldff/leaprc.ff99SB" ;;
        amber03)        amber_ff_leaprc="leaprc.protein.ff03.r1" ;;
        amber96)        amber_ff_leaprc="oldff/leaprc.ff96" ;;
        amber99)        amber_ff_leaprc="oldff/leaprc.ff99" ;;
        ff14SB)         amber_ff_leaprc="leaprc.protein.ff14SB" ;;
        ff19SB)         amber_ff_leaprc="leaprc.protein.ff19SB" ;;
        oldff/*|leaprc.*) amber_ff_leaprc="$FF" ;;
        *)              amber_ff_leaprc="leaprc.protein.${FF%.ff}"
            log "Warning: unrecognised FF '${FF}'; trying ${amber_ff_leaprc}" ;;
    esac

    case "$BOXTYPE" in
        dodecahedron|octahedron) box_cmd="solvateOct" ;;
        *) box_cmd="solvateBox" ;;
    esac

    # tleap uses angstroms; BOXDIST is in nm
    box_dist_a=$(python3 -c "print(${BOXDIST} * 10)")

    log "Building full system with tleap (${FF} → ${amber_ff_leaprc} + GAFF2 + TIP3P)"
    tleap -f - <<EOF_TLEAP
source ${amber_ff_leaprc}
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams ${AC_FRCMOD}
LIG = loadMol2 ${AC_MOL2}
PROT = loadpdb protein.pdb
SYSTEM = combine {PROT LIG}
${box_cmd} SYSTEM TIP3PBOX ${box_dist_a}
addions SYSTEM Na+ 0
addions SYSTEM Cl- 0
saveamberparm SYSTEM ${complex_prmtop} ${complex_inpcrd}
quit
EOF_TLEAP
    [[ -s "$complex_prmtop" ]] || die "tleap did not produce ${complex_prmtop}"
    [[ -s "$complex_inpcrd" ]] || die "tleap did not produce ${complex_inpcrd}"

    log "Converting full system AMBER → GROMACS with parmed"
    python3 - <<EOF_PARMED
from parmed import load_file

amber = load_file('${complex_prmtop}', '${complex_inpcrd}')

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

    # Generate position restraints via make_ndx + genrestr
    log "Generating position restraints"
    gmx make_ndx -f "$parmed_gro" -o "$tmp_ndx" <<EOF_NDX
q
EOF_NDX
    [[ -s "$tmp_ndx" ]] || die "make_ndx did not produce ${tmp_ndx}"

    local prot_grp lig_grp
    prot_grp=$(awk 'BEGIN{i=-1}/^\[/{i++;n=$0;gsub(/^\[ */,"",n);gsub(/ *\]$/,"",n);if(n=="Protein"){print i;exit}}' "$tmp_ndx")
    [[ -n "$prot_grp" ]] || die "Could not find Protein group in parmed index"
    echo "$prot_grp" | gmx genrestr -f "$parmed_gro" -n "$tmp_ndx" -o posre_protein.itp -fc 1000 1000 1000
    [[ -s "posre_protein.itp" ]] || die "posre_protein.itp was not created"

    lig_grp=$(awk -v lig="${LIGID}" 'BEGIN{i=-1}/^\[/{i++;n=$0;gsub(/^\[ */,"",n);gsub(/ *\]$/,"",n);if(n==lig){print i;exit}}' "$tmp_ndx")
    if [[ -z "$lig_grp" ]]; then
        log "Warning: ligand group ${LIGID} not found, trying 'Other'"
        lig_grp=$(awk 'BEGIN{i=-1}/^\[/{i++;n=$0;gsub(/^\[ */,"",n);gsub(/ *\]$/,"",n);if(n=="Other"){print i;exit}}' "$tmp_ndx")
    fi
    [[ -n "$lig_grp" ]] || die "Could not find ligand group in parmed index"

    # Use the ligand mol2 (already available from prepare_ligand) as a standalone
    # single-molecule file so genrestr produces local 1-based indices.
    obabel "$AC_MOL2" -O lig_posre_tmp.pdb 2>/dev/null
    [[ -s "lig_posre_tmp.pdb" ]] || die "obabel could not convert ligand mol2 to PDB for genrestr"
    echo "0" | gmx genrestr -f lig_posre_tmp.pdb -o "${LIG_POSRE_ITP}" -fc 1000 1000 1000
    rm -f lig_posre_tmp.pdb
    [[ -s "${LIG_POSRE_ITP}" ]] || die "${LIG_POSRE_ITP} was not created"
    rm -f "$tmp_ndx"

    # Inject #ifdef POSRES / POSRES_LIG / POSRES_WATER blocks into the topology
    log "Injecting POSRES directives into ${parmed_top}"
    python3 - "$parmed_top" "$LIGID" "posre_protein.itp" "$LIG_POSRE_ITP" <<'EOF_INJECT'
import re, sys

top_file, lig_name, posre_prot, posre_lig = sys.argv[1:5]

with open(top_file) as f:
    lines = f.readlines()

result = []
i = 0
mol_idx = 0       # 1-based count of moleculetype sections seen
current_mol = None

def posres_block(kind, include=None):
    if kind == 'protein':
        return ['\n#ifdef POSRES\n', f'#include "./{posre_prot}"\n', '#endif\n\n']
    elif kind == 'lig':
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
        if mol_idx == 1:
            result.extend(posres_block('protein'))
        elif current_mol == lig_name:
            result.extend(posres_block('lig'))
        elif current_mol == 'SOL':
            result.extend(posres_block('water'))
        current_mol = None

    result.append(line)

    if re.match(r'^\[ *moleculetype', line):
        mol_idx += 1
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
print(f'POSRES injected (mol_idx=1 → protein, {lig_name} → POSRES_LIG, SOL → POSRES_WATER)')
EOF_INJECT
    [[ $? -eq 0 ]] || die "POSRES injection into parmed topology failed"

    FINAL_TOP="$parmed_top"
    COMPLEX_GRO="$parmed_gro"

    generate_tc_index_after_ions "$COMPLEX_GRO"
}

solvate_and_minimize() {
    if [[ "$RUN_SOLVATE_MINIMIZE" -eq 0 ]]; then
        log "Skipping solvation, ionization, and minimization by request"
        return 0
    fi

    local em_coord

    if [[ "$USE_PARMED" -eq 0 ]]; then
        log "Creating box"
        gmx editconf -f "$COMPLEX_GRO" -o newbox.gro -bt "$BOXTYPE" -d "$BOXDIST"

        log "Solvating"
        gmx solvate -cp newbox.gro -cs spc216.gro -p "$FINAL_TOP" -o solv.gro

        log "Preparing ions"
        gmx grompp -f "$IONS_MDP" -c solv.gro -r solv.gro -p "$FINAL_TOP" -o ions.tpr -maxwarn 1

        log "Adding ions"
        echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p "$FINAL_TOP" -pname NA -nname CL -neutral

        generate_tc_index_after_ions "solv_ions.gro"
        em_coord="solv_ions.gro"
    else
        # parmed pathway: tleap already solvated and added ions; index already built
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
Minimization - steepest descent then CG, heavy atoms restrained
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
NVT equilibration - Langevin dynamics, position restraints on (1 ns)
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
NPT equilibration - Langevin dynamics + MC barostat, position restraints on (1 ns)
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
    if [[ "$VERBOSE" -eq 1 ]]; then
        return 0
    fi
    log "Removing intermediates (use --verbose to keep all files)"

    # Raw input / working PDBs
    rm -f "${PROID}.pdb" "$WORK_PDB" "$PROTEIN_PDB" "$LIGAND_PDB" "$COMPLEX_PDB" \
          ligand_sites.tsv

    # Charge-generation intermediates
    rm -f "$LIG_H_MOL2" "$LIG_XYZ" "$LIG_COORD" "$AC_MOL2" "$AC_FRCMOD"
    rm -f ANTECHAMBER_*.AC ANTECHAMBER_*.AC0 ATOMTYPE.INF \
          sqm.in sqm.out sqm.pdb leap.log mdinfo mdout.mdp

    if [[ "$USE_PARMED" -eq 0 ]]; then
        # pdb2gmx-specific intermediates
        rm -f "$AC_PRMTOP" "$AC_INPCRD"
        rm -f "$PROTEIN_GRO" "$PROTEIN_TOP" "$LIG_GMX_TOP" "$LIG_GMX_PDB"
        rm -f newbox.gro solv.gro ions.tpr
        if [[ "$RUN_SOLVATE_MINIMIZE" -eq 1 ]]; then
            rm -f solv_ions.gro "$COMPLEX_GRO" em.tpr em.edr em.log em.trr
        fi
    else
        # parmed-specific intermediates
        if [[ "$RUN_SOLVATE_MINIMIZE" -eq 1 ]]; then
            # em.gro is the final structure; pre-EM GROMACS coords no longer needed
            rm -f "$COMPLEX_GRO" em.tpr em.edr em.log em.trr
        fi
    fi

    # GROMACS backup files (#filename.ext.N#)
    rm -f \#*\#
}

write_summary() {
    log "Writing preparation summary"
    {
        cat <<EOF_SUMMARY
Protein/PDB ID:        ${PROID}
Ligand name:           ${LIGID}
Input mode:            $( [[ "$INPUT_MODE" == "local" ]] && echo "local protein" || echo "RCSB protein" ) + $( (( LOCAL_LIGAND_GIVEN )) && echo "local ligand" || echo "RCSB ligand" )
Topology pathway:      $( [[ "$USE_PARMED" -eq 1 ]] && echo "tleap+parmed" || echo "pdb2gmx+amb2gro" )
Force field:           ${FF}
Water model:           ${WATER}
Charge method:         ${CHARGE_METHOD}
Ligand charge:         ${LIG_CHARGE}
Ligand multiplicity:   ${LIG_MULT}
RESP optimization:     ${RESP_OPT}

--- GROMACS output ---
Final topology:        ${FINAL_TOP}
Complex coordinates:   ${COMPLEX_GRO}
Index file:            ${INDEX_NDX}
Production tc-grps:    Protein_${LIGID} Water_and_ions
EOF_SUMMARY
        if [[ "$USE_PARMED" -eq 1 ]]; then
            cat <<EOF_AMBER

--- AMBER output ---
Topology (prmtop):     complex_solv.prmtop
Coordinates (inpcrd):  complex_solv.inpcrd
Minimization input:    amber_min.in
NVT input:             amber_nvt.in
NPT input:             amber_npt.in
Production MD input:   amber_md.in

Example run (pmemd.cuda):
  pmemd.cuda -O -i amber_min.in -p complex_solv.prmtop -c complex_solv.inpcrd -ref complex_solv.inpcrd -o min.out -r min.rst -x min.nc
  pmemd.cuda -O -i amber_nvt.in -p complex_solv.prmtop -c min.rst              -ref complex_solv.inpcrd -o nvt.out -r nvt.rst -x nvt.nc
  pmemd.cuda -O -i amber_npt.in -p complex_solv.prmtop -c nvt.rst              -ref complex_solv.inpcrd -o npt.out -r npt.rst -x npt.nc
  pmemd.cuda -O -i amber_md.in  -p complex_solv.prmtop -c npt.rst                                       -o md.out  -r md.rst  -x md.nc
EOF_AMBER
        fi
        echo
        echo "Log file:              ${LOGFILE}"
    } > preparation_summary.txt
}

download_pdb
maybe_run_pdbfixer
detect_ligand_sites
extract_coordinates

if [[ "$USE_PARMED" -eq 0 ]]; then
    prepare_protein
    prepare_ligand
    validate_topology_names
    build_final_topology
    merge_gro_files
    generate_index
else
    prepare_ligand           # charges only: antechamber/resp + parmchk2
    prepare_system_parmed    # full system: tleap + parmed + POSRES + index
    write_amber_inputs
fi

solvate_and_minimize
cleanup_files
write_summary

log "All steps completed successfully"
