#!/usr/bin/env bash
set -Eeuo pipefail

###############################################################################
# prepare_complex_gmx_labready.sh
#
# A more robust protein-ligand GROMACS preparation workflow.
#
# Features:
#   - Downloads PDB from RCSB
#   - Detects/selects ligand instance by residue name, chain, residue number
#   - Optional PDBFixer preprocessing if available
#   - Detects protein molecule names from generated protein topology
#   - Parameterizes ligand with GAFF2/AmberTools
#   - Converts ligand Amber topology to GROMACS
#   - Validates [ moleculetype ] vs [ molecules ] naming
#   - Builds final topol.top from scratch, including forcefield, water, ions
#   - Generates index.ndx automatically
#   - Solvates, ionizes, minimizes
###############################################################################

PROID="3HTB"
LIGID="JZ4"
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
SCRIPT_START_DIR="$(pwd)"
MDPDIR="${SCRIPT_START_DIR}/MDP"

usage() {
cat <<EOF
Usage:
  $(basename "$0") [options]

Options:
  -p, --pdb ID             PDB ID from RCSB (default: ${PROID})
  -l, --ligand RES         Ligand residue name in RCSB mode, e.g. JZ4
      --local-protein FILE Local protein PDB file
      --local-ligand FILE  Local ligand PDB file
      --lig-chain CHAIN    Ligand chain ID
      --lig-resi RESI      Ligand residue number
  -f, --ff NAME            Force field (default: ${FF})
  -w, --water NAME         Water model (default: ${WATER})
  -b, --boxtype TYPE       Box type (default: ${BOXTYPE})
  -d, --distance NM        Solvent distance in nm (default: ${BOXDIST})
      --ions-mdp FILE      ions.mdp path (default: $PWD/MDP/ions.mdp)
      --em-mdp FILE        em.mdp path (default: $PWD/MDP/em.mdp)
      --outdir DIR         Output directory (default: prep_<PDB>_<LIG>)
      --cleanup            Remove selected intermediates
      --pdbfixer           Run PDBFixer preprocessing if available
      --ph VALUE           pH for PDBFixer hydrogen addition (default: ${PH})
  -h, --help               Show this help
EOF
}

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
die() { echo "Error: $*" >&2; exit 1; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"; }
optional_cmd() { command -v "$1" >/dev/null 2>&1; }
cleanup_on_error() { echo "Workflow failed." >&2; }
trap cleanup_on_error ERR

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--pdb) PROID="$2"; shift 2 ;;
        -l|--ligand) LIGID="$2"; shift 2 ;;
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
        -h|--help) usage; exit 0 ;;
        *) die "Unknown option: $1" ;;
    esac
done

if [[ "$INPUT_MODE" == "local" ]]; then
    [[ -n "$LOCAL_PROTEIN_PDB" && -n "$LOCAL_LIGAND_PDB" ]] || die "Local mode requires both --local-protein and --local-ligand"
    [[ -f "$LOCAL_PROTEIN_PDB" ]] || die "Local protein PDB not found: $LOCAL_PROTEIN_PDB"
    [[ -f "$LOCAL_LIGAND_PDB"  ]] || die "Local ligand PDB not found: $LOCAL_LIGAND_PDB"
    PROID="PROT"
    LIGID="LIG"
fi

[[ -n "$OUTDIR" ]] || OUTDIR="prep_${PROID}_${LIGID}"
[[ -n "$IONS_MDP" ]] || IONS_MDP="${MDPDIR}/ions.mdp"
[[ -n "$EM_MDP"   ]] || EM_MDP="${MDPDIR}/em.mdp"

for cmd in wget awk grep sed gmx obabel antechamber parmchk2 tleap tr sort uniq head tail wc; do
    require_cmd "$cmd"
done
require_cmd amb2gro_top_gro.py

IONS_MDP="$(realpath "$IONS_MDP")"
EM_MDP="$(realpath "$EM_MDP")"

[[ -f "$IONS_MDP" ]] || die "ions.mdp not found: $IONS_MDP"
[[ -f "$EM_MDP"   ]] || die "em.mdp not found: $EM_MDP"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

LOGFILE="prepare_${PROID}_${LIGID}.log"
exec > >(tee -a "$LOGFILE") 2>&1

RAW_PDB="${PROID}.pdb"
WORK_PDB="input_working.pdb"
PROTEIN_PDB="protein.pdb"
LIGAND_PDB="${LIGID}.pdb"
COMPLEX_PDB="complex.pdb"

PROTEIN_GRO="protein.gro"
PROTEIN_TOP="topol_protein.top"
PROTEIN_ITP="topol_protein.itp"

LIG_H_MOL2="${LIGID}_h.mol2"
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
            if (resn == lig && (alt == "" || alt == "A")) print chain "	" resi
        }
    ' "$WORK_PDB" | sort -u > ligand_sites.tsv

    [[ -s ligand_sites.tsv ]] || die "No ligand ${LIGID} found in ${WORK_PDB}"
    log "Detected ligand site(s):"
    cat ligand_sites.tsv

    if [[ -n "$LIGCHAIN" || -n "$LIGRESI" ]]; then
        [[ -n "$LIGCHAIN" && -n "$LIGRESI" ]] || die "Use both --lig-chain and --lig-resi"
        grep -qP "^${LIGCHAIN}	${LIGRESI}$" ligand_sites.tsv || die "Requested ligand chain/resi not found"
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

prepare_ligand() {
    log "Generating ligand mol2 with hydrogens"
    obabel "$LIGAND_PDB" -O "$LIG_H_MOL2" -h
    [[ -s "$LIG_H_MOL2" ]] || die "Open Babel failed"

    log "Running antechamber"
    antechamber -i "$LIG_H_MOL2" -fi mol2 -o "$AC_MOL2" -fo mol2 -c bcc -s 2

    log "Running parmchk2"
    parmchk2 -i "$AC_MOL2" -f mol2 -o "$AC_FRCMOD"

    log "Running tleap"
    tleap -f - <<EOF
source leaprc.gaff2
lig = loadmol2 ${AC_MOL2}
loadamberparams ${AC_FRCMOD}
saveamberparm lig ${AC_PRMTOP} ${AC_INPCRD}
quit
EOF

    [[ -s "$AC_PRMTOP" ]] || die "Ligand prmtop missing"
    [[ -s "$AC_INPCRD" ]] || die "Ligand inpcrd missing"

    log "Converting ligand Amber files to GROMACS"
    amb2gro_top_gro.py -p "$AC_PRMTOP" -c "$AC_INPCRD" -t "$LIG_GMX_TOP" -g "$LIG_GMX_GRO" -b "$LIG_GMX_PDB"

    [[ -s "$LIG_GMX_TOP" ]] || die "Ligand GROMACS top missing"
    [[ -s "$LIG_GMX_GRO" ]] || die "Ligand GROMACS gro missing"

    log "Splitting ligand topology"
    sed -n '/atomtypes/,/moleculetype/{/moleculetype/b;p}' "$LIG_GMX_TOP" > "$LIG_ATOMTYPES_ITP"
    sed -n '/moleculetype/,/system/{/system/b;p}' "$LIG_GMX_TOP" > "$LIG_MOL_ITP"

    cat >> "$LIG_MOL_ITP" <<EOF

; Include Position restraint file
#ifdef POSRES
#include "posre_${LIGID}.itp"
#endif
EOF

    log "Generating ligand position restraints"
    echo 0 | gmx genrestr -f "$LIG_GMX_GRO" -o "$LIG_POSRE_ITP" -fc 1000 1000 1000

    log "Detecting ligand molecule name from [ moleculetype ]"
    LIG_MOLNAME="$(
        awk '
            /^\[ *moleculetype *\]/ {
                if (getline <= 0) exit
                if (getline <= 0) exit
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
    log "Generating index file"
    gmx make_ndx -f "$COMPLEX_GRO" -o "$INDEX_NDX" <<EOF
q
EOF
    [[ -s "$INDEX_NDX" ]] || die "index.ndx was not created"
}

solvate_and_minimize() {
    log "Creating box"
    gmx editconf -f "$COMPLEX_GRO" -o newbox.gro -bt "$BOXTYPE" -d "$BOXDIST"

    log "Solvating"
    gmx solvate -cp newbox.gro -cs spc216.gro -p "$FINAL_TOP" -o solv.gro

    log "Preparing ions"
    gmx grompp -f "$IONS_MDP" -c solv.gro -p "$FINAL_TOP" -o ions.tpr

    log "Adding ions"
    echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p "$FINAL_TOP" -pname NA -nname CL -neutral

    log "Preparing minimization"
    gmx grompp -f "$EM_MDP" -c solv_ions.gro -p "$FINAL_TOP" -o em.tpr

    log "Running minimization"
    gmx mdrun -v -deffnm em
}

cleanup_files() {
    [[ "$KEEP" -eq 1 ]] && return 0
    log "Cleaning intermediates"
    rm -f "$LIG_H_MOL2" "$AC_MOL2" "$AC_FRCMOD" "$AC_PRMTOP" "$AC_INPCRD" "$LIG_GMX_TOP" "$LIG_GMX_PDB" ligand_sites.tsv "$WORK_PDB"
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

log "All steps completed successfully"
