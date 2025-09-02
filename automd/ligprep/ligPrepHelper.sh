from ase.io import read
from ligPrep import ligPrep
import argparse
import os, sys, shutil
import multiprocessing

nprocs_all = int(multiprocessing.cpu_count())



parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("structure_dir", type=str)
parser.add_argument("add_hydrogen", nargs="?", default="No") # args for bool
parser.add_argument("calculator_type", type=str)
parser.add_argument("optimization_conf", nargs="?", default="No") # args for bool
parser.add_argument("optimization_lig", nargs="?", default="No") # args for bool
parser.add_argument("pre_optimization_lig", nargs="?", default="No") # args for bool
parser.add_argument("genconformer", nargs="?", default="No") # args for bool
parser.add_argument("nprocs", type=int, default=nprocs_all)
parser.add_argument("thr_fmax", type=float, default=0.05)
parser.add_argument("maxiter", type=float, default=500)

parser.add_argument("num_conformers", type=int, default=50)
parser.add_argument("max_attempts", type=int, default=100)
parser.add_argument("prune_rms_thresh", type=float, default=0.2)

args = parser.parse_args()
structure_dir = args.structure_dir
calculator_type = args.calculator_type

def getBoolStr(string):
    string = string.lower()
    if "true" in string or "yes" in string:
        return True
    elif "false" in string or "no" in string:
        return False
    else:
        print("%s is bad input!!! Must be Yes/No or True/False" %string)
        sys.exit(1)


optimization_conf = getBoolStr(args.optimization_conf)
optimization_lig = getBoolStr(args.optimization_lig)
pre_optimization_lig = getBoolStr(args.pre_optimization_lig)
genconformer = getBoolStr(args.genconformer)
add_hydrogen = getBoolStr(args.add_hydrogen)

nprocs = args.nprocs
thr_fmax = args.thr_fmax
maxiter = args.maxiter


#get conformer generator parameters
num_conformers = args.num_conformers
max_attempts = args.max_attempts
prune_rms_thresh = args.prune_rms_thresh

def setG16calculator(lig, file_base, label, WORK_DIR):
    lig.setG16Calculator(
            label="%s/g16_%s/%s"%(WORK_DIR, label, file_base),
            chk="%s.chk"%file_base,
            nprocs=nprocs,
            xc="HF",
            basis="6-31g*",
            scf="maxcycle=100",
            extra="Pop=(MK) IOP(6/50=1)",
            addsec="%s.esp"%file_base,
    )
    return lig


def runLigPrep(file_name):
    "Starting ligand preparetion process... "
    mol_path= "%s/%s"%(structure_dir, file_name)

    file_base = file_name.split(".")[0]
    #create destination directory
    WORK_DIR = file_base
    if os.path.exists(WORK_DIR):
        shutil.rmtree(WORK_DIR)
    os.mkdir(WORK_DIR)

    #Flags
    # default mm calculator set to False
    mmCalculator=False
    # default adding H is False
    addH = False

    # if desire adding H by openbabel
    prefix = ""
    if add_hydrogen:
        addH = True
        prefix += "addH_"
    if optimization_lig or optimization_conf:
        prefix += "opt_"

    # initialize ligPrep
    lig = ligPrep(mol_path, addH, WORK_DIR)
    #  lig.writeRWMol2File("test/test.xyz")

    if "ani2x" in calculator_type.lower():
        lig.setANI2XCalculator()
    elif "g16" in calculator_type.lower():
        lig = setG16calculator(lig, file_base, label="calculation", WORK_DIR=WORK_DIR)
    elif "uff" in calculator_type.lower():
        if optimization_conf:
            print("UFF calculator not support optimization")
            sys.exit(1)
        else:
            mmCalculator=True

    # set optimizetion parameters
    lig.setOptParams(fmax=thr_fmax, maxiter=1000)

    if pre_optimization_lig:
        print("G16 Optimization process.. before generations")
        lig.geomOptimization()

    if genconformer:
        out_file_path="%s/%sminE_conformer.sdf"%(WORK_DIR, prefix)
        lig.genMinEGonformer(
            file_path=out_file_path,
            numConfs=num_conformers,
            maxAttempts=max_attempts,
            pruneRmsThresh=prune_rms_thresh,
            mmCalculator=mmCalculator,
            optimization_conf=optimization_conf,
        )

        print("Conformer generation process is done")
        if not optimization_conf and optimization_lig:
            print("Optimization for minumum energy conformer")
            lig.geomOptimization()
    else:
        out_file_path="%s/%s%s.sdf"%(WORK_DIR, prefix, file_base)
        # geometry optimizaton for ligand
        if  optimization_lig:
            #  ase_atoms = lig.rwMol2AseAtoms()
            lig.geomOptimization()

    # write minimun energy conformer to sdf file
    lig.writeRWMol2File(out_file_path)

    ## g16 calculation for generate ESP chage ###
    #  for the bug of reading sfd file which have charges in ase
    try:
        atoms = read(out_file_path)
    except:
        out_file_path="%s/%s%s.xyz"%(WORK_DIR, prefix, file_base)
        lig.writeRWMol2File(out_file_path)
        atoms = read(out_file_path)

    atoms = lig.rwMol2AseAtoms()

    label="esp_calculation"
    lig = setG16calculator(lig, file_base, label=label, WORK_DIR=WORK_DIR)

    # for the bug of reading gaussian calculation log file in ase
    try:
        lig.calcSPEnergy()
    except:
        if os.path.exists(WORK_DIR + "/g16_" + label + "/" + file_base + ".esp"):
            pass
        else:
            print("Error: Could not generated ESP charges..")
            sys.exit(1)
# Get the directory where the current script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
# Construct the full path to ligPrepHelper.sh
helper_script = os.path.join(script_dir, "ligPrepHelper.sh")

# Call the shell script using the full path


if __name__ == "__main__":
    file_names = [item for item in os.listdir(structure_dir) if not item.startswith(".")]
    failed_csv = open("failed_files.csv", "w")
    failed_csv.write("FileNames,\n")

    for file_name in file_names:
        file_base = file_name.split(".")[0]
        try:
            runLigPrep(file_name)
            os.system(f"bash {helper_script} {file_base}")
        except:
            print("Error for %s file !!! Skipping...")
            failed_csv.write(file_name+",\n")
        #  break
    failed_csv.close()

(gmxMMPBSA) modellab@node3:~/workspace/burak/MNSol_database$ cat run_ligPrep.sh
#! /usr/bin/env bash

ligPrep_DIR="$HOME/mdutils/automd/ligprep"
PYTHON_DIR="$HOME/miniconda3/envs/gmxMMPBSA/bin"

struct_dir=xyz

# adding hydrogen if missing (yes/no)
add_hydrogen=no

# optimization ligand if desired before conformer generation (yes/no)
pre_optimization_lig=no

# generate conformer if desired (yes/no)
genconformer=no

#configuration for conformer generator parameters
num_conformers=1
max_attempts=5000
prune_rms_thresh=0.5

# select caclulator type (ani2x/g16) for optimization conf
caculator_type=g16
#caculator_type="ani2x"

# perform geometry optimization for conformers if desired (yes/no)
optimization_conf=no

# perform geometry optimization for orginal ligand if desired (yes/no)
optimization_lig=no

# set number of procssors for g16 calcultor (default=all cpu)
nprocs=12

# set thrshold fmax for optimization (default=0.01)
thr_fmax=0.2

#maximum iteration for optimization
maxiter=500


$PYTHON_DIR/python $ligPrep_DIR/runLigPrep.py $struct_dir $add_hydrogen $caculator_type\
        $optimization_conf $optimization_lig $pre_optimization_lig $genconformer\
        $nprocs $thr_fmax $maxiter $num_conformers $max_attempts $prune_rms_thresh
