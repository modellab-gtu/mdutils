#!/bin/bash
#SBATCH -p orfoz
#SBATCH -A akocak
#SBATCH -J qmx_test
#SBATCH -N 1
#SBATCH --ntasks=112
#SBATCH --cpus-per-task=1
#SBATCH -C weka
#SBATCH --time=3-00:00:00

echo "SLURM_NODELIST $SLURM_NODELIST"
echo "NUMBER OF CORES $SLURM_NTASKS"

echo "number of CPUs=  $SLURM_NPROCS"



#export GMX_MAXCONSTRWARN=-1
export GMX_MAXBACKUP=-1

workdir=$(pwd)
MDPDIR="$workdir"/MDP
LIGAND=$1
#for lig in comp_3a_26_docked_prep  comp_3b_30_docked_prep  comp_3c_30_docked_prep  comp_3d_30_docked_prep  comp_3e_30_docked_prep
#do
cd "$LIGAND"

for i in {0..4}
do
    gmx grompp -f "$MDPDIR"/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt_"$i".tpr -n index.ndx -maxwarn 1
    mpirun -np 112 gmx_mpi mdrun -v -deffnm nvt_"$i"

    gmx grompp -f "$MDPDIR"/NPT/npt.mdp -c nvt_"$i".gro -r nvt_"$i".gro -p topol.top -o npt_"$i".tpr -n index.ndx
    mpirun -np 112 gmx_mpi mdrun -v -deffnm npt_"$i"

    gmx grompp -f "$MDPDIR"/Production_MD/md.mdp -c npt_"$i".gro -r npt_"$i".gro -p topol.top -o md_"$i".tpr -n index.ndx
    mpirun -np 112 gmx_mpi mdrun -v -deffnm md_"$i"
done
cd $workdir
#done


exit
