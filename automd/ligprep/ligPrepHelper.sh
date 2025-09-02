#! /usr/bin/env bash
ligPrep_DIR="$HOME/mdutils/automd/ligprep"
export GMX_MAXBACKUP=-1
export GMX_MAXWARN=1




file_base=$1
cd $file_base

mkdir MDRUNS MDRUNS/EM_1 MDRUNS/EM_2 MDRUNS/NVT MDRUNS/NPT MDRUNS/Production_MD
lig=MOL
cp g16_esp_calculation/$file_base.esp $lig.esp

# rm *.gro *.top *.itp *.rtp *.ndx
antechamber -i "$lig".esp -fi gesp -o "$lig".mol2 -fo mol2 -c resp -s 2 -rn "$lig" -at gaff2 -nc 0 -pf yes
parmchk2 -i "$lig".mol2 -f mol2 -o "$lig".frcmod

tleap -f- <<EOF
source leaprc.gaff2
loadamberparams $lig.frcmod
$lig = loadmol2 $lig.mol2
saveamberparm $lig $lig.prmtop $lig.inpcrd
quit
EOF

amb2gro_top_gro.py -p "$lig".prmtop -c "$lig".inpcrd -t "$lig"_GMX.top -g "$lig"_GMX.gro -b "$lig"_GMX.pdb

 cp $ligPrep_DIR/templates/topol.ligsol topol.top
 sed -n '/atomtypes/,/moleculetype/{/moleculetype/b;p}' "$lig"_GMX.top > "$lig"_GMX_GAFF.rtp
 sed -n '/moleculetype/,/system/{/system/b;p}' "$lig"_GMX.top > "$lig"_GMX.itp
 sed -i s/LIGAND/"$lig"/g topol.top

 echo 0 | gmx genrestr -f "$lig".gro -o posre_ligand.itp -fc 1000 1000 1000
 gmx editconf -f "$lig"_GMX.gro -o newbox.gro -bt cubic -d 1.6
 gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
 gmx make_ndx -f solv.gro -o index.ndx <<EOF
 q
 EOF


 cd ..
