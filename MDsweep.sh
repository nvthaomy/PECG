#!/bin/bash
# ask for 16 cores on two nodes
#SBATCH --nodes=1 --ntasks-per-node=6
#SBATCH --time=700:00:00
#SBATCH --job-name=2Mnacl_opc_Uext2_NaCl_298K_NVT_fixUwater
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=my@ucsb.edu

module load intel/18

ext=PAANPT_
Name="'PAA'"
ff="'PAA_ff.dat'"

dirNameA=('N12' 'N24' 'N48' 'N60' 'N90')
nPAAs=(15 8 12 10 7)
nNaA=(0 0 0 0 0)
nClA=(0 0 0 0 0)
nHOHA=(4728 4728 15372 15372 15372)
volA=(5270. 5300. 17078. 17119. 17219.)
MolNamesList=["'PAA','HOH'"]

PAAstructureA=(["'A'"]*12 ["'A'"]*24 ["'A'"]*48 ["'A'"]*60 ["'A'"]*90)
#MD
dt=0.05 #0.1
P=8.52
tau=6000.
equilTau=100.
Stride=40
cut=8.
UseOMM=True
UseLammps=False
OMP_NumThread=4
OPENMM_CPU_THREADS=$OMP_NumThread
#FEP

FEPMolNames=["'Na+','Cl-'"]
ThermoSlice=1
TrajSlice=1
nInsert=10
nDelete=10
FEPDir="'${nInsert}insert_NaCl'"

#for converting traj to real unit
#top="'nacl0_initial.pdb'"
#traj="'nacl0_traj.dcd'"
#scale=3.1

# get the length of the arrays
#length=${#nNaA[@]}
length=${#nPAAs[@]}

#echo $length 'concentration values'

# do the loop
for ((i=0;i<$length;i++)); do
    nPAA=${nPAAs[$i]}
    nNa=${nNaA[$i]}
    nCl=${nClA[$i]}
    nHOH=${nHOHA[$i]}
    vol=${volA[$i]}
    PAAstructure=${PAAstructureA[$i]}

    mydir=${dirNameA[$i]}
#    mydir=${nNa}Na_${nCl}Cl_${nHOH}HOH
    mkdir $mydir
    cp * $mydir/.
    mv $mydir/MDmain_template.py $mydir/MDmain.py
    mv $mydir/pod_template.sh $mydir/pod.sh
    echo === ${mydir} ===
    sed -i "s/__Name__/${Name}/g" $mydir/MDmain.py
    sed -i "s/__nNa__/${nNa}/g" $mydir/MDmain.py 
    sed -i "s/__nCl__/${nCl}/g" $mydir/MDmain.py
    sed -i "s/__nHOH__/${nHOH}/g" $mydir/MDmain.py
    sed -i "s/__nPAA__/${nPAA}/g" $mydir/MDmain.py
    sed -i "s/__PAAstructure__/${PAAstructure}/g" $mydir/MDmain.py
    sed -i "s/__MolNamesList__/${MolNamesList}/g" $mydir/MDmain.py

    sed -i "s/__vol__/${vol}/g" $mydir/MDmain.py
    sed -i "s/__dt__/${dt}/g" $mydir/MDmain.py     
    sed -i "s/__P__/${P}/g" $mydir/MDmain.py
    sed -i "s/__tau__/${tau}/g" $mydir/MDmain.py
    sed -i "s/__Stride__/${Stride}/g" $mydir/MDmain.py
    sed -i "s/__ff__/${ff}/g" $mydir/MDmain.py
    sed -i "s/__cut__/${cut}/g" $mydir/MDmain.py
    sed -i "s/__UseOMM__/${UseOMM}/g" $mydir/MDmain.py
    sed -i "s/__UseLammps__/${UseLammps}/g" $mydir/MDmain.py
    sed -i "s/__OMP_NumThread__/${OMP_NumThread}/g" $mydir/MDmain.py
    sed -i "s/__ThermoSlice__/${ThermoSlice}/g" $mydir/MDmain.py
    sed -i "s/__TrajSlice__/${TrajSlice}/g" $mydir/MDmain.py
    sed -i "s/__nInsert__/${nInsert}/g" $mydir/MDmain.py
    sed -i "s/__nDelete__/${nDelete}/g" $mydir/MDmain.py
    sed -i "s/__FEPMolNames__/${FEPMolNames}/g" $mydir/MDmain.py
    sed -i "s/__FEPDir__/${FEPDir}/g" $mydir/MDmain.py
    sed -i "s/__equilTau__/${equilTau}/g" $mydir/MDmain.py

    sed -i "s/__OMP_NumThread__/${OMP_NumThread}/g" $mydir/pod.sh
    sed -i "s/__jobName__/${ext}${mydir}/g" $mydir/pod.sh
    sed -i "s/__pyName__/MDmain.py/g" $mydir/pod.sh    
#    sed -i "s/__top__/${top}/g" $mydir/pod.sh
#    sed -i "s/__traj__/${traj}/g" $mydir/pod.sh
#    sed -i "s/__scale__/${scale}/g" $mydir/pod.sh
    cd $mydir
    qsub pod.sh
    cd ..
done
