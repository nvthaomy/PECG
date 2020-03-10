#!/bin/bash
# ask for 16 cores on two nodes
#SBATCH --nodes=1 --ntasks-per-node=6
#SBATCH --time=700:00:00
#SBATCH --job-name=2Mnacl_opc_Uext2_NaCl_298K_NVT_fixUwater
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=my@ucsb.edu

module load intel/18

cd $SLURM_SUBMIT_DIR
/bin/hostname
#export PATH=/home/mnguyen/bin/lammps/lammps-12Dec18/bin/:$PATH
export PATH=/home/mnguyen/lammps-7Aug19/bin/:$PATH
export PATH="/home/mnguyen/miniconda3/envs/py2/bin/:$PATH"
export PYTHONPATH=/home/mnguyen/bin/sim_git:$PYTHONPATH
#python main.py

ext=NaClUext2_
Name="'nacl'"
ff="'nacl_ff.dat'"

dirNameA=('0.5M' '1M' '3M' '5M')
nPAAs=(0 0 0 0)
nNaA=(61) # 122 368 613)
nClA=(61 122 368 613)
nHOHA=(6708 6616 6247 5880)
volA=(6840. 6840. 6840. 6840. )
nInsertA=(20 50)
nDeleteA=(20 50)
MolNamesList=["'Na+','Cl-','HOH'"]

#MD
dt=0.1
P=8.52
tau=5000.
equilTau=100.
Stride=20
cut=6.
UseOMM=True
UseLammps=False
OMP_NumThread=4
OPENMM_CPU_THREADS=$OMP_NumThread
#FEP
FEPMolNames=["'HOH'"]
ThermoSlice=1
TrajSlice=1
#nInsert=20
#nDelete=20
MolName=HOH
#for converting traj to real unit
#top="'nacl0_initial.pdb'"
#traj="'nacl0_traj.dcd'"
#scale=3.1

# get the length of the arrays
length=${#nNaA[@]}
l2=${#nInsertA[@]}
echo $length 'concentration values'

# do the loop
for ((i=0;i<$length;i++)); do
    nPAA=${nPAAs[$i]}
    nNa=${nNaA[$i]}
    nCl=${nClA[$i]}
    nHOH=${nHOHA[$i]}
    PAAstructure=["'A'"]*12
    vol=${volA[$i]}
    mydir=${dirNameA[$i]}

    for ((j=0;j<$l2;j++)); do
    nInsert=${nInsertA[$j]}
    nDelete=${nDeleteA[$j]}
    FEPDir=${nInsert}insert_$MolName
    cp MDmain_template.py $mydir/FEP_${FEPDir}.py
    cp ~/bin/PECG/FEP_Module_CG.py $mydir/.
    cp pod_template.sh $mydir/pod_${FEPDir}.sh
    echo === ${mydir} ===
    echo $FEPDir
    sed -i "s/__Name__/${Name}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__nPAA__/${nPAA}/g" $mydir/FEP_${FEPDir}.py 
    sed -i "s/__nNa__/${nNa}/g" $mydir/FEP_${FEPDir}.py 
    sed -i "s/__nCl__/${nCl}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__nHOH__/${nHOH}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__vol__/${vol}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__PAAstructure__/${PAAstructure}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__MolNamesList__/${MolNamesList}/g" $mydir/FEP_${FEPDir}.py

    sed -i "s/__dt__/${dt}/g" $mydir/FEP_${FEPDir}.py     
    sed -i "s/__P__/${P}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__tau__/${tau}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__Stride__/${Stride}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__ff__/${ff}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__cut__/${cut}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__UseOMM__/${UseOMM}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__UseLammps__/${UseLammps}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__OMP_NumThread__/${OMP_NumThread}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__ThermoSlice__/${ThermoSlice}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__TrajSlice__/${TrajSlice}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__nInsert__/${nInsert}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__nDelete__/${nDelete}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__FEPDir__/'${FEPDir}'/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__FEPMolNames__/${FEPMolNames}/g" $mydir/FEP_${FEPDir}.py
    sed -i "s/__equilTau__/${equilTau}/g" $mydir/FEP_${FEPDir}.py

    sed -i "s/__OMP_NumThread__/${OMP_NumThread}/g" $mydir/pod_${FEPDir}.sh
    sed -i "s/__jobName__/${ext}${mydir}_${FEPDir}/g" $mydir/pod_${FEPDir}.sh
    sed -i "s/__pyName__/FEP_${FEPDir}.py/g" $mydir/pod_${FEPDir}.sh
    cd $mydir
    qsub pod_${FEPDir}.sh
    cd ..
    done
done
