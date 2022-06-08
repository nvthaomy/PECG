module load intel/18
export PATH=/home/mnguyen/lammps-7Aug19/bin/:$PATH
export PATH="/home/mnguyen/miniconda3/envs/py2/bin/:$PATH"
export PYTHONPATH=/home/mnguyen/bin/sim_nsherck:$PYTHONPATH

dirPrefix=(4IDP-2Yx2A 4IDP-2Yx2A)
dirExt=(0M 0M)
Name="'poly'"
ff=poly_ff.dat

MolNamesList=["'Pol','HOH'"] 
nPols=(10 5)
nNaA=(0 0)
nClA=(0 0)
nHOHA=(25000 25000)
PolstructureA=("['Z', 'Z', 'Z', 'Z-', 'Z', 'Z+'] * 4 + ['X'] * 20" "['Z', 'Z', 'Z', 'Z-', 'Z', 'Z+'] * 4 + ['X'] * 20")

#MD
InitPDB=None
Lxs=(35. 35)
Lys=(35. 35)
Lzs=(35. 35)
dt=0.1 #0.1
P=8.520
PresAx=None #0 1 2
tau=200000.
equilTau=500.
Stride=1000
cut=8.
UseOMM=True
UseLammps=False
UseSim=False
OMP_NumThread=8
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
length=${#nPols[@]}

#echo $length 'concentration values'

# do the loop
for ((i=0;i<$length;i++)); do
    nPol=${nPols[$i]}
    nNa=${nNaA[$i]}
    nCl=${nClA[$i]}
    nHOH=${nHOHA[$i]}
    Lx=${Lxs[$i]}
    Ly=${Lys[$i]}
    Lz=${Lzs[$i]}
    Polstructure=${PolstructureA[$i]}
    mydir=${dirPrefix[$i]}_${nPol}IDP_${nNa}NaCl_${nHOH}HOH_${dirExt[$i]}

    mkdir $mydir
    cp $ff $mydir/.
    cp MDmain_template.py $mydir/MDmain.py
    cp pod_template.sh $mydir/pod.sh
    echo === ${mydir} ===
    sed -i "s/__Name__/${Name}/g" $mydir/MDmain.py
    sed -i "s/__nNa__/${nNa}/g" $mydir/MDmain.py 
    sed -i "s/__nCl__/${nCl}/g" $mydir/MDmain.py
    sed -i "s/__nHOH__/${nHOH}/g" $mydir/MDmain.py
    sed -i "s/__nPol__/${nPol}/g" $mydir/MDmain.py
    sed -i "s/__Polstructure__/${Polstructure}/g" $mydir/MDmain.py
    sed -i "s/__MolNamesList__/${MolNamesList}/g" $mydir/MDmain.py

    sed -i "s/__Lx__/${Lx}/g" $mydir/MDmain.py
    sed -i "s/__Ly__/${Ly}/g" $mydir/MDmain.py
    sed -i "s/__Lz__/${Lz}/g" $mydir/MDmain.py
    sed -i "s/__InitPDB__/${InitPDB}/g" $mydir/MDmain.py
    sed -i "s/__dt__/${dt}/g" $mydir/MDmain.py     
    sed -i "s/__P__/${P}/g" $mydir/MDmain.py
    sed -i "s/__PresAx__/${PresAx}/g" $mydir/MDmain.py
    sed -i "s/__tau__/${tau}/g" $mydir/MDmain.py
    sed -i "s/__Stride__/${Stride}/g" $mydir/MDmain.py
    sed -i "s/__ff__/'${ff}'/g" $mydir/MDmain.py
    sed -i "s/__cut__/${cut}/g" $mydir/MDmain.py
    sed -i "s/__UseOMM__/${UseOMM}/g" $mydir/MDmain.py
    sed -i "s/__UseLammps__/${UseLammps}/g" $mydir/MDmain.py
    sed -i "s/__UseSim__/${UseSim}/g" $mydir/MDmain.py
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
    cd $mydir
#    qsub pod.sh
    python MDmain.py
    python ~/bin/scripts/mdtrajConvert.py ${Name}0_traj.dcd ${Name}0_equilibrated.pdb xyz -s 20
    cd ..
done
