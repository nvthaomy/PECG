#!/bin/bash
# ask for 16 cores on two nodes
#SBATCH --nodes=1 --ntasks-per-node=__OMP_NumThread__
#SBATCH --time=700:00:00
#SBATCH --job-name=__jobName__
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=my@ucsb.edu

module load intel/18

cd $SLURM_SUBMIT_DIR
/bin/hostname
OPENMM_CPU_THREADS=__OMP_NumThread__
top=__top__
traj=__traj__
scale=__scale__

export PATH=/home/mnguyen/lammps-7Aug19/bin/:$PATH
export PATH="/home/mnguyen/miniconda3/envs/py2/bin/:$PATH"
export PYTHONPATH=/home/mnguyen/bin/sim_git:$PYTHONPATH
python __pyName__

#if test -f "AllStats.dat"; then
#    sed -i "s/__traj__/${traj}/g" convertTraj2RealUnit.py 
#    sed -i "s/__top__/${top}/g" convertTraj2RealUnit.py
#    sed -i "s/__scale__/${scale}/g" convertTraj2RealUnit.py
#    python convertTraj2RealUnit.py
#fi

