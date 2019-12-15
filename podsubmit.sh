#!/bin/bash
# ask for 16 cores on two nodes
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --time=700:00:00
#SBATCH --job-name=xp0.1_N12_f0_V157_LJPME_298K_NVT
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=my@ucsb.edu

module load intel/18

cd $SLURM_SUBMIT_DIR
/bin/hostname
OPENMM_CPU_THREADS=THREADS_DUMMY
export PATH="/home/mnguyen/miniconda3/envs/py2/bin/:$PATH"
export PYTHONPATH=/home/mnguyen/bin/sim_git:$PYTHONPATH
export PATH="/home/nsherck/lammps-22Aug18_test/bin:$PATH"
python main.py 
