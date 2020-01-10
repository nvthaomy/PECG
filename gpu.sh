#!/bin/bash
#SBATCH -N 1 --partition=gpu --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --job-name=xp0.1_N12_f0_V157_LJPME_298K_NVT_Spline_offsetBond_testEqomm_6
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=my@ucsb.edu

cd $SLURM_SUBMIT_DIR

/bin/hostname
srun --gres=gpu:1 /usr/bin/nvidia-smi
export PATH=/home/mnguyen/lammps-7Aug19/bin/:$PATH
export PATH="/home/mnguyen/miniconda3/envs/py2/bin/:$PATH"
export PYTHONPATH=/home/mnguyen/bin/sim_git:$PYTHONPATH
python main.py

