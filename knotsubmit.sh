#!/bin/bash
#PBS -l nodes=1:ppn=5
#PBS -l walltime=700:00:00
#PBS -V
#PBS -N xp0.1_N12_f0_V157_LJPME_298K_NVT 
#PBS -M my@ucsb.edu
#PBS -m ae

cd $PBS_O_WORKDIR
OPENMM_CPU_THREADS=
#export PATH=/home/mnguyen/bin/lammps/lammps-12Dec18/bin/:$PATH
export PATH=/home/mnguyen/lammps-7Aug19/bin/:$PATH
export PATH="/home/mnguyen/miniconda3/envs/py2/bin/:$PATH"
export PYTHONPATH=/home/mnguyen/bin/sim_git:$PYTHONPATH
python main.py

