#!/bin/bash -l
#SBATCH -n 1
#SBATCH --job-name="make"
#SBATCH -t 0-0:1:0
#SBATCH --partition=slurm_courtesy,slurm_shortgpu
#SBATCH --output out%j.out
hostname
cd $SLURM_SUBMIT_DIR
module load intel/psxe/2017/compiler
module load intel/psxe/2017/mkl

make
