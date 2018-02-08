#!/bin/bash -l
#SBATCH -n 1
#SBATCH --job-name="extend"
#SBATCH -t 0-0:10:0
#SBATCH --partition=slurm_courtesy,slurm_shortgpu
#SBATCH --output out%j.out
hostname
cd $SLURM_SUBMIT_DIR
./extend
