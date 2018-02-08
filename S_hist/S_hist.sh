#!/bin/bash -l
#SBATCH -n 1 
#SBATCH --job-name="S_hist"
#SBATCH -t 0-0:20:0
#SBATCH --output S_hist_out%j 
#SBATCH --partition=slurm_courtesy,slurm_shortgpu
cd $SLURM_SUBMIT_DIR
hostname
./S_hist
