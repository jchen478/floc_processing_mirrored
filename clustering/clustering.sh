#!/bin/bash -l
#SBATCH -N 1 -n 1
#SBATCH --job-name="clustering"
#SBATCH -t 0-0:20:0
#SBATCH --partition=slurm_courtesy,slurm_shortgpu
#SBATCH --gres=gpu:gtx680:1
#SBATCH --output out%j.out
hostname
cd $SLURM_SUBMIT_DIR
module load cuda/8.0 
./clustering
