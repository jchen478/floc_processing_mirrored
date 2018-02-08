#!/bin/bash -l
#SBATCH -n 1 -N 1
#SBATCH --job-name="sizeChar"
#SBATCH -t 0-0:20:0
#SBATCH --output calcSize%j 
#SBATCH --partition=slurm_courtesy,slurm_shortgpu
#SBATCH --nodelist=euler02
#SBATCH --gres=gpu:gtx1080:1
cd $SLURM_SUBMIT_DIR
hostname
module load cuda/8.0 
./sizeChar
#rm calcSize*
