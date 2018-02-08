#!/bin/bash -l
#SBATCH -N 1 -n 1
#SBATCH --job-name="flocExtract"
#SBATCH -t 0-0:20:0
#SBATCH --output out%j
#SBATCH --partition=slurm_shortgpu
cd $SLURM_SUBMIT_DIR
module load cuda/8.0 
./flocTrajectory
rm out*
