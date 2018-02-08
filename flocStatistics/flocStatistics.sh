#!/bin/bash -l
#SBATCH -N 1 -n 1
#SBATCH --job-name="flocStat"
#SBATCH -t 0-0:20:0
#SBATCH --output flocStat_out%j
#SBATCH --partition=slurm_shortgpu
cd $SLURM_SUBMIT_DIR
module load cuda/8.0 
./flocStatistics
rm flocStat_* 
