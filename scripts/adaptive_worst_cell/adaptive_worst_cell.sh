#!/bin/bash
#SBATCH --job-name=adaptive_worst_cell
#SBATCH --time=06:00:00
#SBATCH --partition=long-cascadelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4                        ## 4 cores for Stan's 4 chains
#SBATCH --mem=8G
#SBATCH --output=/home/viro3403/sampling-paper/logs/adaptive_worst_cell_logs/slurm_%j_%x.out
#SBATCH --error=/home/viro3403/sampling-paper/logs/adaptive_worst_cell_logs/slurm_%j_%x.err

## The log directory must exist before submission:
##   mkdir -p ~/sampling-paper/logs/adaptive_worst_cell_logs

source ~/miniforge3/etc/profile.d/conda.sh
conda activate rstan-clean
cd ~/sampling-paper

Rscript scripts/adaptive_worst_cell/adaptive_worst_cell.R
