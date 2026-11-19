#!/bin/bash
#SBATCH --job-name=question1
#SBATCH --time=8:00:00                           ## Per-task time limit (not total)
#SBATCH --partition=long-cascadelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4                        ## 4 cores for Stan's 4 chains
#SBATCH --mem=8G                                 ## Per-task memory
#SBATCH --array=1-20                             ## One task per l_0
#SBATCH --output=/home/viro3403/sampling-paper/logs/question1_logs/slurm_%A_%a_%x.out   ## %A = array job ID, %a = task ID
#SBATCH --error=/home/viro3403/sampling-paper/logs/question1_logs/slurm_%A_%a_%x.err

## Before submission, from ~/sampling-paper:
##   mkdir -p logs/question1_logs
##   Rscript scripts/question1/setup_param_array.R

source ~/miniforge3/etc/profile.d/conda.sh
conda activate rstan-clean
cd ~/sampling-paper

Rscript scripts/question1/run_condition.R
