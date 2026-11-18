#!/bin/bash
#SBATCH --job-name=nonadaptive_experiment
#SBATCH --time=24:00:00                         
#SBATCH --partition=long-cascadelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4                        ## 4 cores for Stan's 4 chains
#SBATCH --mem=8G
#SBATCH --output=/home/viro3403/sampling-paper/logs/nonadaptive_experiment_logs/slurm_%j_%x.out
#SBATCH --error=/home/viro3403/sampling-paper/logs/nonadaptive_experiment_logs/slurm_%j_%x.err

source ~/miniforge3/etc/profile.d/conda.sh
conda activate rstan-clean
cd ~/sampling-paper

Rscript scripts/nonadaptive_experiment/nonadaptive_experiment.R
