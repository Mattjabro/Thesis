#!/bin/bash
#SBATCH --job-name=soil_r21_clr_avg
#SBATCH --output=logs/soil_r21_cluster_clr_avg.out
#SBATCH --error=logs/soil_r21_cluster_clr_avg.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Load Conda R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run script
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/21_cluster_clr_averages.R
