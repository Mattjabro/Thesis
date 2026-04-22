#!/bin/bash
#SBATCH --job-name=soil_r24_alpha_cluster
#SBATCH --output=logs/soil_r24_alpha_cluster.out
#SBATCH --error=logs/soil_r24_alpha_cluster.err
#SBATCH --time=30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Plot alpha diversity by cluster

# Load Conda R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run script
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/24_plot_alpha_diversity_by_cluster.R
