#!/bin/bash
#SBATCH --job-name=soil_r23_alpha_div
#SBATCH --output=logs/soil_r23_alpha_div.out
#SBATCH --error=logs/soil_r23_alpha_div.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

# Calculate alpha diversity metrics for global soil samples

# Load Conda R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run script
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/23_calculate_alpha_diversity.R
