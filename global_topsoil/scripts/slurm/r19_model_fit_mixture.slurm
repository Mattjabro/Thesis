#!/bin/bash
#SBATCH --job-name=soil_r19_model_fit
#SBATCH --output=logs/soil_r19_model_fit.out
#SBATCH --error=logs/soil_r19_model_fit.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=1028G

# Load Conda R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run script
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/19_model_fit_mixture.R
