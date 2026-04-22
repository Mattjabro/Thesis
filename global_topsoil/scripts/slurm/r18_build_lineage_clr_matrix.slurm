#!/bin/bash
#SBATCH --job-name=soil_r18_build_clr_matrix
#SBATCH --output=logs/soil_r18_build_clr_matrix.out
#SBATCH --error=logs/soil_r18_build_clr_matrix.err
#SBATCH --time=11:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=1028G

# Load Conda R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run script
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/18_build_lineage_clr_matrix.R
