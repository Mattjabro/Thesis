#!/bin/bash
#SBATCH --job-name=emp_r07_collapse_genus_soil
#SBATCH --output=logs/emp_r07_collapse_genus_soil.out
#SBATCH --error=logs/emp_r07_collapse_genus_soil.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=1028G

# Activate Conda environment with DADA2
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run the R script
Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/07_collapse_genus_soil_only.R
