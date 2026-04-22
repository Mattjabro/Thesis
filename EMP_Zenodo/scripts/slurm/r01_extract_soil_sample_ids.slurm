#!/bin/bash
#SBATCH --job-name=emp_r01_extract_soil_ids
#SBATCH --output=logs/emp_r01_extract_soil_ids.out
#SBATCH --error=logs/emp_r01_extract_soil_ids.err
#SBATCH --time=10:10:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=1028G

# Activate Conda environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run the R script
echo "Running soil sample filter..."
Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/01_extract_soil_sample_ids.R
