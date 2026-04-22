#!/bin/bash
#SBATCH --job-name=emp_r14_model_fit
#SBATCH --output=logs/emp_r14_model_fit.out
#SBATCH --error=logs/emp_r14_model_fit.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=1028G

# Load Conda R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run script
Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/14_model_fit_mixture.R
