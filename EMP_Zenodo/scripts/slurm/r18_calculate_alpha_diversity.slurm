#!/bin/bash
#SBATCH --job-name=emp_r18_alpha_div
#SBATCH --output=logs/emp_r18_alpha_div.out
#SBATCH --error=logs/emp_r18_alpha_div.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G

# Calculate alpha diversity metrics for soil samples

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/18_calculate_alpha_diversity.R
