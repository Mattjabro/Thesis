#!/bin/bash
#SBATCH --job-name=emp_r16_clr_avg
#SBATCH --output=logs/emp_r16_cluster_clr_avg.out
#SBATCH --error=logs/emp_r16_cluster_clr_avg.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Load Conda R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run script
Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/16_cluster_clr_averages.R
