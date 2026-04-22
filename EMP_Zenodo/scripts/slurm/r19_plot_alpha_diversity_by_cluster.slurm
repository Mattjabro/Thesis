#!/bin/bash
#SBATCH --job-name=emp_r19_alpha_cluster
#SBATCH --output=logs/emp_r19_alpha_cluster.out
#SBATCH --error=logs/emp_r19_alpha_cluster.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Plot alpha diversity by cluster

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/19_plot_alpha_diversity_by_cluster.R
