#!/bin/bash
#SBATCH --job-name=emp_r20_clr_stats
#SBATCH --output=logs/emp_r20_clr_stats.out
#SBATCH --error=logs/emp_r20_clr_stats.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Calculate cluster CLR statistics with error bars

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/20_cluster_clr_statistics.R
