#!/bin/bash
#SBATCH --job-name=emp_r12_lineage_counts_thresholds
#SBATCH --output=logs/emp_r12_lineage_counts_thresholds.out
#SBATCH --error=logs/emp_r12_lineage_counts_thresholds.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G

# Adjust resources as needed for your cluster.

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/12_lineage_counts_thresholds.R
