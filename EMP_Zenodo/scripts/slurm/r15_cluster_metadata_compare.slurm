#!/bin/bash
#SBATCH --job-name=emp_r15_cluster_metadata
#SBATCH --output=logs/emp_r15_cluster_metadata.out
#SBATCH --error=logs/emp_r15_cluster_metadata.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=1028G

# Adjust resources as needed for your cluster.

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/15_cluster_metadata_compare.R
