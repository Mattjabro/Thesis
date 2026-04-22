#!/bin/bash
#SBATCH --job-name=emp_r03_collapse_genus_all
#SBATCH --output=logs/emp_r03_collapse_genus_all.out
#SBATCH --error=logs/emp_r03_collapse_genus_all.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=1028G

# Adjust resources as needed for your cluster.

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/03_collapse_genus_all_samples.R
