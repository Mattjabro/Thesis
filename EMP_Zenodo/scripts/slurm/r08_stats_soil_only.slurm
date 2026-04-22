#!/bin/bash
#SBATCH --job-name=emp_r08_stats_soil_only
#SBATCH --output=logs/emp_r08_stats_soil_only.out
#SBATCH --error=logs/emp_r08_stats_soil_only.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=1028G

# Adjust resources as needed for your cluster.

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/08_stats_soil_only.R
