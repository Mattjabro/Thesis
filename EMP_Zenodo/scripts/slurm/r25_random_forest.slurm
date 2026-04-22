#!/bin/bash
#SBATCH --job-name=emp_r25_random_forest
#SBATCH --output=logs/emp_r25_random_forest.out
#SBATCH --error=logs/emp_r25_random_forest.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Train Random Forest classification models for biome prediction

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/25_random_forest.R
