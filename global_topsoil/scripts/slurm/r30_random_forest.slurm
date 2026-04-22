#!/bin/bash
#SBATCH --job-name=soil_r30_random_forest
#SBATCH --output=logs/soil_r30_random_forest.out
#SBATCH --error=logs/soil_r30_random_forest.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Train Random Forest classification models for biome prediction

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/30_random_forest.R
