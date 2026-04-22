#!/bin/bash
#SBATCH --job-name=soil_r31_model_eval
#SBATCH --output=logs/soil_r31_model_eval.out
#SBATCH --error=logs/soil_r31_model_eval.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

# Evaluate Random Forest models and extract feature importance

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/31_model_evaluation.R
