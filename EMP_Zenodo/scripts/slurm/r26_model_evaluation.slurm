#!/bin/bash
#SBATCH --job-name=emp_r26_model_eval
#SBATCH --output=logs/emp_r26_model_eval.out
#SBATCH --error=logs/emp_r26_model_eval.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

# Evaluate Random Forest models and extract feature importance

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/26_model_evaluation.R
