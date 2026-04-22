#!/bin/bash
#SBATCH --job-name=emp_r22_tsne
#SBATCH --output=logs/emp_r22_tsne.out
#SBATCH --error=logs/emp_r22_tsne.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# t-SNE dimensionality reduction analysis

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/22_tsne_analysis.R
