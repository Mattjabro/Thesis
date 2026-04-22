#!/bin/bash
#SBATCH --job-name=soil_r27_tsne
#SBATCH --output=logs/soil_r27_tsne.out
#SBATCH --error=logs/soil_r27_tsne.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# t-SNE dimensionality reduction analysis

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/27_tsne_analysis.R
