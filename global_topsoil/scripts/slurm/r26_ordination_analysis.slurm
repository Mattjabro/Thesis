#!/bin/bash
#SBATCH --job-name=soil_r26_ordination
#SBATCH --output=logs/soil_r26_ordination.out
#SBATCH --error=logs/soil_r26_ordination.err
#SBATCH --time=45:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# PCA/PCoA ordination analysis

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/26_ordination_analysis.R
