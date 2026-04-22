#!/bin/bash
#SBATCH --job-name=soil_r14_taxa_thresholds_bact
#SBATCH --output=logs/soil_r14_taxa_thresholds_bact.out
#SBATCH --error=logs/soil_r14_taxa_thresholds_bact.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Load Conda R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/14_taxa_presence_thresholds_bacteria.R
