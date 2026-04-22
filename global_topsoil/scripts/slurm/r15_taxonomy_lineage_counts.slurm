#!/bin/bash
#SBATCH --job-name=soil_r15_lineage_counts
#SBATCH --output=logs/soil_r15_lineage_counts.out
#SBATCH --error=logs/soil_r15_lineage_counts.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Load Conda R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/15_taxonomy_lineage_counts.R
