#!/bin/bash
#SBATCH --job-name=soil_r10_stats
#SBATCH --output=logs/soil_r10_stats.out
#SBATCH --error=logs/soil_r10_stats.err
#SBATCH --time=01:00:00
#SBATCH --mem=1028G
#SBATCH --cpus-per-task=100

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/10_stats.R
