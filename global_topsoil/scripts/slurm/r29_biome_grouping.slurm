#!/bin/bash
#SBATCH --job-name=soil_r29_biome_grouping
#SBATCH --output=logs/soil_r29_biome_grouping.out
#SBATCH --error=logs/soil_r29_biome_grouping.err
#SBATCH --time=0:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Create ecologically meaningful biome groupings

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/29_biome_grouping.R
