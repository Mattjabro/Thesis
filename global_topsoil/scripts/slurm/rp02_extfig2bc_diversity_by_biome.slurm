#!/bin/bash
#SBATCH --job-name=soil_rp02_extfig2bc
#SBATCH --output=logs/soil_rp02_extfig2bc.out
#SBATCH --error=logs/soil_rp02_extfig2bc.err
#SBATCH --time=20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

# Ext Fig 2b/c: Bacterial diversity across biomes (Bahram 2018 replication)
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/p02_extfig2bc_diversity_by_biome.R
