#!/bin/bash
#SBATCH --job-name=soil_r25_clr_stats
#SBATCH --output=logs/soil_r25_clr_stats.out
#SBATCH --error=logs/soil_r25_clr_stats.err
#SBATCH --time=30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Calculate cluster CLR statistics with error bars

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/25_cluster_clr_statistics.R
