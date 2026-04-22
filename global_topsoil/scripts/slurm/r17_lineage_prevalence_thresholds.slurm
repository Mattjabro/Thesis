#!/bin/bash
#SBATCH --job-name=soil_r17_lineage_prevalence
#SBATCH --output=logs/soil_r17_lineage_prevalence.out
#SBATCH --error=logs/soil_r17_lineage_prevalence.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Load Conda R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/17_lineage_prevalence_thresholds.R
