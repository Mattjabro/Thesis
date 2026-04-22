#!/bin/bash
#SBATCH --job-name=soil_r09_classify_filtered
#SBATCH --output=logs/soil_r09_classify_filtered.out
#SBATCH --error=logs/soil_r09_classify_filtered.err
#SBATCH --time=10:15:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=1028G

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

echo " Building outputs for Bacteria-only data..."
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/09_classify_filtered.R
