#!/bin/bash
#SBATCH --job-name=soil_r07_read_tracking
#SBATCH --output=logs/soil_r07_read_tracking.out
#SBATCH --error=logs/soil_r07_read_tracking.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=1028G

# Load Conda environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run tracking summary R script
echo " Running Step 7: Read Tracking Summary"
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/07_read_tracking.R
