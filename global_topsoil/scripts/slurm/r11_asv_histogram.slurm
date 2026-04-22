#!/bin/bash
#SBATCH --job-name=soil_r11_asv_hist
#SBATCH --output=logs/soil_r11_asv_hist.out
#SBATCH --error=logs/soil_r11_asv_hist.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Load Conda R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/11_asv_histogram.R
