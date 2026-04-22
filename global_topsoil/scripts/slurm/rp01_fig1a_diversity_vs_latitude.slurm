#!/bin/bash
#SBATCH --job-name=soil_rp01_fig1a
#SBATCH --output=logs/soil_rp01_fig1a.out
#SBATCH --error=logs/soil_rp01_fig1a.err
#SBATCH --time=20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

# Fig 1a: Bacterial diversity vs. latitude (Bahram 2018 replication)
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/p01_fig1a_diversity_vs_latitude.R
