#!/bin/bash
#SBATCH --job-name=soil_r28_cluster_histograms
#SBATCH --output=logs/soil_r28_cluster_histograms.out
#SBATCH --error=logs/soil_r28_cluster_histograms.err
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Generate histograms of CLR distributions for each cluster

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/28_cluster_histograms.R
