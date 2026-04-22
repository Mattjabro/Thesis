#!/bin/bash
#SBATCH --job-name=emp_r23_cluster_histograms
#SBATCH --output=logs/emp_r23_cluster_histograms.out
#SBATCH --error=logs/emp_r23_cluster_histograms.err
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Generate histograms of CLR distributions for each cluster

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/23_cluster_histograms.R
