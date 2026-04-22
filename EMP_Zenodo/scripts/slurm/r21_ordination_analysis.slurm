#!/bin/bash
#SBATCH --job-name=emp_r21_ordination
#SBATCH --output=logs/emp_r21_ordination.out
#SBATCH --error=logs/emp_r21_ordination.err
#SBATCH --time=1:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# PCA/PCoA ordination analysis

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/21_ordination_analysis.R
