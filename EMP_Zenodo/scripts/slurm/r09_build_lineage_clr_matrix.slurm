#!/bin/bash
#SBATCH --job-name=emp_r09_lineage_clr
#SBATCH --output=logs/emp_r09_lineage_clr.out
#SBATCH --error=logs/emp_r09_lineage_clr.err
#SBATCH --time=11:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem=1028G

# Activate appropriate R environment (avoid qiime2 env)
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2  #  your working R environment

# Run the R script
Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/09_build_lineage_clr_matrix.R
