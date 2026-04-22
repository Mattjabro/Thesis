#!/bin/bash
#SBATCH --job-name=emp_r04_qc_read_counts
#SBATCH --output=logs/emp_r04_qc_read_counts.out
#SBATCH --error=logs/emp_r04_qc_read_counts.err
#SBATCH --time=01:00:00
#SBATCH --mem=1028G
#SBATCH --cpus-per-task=128

module load R

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/04_qc_read_counts_vs_genus.R
