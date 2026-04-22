#!/bin/bash
#SBATCH --job-name=emp_r24_biome_grouping
#SBATCH --output=logs/emp_r24_biome_grouping.out
#SBATCH --error=logs/emp_r24_biome_grouping.err
#SBATCH --time=0:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Create ecologically meaningful biome groupings

eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/24_biome_grouping.R
