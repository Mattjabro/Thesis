#!/bin/bash
#SBATCH --job-name=soil_r08_classify
#SBATCH --output=logs/soil_r08_classify.out
#SBATCH --error=logs/soil_r08_classify.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=1028G

echo " SLURM job started on $(date)"

# Load Conda and activate your R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run your R script
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/08_classify.R

echo " Job completed on $(date)"
