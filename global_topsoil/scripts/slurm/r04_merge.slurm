#!/bin/bash
#SBATCH --job-name=soil_r04_merge
#SBATCH --output=logs/soil_r04_merge.out
#SBATCH --error=logs/soil_r04_merge.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=1028G

# Activate conda environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run R script
echo " Running Step 4: mergePairs()"
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/04_merge.R

# Confirm output
MERGED="/hopper/groups/yoosephlab/datasets/global_topsoil_matthew/data/intermediate/mergers.rds"
if [[ -f "$MERGED" ]]; then
    echo " Created: mergers.rds"
else
    echo " ERROR: mergers.rds was not created."
fi
