#!/bin/bash
#SBATCH --job-name=soil_r01_filter_trim
#SBATCH --output=logs/soil_r01_filter_trim.out
#SBATCH --error=logs/soil_r01_filter_trim.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024G

# Activate your working R + DADA2 env
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run the filtering R script
echo " Running Step 1: Filter and Trim..."
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/01_filter_and_trim.R

# Check output file before printing summary
RDS_PATH="/hopper/groups/yoosephlab/datasets/global_topsoil_matthew/data/intermediate/filtering_summary.rds"
if [[ -f "$RDS_PATH" ]]; then
    echo -e "\n FILTERING SUMMARY:"
    Rscript -e "out <- readRDS('$RDS_PATH'); print(out)"
else
    echo -e "\n ERROR: filtering_summary.rds was not created."
    echo "Check your R script for silent errors."
fi
