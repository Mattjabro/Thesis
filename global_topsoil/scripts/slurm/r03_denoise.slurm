#!/bin/bash
#SBATCH --job-name=soil_r03_denoise
#SBATCH --output=logs/soil_r03_denoise.out
#SBATCH --error=logs/soil_r03_denoise.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=1028G

# Activate conda environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run R script
echo " Running Step 3: Denoising with dada()"
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/03_denoise.R

# Check for outputs
for f in dadaFs.rds dadaRs.rds; do
  path="/hopper/groups/yoosephlab/datasets/global_topsoil_matthew/data/intermediate/$f"
  if [[ -f "$path" ]]; then
    echo " Created: $f"
  else
    echo " Missing: $f"
  fi
done
