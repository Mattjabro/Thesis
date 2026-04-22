#!/bin/bash
#SBATCH --job-name=soil_r02_learn_errors
#SBATCH --output=logs/soil_r02_learn_errors.out
#SBATCH --error=logs/soil_r02_learn_errors.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=1028G

# Activate working R+DADA2 environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run R script
echo " Running Step 2: Learn Error Rates..."
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/02_learn_errors.R

# Show brief confirmation if files were created
for file in errF.rds errR.rds; do
  path="/hopper/groups/yoosephlab/datasets/global_topsoil_matthew/data/intermediate/$file"
  if [[ -f "$path" ]]; then
    echo " Created: $file"
  else
    echo " Missing: $file"
  fi
done
