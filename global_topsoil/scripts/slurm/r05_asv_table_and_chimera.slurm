#!/bin/bash
#SBATCH --job-name=soil_r05_asv_chimera
#SBATCH --output=logs/soil_r05_asv_chimera.out
#SBATCH --error=logs/soil_r05_asv_chimera.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=1028G

# Activate conda env
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run R script
echo " Running Step 5: ASV table + chimera removal..."
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/05_asv_table_and_chimera.R

# Confirm output
for f in seqtab.rds seqtab_nochim.rds; do
  path="/hopper/groups/yoosephlab/datasets/global_topsoil_matthew/data/intermediate/$f"
  if [[ -f "$path" ]]; then
    echo " Created: $f"
  else
    echo " Missing: $f"
  fi
done
