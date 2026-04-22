#!/bin/bash
#SBATCH --job-name=emp_r02_assign_taxonomy
#SBATCH --output=logs/emp_r02_assign_taxonomy.out
#SBATCH --error=logs/emp_r02_assign_taxonomy.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=1028G

# Activate Conda environment with DADA2
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run the R script
echo "Running R pipeline: assign taxonomy + soil ASV table..."
Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/02_assign_taxonomy_and_soil_asv_table.R

# Final message
echo "R script completed."
