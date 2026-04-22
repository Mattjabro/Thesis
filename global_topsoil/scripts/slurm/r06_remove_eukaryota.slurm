#!/bin/bash
#SBATCH --job-name=soil_r06_filter_euk
#SBATCH --output=logs/soil_r06_filter_euk.out
#SBATCH --error=logs/soil_r06_filter_euk.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Activate R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run filtering R script
echo " Filtering out Eukaryota ASVs..."
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/06_remove_eukaryota.R

# Confirm output
OUT="/hopper/groups/yoosephlab/datasets/global_topsoil_matthew/data/intermediate/seqtab_nochim_bacteria.rds"
if [[ -f "$OUT" ]]; then
    echo -e "\n Eukaryota-filtered ASV table saved to: $OUT"
else
    echo -e "\n ERROR: Filtered ASV table not found. Check your R script."
fi
