#!/bin/bash
#SBATCH --job-name=soil_r00_filter_classifier
#SBATCH --output=logs/soil_r00_filter_classifier.out
#SBATCH --error=logs/soil_r00_filter_classifier.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Activate R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run the R script to filter the classifier FASTA
echo " Filtering SILVA classifier to remove Eukaryota..."
Rscript /hopper/groups/yoosephlab/datasets/global_topsoil_matthew/scripts/00_filter_classifier.R

# Confirm output
OUT="/hopper/groups/yoosephlab/datasets/global_topsoil_matthew/data/reference/silva_nr99_v138.2_toGenus_trainset_filtered.fa"
if [[ -f "$OUT" ]]; then
    echo -e "\n Filtered classifier written to: $OUT"
else
    echo -e "\n ERROR: Filtered classifier not found. Check your R script."
fi
