#!/bin/bash
#SBATCH --job-name=emp_r01b_soil_metadata
#SBATCH --output=logs/emp_r01b_soil_metadata.out
#SBATCH --error=logs/emp_r01b_soil_metadata.err
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Activate R environment
eval "$(/hopper/software/miniconda3/bin/conda shell.bash hook)"
conda activate r-dada2

# Run R script inline
Rscript - <<EOF
# Load full metadata
metadata_path <- "/hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/data/raw/emp_metadata_release1_20170912.tsv"
metadata <- read.delim(metadata_path, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

# Load soil sample IDs
soil_ids <- readLines("/hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/data/intermediate/soil_sample_ids.txt")

# Find sample ID column
id_col <- grep("sampleid", names(metadata), ignore.case = TRUE, value = TRUE)[1]

# Filter and write soil-only metadata
soil_metadata <- metadata[metadata[[id_col]] %in% soil_ids, ]
write.table(soil_metadata,
            "/hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/data/intermediate/soil_metadata.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Saved soil-only metadata to: data/intermediate/soil_metadata.tsv\n")
EOF
