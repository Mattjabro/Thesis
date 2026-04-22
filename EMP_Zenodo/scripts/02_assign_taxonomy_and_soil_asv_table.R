#!/usr/bin/env Rscript

timestamp <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  match <- grep(file_arg, cmd_args)
  if (length(match) > 0) {
    return(normalizePath(sub(file_arg, "", cmd_args[match[1]])))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  return(NA_character_)
}

script_path <- get_script_path()
project_root <- if (!is.na(script_path)) {
  normalizePath(file.path(dirname(script_path), ".."))
} else {
  normalizePath(getwd())
}

# =====================
# Load required libs
# =====================
timestamp("Loading packages...")
suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(dplyr)
  library(tidyr)
})

# ===========================
# Define paths
# ===========================
classifier_fp <- file.path(project_root, "data", "reference", "silva_nr99_v138_1_train_set.fasta.gz")
fasta_fp <- file.path(project_root, "data", "raw", "emp_deblur_150bp_min25.fasta")
biom_fp <- file.path(project_root, "data", "raw", "emp_deblur_150bp_release1.biom")
tsv_fp <- file.path(project_root, "data", "intermediate", "asv_table_all_samples.tsv")
taxonomy_out_fp <- file.path(project_root, "data", "intermediate", "taxonomy_all_samples.csv")
asv_table_out_fp <- file.path(project_root, "data", "intermediate", "asv_table_soil_only.csv")
soil_ids_fp <- file.path(project_root, "data", "intermediate", "soil_sample_ids.txt")
biom_exe <- Sys.getenv("BIOM_EXE", "~/.conda/envs/biom-env/bin/biom")

# =============================
# Convert BIOM to TSV if needed
# =============================
if (!file.exists(tsv_fp)) {
timestamp("Converting BIOM to TSV using biom CLI...")
  cmd <- sprintf("%s convert -i %s -o %s --to-tsv", biom_exe, biom_fp, tsv_fp)
  status <- system(cmd, intern = FALSE)
  if (status != 0) stop("Failed to convert biom to tsv.")
} else {
  timestamp("Found existing TSV table.")
}

# ============================
# Read biom table as matrix
# ============================
timestamp("Reading biom TSV...")
lines <- readLines(tsv_fp)
header_line <- grep("^#OTU ID", lines)
header <- strsplit(sub("^#OTU ID\t", "", lines[header_line]), "\t")[[1]]
raw_table <- lines[(header_line + 1):length(lines)]
parsed <- strsplit(raw_table, "\t")
otu_ids <- sapply(parsed, `[`, 1)
counts <- do.call(rbind, lapply(parsed, function(x) as.numeric(x[-1])))
colnames(counts) <- header
rownames(counts) <- otu_ids
biom_table <- as.data.frame(counts)

# =============================
# Filter to soil sample IDs
# =============================
timestamp("Filtering to soil sample IDs...")
soil_ids <- readLines(soil_ids_fp)
soil_ids <- soil_ids[soil_ids %in% colnames(biom_table)]
biom_table <- biom_table[, soil_ids, drop = FALSE]
timestamp(sprintf("Retained %d soil samples.", length(soil_ids)))

# ============================
# Assign taxonomy
# ============================
timestamp("Reading representative sequences...")
rep_seqs <- readDNAStringSet(fasta_fp)
asv_seqs <- as.character(rep_seqs)
names(asv_seqs) <- names(rep_seqs)

timestamp("Assigning taxonomy...")
taxa <- assignTaxonomy(asv_seqs, classifier_fp, multithread = TRUE)
taxa_df <- as.data.frame(taxa)
taxa_df$ASV <- rownames(taxa_df)
taxa_df <- taxa_df %>% relocate(ASV)
write.csv(taxa_df, taxonomy_out_fp, row.names = FALSE)
timestamp(paste("Saved taxonomy to", taxonomy_out_fp))

# ===============================
# Save filtered ASV table
# ===============================
matched_asvs <- intersect(rownames(biom_table), taxa_df$ASV)
biom_table <- biom_table[matched_asvs, , drop = FALSE]
write.csv(biom_table, asv_table_out_fp)
timestamp(paste("Saved ASV table to", asv_table_out_fp))
