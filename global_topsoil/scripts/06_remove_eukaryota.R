#!/usr/bin/env Rscript

cat(" Loading required libraries...\n")
library(dada2)

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

# Set working directory
intermediate_dir <- file.path(project_root, "data", "intermediate")

# Load taxonomy table and ASV table
cat(" Reading in taxonomy and ASV tables...\n")
taxa <- readRDS(file.path(intermediate_dir, "taxa.rds"))
seqtab.nochim <- readRDS(file.path(intermediate_dir, "seqtab_nochim.rds"))

# Identify ASVs classified as Bacteria (not Eukaryota or NA in phylum)
cat(" Filtering out Eukaryota and unclassified ASVs...\n")
is_bacteria <- !is.na(taxa[,1]) & taxa[,1] == "Bacteria"
filtered_seqtab <- seqtab.nochim[, is_bacteria]
filtered_taxa <- taxa[is_bacteria, ]

# Confirm filtering
cat(" Remaining ASVs after filtering: ", ncol(filtered_seqtab), "\n")

# Save filtered outputs
cat(" Saving filtered ASV and taxonomy tables...\n")
saveRDS(filtered_seqtab, file.path(intermediate_dir, "seqtab_nochim_bacteria.rds"))
saveRDS(filtered_taxa, file.path(intermediate_dir, "taxa_bacteria.rds"))
write.csv(filtered_taxa, file.path(intermediate_dir, "taxonomy_table_bacteria.csv"))

cat(" Done! Outputs:\n")
cat(" - seqtab_nochim_bacteria.rds\n")
cat(" - taxa_bacteria.rds\n")
cat(" - taxonomy_table_bacteria.csv\n")
