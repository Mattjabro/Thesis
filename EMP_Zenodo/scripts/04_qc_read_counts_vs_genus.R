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

# Load libraries
timestamp("Loading packages...")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Define paths
genus_fp <- file.path(project_root, "data", "intermediate", "genus_matrix_all_samples.csv")
asv_fp <- file.path(project_root, "data", "intermediate", "asv_table_all_samples.tsv")
taxonomy_fp <- file.path(project_root, "data", "intermediate", "taxonomy_all_samples.csv")
out_comp_fp <- file.path(project_root, "outputs", "qc_tables", "all_samples_read_count_comparison.csv")
out_unclass_fp <- file.path(project_root, "outputs", "qc_tables", "unclassified_reads_per_sample.csv")

# Read genus table
timestamp("Reading genus matrix...")
genus_table <- read_csv(genus_fp, show_col_types = FALSE)
genus_matrix <- as.data.frame(genus_table)
rownames(genus_matrix) <- genus_matrix$Genus
genus_matrix$Genus <- NULL
genus_matrix <- as.data.frame(t(genus_matrix))
genus_sums <- rowSums(genus_matrix)

# Read ASV table
timestamp("Reading ASV table...")
lines <- readLines(asv_fp)
header_line <- grep("^#OTU ID", lines)
header <- strsplit(sub("^#OTU ID\t", "", lines[header_line]), "\t")[[1]]
raw_table <- lines[(header_line + 1):length(lines)]
parsed <- strsplit(raw_table, "\t")
otu_ids <- sapply(parsed, `[`, 1)

timestamp("Parsing ASV matrix...")
num_samples <- length(header)
num_asvs <- length(parsed)
counts <- matrix(0, nrow = num_asvs, ncol = num_samples)
for (i in seq_len(num_asvs)) {
  counts[i, ] <- as.numeric(parsed[[i]][-1])
}
colnames(counts) <- header
rownames(counts) <- otu_ids
asv_matrix <- as.data.frame(counts)
rm(counts)  # free memory
asv_sums <- colSums(asv_matrix)

# Compare totals
shared_samples <- intersect(names(genus_sums), names(asv_sums))
comparison_df <- data.frame(
  SampleID = shared_samples,
  ASV_Reads = asv_sums[shared_samples],
  Genus_Reads = genus_sums[shared_samples]
)
comparison_df$Percent_Retained <- 100 * comparison_df$Genus_Reads / comparison_df$ASV_Reads
write_csv(comparison_df, out_comp_fp)
timestamp(paste("Saved comparison table to", out_comp_fp))

# Optional: check unclassified reads
if (file.exists(taxonomy_fp)) {
  timestamp("Checking unclassified reads...")

  taxonomy <- read_csv(taxonomy_fp, show_col_types = FALSE)
  tax_cols <- colnames(taxonomy)

  if ("Taxon" %in% tax_cols && "ASV" %in% tax_cols) {
    unclass_asvs <- taxonomy %>%
      filter(grepl("Unassigned|unclassified", Taxon, ignore.case = TRUE)) %>%
      pull(ASV)

    shared_asvs <- intersect(unclass_asvs, rownames(asv_matrix))
    if (length(shared_asvs) > 0) {
      unclassified_counts <- colSums(asv_matrix[shared_asvs, , drop = FALSE])
      unclass_df <- data.frame(
        SampleID = names(unclassified_counts),
        Unclassified_Reads = unclassified_counts
      )
      write_csv(unclass_df, out_unclass_fp)
      timestamp(paste("Saved unclassified counts to", out_unclass_fp))
    } else {
      timestamp("No matching unclassified ASVs found in ASV table.")
    }
  } else {
    timestamp("'Taxon' or 'ASV' column not found in taxonomy table; skipping unclassified check.")
  }
} else {
  timestamp("No taxonomy file found; skipping unclassified check.")
}
