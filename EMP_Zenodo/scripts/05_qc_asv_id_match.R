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

# ================
# Load libraries
# ================
timestamp("Loading packages...")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# ================
# File paths
# ================
tsv_fp <- file.path(project_root, "data", "intermediate", "asv_table_all_samples.tsv")
taxa_fp <- file.path(project_root, "data", "intermediate", "asv_taxonomy_matrix_all_samples.csv")

out_missing <- file.path(project_root, "outputs", "qc_tables", "asvs_missing_from_asv_table.csv")
out_extra <- file.path(project_root, "outputs", "qc_tables", "asvs_extra_in_asv_table.csv")
out_summary <- file.path(project_root, "outputs", "qc_tables", "asv_id_match_summary.txt")

# ============================
# Read ASV IDs from ASV table
# ============================
timestamp("Reading ASV table...")
lines <- readLines(tsv_fp)
header_line <- grep("^#OTU ID", lines)
raw_lines <- lines[(header_line + 1):length(lines)]
parsed <- strsplit(raw_lines, "\t")
asv_ids_emp <- sapply(parsed, `[`, 1)

# ============================
# Read ASV IDs from ASV taxonomy matrix
# ============================
timestamp("Reading ASV taxonomy matrix...")
taxa_matrix <- read_csv(taxa_fp, show_col_types = FALSE)
asv_ids_matrix <- taxa_matrix$ASV

# ============================
# Compare ASV sets
# ============================
timestamp("Comparing ASV IDs...")
missing_from_emp <- setdiff(asv_ids_matrix, asv_ids_emp)
extra_in_emp <- setdiff(asv_ids_emp, asv_ids_matrix)

write_csv(data.frame(ASV = missing_from_emp), out_missing)
write_csv(data.frame(ASV = extra_in_emp), out_extra)

summary_lines <- c(
  sprintf("Total ASVs in ASV table: %d", length(asv_ids_emp)),
  sprintf("Total ASVs in ASV taxonomy matrix: %d", length(asv_ids_matrix)),
  sprintf("ASVs in matrix but NOT in ASV table: %d", length(missing_from_emp)),
  sprintf("ASVs in ASV table but NOT in matrix: %d", length(extra_in_emp))
)
writeLines(summary_lines, out_summary)
timestamp("Comparison complete.")
