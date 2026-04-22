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
  library(dplyr)
  library(tidyr)
  library(readr)
})

# =====================================
# Define input/output file paths
# =====================================
taxonomy_fp <- file.path(project_root, "data", "intermediate", "taxonomy_all_samples.csv")
tsv_fp <- file.path(project_root, "data", "intermediate", "asv_table_all_samples.tsv")
output_fp <- file.path(project_root, "data", "intermediate", "genus_matrix_all_samples.csv")

# =====================================
# Read taxonomy table and ASV table
# =====================================
timestamp("Reading taxonomy table and ASV table...")

taxonomy <- read_csv(taxonomy_fp, show_col_types = FALSE)

# Read and clean .tsv table
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

# ===============================
# Merge and match on ASV IDs
# ===============================
timestamp("Merging taxonomy with abundance table...")

if (!"ASV" %in% colnames(taxonomy)) {
  stop("'ASV' column not found in taxonomy table.")
}
matched_asvs <- intersect(rownames(biom_table), taxonomy$ASV)

if (length(matched_asvs) == 0) {
  stop("No matched ASVs between taxonomy and abundance table.")
}

taxa_filtered <- taxonomy[taxonomy$ASV %in% matched_asvs, ]
biom_filtered <- biom_table[matched_asvs, ]

# Ensure the order matches
biom_filtered <- biom_filtered[taxa_filtered$ASV, ]

merged <- cbind(taxa_filtered, biom_filtered)

# ===============================
# Collapse to genus level
# ===============================
timestamp("Collapsing table to genus level...")

if (!"Genus" %in% colnames(merged)) {
  stop("'Genus' column missing in taxonomy table.")
}

genus_table <- merged %>%
  select(Genus, where(is.numeric)) %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum), .groups = "drop")

# ===============================
# Save result
# ===============================
write_csv(genus_table, output_fp)
timestamp(paste("Saved genus-level matrix to", output_fp))
