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

# ===============================
# Define input/output files
# ===============================
taxonomy_fp <- file.path(project_root, "data", "intermediate", "taxonomy_all_samples.csv")
tsv_fp <- file.path(project_root, "data", "intermediate", "asv_table_all_samples.tsv")
soil_meta_fp <- file.path(project_root, "data", "intermediate", "soil_metadata.tsv")
output_matrix_fp <- file.path(project_root, "data", "intermediate", "genus_matrix_soil_only.csv")
output_taxonomy_fp <- file.path(project_root, "data", "intermediate", "taxonomy_soil_only.csv")

# ===============================
# Read input files
# ===============================
timestamp("Reading taxonomy table, ASV table, and soil sample list...")

taxonomy <- read_csv(taxonomy_fp, show_col_types = FALSE)
soil_samples <- read_tsv(soil_meta_fp, col_names = FALSE, show_col_types = FALSE)[[1]]

# === Read and clean .tsv table safely ===
lines <- readLines(tsv_fp, warn = FALSE)

# Drop final line if it's malformed (missing tab character)
if (!grepl("\t", tail(lines, 1))) {
  lines <- lines[-length(lines)]
}

header_line <- grep("^#OTU ID", lines)
header <- strsplit(sub("^#OTU ID\t", "", lines[header_line]), "\t")[[1]]
raw_table <- lines[(header_line + 1):length(lines)]

# Split and filter to well-formed rows
parsed <- strsplit(raw_table, "\t")
parsed <- parsed[sapply(parsed, length) == (length(header) + 1)]

# Convert to matrix
otu_ids <- sapply(parsed, `[`, 1)
counts <- do.call(rbind, lapply(parsed, function(x) as.numeric(x[-1])))
colnames(counts) <- header
rownames(counts) <- otu_ids
biom_table <- as.data.frame(counts)

# ===============================
# Merge taxonomy with counts
# ===============================
timestamp("Merging taxonomy with abundance table...")

matched_asvs <- intersect(rownames(biom_table), taxonomy$ASV)
if (length(matched_asvs) == 0) stop("No matched ASVs between taxonomy and abundance table.")

taxonomy_matched <- taxonomy[taxonomy$ASV %in% matched_asvs, ]
biom_matched <- biom_table[matched_asvs, ]
biom_matched <- biom_matched[taxonomy_matched$ASV, ]  # ensure row order matches taxonomy

merged <- cbind(taxonomy_matched, biom_matched)

# ===============================
# Collapse to genus level
# ===============================
timestamp("Collapsing to genus level...")

if (!"Genus" %in% colnames(merged)) stop("'Genus' column missing in taxonomy table.")

genus_table <- merged %>%
  select(Genus, where(is.numeric)) %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum), .groups = "drop")

# ===============================
# Filter to soil samples with zero-fill
# ===============================
timestamp("Filtering for soil samples (adding zero rows if needed)...")

# Pivot and filter
soil_genus_matrix <- genus_table %>%
  select(Genus, all_of(intersect(soil_samples, colnames(genus_table)))) %>%
  pivot_longer(-Genus, names_to = "SampleID", values_to = "Count") %>%
  pivot_wider(names_from = Genus, values_from = Count, values_fill = 0)

# Add missing samples (zero rows)
missing_samples <- setdiff(soil_samples, soil_genus_matrix$SampleID)
if (length(missing_samples) > 0) {
  zero_row <- as_tibble(matrix(0, nrow = length(missing_samples), ncol = ncol(soil_genus_matrix) - 1))
  colnames(zero_row) <- colnames(soil_genus_matrix)[-1]
  zero_row <- bind_cols(SampleID = missing_samples, zero_row)
  soil_genus_matrix <- bind_rows(soil_genus_matrix, zero_row)
}

# Reorder to match original metadata order
soil_genus_matrix <- soil_genus_matrix %>%
  arrange(match(SampleID, soil_samples))

# ===============================
# Save outputs
# ===============================
write_csv(soil_genus_matrix, output_matrix_fp)
timestamp(paste("Saved sample  genus matrix to", output_matrix_fp))

taxonomy_soil_only <- taxonomy %>%
  filter(ASV %in% matched_asvs)
write_csv(taxonomy_soil_only, output_taxonomy_fp)
timestamp(paste("Saved filtered taxonomy table to", output_taxonomy_fp))
