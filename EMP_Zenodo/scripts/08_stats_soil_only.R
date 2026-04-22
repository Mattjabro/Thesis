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

# ======================
# Load required libs
# ======================
timestamp("Loading packages...")
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# =============================
# Define file paths
# =============================
asv_fp <- file.path(project_root, "data", "intermediate", "asv_table_soil_only.csv")
taxonomy_fp <- file.path(project_root, "data", "intermediate", "taxonomy_all_samples.csv")
soil_meta_fp <- file.path(project_root, "data", "intermediate", "soil_metadata.tsv")
comparison_fp <- file.path(project_root, "outputs", "qc_tables", "all_samples_read_count_comparison.csv")

# =============================
# Load data
# =============================
timestamp("Reading matrices and metadata...")
asv_df <- read_csv(asv_fp, show_col_types = FALSE)
asv_df <- as.data.frame(asv_df)
rownames(asv_df) <- asv_df[[1]]
asv_df <- asv_df[, -1]

taxonomy <- read_csv(taxonomy_fp, show_col_types = FALSE)
soil_samples <- read_tsv(soil_meta_fp, col_names = FALSE, show_col_types = FALSE)[[1]]
comparison <- read_csv(comparison_fp, show_col_types = FALSE)

# =============================
# Validate and filter soil samples
# =============================
valid_soil_samples <- intersect(soil_samples, colnames(asv_df))
missing_soil_samples <- setdiff(soil_samples, colnames(asv_df))

if (length(missing_soil_samples) > 0) {
  warning(sprintf("%d soil samples not found in ASV matrix:\n%s",
                  length(missing_soil_samples),
                  paste(missing_soil_samples, collapse = ", ")))
}

# =============================
# Filter taxonomy to soil-only ASVs
# =============================
timestamp("Filtering taxonomy table to ASVs found in soil samples...")
soil_asv_ids <- rownames(asv_df)[rowSums(asv_df[, valid_soil_samples, drop = FALSE]) > 0]
taxonomy_soil <- taxonomy %>% filter(ASV %in% soil_asv_ids)
write_csv(taxonomy_soil, file.path(project_root, "data", "intermediate", "taxonomy_soil_only.csv"))

# =============================
# Histogram: ASVs per soil sample
# =============================
timestamp("Plotting histogram of ASV richness (non-zero per sample)...")
asv_soil <- asv_df[, valid_soil_samples, drop = FALSE]
asv_counts <- colSums(asv_soil > 0)

png(file.path(project_root, "outputs", "figures", "soil_asv_nonzero_per_sample_hist.png"), width = 1000, height = 700)
hist(asv_counts,
     breaks = 30,
     col = "skyblue",
     border = "white",
     main = "Number of ASVs Detected per Soil Sample",
     xlab = "Non-Zero ASVs per Sample",
     ylab = "Sample Count")
dev.off()

write.table(data.frame(Sample = names(asv_counts), NonzeroASVs = asv_counts),
            file.path(project_root, "outputs", "qc_tables", "soil_asv_nonzero_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(c(Mean = mean(asv_counts),
              Median = median(asv_counts),
              Min = min(asv_counts),
              Max = max(asv_counts),
              SD = sd(asv_counts)),
            file.path(project_root, "outputs", "qc_tables", "soil_asv_nonzero_per_sample_summary.tsv"),
            sep = "\t", col.names = FALSE)

# =============================
# Histogram: Row sums (read depth per sample)
# =============================
timestamp("Plotting read depth per sample...")
row_sums <- colSums(asv_soil)

write.table(data.frame(Sample = names(row_sums), TotalReads = row_sums),
            file.path(project_root, "outputs", "qc_tables", "soil_reads_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(c(Mean = mean(row_sums),
              Median = median(row_sums),
              Min = min(row_sums),
              Max = max(row_sums),
              SD = sd(row_sums)),
            file.path(project_root, "outputs", "qc_tables", "soil_reads_per_sample_summary.tsv"),
            sep = "\t", col.names = FALSE)

png(file.path(project_root, "outputs", "figures", "soil_reads_per_sample_hist.png"), width = 1000, height = 700)
hist(row_sums,
     breaks = 30,
     col = "seagreen",
     border = "white",
     main = "Total Reads per Soil Sample",
     xlab = "Total Reads",
     ylab = "Sample Count")
dev.off()

# =============================
# Histogram: Column D from comparison file
# =============================
if (ncol(comparison) >= 4) {
  col_d <- comparison[[4]]

  png(file.path(project_root, "outputs", "figures", "soil_genus_retention_percent_hist.png"), width = 1000, height = 700)
  hist(col_d,
       breaks = 50,
       col = "coral",
       border = "white",
       main = "% Reads Retained (Simulated)",
       xlab = "% Reads Retained",
       ylab = "Sample Count")
  dev.off()

  write.table(c(Mean = mean(col_d),
                Median = median(col_d),
                Min = min(col_d),
                Max = max(col_d),
                SD = sd(col_d)),
              file.path(project_root, "outputs", "qc_tables", "soil_genus_retention_percent_summary.tsv"),
              sep = "\t", col.names = FALSE)
}

timestamp("Soil-only stats complete.")
