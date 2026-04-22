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

timestamp("Loading EMP matrices...")
taxa_fp <- file.path(project_root, "data", "intermediate", "asv_taxonomy_matrix_all_samples.csv")
genus_fp <- file.path(project_root, "data", "intermediate", "genus_matrix_all_samples.csv")
comparison_fp <- file.path(project_root, "outputs", "qc_tables", "all_samples_read_count_comparison.csv")

taxa_df <- read.csv(taxa_fp, header = TRUE, row.names = 1, check.names = FALSE)
genus_df <- read.csv(genus_fp, header = TRUE, row.names = 1, check.names = FALSE)
comparison <- read.csv(comparison_fp, header = TRUE)

timestamp("Dimensions before cleaning:")
cat("Taxa matrix:", dim(taxa_df)[1], "rows,", dim(taxa_df)[2], "cols\n")

# Convert all columns to numeric and count NAs
convert_safe <- function(x) {
  out <- suppressWarnings(as.numeric(as.character(x)))
  if (any(is.na(out))) warning("NAs introduced in column during coercion.")
  return(out)
}

taxa_df[] <- lapply(taxa_df, convert_safe)
genus_df[] <- lapply(genus_df, convert_safe)

# Drop all-NA rows or columns
taxa_df <- taxa_df[rowSums(is.na(taxa_df)) != ncol(taxa_df), ]
taxa_df <- taxa_df[, colSums(is.na(taxa_df)) != nrow(taxa_df)]

genus_df <- genus_df[rowSums(is.na(genus_df)) != ncol(genus_df), ]
genus_df <- genus_df[, colSums(is.na(genus_df)) != nrow(genus_df)]

timestamp("Dimensions after cleaning:")
cat("Taxa matrix:", dim(taxa_df)[1], "rows,", dim(taxa_df)[2], "cols\n")

taxa_matrix <- as.matrix(taxa_df)
genus_matrix <- as.matrix(genus_df)

summary_stats <- function(x) {
  c(Mean = mean(x), Median = median(x), Min = min(x), Max = max(x), SD = sd(x))
}

if (ncol(taxa_matrix) > 0) {
  col_sums <- colSums(taxa_matrix, na.rm = TRUE)
  write.table(data.frame(Sum = col_sums),
              file.path(project_root, "outputs", "qc_tables", "all_samples_asv_reads_per_sample.tsv"),
              sep = "\t", quote = FALSE)

  png(file.path(project_root, "outputs", "figures", "all_samples_asv_reads_per_sample_hist.png"))
  hist(col_sums, breaks = 50, main = "ASV Reads per Sample (All Samples)", xlab = "Reads")
  dev.off()

  write.table(summary_stats(col_sums),
              file.path(project_root, "outputs", "qc_tables", "all_samples_asv_reads_per_sample_summary.tsv"),
              sep = "\t")
} else {
  warning("No columns in taxa matrix after cleaning.")
}

if (nrow(taxa_matrix) > 0) {
  row_nonzero <- rowSums(taxa_matrix != 0, na.rm = TRUE)
  write.table(data.frame(Nonzero = row_nonzero),
              file.path(project_root, "outputs", "qc_tables", "all_samples_asv_nonzero_samples.tsv"),
              sep = "\t", quote = FALSE)

  png(file.path(project_root, "outputs", "figures", "all_samples_asv_nonzero_samples_hist.png"))
  hist(row_nonzero, breaks = 50, main = "Nonzero Samples per ASV (All Samples)", xlab = "Samples")
  dev.off()

  write.table(summary_stats(row_nonzero),
              file.path(project_root, "outputs", "qc_tables", "all_samples_asv_nonzero_samples_summary.tsv"),
              sep = "\t")
} else {
  warning("No rows in taxa matrix after cleaning.")
}

if (ncol(comparison) >= 4) {
  col_d <- comparison[[4]]
  png(file.path(project_root, "outputs", "figures", "all_samples_genus_retention_percent_hist.png"))
  hist(col_d, breaks = 50, main = "Reads Retained After Genus Collapse", xlab = "Percent Retained")
  dev.off()
  write.table(summary_stats(col_d),
              file.path(project_root, "outputs", "qc_tables", "all_samples_genus_retention_percent_summary.tsv"),
              sep = "\t")
}

if (ncol(genus_matrix) > 0) {
  genus_sums <- colSums(genus_matrix, na.rm = TRUE)
  write.table(data.frame(Sum = genus_sums),
              file.path(project_root, "outputs", "qc_tables", "all_samples_genus_reads_per_sample.tsv"),
              sep = "\t", quote = FALSE)

  png(file.path(project_root, "outputs", "figures", "all_samples_genus_reads_per_sample_hist.png"))
  hist(genus_sums, breaks = 50, main = "Genus Reads per Sample (All Samples)", xlab = "Reads")
  dev.off()

  write.table(summary_stats(genus_sums),
              file.path(project_root, "outputs", "qc_tables", "all_samples_genus_reads_per_sample_summary.tsv"),
              sep = "\t")
}

timestamp("All-sample stats complete.")
