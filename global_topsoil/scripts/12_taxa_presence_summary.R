#!/usr/bin/env Rscript

# === Load libraries ===
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# Resolve project paths
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

intermediate_dir <- file.path(project_root, "data", "intermediate")
processed_dir <- file.path(project_root, "data", "processed")
fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# === File paths ===
input_file <- file.path(intermediate_dir, "sample_taxa_matrix_L6_bacteria.csv")
output_summary <- file.path(processed_dir, "taxa_presence_summary_soil.csv")
output_histogram <- file.path(fig_dir, "taxa_presence_histogram_soil.png")

# === 1. Read data with first column as rownames ===
cat(" Reading matrix with taxa as columns and samples as rows...\n")
matrix <- read_csv(input_file, col_names = TRUE, show_col_types = FALSE)

# Check if first column is actually sample ID
if (!grepl("ERR", matrix[[1]][1])) {
  stop(" First column does not look like sample IDs. Check format.")
}

# Set sample IDs as rownames and remove from data
sample_ids <- matrix[[1]]
matrix <- matrix[, -1]
rownames(matrix) <- sample_ids

# === 2. Binarize: nonzero = 1 ===
cat(" Binarizing matrix...\n")
binary_matrix <- as.data.frame((matrix > 0) * 1)

# === 3. Compute per-taxa prevalence ===
cat(" Summarizing taxa prevalence...\n")
prevalence <- colSums(binary_matrix)
prevalence_df <- tibble(
  Taxa = names(prevalence),
  NumSamplesPresent = as.integer(prevalence),
  PrevalencePercent = round(100 * prevalence / nrow(binary_matrix), 2)
)

# === 4. Save summary table ===
write_csv(prevalence_df, output_summary)

# === 5. Histogram of how many taxa appear in how many samples ===
hist_data <- prevalence_df %>%
  count(NumSamplesPresent, name = "NumTaxa")

ggplot(hist_data, aes(x = NumSamplesPresent, y = NumTaxa)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Taxa Prevalence Across Samples",
    x = "# Samples Where Taxa Appears",
    y = "# Taxa"
  ) +
  theme_minimal()

ggsave(output_histogram, width = 8, height = 5)

cat(" Done!\n")
cat("- Summary CSV: ", output_summary, "\n")
cat("- Histogram PNG: ", output_histogram, "\n")
