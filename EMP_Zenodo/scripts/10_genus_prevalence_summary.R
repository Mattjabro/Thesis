#!/usr/bin/env Rscript

# === Load libraries ===
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

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

# === File paths ===
input_file <- file.path(project_root, "data", "intermediate", "genus_matrix_soil_only.csv")
output_summary <- file.path(project_root, "data", "processed", "genus_prevalence_summary.csv")
output_histogram <- file.path(project_root, "outputs", "figures", "genus_prevalence_histogram.png")

# === 1. Load the sample  genus matrix ===
cat("Reading matrix...\n")
matrix <- read_csv(input_file, show_col_types = FALSE)

# === 2. Binarize: make all nonzero counts into 1 ===
cat("Binarizing counts (nonzero -> 1)...\n")
binary_matrix <- matrix %>%
  select(-SampleID) %>%
  mutate(across(everything(), ~ as.integer(.x > 0)))

# === 3. Compute prevalence per genus ===
cat("Calculating prevalence per genus...\n")
prevalence <- colSums(binary_matrix)
prevalence_df <- tibble(
  Genus = names(prevalence),
  NumSamplesPresent = as.integer(prevalence),
  PrevalencePercent = round(100 * prevalence / nrow(matrix), 2)
)

# === 4. Write prevalence summary CSV ===
cat("Writing summary CSV...\n")
write_csv(prevalence_df, output_summary)

# === 5. Generate histogram ===
cat("Generating histogram...\n")
hist_data <- prevalence_df %>%
  count(NumSamplesPresent, name = "NumGenera")

ggplot(hist_data, aes(x = NumSamplesPresent, y = NumGenera)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Genus Prevalence Across Samples",
    x = "Number of Samples (with Genus Present)",
    y = "Number of Genera"
  ) +
  theme_minimal()

ggsave(output_histogram, width = 8, height = 5)

cat("Done.\n")
cat("- Summary CSV:", output_summary, "\n")
cat("- Histogram PNG:", output_histogram, "\n")
