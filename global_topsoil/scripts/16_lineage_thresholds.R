library(tidyverse)

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

processed_dir <- file.path(project_root, "data", "processed")
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

# Load full lineage counts
lineage_counts <- read_csv(file.path(processed_dir, "taxonomy_full_lineage_counts.csv"))

# Total ASVs before filtering
total_asvs <- sum(lineage_counts$Count)
total_lineages <- nrow(lineage_counts)

# Initialize results
results <- data.frame(
  Threshold = integer(),
  LineagesRemaining = integer(),
  FractionRemaining = numeric(),
  ASVsLost = integer()
)

# Loop through thresholds from 1 to 20
for (t in 1:20) {
  # Keep only lineages that appear in at least t ASVs
  filtered <- lineage_counts %>% filter(Count >= t)

  lineages_remaining <- nrow(filtered)
  fraction_remaining <- lineages_remaining / total_lineages
  asvs_lost <- total_asvs - sum(filtered$Count)

  results <- results %>%
    add_row(
      Threshold = t,
      LineagesRemaining = lineages_remaining,
      FractionRemaining = fraction_remaining,
      ASVsLost = asvs_lost
    )
}

# Save to CSV
write_csv(results, file.path(processed_dir, "soil_fraction_taxa_remaining_1_to_20.csv"))
print(results)
