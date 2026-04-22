# fraction_taxa_remaining_bacteria.R

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

intermediate_dir <- file.path(project_root, "data", "intermediate")
processed_dir <- file.path(project_root, "data", "processed")
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

# Load the taxonomy table
taxonomy <- read_csv(file.path(intermediate_dir, "taxonomy_table_bacteria.csv"))

# Remove duplicate ASVs with identical taxonomy profile
taxonomy_unique <- taxonomy %>% distinct()

# Turn nonzero cells into 1 (presence/absence)
presence_matrix <- taxonomy_unique %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))

# Initialize storage for output
results <- data.frame(threshold = 1:20,
                      taxa_remaining = NA,
                      fraction_lost = NA)

total_taxa <- ncol(presence_matrix)

# For each threshold from 1 to 20
for (t in 1:20) {
  # Count how many samples each taxon appears in
  taxa_prevalence <- colSums(presence_matrix)

  # Keep only taxa that appear in at least t samples
  kept_taxa <- sum(taxa_prevalence >= t)

  # Compute fraction lost
  frac_lost <- 1 - (kept_taxa / total_taxa)

  # Store results
  results[t, "taxa_remaining"] <- kept_taxa
  results[t, "fraction_lost"] <- frac_lost
}

# Save results
write_csv(results, file.path(processed_dir, "fraction_taxa_remaining_bacteria.csv"))
