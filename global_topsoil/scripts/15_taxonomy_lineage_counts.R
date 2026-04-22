# Load required library
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

# Load taxonomy table
taxonomy <- read_csv(file.path(intermediate_dir, "taxonomy_table_bacteria.csv"))

# Define taxonomic columns to use
taxonomic_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Unite all taxonomy levels into one lineage string
taxonomy_lineages <- taxonomy %>%
  unite("FullLineage", all_of(taxonomic_ranks), sep = ";", remove = FALSE) %>%
  count(FullLineage, name = "Count") %>%
  arrange(desc(Count))

# Save to CSV
write_csv(taxonomy_lineages, file.path(processed_dir, "taxonomy_full_lineage_counts.csv"))
