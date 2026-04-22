library(tidyverse)

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

# === Step 1: Load the correct taxonomy table ===
taxonomy <- read_csv(file.path(project_root, "data", "intermediate", "taxonomy_soil_only.csv"))

# Define the taxonomic levels
tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Create full lineage string (KPCOFG)
lineage_counts <- taxonomy %>%
  unite("FullLineage", all_of(tax_levels), sep = ";", remove = FALSE) %>%
  count(FullLineage, name = "Count") %>%
  arrange(desc(Count))

# Save lineage count table
write_csv(lineage_counts, file.path(project_root, "data", "processed", "lineage_counts_all.csv"))

# === Step 2: Generate lineage filtering results for thresholds 120 ===
total_asvs <- sum(lineage_counts$Count)
total_lineages <- nrow(lineage_counts)

# Initialize result table
results <- data.frame(
  Threshold = integer(),
  LineagesRemaining = integer(),
  FractionRemaining = numeric(),
  ASVsLost = integer()
)

for (t in 1:20) {
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

# Save output
write_csv(results, file.path(project_root, "data", "processed", "lineage_counts_thresholds_1_to_20.csv"))
print(results)
