# Load required packages
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

# Load the soil taxa matrix
taxa_matrix <- read.csv(
  file.path(intermediate_dir, "sample_taxa_matrix_L6_bacteria.csv"),
  row.names = 1,
  check.names = FALSE
)

# Convert to presence/absence matrix
presence_absence <- taxa_matrix > 0

# Count in how many samples each taxon appears
taxa_prevalence <- colSums(presence_absence)

# Total number of taxa
total_taxa <- length(taxa_prevalence)

# Thresholds from 1 to 20
thresholds <- 1:20

# Calculate fraction of taxa remaining at each threshold
fraction_remaining <- sapply(thresholds, function(thresh) {
  sum(taxa_prevalence >= thresh) / total_taxa
})

# Compile results into a data frame
result <- data.frame(
  Threshold = thresholds,
  FractionRemaining = fraction_remaining
)

# Print results
print(result)

# Save to CSV
write.csv(
  result,
  file.path(processed_dir, "soil_fraction_taxa_remaining_1_to_20.csv"),
  row.names = FALSE
)

# Optional: Plot
ggplot(result, aes(x = Threshold, y = FractionRemaining)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Fraction of Soil Taxa Remaining vs. Sample Presence Threshold",
    x = "Minimum # Samples a Taxon Must Appear In",
    y = "Fraction of Taxa Remaining"
  ) +
  theme_minimal()
