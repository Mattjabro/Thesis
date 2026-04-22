# Load necessary packages
library(tidyverse)

# Resolve project root (works with Rscript and interactive runs)
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

# Load the EMP taxa matrix
taxa_matrix <- read.csv(
  file.path(project_root, "data", "intermediate", "genus_matrix_soil_only.csv"),
  row.names = 1, check.names = FALSE
)

# Convert to presence/absence (1 = present, 0 = absent)
presence_absence <- taxa_matrix > 0

# Count how many samples each taxon appears in
taxa_prevalence <- colSums(presence_absence)

# Total number of taxa
total_taxa <- length(taxa_prevalence)

# Generate threshold values from 1 to 20
thresholds <- 1:20

# Calculate fraction of taxa remaining for each threshold
fraction_remaining <- sapply(thresholds, function(thresh) {
  sum(taxa_prevalence >= thresh) / total_taxa
})

# Combine into a data frame
result <- data.frame(
  Threshold = thresholds,
  FractionRemaining = fraction_remaining
)

# Print the result
print(result)

# Save to CSV
write.csv(
  result,
  file.path(project_root, "data", "processed", "genus_prevalence_thresholds_1_to_20.csv"),
  row.names = FALSE
)

# Optional: plot the result
ggplot(result, aes(x = Threshold, y = FractionRemaining)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Fraction of Taxa Remaining vs. Sample Presence Threshold",
    x = "Minimum # Samples a Taxon Must Appear In",
    y = "Fraction of Taxa Remaining"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(project_root, "outputs", "figures", "genus_prevalence_thresholds_1_to_20.png"),
  width = 6, height = 4, dpi = 300
)
