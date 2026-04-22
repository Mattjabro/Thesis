library(tidyverse)

# Resolve project root
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

# === Step 1: Load taxonomy and taxa matrix ===

taxonomy <- read_csv(file.path(project_root, "data", "intermediate", "taxonomy_soil_only.csv")) %>%
  unite("Lineage", Kingdom:Genus, sep = ";", remove = FALSE)


taxa_matrix <- read.csv(
  file.path(project_root, "data", "intermediate", "genus_matrix_soil_only.csv"),
  row.names = 1, check.names = FALSE
)

# === Step 2: Lineage prevalence across samples ===
presence_matrix <- taxa_matrix > 0
lineage_prevalence <- colSums(presence_matrix)  # named vector of lineages

# === Step 3: Match each ASV to its lineage prevalence ===
taxonomy <- taxonomy %>%
  mutate(LineagePrevalence = lineage_prevalence[Lineage]) %>%
  replace_na(list(LineagePrevalence = 0))

# === Step 4: Setup for counts ===
total_asvs <- nrow(taxonomy)
total_lineages <- length(lineage_prevalence)
total_matrix_sum <- sum(taxa_matrix)

# === Step 5: Threshold loop ===
thresholds <- 1:250

results <- data.frame(
  Threshold = thresholds,
  LineagesRemaining = NA,
  FractionLineagesRemaining = NA,
  FractionCountsRemaining = NA,
  FractionMatrixRemaining = NA
)

for (t in thresholds) {
  # Lineages retained
  kept_lineages <- names(lineage_prevalence[lineage_prevalence >= t])
  n_lineages <- length(kept_lineages)
  frac_lineages <- n_lineages / total_lineages
  
  # ASVs retained
  retained_asvs <- taxonomy %>% filter(Lineage %in% kept_lineages)
  frac_asvs <- nrow(retained_asvs) / total_asvs
  
  # Matrix signal retained
  matrix_kept <- taxa_matrix[, colnames(taxa_matrix) %in% kept_lineages]
  frac_matrix <- sum(matrix_kept) / total_matrix_sum
  
  results[results$Threshold == t, ] <- list(
    Threshold = t,
    LineagesRemaining = n_lineages,
    FractionLineagesRemaining = frac_lineages,
    FractionCountsRemaining = frac_asvs,
    FractionMatrixRemaining = frac_matrix
  )
}

# === Step 6: Save table ===
write_csv(results, file.path(project_root, "data", "processed", "lineage_prevalence_thresholds.csv"))

# === Step 7: Plot 1  Fraction of lineages retained ===
p1 <- ggplot(results, aes(x = Threshold, y = FractionLineagesRemaining)) +
  geom_line() + geom_point() +
  labs(
    title = "Fraction of Lineages Remaining vs Threshold",
    x = "Minimum Sample Prevalence",
    y = "Fraction of Lineages Remaining"
  ) +
  theme_minimal()

ggsave(
  file.path(project_root, "outputs", "figures", "lineage_fraction_remaining_plot.png"),
  p1, width = 6, height = 4, dpi = 300
)

# === Step 8: Plot 2  Fraction of matrix signal retained ===
p2 <- ggplot(results, aes(x = Threshold, y = FractionMatrixRemaining)) +
  geom_line() + geom_point() +
  labs(
    title = "Fraction of Total Matrix Signal Remaining vs Threshold",
    x = "Minimum Sample Prevalence",
    y = "Fraction of Relative Abundance Retained"
  ) +
  theme_minimal()

ggsave(
  file.path(project_root, "outputs", "figures", "lineage_fraction_matrix_remaining_plot.png"),
  p2, width = 6, height = 4, dpi = 300
)
