#!/usr/bin/env Rscript
# ============================================================================
# 23_calculate_alpha_diversity.R
#
# Calculates alpha diversity metrics (Observed OTUs, Chao1, Shannon) for
# global soil samples from the DADA2-derived ASV table (bacteria only).
# Outputs a summary table and distribution plots for quality control.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

timestamp <- function(msg) {
  cat(sprintf("[%s]  %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
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

# ============================================================================
# Define file paths
# ============================================================================
seqtab_fp <- file.path(project_root, "data", "intermediate", "seqtab_nochim_bacteria.rds")
out_dir_processed <- file.path(project_root, "data", "processed")
out_dir_figures <- file.path(project_root, "outputs", "figures")

if (!dir.exists(out_dir_processed)) {
  dir.create(out_dir_processed, recursive = TRUE)
}
if (!dir.exists(out_dir_figures)) {
  dir.create(out_dir_figures, recursive = TRUE)
}

alpha_div_fp <- file.path(out_dir_processed, "alpha_diversity_soil.csv")
fig_fp <- file.path(out_dir_figures, "alpha_diversity_distributions.png")

cat("============================================================\n")
cat("  Calculating Alpha Diversity Metrics (Global Soil)\n")
cat("============================================================\n")

# ============================================================================
# Step 1: Load ASV Table
# ============================================================================
timestamp("Loading ASV table (seqtab_nochim_bacteria.rds)...")
cat("  Input: ", seqtab_fp, "\n")

seqtab <- readRDS(seqtab_fp)
cat("  ASV table dimensions: ", nrow(seqtab), " samples x ", ncol(seqtab), " ASVs\n")

# Matrix is already in samples (rows) x ASVs (columns) format for vegan

# ============================================================================
# Step 2: Calculate Alpha Diversity Metrics
# ============================================================================
timestamp("Calculating alpha diversity metrics...")

# Observed OTUs (ASVs)
cat("  - Observed OTUs/ASVs...\n")
observed <- specnumber(seqtab)

# Chao1 estimator
cat("  - Chao1 estimator...\n")
richness_estimates <- estimateR(seqtab)
chao1 <- richness_estimates[2, ]  # Row 2 is Chao1

# Shannon diversity
cat("  - Shannon diversity index...\n")
shannon <- diversity(seqtab, index = "shannon")

# ============================================================================
# Step 3: Build Output Data Frame
# ============================================================================
timestamp("Building output table...")

# Extract sample IDs from seqtab rownames
sample_ids <- rownames(seqtab)

alpha_div <- data.frame(
  SampleID = sample_ids,
  adiv_observed_otus = observed,
  adiv_chao1 = chao1,
  adiv_shannon = shannon,
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Sort by sample ID for consistency
alpha_div <- alpha_div %>% arrange(SampleID)

cat("  Summary statistics:\n")
cat("    Observed OTUs - Mean: ", round(mean(alpha_div$adiv_observed_otus), 2),
    ", Median: ", round(median(alpha_div$adiv_observed_otus), 2), "\n")
cat("    Chao1         - Mean: ", round(mean(alpha_div$adiv_chao1), 2),
    ", Median: ", round(median(alpha_div$adiv_chao1), 2), "\n")
cat("    Shannon       - Mean: ", round(mean(alpha_div$adiv_shannon), 2),
    ", Median: ", round(median(alpha_div$adiv_shannon), 2), "\n")

# ============================================================================
# Step 4: Write Output
# ============================================================================
timestamp("Writing alpha diversity table...")
write_csv(alpha_div, alpha_div_fp)
cat("   Written: ", alpha_div_fp, "\n")

# ============================================================================
# Step 5: Create Distribution Plots
# ============================================================================
timestamp("Creating distribution plots...")

png(fig_fp, width = 1400, height = 500)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

# Observed OTUs histogram
hist(alpha_div$adiv_observed_otus,
     breaks = 30,
     col = "steelblue",
     border = "white",
     main = "Observed ASVs per Sample",
     xlab = "Number of Observed ASVs",
     ylab = "Sample Count",
     cex.main = 1.3,
     cex.lab = 1.2)
abline(v = median(alpha_div$adiv_observed_otus), col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("Median =", round(median(alpha_div$adiv_observed_otus), 0)),
       col = "red", lty = 2, lwd = 2, cex = 0.9)

# Chao1 histogram
hist(alpha_div$adiv_chao1,
     breaks = 30,
     col = "seagreen",
     border = "white",
     main = "Chao1 Estimated Richness",
     xlab = "Chao1 Estimate",
     ylab = "Sample Count",
     cex.main = 1.3,
     cex.lab = 1.2)
abline(v = median(alpha_div$adiv_chao1), col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("Median =", round(median(alpha_div$adiv_chao1), 0)),
       col = "red", lty = 2, lwd = 2, cex = 0.9)

# Shannon histogram
hist(alpha_div$adiv_shannon,
     breaks = 30,
     col = "coral",
     border = "white",
     main = "Shannon Diversity Index",
     xlab = "Shannon Index",
     ylab = "Sample Count",
     cex.main = 1.3,
     cex.lab = 1.2)
abline(v = median(alpha_div$adiv_shannon), col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("Median =", round(median(alpha_div$adiv_shannon), 2)),
       col = "red", lty = 2, lwd = 2, cex = 0.9)

dev.off()
cat("   Written: ", fig_fp, "\n")

cat("\n============================================================\n")
cat("   Alpha diversity calculation complete!\n")
cat("============================================================\n")
