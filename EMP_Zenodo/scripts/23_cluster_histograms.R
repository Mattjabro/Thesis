#!/usr/bin/env Rscript
# ============================================================================
# 23_cluster_histograms.R
#
# Generates one histogram per cluster showing the distribution of CLR values
# across all taxa for that cluster. Provides visual summary of cluster-specific
# taxonomic profiles.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
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

# ============ Parse arguments or use defaults ============
args <- commandArgs(trailingOnly = TRUE)
cluster_avg_fp <- if (length(args) >= 1) args[1] else {
  file.path(project_root, "outputs", "cluster_metadata", "cluster_avg_clr_matrix.csv")
}
out_dir <- file.path(project_root, "outputs", "cluster_histograms")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("============================================================\n")
cat("  Cluster Average CLR Histograms (EMP)\n")
cat("============================================================\n")
cat("Reading cluster averages from:", cluster_avg_fp, "\n")
cat("Output directory:", out_dir, "\n\n")

# ============ Read cluster average matrix ============
cluster_avg <- read.csv(cluster_avg_fp, row.names = 1, check.names = FALSE)
cat("Cluster average matrix dimensions:", nrow(cluster_avg), "taxa x", ncol(cluster_avg), "clusters\n\n")

# ============ Generate histogram for each cluster ============
summary_stats <- data.frame(
  cluster = character(),
  n_taxa = integer(),
  mean_clr = numeric(),
  median_clr = numeric(),
  sd_clr = numeric(),
  min_clr = numeric(),
  max_clr = numeric(),
  q25_clr = numeric(),
  q75_clr = numeric(),
  stringsAsFactors = FALSE
)

cat("Generating bar plots for each cluster (each bar = one taxon)...\n")
for (cluster_col in colnames(cluster_avg)) {
  clr_values <- cluster_avg[[cluster_col]]
  taxa_names <- rownames(cluster_avg)

  # Extract cluster number for labeling
  cluster_num <- as.integer(gsub("cluster_", "", cluster_col))

  # Sort taxa by CLR value (descending) for better visualization
  sorted_idx <- order(clr_values, decreasing = TRUE)
  clr_sorted <- clr_values[sorted_idx]
  taxa_sorted <- taxa_names[sorted_idx]

  # Shorten taxa names for plotting (use Phylum;Genus format)
  taxa_short <- sapply(strsplit(taxa_sorted, ";"), function(x) {
    if (length(x) >= 6) {
      paste(x[2], x[6], sep = ";")
    } else if (length(x) >= 2) {
      x[2]
    } else {
      x[1]
    }
  })

  # Create PDF for this cluster
  pdf_file <- file.path(out_dir, sprintf("cluster_%02d_histogram.pdf", cluster_num))
  pdf(pdf_file, width = 16, height = 10)

  # Create bar plot
  par(mar = c(12, 4, 4, 2))
  barplot(
    clr_sorted,
    names.arg = taxa_short,
    main = sprintf("Cluster %d: Taxa CLR Values (n=%d taxa)",
                   cluster_num, length(clr_values)),
    xlab = "",
    ylab = "Mean CLR Value",
    col = ifelse(clr_sorted > 0, "steelblue", "coral"),
    border = NA,
    las = 2,
    cex.names = 0.35,
    cex.axis = 0.9
  )
  mtext("Taxa (Phylum;Genus)", side = 1, line = 10.5, cex = 1)

  # Add horizontal line at zero
  abline(h = 0, col = "black", lwd = 1.5, lty = 2)

  # Add text box with summary stats
  legend(
    "topright",
    legend = c(
      sprintf("Mean: %.2f", mean(clr_values)),
      sprintf("Median: %.2f", median(clr_values)),
      sprintf("SD: %.2f", sd(clr_values)),
      sprintf("Range: [%.2f, %.2f]", min(clr_values), max(clr_values))
    ),
    bty = "n",
    cex = 0.9
  )

  dev.off()

  # Calculate summary statistics
  stats <- data.frame(
    cluster = cluster_col,
    n_taxa = length(clr_values),
    mean_clr = mean(clr_values),
    median_clr = median(clr_values),
    sd_clr = sd(clr_values),
    min_clr = min(clr_values),
    max_clr = max(clr_values),
    q25_clr = quantile(clr_values, 0.25),
    q75_clr = quantile(clr_values, 0.75),
    stringsAsFactors = FALSE
  )

  summary_stats <- rbind(summary_stats, stats)

  cat(sprintf("  Cluster %02d: %.2f  %.2f (median: %.2f, range: [%.2f, %.2f])\n",
              cluster_num,
              mean(clr_values),
              sd(clr_values),
              median(clr_values),
              min(clr_values),
              max(clr_values)))
}

# ============ Save summary statistics ============
write.csv(
  summary_stats,
  file.path(out_dir, "histogram_summary_stats.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\n============================================================\n")
cat("  Summary\n")
cat("============================================================\n")
cat("Generated", nrow(summary_stats), "histograms\n")
cat("Summary statistics saved to: histogram_summary_stats.csv\n")
cat("\nCluster histograms complete!\n")
cat("============================================================\n")
