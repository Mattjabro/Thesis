#!/usr/bin/env Rscript
# ============================================================================
# 20_cluster_clr_statistics.R
#
# Calculates enhanced cluster statistics with error bars for CLR values.
# Computes mean, SD, SE, and 95% CI for each lineage within each cluster.
# Generates error bar visualizations and variance heatmaps.
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
clr_fp <- if (length(args) >= 1) args[1] else {
  file.path(project_root, "data", "processed", "lineage_clr_matrix_soil_filtered.csv")
}
cluster_fp <- if (length(args) >= 2) args[2] else {
  file.path(project_root, "outputs", "cluster_metadata", "cluster_assignments.tsv")
}
out_dir_meta <- file.path(project_root, "outputs", "cluster_metadata")
out_dir_ord <- file.path(project_root, "outputs", "ordination")

if (!dir.exists(out_dir_meta)) {
  dir.create(out_dir_meta, recursive = TRUE)
}
if (!dir.exists(out_dir_ord)) {
  dir.create(out_dir_ord, recursive = TRUE)
}

cat("============================================================\n")
cat("  Computing Cluster CLR Statistics with Error Bars (EMP)\n")
cat("============================================================\n")

# ============ Step 1: Load CLR Matrix ============
cat("\n[Step 1] Loading CLR matrix from:\n  ", clr_fp, "\n")
clr_matrix <- as.matrix(read.csv(clr_fp, row.names = 1, check.names = FALSE))
cat("  CLR matrix dimensions: ", nrow(clr_matrix), " taxa x ", ncol(clr_matrix), " samples\n")

# ============ Step 2: Load Cluster Assignments ============
cat("\n[Step 2] Loading cluster assignments from:\n  ", cluster_fp, "\n")
clusters <- read_tsv(cluster_fp, show_col_types = FALSE)
cat("  Total samples with cluster assignments: ", nrow(clusters), "\n")
cat("  Unique clusters: ", n_distinct(clusters$cluster), "\n")

# ============ Step 3: Match Samples ============
cat("\n[Step 3] Matching samples...\n")
clr_samples <- colnames(clr_matrix)
cluster_samples <- clusters$X.SampleID

shared_samples <- intersect(clr_samples, cluster_samples)
cat("  Shared samples: ", length(shared_samples), "\n")

if (length(shared_samples) == 0) {
  stop("ERROR: No shared samples!")
}

# Filter to shared samples
clr_matrix <- clr_matrix[, shared_samples, drop = FALSE]
clusters <- clusters %>% filter(X.SampleID %in% shared_samples)

# ============ Step 4: Convert to Long Format ============
cat("\n[Step 4] Converting matrix to long format...\n")

clr_long <- clr_matrix %>%
  as.data.frame() %>%
  rownames_to_column("lineage") %>%
  pivot_longer(cols = -lineage, names_to = "SampleID", values_to = "clr_value") %>%
  left_join(clusters %>% select(X.SampleID, cluster), by = c("SampleID" = "X.SampleID"))

cat("  Long format: ", nrow(clr_long), " rows\n")

# ============ Step 5: Calculate Statistics per Lineage per Cluster ============
cat("\n[Step 5] Calculating statistics (mean, SD, SE, CI)...\n")

cluster_stats <- clr_long %>%
  group_by(lineage, cluster) %>%
  summarise(
    n_samples = n(),
    mean_clr = mean(clr_value, na.rm = TRUE),
    sd_clr = sd(clr_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_clr = sd_clr / sqrt(n_samples),
    ci_lower = mean_clr - 1.96 * se_clr,
    ci_upper = mean_clr + 1.96 * se_clr
  )

cat("  Statistics calculated for ", nrow(cluster_stats), " lineage-cluster combinations\n")

# ============ Step 6: Write Statistics Table ============
cat("\n[Step 6] Writing statistics table...\n")
stats_fp <- file.path(out_dir_meta, "cluster_clr_statistics.csv")
write_csv(cluster_stats, stats_fp)
cat("  Written: ", stats_fp, "\n")

# Summary by cluster
cluster_summary <- cluster_stats %>%
  group_by(cluster) %>%
  summarise(
    n_lineages = n(),
    mean_sd = mean(sd_clr, na.rm = TRUE),
    median_sd = median(sd_clr, na.rm = TRUE),
    max_sd = max(sd_clr, na.rm = TRUE),
    .groups = "drop"
  )

summary_fp <- file.path(out_dir_meta, "cluster_clr_statistics_summary.tsv")
write_tsv(cluster_summary, summary_fp)
cat("  Written: ", summary_fp, "\n")

# ============ Step 7: Identify Top Variable Taxa ============
cat("\n[Step 7] Identifying top variable taxa across clusters...\n")

# Calculate variance across clusters for each lineage
taxa_variance <- cluster_stats %>%
  group_by(lineage) %>%
  summarise(
    mean_across_clusters = mean(mean_clr),
    var_across_clusters = var(mean_clr),
    .groups = "drop"
  ) %>%
  arrange(desc(var_across_clusters))

top20_taxa <- head(taxa_variance$lineage, 20)
cat("  Top 20 most variable taxa identified\n")

# ============ Step 8: Create Error Bar Plot (Top 20 Taxa) ============
cat("\n[Step 8] Creating error bar plot for top 20 taxa...\n")

plot_data_top20 <- cluster_stats %>%
  filter(lineage %in% top20_taxa) %>%
  mutate(lineage_short = substr(lineage, 1, 60))

pdf(file.path(out_dir_ord, "cluster_errorbar_top20_taxa.pdf"), width = 16, height = 12)
p <- ggplot(plot_data_top20, aes(x = factor(cluster), y = mean_clr, color = lineage_short)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.3, alpha = 0.7) +
  facet_wrap(~lineage_short, scales = "free_y", ncol = 4) +
  theme_bw() +
  labs(
    title = "Mean CLR  95% CI by Cluster (Top 20 Most Variable Taxa)",
    x = "Cluster",
    y = "CLR Value"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 7)
  )
print(p)
dev.off()
cat("  Written: ", file.path(out_dir_ord, "cluster_errorbar_top20_taxa.pdf"), "\n")

# ============ Step 9: Create Error Bar Plot for Selected Phyla ============
cat("\n[Step 9] Creating error bar plot for key phyla...\n")

# Identify major phyla
key_phyla_patterns <- c("Proteobacteria", "Actinobacteria", "Acidobacteria",
                       "Bacteroidetes", "Firmicutes", "Verrucomicrobia")

selected_phyla_taxa <- cluster_stats %>%
  filter(str_detect(lineage, paste(key_phyla_patterns, collapse = "|"))) %>%
  group_by(lineage) %>%
  summarise(var_clr = var(mean_clr), .groups = "drop") %>%
  arrange(desc(var_clr)) %>%
  head(12) %>%
  pull(lineage)

if (length(selected_phyla_taxa) > 0) {
  plot_data_phyla <- cluster_stats %>%
    filter(lineage %in% selected_phyla_taxa) %>%
    mutate(lineage_short = substr(lineage, 1, 60))

  pdf(file.path(out_dir_ord, "cluster_errorbar_selected_phyla.pdf"), width = 14, height = 10)
  p <- ggplot(plot_data_phyla, aes(x = factor(cluster), y = mean_clr, color = lineage_short)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.3, alpha = 0.7) +
    facet_wrap(~lineage_short, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(
      title = "Mean CLR  95% CI by Cluster (Key Phyla)",
      x = "Cluster",
      y = "CLR Value"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 8)
    )
  print(p)
  dev.off()
  cat("  Written: ", file.path(out_dir_ord, "cluster_errorbar_selected_phyla.pdf"), "\n")
} else {
  cat("  No key phyla taxa found, skipping phyla plot\n")
}

# ============ Step 10: Create Variance Heatmap ============
cat("\n[Step 10] Creating variance heatmap...\n")

# Create matrix of SD values: taxa  clusters
sd_matrix <- cluster_stats %>%
  select(lineage, cluster, sd_clr) %>%
  pivot_wider(names_from = cluster, values_from = sd_clr, names_prefix = "cluster_") %>%
  column_to_rownames("lineage") %>%
  as.matrix()

# Select top 50 most variable taxa for heatmap
top50_for_heatmap <- head(taxa_variance$lineage, 50)
sd_matrix_top50 <- sd_matrix[top50_for_heatmap, , drop = FALSE]

# Truncate row names for readability
rownames(sd_matrix_top50) <- substr(rownames(sd_matrix_top50), 1, 60)

pdf(file.path(out_dir_ord, "cluster_variance_heatmap.pdf"), width = 10, height = 14)
heatmap(sd_matrix_top50,
        Rowv = NA,
        Colv = NA,
        col = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
        scale = "none",
        margins = c(10, 15),
        main = "Standard Deviation of CLR Values (Top 50 Taxa)",
        cex.main = 1.2,
        cex.axis = 0.8)
dev.off()
cat("  Written: ", file.path(out_dir_ord, "cluster_variance_heatmap.pdf"), "\n")

cat("\n============================================================\n")
cat("  Cluster CLR statistics calculation complete!\n")
cat("  Created 3-4 PDF visualizations and 2 data tables\n")
cat("============================================================\n")
