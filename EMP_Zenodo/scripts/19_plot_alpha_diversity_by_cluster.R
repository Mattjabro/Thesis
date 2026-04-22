#!/usr/bin/env Rscript
# ============================================================================
# 19_plot_alpha_diversity_by_cluster.R
#
# Joins alpha diversity metrics with cluster assignments and creates
# visualizations showing how alpha diversity varies across clusters.
# Generates PDF plots similar to the cluster-metadata analysis.
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

# ============================================================================
# Define file paths
# ============================================================================
alpha_div_fp <- file.path(project_root, "data", "processed", "alpha_diversity_soil.csv")
cluster_fp <- file.path(project_root, "outputs", "cluster_metadata", "cluster_assignments.tsv")
out_dir <- file.path(project_root, "outputs", "cluster_metadata")
plot_dir <- file.path(out_dir, "plots")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

cluster_alpha_fp <- file.path(out_dir, "cluster_alpha_diversity.tsv")

cat("============================================================\n")
cat("  Plotting Alpha Diversity by Cluster (EMP)\n")
cat("============================================================\n")

# ============================================================================
# Step 1: Load Data
# ============================================================================
cat("\n[Step 1] Loading alpha diversity and cluster assignments...\n")
alpha_div <- read_csv(alpha_div_fp, show_col_types = FALSE)
clusters <- read_tsv(cluster_fp, show_col_types = FALSE)

cat("  Alpha diversity: ", nrow(alpha_div), " samples\n")
cat("  Cluster assignments: ", nrow(clusters), " samples\n")

# ============================================================================
# Step 2: Join Data
# ============================================================================
cat("\n[Step 2] Joining alpha diversity with clusters...\n")

# Standardize sample ID column names
if ("X.SampleID" %in% colnames(clusters)) {
  clusters <- clusters %>% rename(SampleID = X.SampleID)
}

# Join
cluster_alpha <- alpha_div %>%
  inner_join(clusters, by = "SampleID")

cat("  Joined data: ", nrow(cluster_alpha), " samples\n")
cat("  Clusters present: ", n_distinct(cluster_alpha$cluster), "\n")

# ============================================================================
# Step 3: Write Joined Table
# ============================================================================
cat("\n[Step 3] Writing joined table...\n")
write_tsv(cluster_alpha, cluster_alpha_fp)
cat("  Written: ", cluster_alpha_fp, "\n")

# ============================================================================
# Step 4: Summary Statistics by Cluster
# ============================================================================
cat("\n[Step 4] Computing summary statistics by cluster...\n")
cluster_summary <- cluster_alpha %>%
  group_by(cluster) %>%
  summarise(
    n_samples = n(),
    mean_observed = mean(adiv_observed_otus, na.rm = TRUE),
    median_observed = median(adiv_observed_otus, na.rm = TRUE),
    sd_observed = sd(adiv_observed_otus, na.rm = TRUE),
    mean_chao1 = mean(adiv_chao1, na.rm = TRUE),
    median_chao1 = median(adiv_chao1, na.rm = TRUE),
    sd_chao1 = sd(adiv_chao1, na.rm = TRUE),
    mean_shannon = mean(adiv_shannon, na.rm = TRUE),
    median_shannon = median(adiv_shannon, na.rm = TRUE),
    sd_shannon = sd(adiv_shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(cluster)

summary_fp <- file.path(out_dir, "cluster_alpha_diversity_summary.tsv")
write_tsv(cluster_summary, summary_fp)
cat("  Written: ", summary_fp, "\n")

print(cluster_summary)

# ============================================================================
# Step 5: Create Boxplot Visualizations
# ============================================================================
cat("\n[Step 5] Creating boxplot visualizations...\n")

# Function to create a boxplot PDF
create_boxplot <- function(data, metric_col, metric_name, filename) {
  pdf(file.path(plot_dir, filename), width = 12, height = 6)

  par(mar = c(5, 5, 4, 2))

  # Create boxplot
  boxplot(as.formula(paste(metric_col, "~ cluster")),
          data = data,
          col = "lightblue",
          border = "darkblue",
          main = paste(metric_name, "by Cluster"),
          xlab = "Cluster",
          ylab = metric_name,
          cex.main = 1.5,
          cex.lab = 1.3,
          cex.axis = 1.1,
          las = 1)

  # Add median line
  abline(h = median(data[[metric_col]], na.rm = TRUE), col = "red", lty = 2, lwd = 2)

  # Add sample counts on top
  cluster_counts <- table(data$cluster)
  axis(3, at = 1:length(cluster_counts), labels = paste0("n=", cluster_counts),
       cex.axis = 0.9, col.axis = "darkgray")

  dev.off()
  cat("  Written: ", file.path(plot_dir, filename), "\n")
}

# Create boxplots for each metric
create_boxplot(cluster_alpha, "adiv_observed_otus",
               "Observed ASVs", "alpha_div_observed_by_cluster.pdf")
create_boxplot(cluster_alpha, "adiv_chao1",
               "Chao1 Estimate", "alpha_div_chao1_by_cluster.pdf")
create_boxplot(cluster_alpha, "adiv_shannon",
               "Shannon Index", "alpha_div_shannon_by_cluster.pdf")

# ============================================================================
# Step 6: Create Violin Plot with Points
# ============================================================================
cat("\n[Step 6] Creating detailed violin plots...\n")

create_violin_plot <- function(data, metric_col, metric_name, filename) {
  pdf(file.path(plot_dir, filename), width = 14, height = 7)

  # Base plot
  plot.new()
  plot.window(xlim = c(0.5, max(data$cluster) + 0.5),
              ylim = range(data[[metric_col]], na.rm = TRUE))

  # Add axes
  axis(1, at = sort(unique(data$cluster)), labels = sort(unique(data$cluster)))
  axis(2, las = 1)
  box()

  # Title and labels
  title(main = paste(metric_name, "Distribution by Cluster"),
        xlab = "Cluster",
        ylab = metric_name,
        cex.main = 1.5,
        cex.lab = 1.3)

  # Add points for each cluster
  for (cl in sort(unique(data$cluster))) {
    cluster_data <- data %>% filter(cluster == cl)
    points(jitter(rep(cl, nrow(cluster_data)), amount = 0.1),
           cluster_data[[metric_col]],
           col = rgb(0, 0, 1, 0.3),
           pch = 16,
           cex = 0.5)

    # Add median line
    median_val <- median(cluster_data[[metric_col]], na.rm = TRUE)
    segments(cl - 0.3, median_val, cl + 0.3, median_val,
             col = "red", lwd = 3)
  }

  dev.off()
  cat("  Written: ", file.path(plot_dir, filename), "\n")
}

create_violin_plot(cluster_alpha, "adiv_observed_otus",
                   "Observed ASVs", "alpha_div_observed_scatter_by_cluster.pdf")
create_violin_plot(cluster_alpha, "adiv_chao1",
                   "Chao1 Estimate", "alpha_div_chao1_scatter_by_cluster.pdf")
create_violin_plot(cluster_alpha, "adiv_shannon",
                   "Shannon Index", "alpha_div_shannon_scatter_by_cluster.pdf")

# ============================================================================
# Step 7: Create Heatmap of Mean Alpha Diversity
# ============================================================================
cat("\n[Step 7] Creating heatmap summary...\n")

pdf(file.path(plot_dir, "alpha_div_heatmap_by_cluster.pdf"), width = 10, height = 6)

# Normalize each metric to 0-1 scale for comparison
heatmap_data <- cluster_summary %>%
  select(cluster, mean_observed, mean_chao1, mean_shannon) %>%
  mutate(
    norm_observed = (mean_observed - min(mean_observed)) / (max(mean_observed) - min(mean_observed)),
    norm_chao1 = (mean_chao1 - min(mean_chao1)) / (max(mean_chao1) - min(mean_chao1)),
    norm_shannon = (mean_shannon - min(mean_shannon)) / (max(mean_shannon) - min(mean_shannon))
  ) %>%
  select(cluster, norm_observed, norm_chao1, norm_shannon) %>%
  as.data.frame()

rownames(heatmap_data) <- paste0("Cluster_", heatmap_data$cluster)
heatmap_matrix <- as.matrix(heatmap_data[, -1])
colnames(heatmap_matrix) <- c("Observed ASVs", "Chao1", "Shannon")

# Create heatmap
heatmap(heatmap_matrix,
        Rowv = NA,
        Colv = NA,
        col = colorRampPalette(c("blue", "white", "red"))(100),
        scale = "none",
        margins = c(10, 10),
        main = "Normalized Mean Alpha Diversity by Cluster",
        cex.main = 1.3)

dev.off()
cat("  Written: ", file.path(plot_dir, "alpha_div_heatmap_by_cluster.pdf"), "\n")

cat("\n============================================================\n")
cat("  Alpha diversity cluster visualization complete!\n")
cat("  Created 7 PDF plots in: ", plot_dir, "\n")
cat("============================================================\n")
