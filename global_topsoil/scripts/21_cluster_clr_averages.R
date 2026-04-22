#!/usr/bin/env Rscript
# ============================================================================
# 21_cluster_clr_averages.R
# 
# Computes the average CLR vector for each cluster.
# Reads the CLR matrix and cluster assignments, calculates mean CLR values
# per lineage for each cluster, and writes the output.
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
out_dir <- if (length(args) >= 3) args[3] else {
  file.path(project_root, "outputs", "cluster_metadata")
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("============================================================\n")
cat("  Computing Average CLR Vectors per Cluster (Global Topsoil)\n")
cat("============================================================\n")

# ============ Step 1: Load CLR Matrix ============
cat("\n[Step 1] Loading CLR matrix from:\n  ", clr_fp, "\n")
clr_matrix <- as.matrix(read.csv(clr_fp, row.names = 1, check.names = FALSE))
# Matrix is taxa (rows) x samples (columns)
cat("  CLR matrix dimensions: ", nrow(clr_matrix), " taxa x ", ncol(clr_matrix), " samples\n")

# ============ Step 2: Load Cluster Assignments ============
cat("\n[Step 2] Loading cluster assignments from:\n  ", cluster_fp, "\n")
clusters <- read_tsv(cluster_fp, show_col_types = FALSE)
cat("  Total samples with cluster assignments: ", nrow(clusters), "\n")
cat("  Unique clusters: ", n_distinct(clusters$cluster), "\n")

# ============ Step 3: Match Samples ============
cat("\n[Step 3] Matching samples between CLR matrix and cluster assignments...\n")
clr_samples <- colnames(clr_matrix)
cluster_samples <- clusters$X.SampleID

shared_samples <- intersect(clr_samples, cluster_samples)
cat("  Samples in CLR matrix: ", length(clr_samples), "\n")
cat("  Samples in cluster file: ", length(cluster_samples), "\n")
cat("  Shared samples: ", length(shared_samples), "\n")

if (length(shared_samples) == 0) {
  stop("ERROR: No shared samples between CLR matrix and cluster assignments!")
}

# Filter to shared samples
clr_matrix <- clr_matrix[, shared_samples, drop = FALSE]
clusters <- clusters %>% filter(X.SampleID %in% shared_samples)

# ============ Step 4: Compute Average CLR per Cluster ============
cat("\n[Step 4] Computing average CLR vector for each cluster...\n")

# Create a named vector of cluster assignments
cluster_vec <- setNames(clusters$cluster, clusters$X.SampleID)

# Get unique clusters sorted
unique_clusters <- sort(unique(cluster_vec))
cat("  Clusters found: ", paste(unique_clusters, collapse = ", "), "\n")

# Initialize result matrix
avg_clr_matrix <- matrix(
  NA_real_,
  nrow = nrow(clr_matrix),
  ncol = length(unique_clusters),
  dimnames = list(rownames(clr_matrix), paste0("cluster_", unique_clusters))
)

# Compute mean for each cluster
for (cl in unique_clusters) {
  samples_in_cluster <- names(cluster_vec[cluster_vec == cl])
  n_samples <- length(samples_in_cluster)
  
  if (n_samples == 1) {
    avg_clr_matrix[, paste0("cluster_", cl)] <- clr_matrix[, samples_in_cluster]
  } else {
    avg_clr_matrix[, paste0("cluster_", cl)] <- rowMeans(clr_matrix[, samples_in_cluster, drop = FALSE])
  }
  
  cat("  Cluster ", cl, ": ", n_samples, " samples\n")
}

# ============ Step 5: Write Outputs ============
cat("\n[Step 5] Writing output files...\n")

# Output 1: Full matrix (taxa x clusters)
out_matrix_fp <- file.path(out_dir, "cluster_avg_clr_matrix.csv")
avg_clr_df <- as.data.frame(avg_clr_matrix)
avg_clr_df$lineage <- rownames(avg_clr_matrix)
avg_clr_df <- avg_clr_df[, c("lineage", paste0("cluster_", unique_clusters))]
write_csv(avg_clr_df, out_matrix_fp)
cat("  Written: ", out_matrix_fp, "\n")

# Output 2: Long format (for easier plotting)
out_long_fp <- file.path(out_dir, "cluster_avg_clr_long.tsv")
long_df <- avg_clr_df %>%
  pivot_longer(
    cols = starts_with("cluster_"),
    names_to = "cluster",
    names_prefix = "cluster_",
    values_to = "avg_clr"
  ) %>%
  mutate(cluster = as.integer(cluster)) %>%
  arrange(cluster, lineage)
write_tsv(long_df, out_long_fp)
cat("  Written: ", out_long_fp, "\n")

# Output 3: Summary statistics
out_summary_fp <- file.path(out_dir, "cluster_avg_clr_summary.tsv")
summary_df <- long_df %>%
  group_by(cluster) %>%
  summarise(
    n_lineages = n(),
    mean_clr = mean(avg_clr, na.rm = TRUE),
    sd_clr = sd(avg_clr, na.rm = TRUE),
    min_clr = min(avg_clr, na.rm = TRUE),
    max_clr = max(avg_clr, na.rm = TRUE),
    .groups = "drop"
  )
write_tsv(summary_df, out_summary_fp)
cat("  Written: ", out_summary_fp, "\n")

cat("\n============================================================\n")
cat("  Done! Average CLR computation complete.\n")
cat("============================================================\n")
