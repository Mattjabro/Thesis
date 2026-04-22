#!/usr/bin/env Rscript
# ============================================================================
# 22_tsne_analysis.R
#
# Performs t-SNE dimensionality reduction on CLR-transformed lineage data.
# Generates 2D embeddings colored by cluster, original biome, and grouped biome.
# Complements existing PCA/PCoA ordination analyses.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Rtsne)
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
metadata_fp <- if (length(args) >= 3) args[3] else {
  file.path(project_root, "data", "intermediate", "soil_metadata.tsv")
}
biome_mapping_fp <- if (length(args) >= 4) args[4] else {
  file.path(project_root, "outputs", "biome_groupings", "biome_mapping.tsv")
}
out_dir <- file.path(project_root, "outputs", "tsne")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("============================================================\n")
cat("  t-SNE Analysis (EMP)\n")
cat("============================================================\n")
cat("CLR matrix:", clr_fp, "\n")
cat("Cluster assignments:", cluster_fp, "\n")
cat("Metadata:", metadata_fp, "\n")
cat("Biome mapping:", biome_mapping_fp, "\n")
cat("Output directory:", out_dir, "\n\n")

# ============ Read data ============
cat("Reading CLR matrix...\n")
clr_matrix <- read.csv(clr_fp, row.names = 1, check.names = FALSE)
cat("  Matrix dimensions:", nrow(clr_matrix), "lineages x", ncol(clr_matrix), "samples\n")

cat("Reading cluster assignments...\n")
clusters <- read_tsv(cluster_fp, show_col_types = FALSE)
cat("  Clusters:", length(unique(clusters$cluster)), "\n")

cat("Reading metadata...\n")
metadata <- read_tsv(metadata_fp, show_col_types = FALSE)
cat("  Metadata rows:", nrow(metadata), "\n")

cat("Reading biome mapping...\n")
biome_mapping <- read_tsv(biome_mapping_fp, show_col_types = FALSE)
cat("  Biome groups:", length(unique(biome_mapping$group_biome)), "\n\n")

# ============ Prepare data for t-SNE ============
# Transpose: samples in rows, taxa in columns
cat("Transposing CLR matrix for t-SNE...\n")
clr_t <- t(clr_matrix)
cat("  Transposed dimensions:", nrow(clr_t), "samples x", ncol(clr_t), "lineages\n\n")

# ============ Run t-SNE ============
cat("Running t-SNE (this may take a while)...\n")
cat("  Parameters: perplexity=30, max_iter=1000\n")

set.seed(42)
tsne_res <- Rtsne(
  clr_t,
  dims = 2,
  perplexity = 30,
  theta = 0.5,
  check_duplicates = FALSE,
  pca = TRUE,
  max_iter = 1000,
  verbose = TRUE
)

cat("t-SNE complete!\n\n")

# ============ Extract coordinates and merge with metadata ============
cat("Merging t-SNE coordinates with metadata...\n")
tsne_coords <- data.frame(
  SampleID = rownames(clr_t),
  tSNE1 = tsne_res$Y[, 1],
  tSNE2 = tsne_res$Y[, 2],
  stringsAsFactors = FALSE
) %>%
  left_join(clusters, by = c("SampleID" = "X.SampleID")) %>%
  left_join(metadata %>% select(X.SampleID, env_biome), by = c("SampleID" = "X.SampleID")) %>%
  left_join(biome_mapping, by = "env_biome") %>%
  mutate(
    cluster = as.factor(cluster),
    group_biome = if_else(is.na(group_biome), "UNKNOWN", group_biome)
  )

cat("  Merged data rows:", nrow(tsne_coords), "\n\n")

# ============ Save outputs ============
cat("Saving outputs...\n")
saveRDS(tsne_res, file.path(out_dir, "tsne_results.rds"))
cat("  Saved: tsne_results.rds\n")

write.csv(
  tsne_coords,
  file.path(out_dir, "tsne_coordinates.csv"),
  row.names = FALSE,
  quote = FALSE
)
cat("  Saved: tsne_coordinates.csv\n\n")

# ============ Create plots ============
cat("Generating t-SNE plots...\n")

# Plot 1: Numbered by cluster
pdf(file.path(out_dir, "tsne_by_cluster.pdf"), width = 12, height = 8)
p1 <- ggplot(tsne_coords, aes(x = tSNE1, y = tSNE2, label = cluster)) +
  geom_text(alpha = 0.7, size = 2.5) +
  theme_bw() +
  labs(
    title = "t-SNE: Samples Numbered by Cluster",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2",
    caption = paste("Cluster labels: 1 -", length(unique(tsne_coords$cluster)))
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0, size = 8)
  )
print(p1)
dev.off()
cat("  Saved: tsne_by_cluster.pdf\n")

# Plot 2: Numbered by original biome
n_biomes <- length(unique(tsne_coords$env_biome[!is.na(tsne_coords$env_biome)]))

# Create numeric codes for original biomes
original_biome_levels <- sort(unique(tsne_coords$env_biome[!is.na(tsne_coords$env_biome)]))
original_biome_code_map <- data.frame(
  env_biome = original_biome_levels,
  biome_code = 1:length(original_biome_levels),
  stringsAsFactors = FALSE
)

tsne_coords_biome_numbered <- tsne_coords %>%
  left_join(original_biome_code_map, by = "env_biome")

pdf(file.path(out_dir, "tsne_by_biome.pdf"), width = 12, height = 8)
p2 <- ggplot(tsne_coords_biome_numbered, aes(x = tSNE1, y = tSNE2, label = biome_code)) +
  geom_text(alpha = 0.7, size = 2.5) +
  theme_bw() +
  labs(
    title = "t-SNE: Samples Numbered by Original Biome",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2",
    caption = paste(
      "Biome codes:\n",
      paste(original_biome_code_map$biome_code, original_biome_code_map$env_biome, sep = " = ", collapse = "\n")
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0, size = 7)
  )
print(p2)
dev.off()
cat("  Saved: tsne_by_biome.pdf\n")

# Save original biome code legend
write.csv(
  original_biome_code_map,
  file.path(out_dir, "original_biome_code_legend.csv"),
  row.names = FALSE,
  quote = FALSE
)
cat("  Saved: original_biome_code_legend.csv\n")

# Plot 3: Numbered by grouped biome
# Create numeric codes for each biome group
biome_levels <- sort(unique(tsne_coords$group_biome))
biome_code_map <- data.frame(
  group_biome = biome_levels,
  biome_code = 1:length(biome_levels),
  stringsAsFactors = FALSE
)

tsne_coords_numbered <- tsne_coords %>%
  left_join(biome_code_map, by = "group_biome")

pdf(file.path(out_dir, "tsne_by_biome_grouped.pdf"), width = 12, height = 8)
p3 <- ggplot(tsne_coords_numbered, aes(x = tSNE1, y = tSNE2, label = biome_code)) +
  geom_text(alpha = 0.7, size = 2.5) +
  theme_bw() +
  labs(
    title = "t-SNE: Samples Numbered by Grouped Biome",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2",
    caption = paste(
      "Biome codes:\n",
      paste(biome_code_map$biome_code, biome_code_map$group_biome, sep = " = ", collapse = "\n")
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0, size = 8)
  )
print(p3)
dev.off()
cat("  Saved: tsne_by_biome_grouped.pdf\n")

# Save biome code legend
write.csv(
  biome_code_map,
  file.path(out_dir, "biome_code_legend.csv"),
  row.names = FALSE,
  quote = FALSE
)
cat("  Saved: biome_code_legend.csv\n")

cat("\n============================================================\n")
cat("  Summary\n")
cat("============================================================\n")
cat("t-SNE completed successfully\n")
cat("Samples embedded:", nrow(tsne_coords), "\n")
cat("Clusters:", length(unique(tsne_coords$cluster)), "\n")
cat("Original biomes:", n_biomes, "\n")
cat("Grouped biomes:", length(unique(tsne_coords$group_biome)), "\n")
cat("\nt-SNE analysis complete!\n")
cat("============================================================\n")
