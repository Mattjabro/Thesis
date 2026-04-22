#!/usr/bin/env Rscript
# ============================================================================
# 26_ordination_analysis.R
#
# Performs PCA and PCoA ordination analysis on CLR-transformed lineage data.
# Generates biplots colored by cluster, scree plots, and loadings analysis.
# Identifies taxa driving principal components.
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
out_dir <- file.path(project_root, "outputs", "ordination")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("============================================================\n")
cat("  PCA/PCoA Ordination Analysis (Global Soil)\n")
cat("============================================================\n")

# ============ Step 1: Load CLR Matrix ============
cat("\n[Step 1] Loading CLR matrix from:\n  ", clr_fp, "\n")
clr_matrix <- as.matrix(read.csv(clr_fp, row.names = 1, check.names = FALSE))
cat("  CLR matrix dimensions: ", nrow(clr_matrix), " taxa x ", ncol(clr_matrix), " samples\n")

# ============ Step 2: Load Cluster Assignments ============
cat("\n[Step 2] Loading cluster assignments from:\n  ", cluster_fp, "\n")
clusters <- read_tsv(cluster_fp, show_col_types = FALSE)
cat("  Clusters: ", n_distinct(clusters$cluster), "\n")

# ============ Step 3: Transpose Matrix for PCA ============
cat("\n[Step 3] Transposing matrix (taxa  samples  samples  taxa)...\n")
clr_t <- t(clr_matrix)
cat("  Transposed dimensions: ", nrow(clr_t), " samples x ", ncol(clr_t), " taxa\n")

# ============ Step 4: Run PCA ============
cat("\n[Step 4] Running PCA using prcomp()...\n")
pca_res <- prcomp(clr_t, center = TRUE, scale. = FALSE)
cat("  PCA complete. Calculating variance explained...\n")

# Variance explained
variance_explained <- data.frame(
  PC = paste0("PC", 1:min(20, length(pca_res$sdev))),
  sdev = pca_res$sdev[1:min(20, length(pca_res$sdev))],
  variance = pca_res$sdev[1:min(20, length(pca_res$sdev))]^2
) %>%
  mutate(
    prop_variance = variance / sum(pca_res$sdev^2),
    cumulative_variance = cumsum(prop_variance),
    pc_index = row_number()
  )

cat("  PC1 explains ", round(variance_explained$prop_variance[1] * 100, 2), "% variance\n")
cat("  PC2 explains ", round(variance_explained$prop_variance[2] * 100, 2), "% variance\n")
cat("  First 10 PCs explain ", round(sum(variance_explained$prop_variance[1:10]) * 100, 2), "% cumulative variance\n")

# Save PCA results
saveRDS(pca_res, file.path(out_dir, "pca_results.rds"))
write_csv(variance_explained, file.path(out_dir, "pca_variance_explained.csv"))

# ============ Step 5: Extract Sample Scores ============
cat("\n[Step 5] Extracting sample scores and merging with clusters...\n")
pca_scores <- as.data.frame(pca_res$x[, 1:min(10, ncol(pca_res$x))]) %>%
  rownames_to_column("SampleID") %>%
  left_join(clusters, by = c("SampleID" = "X.SampleID"))

write_csv(pca_scores, file.path(out_dir, "pca_sample_scores.csv"))
cat("  Written: pca_sample_scores.csv\n")

# For plotting, drop samples without cluster labels
pca_scores_plot <- pca_scores %>% filter(!is.na(cluster))

# ============ Step 6: Extract Loadings ============
cat("\n[Step 6] Extracting taxa loadings...\n")
loadings <- as.data.frame(pca_res$rotation[, 1:min(10, ncol(pca_res$rotation))]) %>%
  rownames_to_column("lineage")

write_csv(loadings, file.path(out_dir, "pca_loadings.csv"))

# Top 50 for PC1 and PC2
top50_pc1 <- loadings %>%
  arrange(desc(abs(PC1))) %>%
  head(50)
write_csv(top50_pc1, file.path(out_dir, "pca_loadings_top50_pc1.csv"))

top50_pc2 <- loadings %>%
  arrange(desc(abs(PC2))) %>%
  head(50)
write_csv(top50_pc2, file.path(out_dir, "pca_loadings_top50_pc2.csv"))
cat("  Written: loadings files\n")

# ============ Step 7: Scree Plot ============
cat("\n[Step 7] Creating scree plot...\n")
pdf(file.path(out_dir, "pca_scree_plot.pdf"), width = 10, height = 6)
p <- ggplot(variance_explained, aes(x = pc_index, y = prop_variance * 100)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_line(aes(x = pc_index, y = cumulative_variance * 100, group = 1), color = "red", size = 1) +
  geom_point(aes(x = pc_index, y = cumulative_variance * 100), color = "red", size = 3) +
  theme_bw() +
  labs(
    title = "PCA Scree Plot - Variance Explained",
    x = "Principal Component",
    y = "Variance Explained (%)",
    subtitle = "Blue bars: individual variance | Red line: cumulative variance"
  ) +
  scale_x_continuous(breaks = variance_explained$pc_index, labels = variance_explained$PC) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()
cat("  Written: pca_scree_plot.pdf\n")

# ============ Step 8: PCA Biplot PC1 vs PC2 ============
cat("\n[Step 8] Creating PCA biplot (PC1 vs PC2)...\n")
pdf(file.path(out_dir, "pca_biplot_pc1_pc2.pdf"), width = 10, height = 8)
p <- ggplot(pca_scores_plot, aes(x = PC1, y = PC2)) +
  geom_text(aes(label = cluster), size = 2, alpha = 0.7, color = "black") +
  stat_ellipse(aes(color = factor(cluster)), level = 0.95, linetype = 2, size = 0.8) +
  scale_color_manual(
    values = c(
      "#000000", "#F3C300", "#875692", "#F38400", "#A1CAF1",
      "#BE0032", "#C2B280", "#848482", "#008856", "#E68FAC",
      "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C",
      "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822"
    ),
    name = "Cluster"
  ) +
  theme_bw() +
  labs(
    title = sprintf("PCA Biplot: PC1 (%.1f%%) vs PC2 (%.1f%%)",
                    variance_explained$prop_variance[1] * 100,
                    variance_explained$prop_variance[2] * 100),
    x = sprintf("PC1 (%.1f%% variance)", variance_explained$prop_variance[1] * 100),
    y = sprintf("PC2 (%.1f%% variance)", variance_explained$prop_variance[2] * 100)
  ) +
  theme(legend.position = "right")
print(p)
dev.off()
cat("  Written: pca_biplot_pc1_pc2.pdf\n")

# ============ Step 9: Additional Biplots ============
cat("\n[Step 9] Creating additional biplots...\n")

# PC1 vs PC3
pdf(file.path(out_dir, "pca_biplot_pc1_pc3.pdf"), width = 10, height = 8)
p <- ggplot(pca_scores_plot, aes(x = PC1, y = PC3)) +
  geom_text(aes(label = cluster), size = 2, alpha = 0.7, color = "black") +
  stat_ellipse(aes(color = factor(cluster)), level = 0.95, linetype = 2, size = 0.8) +
  scale_color_manual(
    values = c(
      "#000000", "#F3C300", "#875692", "#F38400", "#A1CAF1",
      "#BE0032", "#C2B280", "#848482", "#008856", "#E68FAC",
      "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C",
      "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822"
    ),
    name = "Cluster"
  ) +
  theme_bw() +
  labs(
    title = sprintf("PCA Biplot: PC1 (%.1f%%) vs PC3 (%.1f%%)",
                    variance_explained$prop_variance[1] * 100,
                    variance_explained$prop_variance[3] * 100),
    x = sprintf("PC1 (%.1f%% variance)", variance_explained$prop_variance[1] * 100),
    y = sprintf("PC3 (%.1f%% variance)", variance_explained$prop_variance[3] * 100)
  )
print(p)
dev.off()

# PC2 vs PC3
pdf(file.path(out_dir, "pca_biplot_pc2_pc3.pdf"), width = 10, height = 8)
p <- ggplot(pca_scores_plot, aes(x = PC2, y = PC3)) +
  geom_text(aes(label = cluster), size = 2, alpha = 0.7, color = "black") +
  stat_ellipse(aes(color = factor(cluster)), level = 0.95, linetype = 2, size = 0.8) +
  scale_color_manual(
    values = c(
      "#000000", "#F3C300", "#875692", "#F38400", "#A1CAF1",
      "#BE0032", "#C2B280", "#848482", "#008856", "#E68FAC",
      "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C",
      "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822"
    ),
    name = "Cluster"
  ) +
  theme_bw() +
  labs(
    title = sprintf("PCA Biplot: PC2 (%.1f%%) vs PC3 (%.1f%%)",
                    variance_explained$prop_variance[2] * 100,
                    variance_explained$prop_variance[3] * 100),
    x = sprintf("PC2 (%.1f%% variance)", variance_explained$prop_variance[2] * 100),
    y = sprintf("PC3 (%.1f%% variance)", variance_explained$prop_variance[3] * 100)
  )
print(p)
dev.off()
cat("  Written: additional biplots\n")

# ============ Step 10: Loadings Bar Plots ============
cat("\n[Step 10] Creating loadings bar plots...\n")

# Top 30 taxa on PC1
top30_pc1 <- loadings %>%
  arrange(desc(abs(PC1))) %>%
  head(30) %>%
  mutate(
    lineage_short = if_else(
      nchar(lineage) > 60,
      paste0(substr(lineage, 1, 55), "", sprintf("%02d", row_number())),
      lineage
    )
  )

pdf(file.path(out_dir, "pca_loadings_pc1_top30.pdf"), width = 10, height = 12)
p <- ggplot(top30_pc1, aes(x = reorder(lineage_short, PC1), y = PC1)) +
  geom_col(aes(fill = PC1 > 0)) +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "coral"), guide = "none") +
  theme_bw() +
  labs(
    title = "Top 30 Taxa Loadings on PC1",
    subtitle = "Blue = negative loading, Red = positive loading",
    x = "Taxon",
    y = "Loading on PC1"
  ) +
  theme(axis.text.y = element_text(size = 8))
print(p)
dev.off()

# Top 30 taxa on PC2
top30_pc2 <- loadings %>%
  arrange(desc(abs(PC2))) %>%
  head(30) %>%
  mutate(
    lineage_short = if_else(
      nchar(lineage) > 60,
      paste0(substr(lineage, 1, 55), "", sprintf("%02d", row_number())),
      lineage
    )
  )

pdf(file.path(out_dir, "pca_loadings_pc2_top30.pdf"), width = 10, height = 12)
p <- ggplot(top30_pc2, aes(x = reorder(lineage_short, PC2), y = PC2)) +
  geom_col(aes(fill = PC2 > 0)) +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "coral"), guide = "none") +
  theme_bw() +
  labs(
    title = "Top 30 Taxa Loadings on PC2",
    subtitle = "Blue = negative loading, Red = positive loading",
    x = "Taxon",
    y = "Loading on PC2"
  ) +
  theme(axis.text.y = element_text(size = 8))
print(p)
dev.off()
cat("  Written: loadings bar plots\n")

# ============ Step 11: Loadings Biplot with Arrows ============
cat("\n[Step 11] Creating loadings biplot with taxa arrows...\n")

# Select top 20 most influential taxa (based on length of loading vector)
top20_influential <- loadings %>%
  mutate(loading_magnitude = sqrt(PC1^2 + PC2^2)) %>%
  arrange(desc(loading_magnitude)) %>%
  head(20) %>%
  mutate(
    lineage_short = if_else(
      nchar(lineage) > 40,
      paste0(substr(lineage, 1, 35), "", sprintf("%02d", row_number())),
      lineage
    )
  )

pdf(file.path(out_dir, "pca_loadings_biplot.pdf"), width = 12, height = 10)
p <- ggplot(pca_scores_plot, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = factor(cluster)), size = 1.5, alpha = 0.3) +
  scale_color_brewer(palette = "Set3", name = "Cluster") +
  geom_segment(
    data = top20_influential,
    aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5),
    arrow = arrow(length = unit(0.3, "cm")),
    color = "black",
    size = 0.5,
    alpha = 0.7
  ) +
  geom_text(
    data = top20_influential,
    aes(x = PC1 * 5.5, y = PC2 * 5.5, label = lineage_short),
    size = 2.5,
    check_overlap = TRUE
  ) +
  theme_bw() +
  labs(
    title = "PCA Biplot with Taxa Loadings (Top 20)",
    x = sprintf("PC1 (%.1f%% variance)", variance_explained$prop_variance[1] * 100),
    y = sprintf("PC2 (%.1f%% variance)", variance_explained$prop_variance[2] * 100)
  )
print(p)
dev.off()
cat("  Written: pca_loadings_biplot.pdf\n")

# ============ Step 12: PCoA using cmdscale ============
cat("\n[Step 12] Running PCoA (classical MDS)...\n")
dist_matrix <- dist(clr_t, method = "euclidean")
pcoa_res <- cmdscale(dist_matrix, k = 10, eig = TRUE)
cat("  PCoA complete\n")

# Calculate variance explained for PCoA
pcoa_eigenvalues <- pcoa_res$eig[pcoa_res$eig > 0]
pcoa_var_explained <- pcoa_eigenvalues / sum(pcoa_eigenvalues) * 100

saveRDS(pcoa_res, file.path(out_dir, "pcoa_results.rds"))

# PCoA scores
pcoa_scores <- as.data.frame(pcoa_res$points) %>%
  setNames(paste0("Axis", 1:ncol(.))) %>%
  rownames_to_column("SampleID") %>%
  left_join(clusters, by = c("SampleID" = "X.SampleID"))

write_csv(pcoa_scores, file.path(out_dir, "pcoa_sample_scores.csv"))
cat("  Written: pcoa_sample_scores.csv\n")

# For plotting, drop samples without cluster labels
pcoa_scores_plot <- pcoa_scores %>% filter(!is.na(cluster))

# ============ Step 13: PCoA Biplot ============
cat("\n[Step 13] Creating PCoA biplot...\n")
pdf(file.path(out_dir, "pcoa_biplot_axis1_axis2.pdf"), width = 10, height = 8)
p <- ggplot(pcoa_scores_plot, aes(x = Axis1, y = Axis2)) +
  geom_text(aes(label = cluster), size = 2, alpha = 0.7, color = "black") +
  stat_ellipse(aes(color = factor(cluster)), level = 0.95, linetype = 2, size = 0.8) +
  scale_color_manual(
    values = c(
      "#000000", "#F3C300", "#875692", "#F38400", "#A1CAF1",
      "#BE0032", "#C2B280", "#848482", "#008856", "#E68FAC",
      "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C",
      "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822"
    ),
    name = "Cluster"
  ) +
  theme_bw() +
  labs(
    title = sprintf("PCoA Biplot: Axis1 (%.1f%%) vs Axis2 (%.1f%%)",
                    pcoa_var_explained[1], pcoa_var_explained[2]),
    x = sprintf("Axis 1 (%.1f%% variance)", pcoa_var_explained[1]),
    y = sprintf("Axis 2 (%.1f%% variance)", pcoa_var_explained[2]),
    subtitle = "Euclidean distance on CLR-transformed data"
  )
print(p)
dev.off()
cat("  Written: pcoa_biplot_axis1_axis2.pdf\n")

cat("\n============================================================\n")
cat("  Ordination analysis complete!\n")
cat("  Created 8 PDF plots and 6 data files\n")
cat("============================================================\n")
