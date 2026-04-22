#!/usr/bin/env Rscript

library(dada2)
library(ggplot2)
library(reshape2)

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

# Set paths
intermediate_dir <- file.path(project_root, "data", "intermediate")
fig_dir <- file.path(project_root, "outputs", "figures")

# Load filtered ASV table and taxonomy
cat(" Loading filtered ASV table and taxonomy...\n")
seqtab.nochim <- readRDS(file.path(intermediate_dir, "seqtab_nochim_bacteria.rds"))
taxa <- readRDS(file.path(intermediate_dir, "taxa_bacteria.rds"))

# Save taxonomy table
write.csv(taxa, file.path(intermediate_dir, "taxonomy_table_bacteria.csv"))

# -------------------------
# Barplot & Sample-Taxa Matrix at Genus Level (Level 6)
# -------------------------

cat(" Building barplot and taxa matrix at Genus level (Bacteria only)...\n")

# Extract Genus-level assignments
tax_level6 <- taxa[,6]
tax_level6[is.na(tax_level6)] <- "Unclassified"
asv_counts <- as.data.frame(seqtab.nochim)
colnames(asv_counts) <- tax_level6

# Collapse identical taxa by summing ASVs
taxa_genus <- as.data.frame(t(sapply(split.default(asv_counts, colnames(asv_counts)), rowSums)))
taxa_genus <- as.data.frame(t(taxa_genus))  # Samples as rows again

# Normalize to relative abundance
rel_abund <- sweep(taxa_genus, 1, rowSums(taxa_genus), "/")

# Save sample-by-genus matrix
write.csv(rel_abund, file.path(intermediate_dir, "sample_taxa_matrix_L6_bacteria.csv"))

# Melt for ggplot
rel_abund$SampleID <- rownames(rel_abund)
df_melt <- melt(rel_abund, id.vars = "SampleID")
colnames(df_melt) <- c("Sample", "Genus", "RelativeAbundance")

# Plot
p <- ggplot(df_melt, aes(x = Sample, y = RelativeAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Genus-level Relative Abundance (Bacteria only)", y = "Proportion", x = "Sample")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(fig_dir, "genus_level_barplot_bacteria.png"), plot = p, width = 12, height = 6)

cat(" Done! Outputs:\n")
cat(" - Taxonomy assignments: taxonomy_table_bacteria.csv\n")
cat(" - Genus-level matrix: sample_taxa_matrix_L6_bacteria.csv\n")
cat(" - Barplot: genus_level_barplot_bacteria.png\n")
