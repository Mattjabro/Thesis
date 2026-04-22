#!/usr/bin/env Rscript
# ============================================================================
# p02_extfig2bc_diversity_by_biome.R
#
# Recreates Extended Data Figure 2b/c from Bahram et al. 2018 (Nature 560):
#   Bacterial OTU richness (b) and taxonomic diversity/inverse Simpson (c)
#   across major biome categories: Tropical / Temperate / Boreal-Arctic.
#
#   Shows mean  SD bars with pairwise significance letters (Kruskal-Wallis
#   + Wilcoxon, Benjamini-Hochberg corrected), matching the paper's style.
#
# Input:
#   data/intermediate/seqtab_nochim_bacteria.rds
#   data/intermediate/soil_metadata.tsv
#
# Output:
#   outputs/paper_figures/extfig2bc_diversity_by_biome.pdf
#   outputs/paper_figures/extfig2bc_diversity_by_biome.png
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg  <- "--file="
  match     <- grep(file_arg, cmd_args)
  if (length(match) > 0)
    return(normalizePath(sub(file_arg, "", cmd_args[match[1]])))
  if (!is.null(sys.frames()[[1]]$ofile))
    return(normalizePath(sys.frames()[[1]]$ofile))
  return(NA_character_)
}

script_path  <- get_script_path()
project_root <- if (!is.na(script_path)) {
  normalizePath(file.path(dirname(script_path), ".."))
} else {
  normalizePath(getwd())
}

# ---- Paths ------------------------------------------------------------------
seqtab_fp <- file.path(project_root, "data", "intermediate", "seqtab_nochim_bacteria.rds")
meta_fp   <- file.path(project_root, "data", "intermediate", "soil_metadata.tsv")
out_dir   <- file.path(project_root, "outputs", "paper_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("  Ext Fig 2b/c: Diversity by Biome (Bahram 2018)\n")
cat("============================================================\n\n")

# ---- Step 1: Load & rarefy --------------------------------------------------
cat("[1] Loading ASV table...\n")
seqtab <- readRDS(seqtab_fp)
cat("    Dimensions:", nrow(seqtab), "samples x", ncol(seqtab), "ASVs\n")

cat("[2] Loading metadata...\n")
biome_to_hl <- c(
  Moist_tropical_forests       = "Tropical",
  Tropical_montane_forests     = "Tropical",
  Dry_tropical_forests         = "Tropical",
  Savannas                     = "Tropical",
  Temperate_coniferous_forests = "Temperate",
  Grasslands_and_shrublands    = "Temperate",
  Southern_temperate_forests   = "Temperate",
  Temperate_deciduous_forests  = "Temperate",
  Mediterrean                  = "Temperate",
  Boreal_forests               = "Boreal-Arctic",
  Arctic_tundra                = "Boreal-Arctic"
)

meta <- read_tsv(meta_fp, show_col_types = FALSE) %>%
  rename(SampleID = X.SampleID, biome = env_biome) %>%
  filter(biome != "-", !is.na(biome)) %>%
  mutate(
    hl_biome = biome_to_hl[biome],
    hl_biome = factor(hl_biome, levels = c("Tropical", "Temperate", "Boreal-Arctic"))
  ) %>%
  filter(!is.na(hl_biome))

cat("    Samples with valid biome:", nrow(meta), "\n")
print(table(meta$hl_biome))
cat("\n")

shared_ids  <- intersect(rownames(seqtab), meta$SampleID)
seqtab_sub  <- seqtab[shared_ids, ]
cat("[3] Rarefying", length(shared_ids), "shared samples...\n")
min_depth   <- min(rowSums(seqtab_sub))
cat("    Rarefaction depth:", min_depth, "\n")

set.seed(42)
seqtab_rare <- rrarefy(seqtab_sub, sample = min_depth)

# ---- Step 2: Compute diversity metrics --------------------------------------
cat("[4] Computing observed OTUs and inverse Simpson...\n")

div_df <- tibble(
  SampleID    = rownames(seqtab_rare),
  obs_otus    = specnumber(seqtab_rare),
  inv_simpson = diversity(seqtab_rare, index = "invsimpson")
) %>%
  inner_join(meta %>% select(SampleID, biome, hl_biome), by = "SampleID")

cat("    n per high-level biome:\n")
print(div_df %>% count(hl_biome))
cat("\n")

# ---- Step 3: Scale diversity metrics (matching paper's "Scaled" axes) -------
# The paper scales diversity values so different groups are comparable on the
# same axis. We z-score within each metric.
div_df <- div_df %>%
  mutate(
    obs_otus_scaled    = scale(obs_otus)[,1],
    inv_simpson_scaled = scale(inv_simpson)[,1]
  )

# ---- Step 4: Statistical tests ----------------------------------------------
cat("[5] Running Kruskal-Wallis + pairwise Wilcoxon tests...\n")

kw_otus <- kruskal.test(obs_otus ~ hl_biome, data = div_df)
kw_div  <- kruskal.test(inv_simpson ~ hl_biome, data = div_df)
cat(sprintf("    KW OTU richness : chi=%.2f, df=%d, P=%.4g\n",
            kw_otus$statistic, kw_otus$parameter, kw_otus$p.value))
cat(sprintf("    KW inv. Simpson : chi=%.2f, df=%d, P=%.4g\n",
            kw_div$statistic,  kw_div$parameter,  kw_div$p.value))

pw_otus <- pairwise.wilcox.test(div_df$obs_otus,    div_df$hl_biome, p.adjust.method = "BH")
pw_div  <- pairwise.wilcox.test(div_df$inv_simpson, div_df$hl_biome, p.adjust.method = "BH")

# Assign compact letter display (CLD) matching the paper's "a, b, c" notation
# We use a simple custom function since multcompView may not be installed
assign_letters <- function(pw_result) {
  groups  <- levels(div_df$hl_biome)
  n_grp   <- length(groups)
  # Start all groups with same letter; split when p < 0.05
  letters_out <- rep("a", n_grp)
  names(letters_out) <- groups
  pmat <- pw_result$p.value
  letter_pool <- letters

  # Build adjacency: groups that are NOT significantly different share a letter
  used_letter <- 1L
  for (i in seq_len(n_grp - 1)) {
    for (j in (i + 1):n_grp) {
      gi <- groups[i]; gj <- groups[j]
      p_ij <- if (gi %in% rownames(pmat) && gj %in% colnames(pmat)) {
        pmat[gi, gj]
      } else if (gj %in% rownames(pmat) && gi %in% colnames(pmat)) {
        pmat[gj, gi]
      } else NA_real_
      if (!is.na(p_ij) && p_ij < 0.05) {
        # Significant: ensure different letters
        if (letters_out[gi] == letters_out[gj]) {
          used_letter <- used_letter + 1L
          letters_out[gj] <- letter_pool[used_letter]
        }
      }
    }
  }
  letters_out
}

cld_otus <- assign_letters(pw_otus)
cld_div  <- assign_letters(pw_div)
cat("    OTU letter groups: "); print(cld_otus)
cat("    Div letter groups: "); print(cld_div)
cat("\n")

# ---- Step 5: Summary stats for bar plot -------------------------------------
summ_otus <- div_df %>%
  group_by(hl_biome) %>%
  summarise(
    mean_val = mean(obs_otus_scaled),
    sd_val   = sd(obs_otus_scaled),
    se_val   = sd_val / sqrt(n()),
    n        = n(),
    .groups  = "drop"
  ) %>%
  mutate(
    letter   = cld_otus[as.character(hl_biome)],
    metric   = "Bacterial OTU"
  )

summ_div <- div_df %>%
  group_by(hl_biome) %>%
  summarise(
    mean_val = mean(inv_simpson_scaled),
    sd_val   = sd(inv_simpson_scaled),
    se_val   = sd_val / sqrt(n()),
    n        = n(),
    .groups  = "drop"
  ) %>%
  mutate(
    letter   = cld_div[as.character(hl_biome)],
    metric   = "Bacterial Diversity\n(inverse Simpson)"
  )

summ_all <- bind_rows(summ_otus, summ_div) %>%
  mutate(metric = factor(metric, levels = c("Bacterial OTU",
                                            "Bacterial Diversity\n(inverse Simpson)")))

# ---- Step 6: Plot -----------------------------------------------------------
cat("[6] Building plots...\n")

# Biome colours matching paper (blue=tropical, orange=temperate, red=boreal)
hl_colours <- c(
  Tropical      = "#1565C0",
  Temperate     = "#FB8C00",
  `Boreal-Arctic` = "#B71C1C"
)

# Shared theme
base_theme <- theme_classic(base_size = 11) +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 10),
    axis.title        = element_text(size = 10),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    legend.position   = "none",
    plot.title        = element_text(face = "bold", size = 10)
  )

# --- Panel b: OTU richness ---
p_b <- ggplot(
    summ_otus,
    aes(x = hl_biome, y = mean_val, fill = hl_biome)
  ) +
  geom_col(width = 0.55, colour = "grey30", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.2, linewidth = 0.5
  ) +
  geom_text(
    aes(y = mean_val + sd_val + 0.05, label = letter),
    vjust = 0, size = 3.5, fontface = "bold"
  ) +
  geom_text(
    aes(y = -Inf, label = paste0("n=", n)),
    vjust = -0.4, size = 2.8, colour = "grey40"
  ) +
  scale_fill_manual(values = hl_colours) +
  scale_x_discrete(labels = c("Tropical", "Temperate", "Boreal-\nArctic")) +
  labs(
    x     = NULL,
    y     = "Scaled richness",
    title = sprintf("b  Richness (KW P = %.3g)", kw_otus$p.value)
  ) +
  base_theme

# --- Panel c: Inverse Simpson diversity ---
p_c <- ggplot(
    summ_div,
    aes(x = hl_biome, y = mean_val, fill = hl_biome)
  ) +
  geom_col(width = 0.55, colour = "grey30", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.2, linewidth = 0.5
  ) +
  geom_text(
    aes(y = mean_val + sd_val + 0.05, label = letter),
    vjust = 0, size = 3.5, fontface = "bold"
  ) +
  geom_text(
    aes(y = -Inf, label = paste0("n=", n)),
    vjust = -0.4, size = 2.8, colour = "grey40"
  ) +
  scale_fill_manual(values = hl_colours) +
  scale_x_discrete(labels = c("Tropical", "Temperate", "Boreal-\nArctic")) +
  labs(
    x     = NULL,
    y     = "Scaled diversity",
    title = sprintf("c  Diversity (KW P = %.3g)", kw_div$p.value)
  ) +
  base_theme

# --- Individual scatter panels (like paper's jitter) ---
scatter_theme <- base_theme +
  theme(legend.position = "none")

scatter_otus <- ggplot(div_df, aes(x = hl_biome, y = obs_otus_scaled, colour = hl_biome)) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.6) +
  geom_crossbar(
    data = summ_otus,
    aes(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.4, fatten = 2, linewidth = 0.5
  ) +
  geom_text(
    data = summ_otus,
    aes(y = mean_val + sd_val + 0.15, label = letter),
    vjust = 0, size = 3.5, fontface = "bold"
  ) +
  scale_colour_manual(values = hl_colours) +
  scale_x_discrete(labels = c("Tropical", "Temperate", "Boreal-\nArctic")) +
  labs(x = NULL, y = "Scaled richness (OTUs)",
       title = sprintf("b  Bacterial OTU Richness  (KW P = %.3g)", kw_otus$p.value)) +
  scatter_theme

scatter_div <- ggplot(div_df, aes(x = hl_biome, y = inv_simpson_scaled, colour = hl_biome)) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.6) +
  geom_crossbar(
    data = summ_div,
    aes(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    width = 0.4, fatten = 2, linewidth = 0.5
  ) +
  geom_text(
    data = summ_div,
    aes(y = mean_val + sd_val + 0.15, label = letter),
    vjust = 0, size = 3.5, fontface = "bold"
  ) +
  scale_colour_manual(values = hl_colours) +
  scale_x_discrete(labels = c("Tropical", "Temperate", "Boreal-\nArctic")) +
  labs(x = NULL, y = "Scaled diversity (inv. Simpson)",
       title = sprintf("c  Bacterial Diversity  (KW P = %.3g)", kw_div$p.value)) +
  scatter_theme

# Combine b + c side by side
# Use patchwork if available, otherwise cowplot, otherwise save separately
combine_plots <- function(p1, p2, fp_pdf, fp_png) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    combined <- p1 + p2 + plot_layout(ncol = 2) +
      plot_annotation(
        title    = "Extended Data Fig. 2b/c  Bacterial diversity across biomes",
        subtitle = sprintf(
          "n = %d samples | Mean  SD | Letters = pairwise Wilcoxon BH-corrected (P<0.05)",
          nrow(div_df)
        ),
        theme = theme(
          plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, colour = "grey40")
        )
      )
    pdf(fp_pdf, width = 9, height = 5)
    print(combined)
    dev.off()
    png(fp_png, width = 1800, height = 1000, res = 200)
    print(combined)
    dev.off()
  } else {
    # Save separately
    pdf(sub(".pdf", "_b.pdf", fp_pdf), width = 5, height = 5); print(p1); dev.off()
    pdf(sub(".pdf", "_c.pdf", fp_pdf), width = 5, height = 5); print(p2); dev.off()
    png(sub(".png", "_b.png", fp_png), width = 900, height = 1000, res = 200); print(p1); dev.off()
    png(sub(".png", "_c.png", fp_png), width = 900, height = 1000, res = 200); print(p2); dev.off()
    cat("    (patchwork not available; saved panels separately)\n")
  }
}

pdf_fp <- file.path(out_dir, "extfig2bc_diversity_by_biome.pdf")
png_fp <- file.path(out_dir, "extfig2bc_diversity_by_biome.png")

cat("[7] Saving outputs...\n")
combine_plots(scatter_otus, scatter_div, pdf_fp, png_fp)
cat("    Written:", pdf_fp, "\n")
cat("    Written:", png_fp, "\n")

# ---- Summary ----------------------------------------------------------------
cat("\n============================================================\n")
cat("  Summary\n")
cat("============================================================\n")
cat(sprintf("  n samples plotted : %d\n", nrow(div_df)))
cat(sprintf("  Rarefaction depth : %d reads\n", min_depth))
cat(sprintf("  KW OTU richness   : P = %.4g\n", kw_otus$p.value))
cat(sprintf("  KW inv. Simpson   : P = %.4g\n", kw_div$p.value))
cat("  Pairwise p-values (OTU richness, BH corrected):\n")
print(pw_otus$p.value)
cat("  Pairwise p-values (inv. Simpson, BH corrected):\n")
print(pw_div$p.value)
cat("============================================================\n")
