#!/usr/bin/env Rscript
# ============================================================================
# p01_fig1a_diversity_vs_latitude.R
#
# Recreates Figure 1a from Bahram et al. 2018 (Nature 560:233-237):
#   Bacterial taxonomic diversity (inverse Simpson) vs. absolute latitude,
#   with 1st- and 2nd-order polynomial regression fits.
#
# Input:
#   data/intermediate/seqtab_nochim_bacteria.rds  (samples x ASVs count matrix)
#   data/intermediate/soil_metadata.tsv           (SampleID, env_biome, lat, lon)
#
# Output:
#   outputs/paper_figures/fig1a_diversity_vs_latitude.pdf
#   outputs/paper_figures/fig1a_diversity_vs_latitude.png
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
seqtab_fp   <- file.path(project_root, "data", "intermediate", "seqtab_nochim_bacteria.rds")
meta_fp     <- file.path(project_root, "data", "intermediate", "soil_metadata.tsv")
out_dir     <- file.path(project_root, "outputs", "paper_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("  Fig 1a: Bacterial Diversity vs. Latitude (Bahram 2018)\n")
cat("============================================================\n\n")

# ---- Step 1: Load data ------------------------------------------------------
cat("[1] Loading ASV table...\n")
seqtab <- readRDS(seqtab_fp)
cat("    Dimensions:", nrow(seqtab), "samples x", ncol(seqtab), "ASVs\n")

cat("[2] Loading metadata...\n")
meta <- read_tsv(meta_fp, show_col_types = FALSE) %>%
  rename(SampleID = X.SampleID, biome = env_biome,
         lat = latitude_deg, lon = longitude_deg) %>%
  filter(biome != "-", !is.na(lat)) %>%
  mutate(abs_lat = abs(lat))
cat("    Samples with valid biome + lat:", nrow(meta), "\n")

# ---- Step 2: Rarefy & compute inverse Simpson --------------------------------
cat("[3] Rarefying to minimum depth and computing inverse Simpson...\n")

# Keep only samples present in both seqtab and metadata
shared_ids <- intersect(rownames(seqtab), meta$SampleID)
seqtab_sub <- seqtab[shared_ids, ]
cat("    Shared samples:", length(shared_ids), "\n")

min_depth <- min(rowSums(seqtab_sub))
cat("    Rarefaction depth:", min_depth, "reads\n")

set.seed(42)
seqtab_rare <- rrarefy(seqtab_sub, sample = min_depth)

inv_simp <- diversity(seqtab_rare, index = "invsimpson")
obs_otus <- specnumber(seqtab_rare)

div_df <- tibble(
  SampleID = names(inv_simp),
  inv_simpson = inv_simp,
  obs_otus    = obs_otus
) %>%
  inner_join(meta %>% select(SampleID, biome, abs_lat, lat, lon), by = "SampleID")

cat("    Final n for plotting:", nrow(div_df), "\n\n")

# ---- Step 3: AICc for polynomial fits ---------------------------------------
cat("[4] Fitting 1st- and 2nd-order polynomial models...\n")

aic_c <- function(fit, n) {
  k   <- length(coef(fit)) + 1   # +1 for sigma
  aic <- AIC(fit)
  aic + (2 * k * (k + 1)) / (n - k - 1)
}

n   <- nrow(div_df)
lm1 <- lm(inv_simpson ~ abs_lat,              data = div_df)
lm2 <- lm(inv_simpson ~ abs_lat + I(abs_lat^2), data = div_df)

r2_1  <- summary(lm1)$r.squared
r2_2  <- summary(lm2)$r.squared
p1    <- summary(lm1)$coefficients["abs_lat", "Pr(>|t|)"]
p2    <- anova(lm1, lm2)$`Pr(>F)`[2]   # test if quadratic term improves fit
aicc1 <- aic_c(lm1, n)
aicc2 <- aic_c(lm2, n)

cat(sprintf("    Linear:    r=%.3f  AICc=%.1f\n", r2_1, aicc1))
cat(sprintf("    Quadratic: r=%.3f  AICc=%.1f\n\n", r2_2, aicc2))

best_model <- if (aicc2 < aicc1) "quadratic" else "linear"
cat("    Best model (lower AICc):", best_model, "\n\n")

# ---- Step 4: Biome colours (matches paper legend) ---------------------------
# Biome  high-level category for the top strip
biome_levels <- c(
  "Moist_tropical_forests", "Tropical_montane_forests",
  "Dry_tropical_forests", "Savannas",
  "Temperate_coniferous_forests", "Grasslands_and_shrublands",
  "Southern_temperate_forests", "Temperate_deciduous_forests",
  "Mediterrean",
  "Boreal_forests", "Arctic_tundra"
)

biome_colours <- c(
  Moist_tropical_forests       = "#1565C0",
  Tropical_montane_forests     = "#00897B",
  Dry_tropical_forests         = "#80DEEA",
  Savannas                     = "#FF8F00",
  Temperate_coniferous_forests = "#FDD835",
  Grasslands_and_shrublands    = "#F48FB1",
  Southern_temperate_forests   = "#FB8C00",
  Temperate_deciduous_forests  = "#BF360C",
  Mediterrean                  = "#6A1A4C",
  Boreal_forests               = "#B71C1C",
  Arctic_tundra                = "#37474F"
)

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

hl_colours <- c(
  Tropical        = "#1565C0",
  Temperate       = "#FB8C00",
  `Boreal-Arctic` = "#B71C1C"
)

div_df <- div_df %>%
  mutate(
    biome    = factor(biome, levels = biome_levels),
    hl_biome = factor(biome_to_hl[as.character(biome)],
                      levels = c("Tropical", "Temperate", "Boreal-Arctic"))
  )

# ---- Step 5: Regression curve data -----------------------------------------
lat_seq  <- seq(0, 80, length.out = 300)
pred_df  <- tibble(abs_lat = lat_seq)
pred_df$fit1 <- predict(lm1, newdata = pred_df)
pred_df$fit2 <- predict(lm2, newdata = pred_df)

# ---- Step 6: Build plot -----------------------------------------------------
cat("[5] Building plot...\n")

# Format p-values for annotations (match paper style)
fmt_p <- function(p) {
  if (p < 1e-10) "P < 1e-10"
  else if (p < 0.001) sprintf("P = %.2e", p)
  else sprintf("P = %.3f", p)
}

# AICc: underline best model by making it bold in the label
ann1 <- sprintf("r = %.3f\n%s\nAICc = %.0f", r2_1, fmt_p(p1), aicc1)
ann2 <- sprintf("r = %.3f\nAICc = %.0f", r2_2, aicc2)

# Top annotation: show model stats in two clusters like the paper
ann_x <- 2
ann_y_top <- max(div_df$inv_simpson, na.rm = TRUE) * 0.98

p_main <- ggplot(div_df, aes(x = abs_lat, y = inv_simpson)) +

  # Top biome-zone bands (Tropical / Temperate / Boreal-Arctic), like paper header
  annotate("rect", xmin =  0, xmax = 23, ymin = 420, ymax = 450,
           fill = hl_colours["Tropical"],   alpha = 0.25) +
  annotate("rect", xmin = 23, xmax = 55, ymin = 420, ymax = 450,
           fill = hl_colours["Temperate"],  alpha = 0.25) +
  annotate("rect", xmin = 55, xmax = 80, ymin = 420, ymax = 450,
           fill = hl_colours[["Boreal-Arctic"]], alpha = 0.25) +
  annotate("text", x = 11.5, y = 435, label = "Tropical",
           size = 2.5, fontface = "bold", colour = hl_colours["Tropical"]) +
  annotate("text", x = 39,   y = 435, label = "Temperate",
           size = 2.5, fontface = "bold", colour = hl_colours["Temperate"]) +
  annotate("text", x = 67.5, y = 435, label = "BorealArctic",
           size = 2.5, fontface = "bold", colour = hl_colours[["Boreal-Arctic"]]) +

  # Points coloured by high-level biome (3 categories)
  geom_point(aes(colour = hl_biome), size = 2.2, alpha = 0.85) +

  # 1st-order fit (grey dashed, like paper)
  geom_line(data = pred_df, aes(x = abs_lat, y = fit1),
            colour = "grey60", linetype = "dashed", linewidth = 0.8) +

  # 2nd-order fit (black solid, like paper)
  geom_line(data = pred_df, aes(x = abs_lat, y = fit2),
            colour = "black", linetype = "solid", linewidth = 0.9) +

  scale_colour_manual(
    values = hl_colours,
    name   = NULL
  ) +

  # Stats annotations: linear (upper-left) and quadratic (beside it)
  annotate("text", x = 2, y = 390,
           label = sprintf("r = %.3f\n%s\nAICc = %.0f", r2_1, fmt_p(p1), aicc1),
           hjust = 0, vjust = 1, size = 2.8, colour = "grey50",
           fontface = "plain") +
  annotate("text", x = 22, y = 390,
           label = sprintf("r = %.3f\n%s\nAICc = %.0f", r2_2,
                           fmt_p(summary(lm2)$coefficients[3, 4]), aicc2),
           hjust = 0, vjust = 1, size = 2.8, colour = "black",
           fontface = "bold") +

  scale_x_continuous(breaks = seq(0, 80, 20), limits = c(0, 80)) +
  scale_y_continuous(limits = c(0, 450)) +
  labs(
    x = "Absolute latitude (degrees)",
    y = "Taxonomic diversity\n(inverse Simpson)",
    title = "Fig. 1a  Bacterial taxonomic diversity vs. latitude",
    subtitle = sprintf(
      "n = %d samples | Grey dashed = linear fit | Black solid = quadratic fit",
      nrow(div_df)
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.text       = element_text(size = 9),
    legend.key.size   = unit(0.4, "cm"),
    legend.title      = element_blank(),
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(size = 9, colour = "grey40"),
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
    axis.title        = element_text(size = 10)
  ) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3)))

# ---- Step 7: Save -----------------------------------------------------------
cat("[6] Saving outputs...\n")

pdf_fp <- file.path(out_dir, "fig1a_diversity_vs_latitude.pdf")
png_fp <- file.path(out_dir, "fig1a_diversity_vs_latitude.png")

pdf(pdf_fp, width = 7, height = 5.5)
print(p_main)
dev.off()
cat("    Written:", pdf_fp, "\n")

png(png_fp, width = 1400, height = 1100, res = 200)
print(p_main)
dev.off()
cat("    Written:", png_fp, "\n")

# ---- Summary ----------------------------------------------------------------
cat("\n============================================================\n")
cat("  Summary\n")
cat("============================================================\n")
cat(sprintf("  n samples plotted : %d\n", nrow(div_df)))
cat(sprintf("  Rarefaction depth : %d reads\n", min_depth))
cat(sprintf("  Inv. Simpson range: %.1f  %.1f\n",
            min(div_df$inv_simpson), max(div_df$inv_simpson)))
cat(sprintf("  Linear    r=%.3f  AICc=%.1f\n", r2_1, aicc1))
cat(sprintf("  Quadratic r=%.3f  AICc=%.1f  (paper reports r=0.160, AICc=2180)\n",
            r2_2, aicc2))
cat("============================================================\n")
