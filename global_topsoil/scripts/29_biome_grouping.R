#!/usr/bin/env Rscript
# ============================================================================
# 29_biome_grouping.R
#
# Creates ecologically meaningful biome groupings based on environmental
# similarity. Maps original env_biome categories to broader ecological groups.
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
metadata_fp <- if (length(args) >= 1) args[1] else {
  file.path(project_root, "data", "intermediate", "soil_metadata.tsv")
}
out_dir <- file.path(project_root, "outputs", "biome_groupings")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cat("============================================================\n")
cat("  Biome Grouping Definitions (Global Topsoil)\n")
cat("============================================================\n")
cat("Reading metadata from:", metadata_fp, "\n")
cat("Output directory:", out_dir, "\n\n")

# ============ Read metadata ============
metadata <- read_tsv(metadata_fp, show_col_types = FALSE)
cat("Total samples in metadata:", nrow(metadata), "\n")
cat("Unique biomes found:", length(unique(metadata$env_biome)), "\n\n")

# ============ Define biome groupings ============
# Based on ecological similarity: forest types, grasslands, thermal regimes
biome_groups_topsoil <- list(
  TROPICAL_FOREST = c("Moist_tropical_forests", "Dry_tropical_forests",
                      "Tropical_montane_forests"),

  TEMPERATE_FOREST = c("Temperate_deciduous_forests", "Temperate_coniferous_forests",
                       "Southern_temperate_forests"),

  BOREAL_FOREST = c("Boreal_forests"),

  GRASSLAND_SHRUBLAND = c("Savannas", "Grasslands_and_shrublands"),

  MEDITERRANEAN = c("Mediterrean"),

  ARCTIC_TUNDRA = c("Arctic_tundra"),

  UNKNOWN = c("-", "", "missing: Not provided", "not provided")
)

# ============ Create mapping dataframe ============
cat("Creating biome group mapping...\n")
mapping <- map_dfr(names(biome_groups_topsoil), function(group) {
  data.frame(
    env_biome = biome_groups_topsoil[[group]],
    group_biome = group,
    stringsAsFactors = FALSE
  )
})

# ============ Merge with metadata to get actual counts ============
metadata_with_groups <- metadata %>%
  left_join(mapping, by = "env_biome") %>%
  mutate(group_biome = if_else(is.na(group_biome), "UNKNOWN", group_biome))

# ============ Count samples per group ============
group_counts <- metadata_with_groups %>%
  group_by(group_biome) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  arrange(desc(n_samples))

cat("\nSample counts per biome group:\n")
print(group_counts, n = Inf)

# ============ Count original biomes per group ============
biome_breakdown <- metadata_with_groups %>%
  group_by(group_biome, env_biome) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  arrange(group_biome, desc(n_samples))

# ============ Save outputs ============
write.csv(
  mapping,
  file.path(out_dir, "biome_group_definitions.csv"),
  row.names = FALSE,
  quote = FALSE
)
cat("\nSaved biome group definitions to: biome_group_definitions.csv\n")

write.csv(
  group_counts,
  file.path(out_dir, "biome_group_counts.csv"),
  row.names = FALSE,
  quote = FALSE
)
cat("Saved biome group counts to: biome_group_counts.csv\n")

write_tsv(
  mapping,
  file.path(out_dir, "biome_mapping.tsv")
)
cat("Saved biome mapping to: biome_mapping.tsv\n")

write.csv(
  biome_breakdown,
  file.path(out_dir, "biome_breakdown_by_group.csv"),
  row.names = FALSE,
  quote = FALSE
)
cat("Saved detailed breakdown to: biome_breakdown_by_group.csv\n")

# ============ Summary statistics ============
cat("\n============================================================\n")
cat("  Summary\n")
cat("============================================================\n")
cat("Original biome categories:", length(unique(metadata$env_biome)), "\n")
cat("Grouped biome categories:", nrow(group_counts), "\n")
cat("Reduction ratio:", round(length(unique(metadata$env_biome)) / nrow(group_counts), 2), "x\n")
cat("Samples with unknown mapping:", sum(metadata_with_groups$group_biome == "UNKNOWN"), "\n")
cat("\nBiome grouping complete!\n")
cat("============================================================\n")
