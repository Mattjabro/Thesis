#!/usr/bin/env Rscript
# ============================================================================
# 25_random_forest.R
#
# Trains Random Forest classification models to predict biome categories
# from CLR-transformed taxa abundances. Creates models for both original
# and grouped biome categories. Uses 70/30 train/test split.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(randomForest)
  library(caret)
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
metadata_fp <- if (length(args) >= 2) args[2] else {
  file.path(project_root, "data", "intermediate", "soil_metadata.tsv")
}
biome_mapping_fp <- if (length(args) >= 3) args[3] else {
  file.path(project_root, "outputs", "biome_groupings", "biome_mapping.tsv")
}
out_dir <- file.path(project_root, "outputs", "ml_models")

# Create output directories
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "rf_original_biomes"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "rf_grouped_biomes"), recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("  Random Forest Classification (EMP)\n")
cat("============================================================\n")
cat("CLR matrix:", clr_fp, "\n")
cat("Metadata:", metadata_fp, "\n")
cat("Biome mapping:", biome_mapping_fp, "\n")
cat("Output directory:", out_dir, "\n\n")

# ============ Read data ============
cat("Reading CLR matrix...\n")
clr_matrix <- read.csv(clr_fp, row.names = 1, check.names = FALSE)
cat("  Matrix dimensions:", nrow(clr_matrix), "lineages x", ncol(clr_matrix), "samples\n")

cat("Reading metadata...\n")
metadata <- read_tsv(metadata_fp, show_col_types = FALSE)
cat("  Metadata rows:", nrow(metadata), "\n")

cat("Reading biome mapping...\n")
biome_mapping <- read_tsv(biome_mapping_fp, show_col_types = FALSE)
cat("  Biome groups:", length(unique(biome_mapping$group_biome)), "\n\n")

# ============ Prepare data ============
cat("Preparing data for classification...\n")

# Transpose CLR matrix: samples  taxa
clr_t <- as.data.frame(t(clr_matrix))
clr_t$SampleID <- rownames(clr_t)

# Merge with metadata and biome grouping
df <- clr_t %>%
  left_join(metadata %>% select(X.SampleID, env_biome), by = c("SampleID" = "X.SampleID")) %>%
  left_join(biome_mapping, by = "env_biome") %>%
  filter(!is.na(env_biome)) %>%
  mutate(group_biome = if_else(is.na(group_biome), "UNKNOWN", group_biome))

cat("  Total samples with biome labels:", nrow(df), "\n")

# Filter rare classes (< 10 samples per class)
cat("\nFiltering rare biome classes (< 10 samples)...\n")
biome_counts_orig <- table(df$env_biome)
valid_biomes_orig <- names(biome_counts_orig[biome_counts_orig >= 10])
df_original <- df %>% filter(env_biome %in% valid_biomes_orig)
cat("  Original biomes retained:", length(valid_biomes_orig), "out of", length(unique(df$env_biome)), "\n")
cat("  Samples for original biomes model:", nrow(df_original), "\n")

biome_counts_grouped <- table(df$group_biome)
valid_groups <- names(biome_counts_grouped[biome_counts_grouped >= 10])
df_grouped <- df %>% filter(group_biome %in% valid_groups)
cat("  Grouped biomes retained:", length(valid_groups), "out of", length(unique(df$group_biome)), "\n")
cat("  Samples for grouped biomes model:", nrow(df_grouped), "\n\n")

# ============ Train/Test Split and Model Training ============
set.seed(42)

# ============ ORIGINAL BIOMES MODEL ============
cat("============================================================\n")
cat("  Training Random Forest: ORIGINAL BIOMES\n")
cat("============================================================\n")

if (nrow(df_original) > 20 && length(valid_biomes_orig) >= 2) {
  # Create train/test split
  train_indices_original <- createDataPartition(
    df_original$env_biome,
    p = 0.7,
    list = FALSE
  )

  train_original <- df_original[train_indices_original, ]
  test_original <- df_original[-train_indices_original, ]

  cat("Training set:", nrow(train_original), "samples\n")
  cat("Test set:", nrow(test_original), "samples\n")
  cat("Number of classes:", length(unique(train_original$env_biome)), "\n")
  cat("Number of features:", ncol(train_original) - 3, "\n\n")

  # Train Random Forest
  cat("Training Random Forest...\n")
  rf_original <- randomForest(
    x = train_original %>% select(-SampleID, -env_biome, -group_biome),
    y = as.factor(train_original$env_biome),
    ntree = 500,
    mtry = floor(sqrt(ncol(train_original) - 3)),
    importance = TRUE,
    proximity = FALSE,
    keep.forest = TRUE
  )

  cat("\nModel OOB error rate:", round(tail(rf_original$err.rate[, 1], 1) * 100, 2), "%\n")

  # Predict on test set
  cat("Predicting on test set...\n")
  pred_original <- predict(rf_original, test_original %>% select(-SampleID, -env_biome, -group_biome))

  # Save model
  saveRDS(rf_original, file.path(out_dir, "rf_original_biomes", "model.rds"))
  cat("Saved model to: rf_original_biomes/model.rds\n")

  # Save predictions
  pred_df_original <- data.frame(
    SampleID = test_original$SampleID,
    true_biome = test_original$env_biome,
    predicted_biome = pred_original,
    stringsAsFactors = FALSE
  )
  write.csv(
    pred_df_original,
    file.path(out_dir, "rf_original_biomes", "predictions.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  cat("Saved predictions to: rf_original_biomes/predictions.csv\n")

  # Calculate and save train set info
  train_info_original <- data.frame(
    n_samples = nrow(train_original),
    n_features = ncol(train_original) - 3,
    n_classes = length(unique(train_original$env_biome)),
    ntree = 500,
    mtry = floor(sqrt(ncol(train_original) - 3)),
    oob_error = tail(rf_original$err.rate[, 1], 1)
  )
  write.csv(
    train_info_original,
    file.path(out_dir, "rf_original_biomes", "train_info.csv"),
    row.names = FALSE
  )

} else {
  cat("Insufficient data for original biomes model (need > 20 samples and >= 2 classes)\n")
}

# ============ GROUPED BIOMES MODEL ============
cat("\n============================================================\n")
cat("  Training Random Forest: GROUPED BIOMES\n")
cat("============================================================\n")

if (nrow(df_grouped) > 20 && length(valid_groups) >= 2) {
  # Create train/test split
  train_indices_grouped <- createDataPartition(
    df_grouped$group_biome,
    p = 0.7,
    list = FALSE
  )

  train_grouped <- df_grouped[train_indices_grouped, ]
  test_grouped <- df_grouped[-train_indices_grouped, ]

  cat("Training set:", nrow(train_grouped), "samples\n")
  cat("Test set:", nrow(test_grouped), "samples\n")
  cat("Number of classes:", length(unique(train_grouped$group_biome)), "\n")
  cat("Number of features:", ncol(train_grouped) - 3, "\n\n")

  # Train Random Forest
  cat("Training Random Forest...\n")
  rf_grouped <- randomForest(
    x = train_grouped %>% select(-SampleID, -env_biome, -group_biome),
    y = as.factor(train_grouped$group_biome),
    ntree = 500,
    mtry = floor(sqrt(ncol(train_grouped) - 3)),
    importance = TRUE,
    proximity = FALSE,
    keep.forest = TRUE
  )

  cat("\nModel OOB error rate:", round(tail(rf_grouped$err.rate[, 1], 1) * 100, 2), "%\n")

  # Predict on test set
  cat("Predicting on test set...\n")
  pred_grouped <- predict(rf_grouped, test_grouped %>% select(-SampleID, -env_biome, -group_biome))

  # Save model
  saveRDS(rf_grouped, file.path(out_dir, "rf_grouped_biomes", "model.rds"))
  cat("Saved model to: rf_grouped_biomes/model.rds\n")

  # Save predictions
  pred_df_grouped <- data.frame(
    SampleID = test_grouped$SampleID,
    true_biome = test_grouped$group_biome,
    predicted_biome = pred_grouped,
    stringsAsFactors = FALSE
  )
  write.csv(
    pred_df_grouped,
    file.path(out_dir, "rf_grouped_biomes", "predictions.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  cat("Saved predictions to: rf_grouped_biomes/predictions.csv\n")

  # Calculate and save train set info
  train_info_grouped <- data.frame(
    n_samples = nrow(train_grouped),
    n_features = ncol(train_grouped) - 3,
    n_classes = length(unique(train_grouped$group_biome)),
    ntree = 500,
    mtry = floor(sqrt(ncol(train_grouped) - 3)),
    oob_error = tail(rf_grouped$err.rate[, 1], 1)
  )
  write.csv(
    train_info_grouped,
    file.path(out_dir, "rf_grouped_biomes", "train_info.csv"),
    row.names = FALSE
  )

} else {
  cat("Insufficient data for grouped biomes model (need > 20 samples and >= 2 classes)\n")
}

cat("\n============================================================\n")
cat("  Random Forest Training Complete\n")
cat("============================================================\n")
cat("Models and predictions saved to:", out_dir, "\n")
cat("\nNext step: Run model evaluation script (26) to assess performance\n")
cat("============================================================\n")
