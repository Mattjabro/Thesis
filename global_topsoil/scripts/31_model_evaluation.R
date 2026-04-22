#!/usr/bin/env Rscript
# ============================================================================
# 31_model_evaluation.R
#
# Evaluates Random Forest model performance and extracts feature importance.
# Generates confusion matrices, performance metrics, and visualizations of
# the most important taxa driving biome classification.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(randomForest)
  library(caret)
  library(pheatmap)
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

# ============ Setup ============
ml_dir <- file.path(project_root, "outputs", "ml_models")
plots_dir <- file.path(ml_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("  Random Forest Model Evaluation (Global Topsoil)\n")
cat("============================================================\n")
cat("ML models directory:", ml_dir, "\n")
cat("Plots directory:", plots_dir, "\n\n")

# ============ Helper Functions ============
evaluate_model <- function(model_name, model_path, pred_path, out_subdir, plots_dir) {
  cat("\n============================================================\n")
  cat(" ", toupper(model_name), "\n")
  cat("============================================================\n")

  # Check if files exist
  if (!file.exists(model_path)) {
    cat("Model file not found:", model_path, "\n")
    return(NULL)
  }
  if (!file.exists(pred_path)) {
    cat("Predictions file not found:", pred_path, "\n")
    return(NULL)
  }

  # Load model and predictions
  cat("Loading model...\n")
  rf_model <- readRDS(model_path)

  cat("Loading predictions...\n")
  predictions <- read.csv(pred_path, stringsAsFactors = FALSE)

  cat("Samples in test set:", nrow(predictions), "\n")
  cat("Classes:", length(unique(predictions$true_biome)), "\n\n")

  # ============ Confusion Matrix ============
  cat("Calculating confusion matrix...\n")
  cm <- confusionMatrix(
    as.factor(predictions$predicted_biome),
    as.factor(predictions$true_biome)
  )

  write.csv(
    cm$table,
    file.path(out_subdir, "confusion_matrix.csv"),
    quote = FALSE
  )
  cat("  Saved: confusion_matrix.csv\n")

  # ============ Performance Metrics ============
  cat("Calculating performance metrics...\n")
  metrics <- data.frame(
    Accuracy = cm$overall["Accuracy"],
    Kappa = cm$overall["Kappa"],
    Precision = mean(cm$byClass[, "Precision"], na.rm = TRUE),
    Recall = mean(cm$byClass[, "Recall"], na.rm = TRUE),
    F1 = mean(cm$byClass[, "F1"], na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  write.csv(
    metrics,
    file.path(out_subdir, "metrics.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  cat("  Saved: metrics.csv\n")
  cat("  Accuracy:", round(metrics$Accuracy * 100, 2), "%\n")
  cat("  Kappa:", round(metrics$Kappa, 3), "\n")
  cat("  F1:", round(metrics$F1, 3), "\n\n")

  # ============ Feature Importance ============
  cat("Extracting feature importance...\n")
  importance_matrix <- importance(rf_model)

  importance_df <- data.frame(
    lineage = rownames(importance_matrix),
    MeanDecreaseAccuracy = importance_matrix[, "MeanDecreaseAccuracy"],
    MeanDecreaseGini = importance_matrix[, "MeanDecreaseGini"],
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(MeanDecreaseAccuracy))

  write.csv(
    importance_df,
    file.path(out_subdir, "feature_importance.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  cat("  Saved: feature_importance.csv\n")

  # ============ Taxonomy Breakdown ============
  cat("Analyzing taxonomy of important features...\n")
  top50 <- importance_df %>% head(50)

  # Extract phylum from lineage (assuming format: Bacteria;Phylum;...)
  top50 <- top50 %>%
    mutate(
      phylum = sapply(strsplit(lineage, ";"), function(x) {
        if (length(x) >= 2) x[2] else "Unknown"
      })
    )

  taxonomy_summary <- top50 %>%
    group_by(phylum) %>%
    summarise(
      count = n(),
      avg_importance = mean(MeanDecreaseAccuracy),
      .groups = "drop"
    ) %>%
    arrange(desc(count))

  write.csv(
    taxonomy_summary,
    file.path(out_subdir, "feature_taxonomy_summary.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  cat("  Saved: feature_taxonomy_summary.csv\n\n")

  # ============ Plots ============
  cat("Generating plots...\n")

  # Confusion matrix heatmap
  pdf_name <- paste0("confusion_matrix_", tolower(gsub(" ", "_", model_name)), ".pdf")
  pdf(file.path(plots_dir, pdf_name), width = 12, height = 10)
  pheatmap(
    cm$table,
    display_numbers = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = paste("Confusion Matrix:", model_name),
    angle_col = 45,
    fontsize_number = 8
  )
  dev.off()
  cat("  Saved:", pdf_name, "\n")

  # Feature importance plot (top 50)
  pdf_name <- paste0("feature_importance_", tolower(gsub(" ", "_", model_name)), "_top50.pdf")
  pdf(file.path(plots_dir, pdf_name), width = 16, height = 14)

  p <- ggplot(top50, aes(x = reorder(lineage, MeanDecreaseAccuracy),
                          y = MeanDecreaseAccuracy)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    theme_bw() +
    labs(
      title = paste("Top 50 Important Taxa:", model_name),
      x = "Taxon (Full Lineage)",
      y = "Mean Decrease Accuracy"
    ) +
    theme(
      axis.text.y = element_text(size = 5),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  print(p)
  dev.off()
  cat("  Saved:", pdf_name, "\n")

  cat("\n", toupper(model_name), "evaluation complete!\n")

  return(list(
    metrics = metrics,
    confusion_matrix = cm$table,
    importance = importance_df,
    taxonomy = taxonomy_summary
  ))
}

# ============ Evaluate Original Biomes Model ============
original_results <- evaluate_model(
  model_name = "Original Biomes",
  model_path = file.path(ml_dir, "rf_original_biomes", "model.rds"),
  pred_path = file.path(ml_dir, "rf_original_biomes", "predictions.csv"),
  out_subdir = file.path(ml_dir, "rf_original_biomes"),
  plots_dir = plots_dir
)

# ============ Evaluate Grouped Biomes Model ============
grouped_results <- evaluate_model(
  model_name = "Grouped Biomes",
  model_path = file.path(ml_dir, "rf_grouped_biomes", "model.rds"),
  pred_path = file.path(ml_dir, "rf_grouped_biomes", "predictions.csv"),
  out_subdir = file.path(ml_dir, "rf_grouped_biomes"),
  plots_dir = plots_dir
)

# ============ Summary ============
cat("\n============================================================\n")
cat("  Model Evaluation Complete\n")
cat("============================================================\n")

if (!is.null(original_results)) {
  cat("\nOriginal Biomes Model:\n")
  cat("  Accuracy:", round(original_results$metrics$Accuracy * 100, 2), "%\n")
  cat("  Top 5 important taxa:\n")
  for (i in 1:min(5, nrow(original_results$importance))) {
    cat("   ", i, ".", original_results$importance$lineage[i], "\n")
  }
}

if (!is.null(grouped_results)) {
  cat("\nGrouped Biomes Model:\n")
  cat("  Accuracy:", round(grouped_results$metrics$Accuracy * 100, 2), "%\n")
  cat("  Top 5 important taxa:\n")
  for (i in 1:min(5, nrow(grouped_results$importance))) {
    cat("   ", i, ".", grouped_results$importance$lineage[i], "\n")
  }
}

cat("\nAll results saved to:", ml_dir, "\n")
cat("============================================================\n")
