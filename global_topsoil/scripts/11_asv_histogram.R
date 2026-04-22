#!/usr/bin/env Rscript

timestamp <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

# =====================
#  Load required libs
# =====================
timestamp(" Loading packages...")
suppressPackageStartupMessages({
  library(ggplot2)
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

intermediate_dir <- file.path(project_root, "data", "intermediate")
fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================
#  Define file paths
# ============================
seqtab_fp <- file.path(intermediate_dir, "seqtab_nochim.rds")
hist_out_fp <- file.path(fig_dir, "soil_hist_asv_counts_per_sample.png")

# ============================
#  Read in ASV table
# ============================
timestamp(" Reading ASV table (post-chimera)...")
seqtab <- readRDS(seqtab_fp)

# ============================
#  Count ASVs per sample
# ============================
timestamp(" Counting non-zero ASVs per sample...")
asv_counts <- rowSums(seqtab > 0)

# ============================
#  Plot histogram
# ============================
timestamp(" Plotting histogram...")

png(hist_out_fp, width = 1000, height = 700)
hist(asv_counts,
     breaks = 30,
     col = "skyblue",
     border = "white",
     main = "Number of ASVs Detected per Sample",
     xlab = "Number of Non-Zero ASVs",
     ylab = "Number of Samples")
dev.off()

timestamp(paste(" Histogram saved to", hist_out_fp))
