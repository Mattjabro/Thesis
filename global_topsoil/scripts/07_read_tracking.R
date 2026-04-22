#!/usr/bin/env Rscript

library(dada2)

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

# Define paths
intermediate_dir <- file.path(project_root, "data", "intermediate")
qc_dir <- file.path(project_root, "outputs", "qc_tables")

# Load saved .rds objects from previous steps
filt_summary <- readRDS(file.path(intermediate_dir, "filtering_summary.rds"))
dadaFs <- readRDS(file.path(intermediate_dir, "dadaFs.rds"))
dadaRs <- readRDS(file.path(intermediate_dir, "dadaRs.rds"))
mergers <- readRDS(file.path(intermediate_dir, "mergers.rds"))
seqtab.nochim <- readRDS(file.path(intermediate_dir, "seqtab_nochim.rds"))

# Derive sample names automatically from dadaFs
samples <- names(dadaFs)

# Helper to count number of reads per step
getN <- function(x) sum(getUniques(x))

# Create tracking table
track <- data.frame(
  sample_id = samples,
  input = filt_summary[,1],
  filtered = filt_summary[,2],
  `percentage of input passed filter` = round(100 * filt_summary[,2] / filt_summary[,1], 2),
  denoised = sapply(dadaFs, getN),
  merged = sapply(mergers, getN),
  `percentage of input merged` = round(100 * sapply(mergers, getN) / filt_summary[,1], 2),
  non_chimeric = rowSums(seqtab.nochim),
  `percentage of input non-chimeric` = round(100 * rowSums(seqtab.nochim) / filt_summary[,1], 2),
  `percentage filtered that merged` = round(100 * sapply(mergers, getN) / filt_summary[,2], 2),
  `percentage merged that non-chimeric` = round(100 * rowSums(seqtab.nochim) / sapply(mergers, getN), 2),
  num_ASVs = rowSums(seqtab.nochim > 0)
)

# Save output
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(track, file.path(qc_dir, "dada2_read_tracking.csv"), row.names = FALSE)

# Print preview
print(track)
