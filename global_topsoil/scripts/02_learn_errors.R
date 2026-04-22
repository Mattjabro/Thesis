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

# Define input/output paths
intermediate_dir <- file.path(project_root, "data", "intermediate")
filt_path <- file.path(intermediate_dir, "filtered_fastq")

# Automatically detect sample names based on _F_filt.fastq.gz files
filtFs <- sort(list.files(filt_path, pattern = "_F_filt.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(filt_path, pattern = "_R_filt.fastq.gz", full.names = TRUE))

# Optional: sanity check that file pairs match
samples <- gsub("_F_filt.fastq.gz", "", basename(filtFs))
stopifnot(all(samples == gsub("_R_filt.fastq.gz", "", basename(filtRs))))

cat(" Forward filtered files:\n"); print(filtFs)
cat(" Reverse filtered files:\n"); print(filtRs)

# Learn error rates
cat(" Learning error rates (forward)...\n")
errF <- learnErrors(filtFs, multithread=TRUE)

cat(" Learning error rates (reverse)...\n")
errR <- learnErrors(filtRs, multithread=TRUE)

# Save error models
cat(" Saving error models...\n")
saveRDS(errF, file.path(intermediate_dir, "errF.rds"))
saveRDS(errR, file.path(intermediate_dir, "errR.rds"))

cat(" Done learning error rates.\n")
