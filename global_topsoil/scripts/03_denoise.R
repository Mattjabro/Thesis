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

# Set paths
intermediate_dir <- file.path(project_root, "data", "intermediate")
filt_path <- file.path(intermediate_dir, "filtered_fastq")

# Automatically detect forward and reverse filtered reads
filtFs <- sort(list.files(filt_path, pattern = "_F_filt.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(filt_path, pattern = "_R_filt.fastq.gz", full.names = TRUE))

# Ensure forward and reverse files match
samples <- gsub("_F_filt.fastq.gz", "", basename(filtFs))
stopifnot(all(samples == gsub("_R_filt.fastq.gz", "", basename(filtRs))))

cat(" Reading error models...\n")
errF <- readRDS(file.path(intermediate_dir, "errF.rds"))
errR <- readRDS(file.path(intermediate_dir, "errR.rds"))

cat(" Running dada() on forward reads...\n")
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

cat(" Running dada() on reverse reads...\n")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

cat(" Saving denoised outputs...\n")
saveRDS(dadaFs, file.path(intermediate_dir, "dadaFs.rds"))
saveRDS(dadaRs, file.path(intermediate_dir, "dadaRs.rds"))

cat(" Denoising complete.\n")
