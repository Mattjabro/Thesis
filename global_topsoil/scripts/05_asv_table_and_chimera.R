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

cat(" Reading merged pairs...\n")
mergers <- readRDS(file.path(intermediate_dir, "mergers.rds"))

cat(" Creating ASV sequence table...\n")
seqtab <- makeSequenceTable(mergers)
cat(" Sequence table dimensions: ", dim(seqtab)[1], "samples x", dim(seqtab)[2], "ASVs\n")
saveRDS(seqtab, file.path(intermediate_dir, "seqtab.rds"))

cat(" Removing chimeras...\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
cat(" Chimera-free table dimensions: ", dim(seqtab.nochim)[1], "samples x", dim(seqtab.nochim)[2], "ASVs\n")

saveRDS(seqtab.nochim, file.path(intermediate_dir, "seqtab_nochim.rds"))

cat(" ASV table created and chimeras removed.\n")
