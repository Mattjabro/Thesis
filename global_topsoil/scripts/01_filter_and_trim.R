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

fastq_path <- file.path(project_root, "data", "raw", "fastq_files")
intermediate_dir <- file.path(project_root, "data", "intermediate")
filtered_dir <- file.path(intermediate_dir, "filtered_fastq")

# List forward and reverse reads
fnFs_all <- sort(list.files(fastq_path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs_all <- sort(list.files(fastq_path, pattern = "_2.fastq.gz", full.names = TRUE))

# Extract sample names (before "_1.fastq.gz")
samples <- gsub("_1.fastq.gz", "", basename(fnFs_all))

# Set filtered file paths
filtFs <- file.path(filtered_dir, paste0(samples, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(samples, "_R_filt.fastq.gz"))
names(filtFs) <- samples
names(filtRs) <- samples

cat(" Creating output directory (if needed)...\n")
dir.create(filtered_dir, showWarnings = FALSE, recursive = TRUE)

cat(" Running filterAndTrim...\n")
out <- filterAndTrim(fnFs_all, filtFs, fnRs_all, filtRs,
                     truncLen = c(250,179),
                     maxN = 0,
                     maxEE = c(2,2),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE)

cat(" Filtering complete. Output:\n")
print(out)

saveRDS(out, file.path(intermediate_dir, "filtering_summary.rds"))

