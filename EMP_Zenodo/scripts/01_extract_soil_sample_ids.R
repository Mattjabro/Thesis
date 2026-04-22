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

metadata_path <- file.path(project_root, "data", "raw", "emp_metadata_release1_20170912.tsv")
metadata <- read.delim(metadata_path, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

# Normalize and match
env_clean <- tolower(metadata$env_material)
env_clean[is.na(env_clean)] <- ""
is_soil <- grepl("soil|sediment|rhizosphere", env_clean)

cat("Number of samples matching soil terms:", sum(is_soil), "\n")

soil_samples <- metadata[is_soil, ]

# Extract the SampleID column safely
sample_col <- grep("sampleid", names(soil_samples), ignore.case = TRUE, value = TRUE)[1]
cat("Using sample ID column:", sample_col, "\n")
soil_sample_ids <- soil_samples[[sample_col]]

# Write to file
out_path <- file.path(project_root, "data", "intermediate", "soil_sample_ids.txt")
write.table(soil_sample_ids, file = out_path, quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(length(soil_sample_ids), "sample IDs saved to:", out_path, "\n")
