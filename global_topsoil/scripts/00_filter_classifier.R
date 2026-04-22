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

# Load and filter a FASTA classifier to exclude unwanted taxa
ref_dir <- file.path(project_root, "data", "reference")
input_file <- file.path(ref_dir, "silva_nr99_v138.2_toGenus_trainset.fa")
output_file <- file.path(ref_dir, "silva_nr99_v138.2_toGenus_trainset_filtered.fa")

# Read all lines
lines <- readLines(input_file)

# Identify header lines (taxa) and match indices
header_idx <- grep("^>", lines)
keep <- sapply(header_idx, function(i) {
  tax <- lines[i]
  !grepl("Eukaryota|Chloroplast|Mitochondria", tax, ignore.case = TRUE)
})

# Flatten: keep header + sequence line after it
filtered_lines <- unlist(mapply(function(i, k) {
  if (k) return(lines[i:(i+1)])
  else return(NULL)
}, header_idx, keep, SIMPLIFY = FALSE))

# Write filtered FASTA
writeLines(filtered_lines, output_file)
cat("Filtered classifier written to:", output_file, "\n")
