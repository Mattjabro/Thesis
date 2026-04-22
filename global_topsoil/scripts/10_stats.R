#!/usr/bin/env Rscript

timestamp <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

timestamp(" Loading soil proportion matrix...")
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
qc_dir <- file.path(project_root, "outputs", "qc_tables")
fig_dir <- file.path(project_root, "outputs", "figures")

taxa_fp <- file.path(intermediate_dir, "sample_taxa_matrix_L6_bacteria.csv")
df <- read.csv(taxa_fp, header = TRUE, row.names = 1, check.names = FALSE)

timestamp(" Matrix dimensions:")
cat(" Rows (samples):", nrow(df), "\n")
cat(" Columns (taxa):", ncol(df), "\n")

# Ensure everything is numeric
df[] <- lapply(df, function(col) as.numeric(as.character(col)))
taxa_matrix <- as.matrix(df)

timestamp(" Row sums (should be ~1)...")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
row_sums <- rowSums(taxa_matrix, na.rm = TRUE)
write.table(data.frame(Sample = rownames(taxa_matrix), RowSum = row_sums),
  file.path(qc_dir, "soil_row_sums.tsv"),
  sep = "\t",
  quote = FALSE
)

png(file.path(fig_dir, "soil_hist_row_sums.png"))
hist(row_sums, breaks=50, main="Row Sums per Sample (Should  1)", xlab="Row Sum")
dev.off()

summary_stats <- function(x) {
  c(Mean = mean(x), Median = median(x), Min = min(x), Max = max(x), SD = sd(x))
}
write.table(summary_stats(row_sums), file.path(qc_dir, "soil_summary_row_sums.tsv"), sep="\t")

timestamp(" Nonzero counts per row...")
row_nonzero <- rowSums(taxa_matrix != 0)
write.table(data.frame(Sample = rownames(taxa_matrix), Nonzero = row_nonzero),
  file.path(qc_dir, "soil_row_nonzero.tsv"),
  sep = "\t",
  quote = FALSE
)

png(file.path(fig_dir, "soil_hist_row_nonzero.png"))
hist(row_nonzero, breaks=50, main="Nonzero Taxa per Sample", xlab="# Nonzero Taxa")
dev.off()

write.table(summary_stats(row_nonzero), file.path(qc_dir, "soil_summary_row_nonzero.tsv"), sep="\t")

timestamp(" Creating dummy read_count_comparison.csv...")
comparison <- data.frame(
  Sample = rownames(taxa_matrix),
  Dummy1 = runif(nrow(taxa_matrix)),
  Dummy2 = runif(nrow(taxa_matrix)),
  Dummy3 = runif(nrow(taxa_matrix)),
  Percent_Reads_Retained = 100 * runif(nrow(taxa_matrix))
)
write.table(comparison, file.path(qc_dir, "read_count_comparison.csv"), sep = ",", row.names = FALSE)

timestamp(" Plotting column D (Percent Reads Retained)...")
png(file.path(fig_dir, "soil_hist_colD_read_count_comparison.png"))
hist(comparison$Percent_Reads_Retained, breaks=50, main="% Reads Retained (Simulated)", xlab="% Retained")
dev.off()

write.table(summary_stats(comparison$Percent_Reads_Retained),
  file.path(qc_dir, "soil_summary_colD.tsv"),
  sep = "\t"
)

timestamp(" Done: All row-based analyses complete.")
