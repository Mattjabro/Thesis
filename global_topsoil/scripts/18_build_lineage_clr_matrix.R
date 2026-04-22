# Load required libraries
library(tidyverse)
library(phyloseq)

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
processed_dir <- file.path(project_root, "data", "processed")
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

# Load ASV counts and taxonomy table
asv_counts <- readRDS(file.path(intermediate_dir, "seqtab_nochim_bacteria.rds"))
tax_table_raw <- readRDS(file.path(intermediate_dir, "taxa_bacteria.rds"))

# Check dimensions match
stopifnot(identical(colnames(asv_counts), rownames(tax_table_raw)))

# Create phyloseq object
OTU <- otu_table(asv_counts, taxa_are_rows = FALSE)
TAX <- tax_table(tax_table_raw)
ps <- phyloseq(OTU, TAX)

# Collapse to Genus level (L6)
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)

# Extract full taxonomy lineage for each Genus (taxa names)
lineage_strings <- apply(tax_table(ps_genus), 1, function(x) {
  paste(na.omit(x), collapse = ";")
})
# Assign to matrix column names
collapsed_matrix <- as(otu_table(ps_genus), "matrix")
colnames(collapsed_matrix) <- make.unique(lineage_strings)

# Step 1: Filter samples with <5000 total counts
sample_sums <- rowSums(collapsed_matrix)
filtered <- collapsed_matrix[sample_sums >= 5000, ]

# Step 2: Filter taxa present in <20 samples
presence <- filtered > 0
prevalence <- colSums(presence)
filtered <- filtered[, prevalence >= 20]

# Step 3: CLR transform
clr_input <- filtered + 1
clr_matrix <- log(clr_input)
clr_matrix <- clr_matrix - rowMeans(clr_matrix)

# Step 4: Transpose the matrix (taxa become rows)
clr_matrix_transposed <- t(clr_matrix)

# Save output
write.csv(clr_matrix_transposed, file.path(processed_dir, "lineage_clr_matrix_soil_filtered.csv"))

# Optional: Print preview
cat("Taxa remaining:", nrow(clr_matrix_transposed), "\n")
cat("Samples remaining:", ncol(clr_matrix_transposed), "\n")
print(head(clr_matrix_transposed[, 1:min(5, ncol(clr_matrix_transposed))]))
