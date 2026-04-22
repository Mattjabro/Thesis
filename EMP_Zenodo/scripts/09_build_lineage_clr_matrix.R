# Load libraries
library(tidyverse)

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

cat("Starting EMP soil-only lineage matrix build...\n")
Sys.setenv("VROOM_CONNECTION_SIZE" = 5000000)

# ============ Step 1: Load Inputs ============

cat("Reading full EMP ASV table...\n")
asv_counts <- read_tsv(file.path(project_root, "data", "intermediate", "asv_table_all_samples.tsv"),
                       skip = 1, show_col_types = FALSE) %>%
  rename(ASV = `#OTU ID`) %>%
  column_to_rownames("ASV") %>%
  as.data.frame(check.names = FALSE)

cat("Reading list of soil sample IDs...\n")
soil_samples <- read_lines(file.path(project_root, "data", "intermediate", "soil_sample_ids.txt"))

cat("Reading taxonomy for soil-only ASVs...\n")
taxonomy <- read_csv(file.path(project_root, "data", "intermediate", "taxonomy_soil_only.csv"),
                     show_col_types = FALSE)

# ============ Step 2: Filter and Match ============

cat("Subsetting ASV table to soil samples...\n")
soil_samples <- intersect(soil_samples, colnames(asv_counts))
asv_counts_soil <- asv_counts[, soil_samples]

cat("Matching ASVs to taxonomy...\n")
shared_asvs <- intersect(rownames(asv_counts_soil), taxonomy$ASV)
asv_counts_soil <- asv_counts_soil[shared_asvs, ]
taxonomy <- taxonomy %>% filter(ASV %in% shared_asvs)

cat("Matched ASVs:", length(shared_asvs), "\n")

# ============ Step 3: Build lineage strings only (no ASV sequences) ============

cat("Building full lineage strings...\n")
taxonomy$Lineage <- apply(taxonomy[, -which(names(taxonomy) == "ASV")], 1, function(x) paste(na.omit(x), collapse = ";"))

asv_to_lineage <- taxonomy %>% select(ASV, Lineage)

# ============ Step 4: Collapse ASV counts to Lineage ============

cat("Collapsing ASV counts to full lineage level...\n")
asv_counts_soil$ASV <- rownames(asv_counts_soil)
long_counts <- pivot_longer(asv_counts_soil, -ASV, names_to = "SampleID", values_to = "Count")

# Merge with lineage
long_lineage <- long_counts %>%
  left_join(asv_to_lineage, by = "ASV") %>%
  filter(!is.na(Lineage)) %>%
  group_by(Lineage, SampleID) %>%
  summarise(Count = sum(Count), .groups = "drop")

# Pivot to sample  lineage matrix
wide_matrix <- long_lineage %>%
  pivot_wider(names_from = SampleID, values_from = Count, values_fill = 0) %>%
  column_to_rownames("Lineage")

# ============ Step 5: Filter Samples and Taxa ============

cat("Filtering samples with <5000 reads...\n")
sample_sums <- colSums(wide_matrix)
wide_matrix <- wide_matrix[, sample_sums >= 5000]
cat("Remaining samples:", ncol(wide_matrix), "\n")

cat("Filtering taxa found in <50 samples...\n")
taxa_prevalence <- rowSums(wide_matrix > 0)
wide_matrix <- wide_matrix[taxa_prevalence >= 50, ]
cat("Remaining taxa:", nrow(wide_matrix), "\n")

# ============ Step 6: CLR Transform ============

cat("Performing CLR transformation...\n")
clr_input <- wide_matrix + 1
clr_matrix <- log(clr_input)
clr_matrix <- clr_matrix - rowMeans(clr_matrix)
cat("CLR transformation complete.\n")

# ============ Step 7: Save Output ============

output_fp <- file.path(project_root, "data", "processed", "lineage_clr_matrix_soil_filtered.csv")
cat("Saving final matrix to", output_fp, "...\n")
write.csv(clr_matrix, output_fp)
cat("Done! Matrix has", nrow(clr_matrix), "lineages ", ncol(clr_matrix), "samples.\n")
