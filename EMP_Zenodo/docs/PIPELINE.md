# EMP R-only pipeline

This document describes what the R pipeline did, what was run, and what outputs
were produced. QIIME artifacts are legacy and live under `legacy/qiime`.

## Step-by-step
All R steps are executed through SLURM wrappers in `scripts/slurm`.

1) Extract soil sample IDs
- Script: `scripts/01_extract_soil_sample_ids.R`
- SLURM: `scripts/slurm/r01_extract_soil_sample_ids.slurm`
- Inputs: `data/raw/emp_metadata_release1_20170912.tsv`
- Outputs: `data/intermediate/soil_sample_ids.txt`

2) Extract soil metadata (helper step)
- Script: `scripts/slurm/r01b_extract_soil_metadata.slurm` (inline R)
- SLURM: `scripts/slurm/r01b_extract_soil_metadata.slurm`
- Inputs: `data/raw/emp_metadata_release1_20170912.tsv`, `data/intermediate/soil_sample_ids.txt`
- Outputs: `data/intermediate/soil_metadata.tsv`

3) Assign taxonomy + build soil ASV table
- Script: `scripts/02_assign_taxonomy_and_soil_asv_table.R`
- SLURM: `scripts/slurm/r02_assign_taxonomy_and_soil_asv_table.slurm`
- Inputs:
  - `data/raw/emp_deblur_150bp_release1.biom`
  - `data/raw/emp_deblur_150bp_min25.fasta`
  - `data/reference/silva_nr99_v138_1_train_set.fasta.gz`
  - `data/intermediate/soil_sample_ids.txt`
- Outputs:
  - `data/intermediate/asv_table_all_samples.tsv` (converted from BIOM if needed)
  - `data/intermediate/taxonomy_all_samples.csv`
  - `data/intermediate/asv_table_soil_only.csv`

4) Collapse to genus (all samples)
- Script: `scripts/03_collapse_genus_all_samples.R`
- SLURM: `scripts/slurm/r03_collapse_genus_all_samples.slurm`
- Inputs: `data/intermediate/taxonomy_all_samples.csv`, `data/intermediate/asv_table_all_samples.tsv`
- Outputs: `data/intermediate/genus_matrix_all_samples.csv`

5) Collapse to genus (soil only) + soil taxonomy
- Script: `scripts/07_collapse_genus_soil_only.R`
- SLURM: `scripts/slurm/r07_collapse_genus_soil_only.slurm`
- Inputs: `data/intermediate/taxonomy_all_samples.csv`, `data/intermediate/asv_table_all_samples.tsv`, `data/intermediate/soil_metadata.tsv`
- Outputs:
  - `data/intermediate/genus_matrix_soil_only.csv`
  - `data/intermediate/taxonomy_soil_only.csv`

6) QC and summaries (optional but recommended)
- `scripts/04_qc_read_counts_vs_genus.R`
  - SLURM: `scripts/slurm/r04_qc_read_counts_vs_genus.slurm`
  - Inputs: `data/intermediate/genus_matrix_all_samples.csv`, `data/intermediate/asv_table_all_samples.tsv`, `data/intermediate/taxonomy_all_samples.csv`
  - Outputs: `outputs/qc_tables/all_samples_read_count_comparison.csv`, `outputs/qc_tables/unclassified_reads_per_sample.csv` (optional)
- `scripts/05_qc_asv_id_match.R`
  - SLURM: `scripts/slurm/r05_qc_asv_id_match.slurm`
  - Inputs: `data/intermediate/asv_table_all_samples.tsv`, `data/intermediate/asv_taxonomy_matrix_all_samples.csv`
  - Outputs: `outputs/qc_tables/asv_id_match_summary.txt`, `outputs/qc_tables/asvs_missing_from_asv_table.csv`, `outputs/qc_tables/asvs_extra_in_asv_table.csv`
- `scripts/06_stats_all_samples.R`
  - SLURM: `scripts/slurm/r06_stats_all_samples.slurm`
  - Outputs: `outputs/qc_tables/*` and `outputs/figures/*` for all-sample summaries
- `scripts/08_stats_soil_only.R`
  - SLURM: `scripts/slurm/r08_stats_soil_only.slurm`
  - Outputs: `outputs/qc_tables/*` and `outputs/figures/*` for soil-only summaries

7) Prevalence and threshold analyses
- `scripts/10_genus_prevalence_summary.R`
  - SLURM: `scripts/slurm/r10_genus_prevalence_summary.slurm`
  - Inputs: `data/intermediate/genus_matrix_soil_only.csv`
  - Outputs: `data/processed/genus_prevalence_summary.csv`, `outputs/figures/genus_prevalence_histogram.png`
- `scripts/11_genus_prevalence_thresholds.R`
  - SLURM: `scripts/slurm/r11_genus_prevalence_thresholds.slurm`
  - Inputs: `data/intermediate/genus_matrix_soil_only.csv`
  - Outputs: `data/processed/genus_prevalence_thresholds_1_to_20.csv`, `outputs/figures/genus_prevalence_thresholds_1_to_20.png`
- `scripts/12_lineage_counts_thresholds.R`
  - SLURM: `scripts/slurm/r12_lineage_counts_thresholds.slurm`
  - Inputs: `data/intermediate/taxonomy_soil_only.csv`
  - Outputs: `data/processed/lineage_counts_all.csv`, `data/processed/lineage_counts_thresholds_1_to_20.csv`
- `scripts/13_lineage_prevalence_thresholds.R`
  - SLURM: `scripts/slurm/r13_lineage_prevalence_thresholds.slurm`
  - Inputs: `data/intermediate/taxonomy_soil_only.csv`, `data/intermediate/genus_matrix_soil_only.csv`
  - Outputs: `data/processed/lineage_prevalence_thresholds.csv`,
    `outputs/figures/lineage_fraction_remaining_plot.png`,
    `outputs/figures/lineage_fraction_matrix_remaining_plot.png`

8) Build lineage CLR matrix (analysis-ready)
- Script: `scripts/09_build_lineage_clr_matrix.R`
- SLURM: `scripts/slurm/r09_build_lineage_clr_matrix.slurm`
- Inputs: `data/intermediate/asv_table_all_samples.tsv`, `data/intermediate/soil_sample_ids.txt`, `data/intermediate/taxonomy_soil_only.csv`
- Outputs: `data/processed/lineage_clr_matrix_soil_filtered.csv`

9) Modeling
- Script: `scripts/14_model_fit_mixture.R`
- SLURM: `scripts/slurm/r14_model_fit_mixture.slurm`
- Inputs: `data/processed/lineage_clr_matrix_soil_filtered.csv`
- Outputs: `models/best_model_emp_k20.RData`

10) Cluster vs metadata comparison
- Script: `scripts/15_cluster_metadata_compare.R`
- SLURM: `scripts/slurm/r15_cluster_metadata_compare.slurm`
- Inputs:
  - `models/best_model_emp_k20.RData`
  - `data/intermediate/soil_metadata.tsv`
  - `data/processed/lineage_clr_matrix_soil_filtered.csv` (for sample IDs)
- Outputs:
  - `outputs/cluster_metadata/cluster_assignments.tsv`
  - `outputs/cluster_metadata/cluster_metadata_join.tsv`
  - `outputs/cluster_metadata/cluster_sizes.tsv`
  - `outputs/cluster_metadata/cluster_metadata_tests.tsv`

11) Cluster CLR averages
- Script: `scripts/16_cluster_clr_averages.R`
- SLURM: `scripts/slurm/r16_cluster_clr_averages.slurm`
- Inputs:
  - `data/processed/lineage_clr_matrix_soil_filtered.csv`
  - `outputs/cluster_metadata/cluster_assignments.tsv`
- Outputs:
  - `outputs/cluster_metadata/cluster_avg_clr_matrix.csv`
  - `outputs/cluster_metadata/cluster_avg_clr_long.tsv`
  - `outputs/cluster_metadata/cluster_avg_clr_summary.tsv`

12) Plot cluster metadata
- Script: `scripts/17_plot_cluster_metadata.py`
- SLURM: `scripts/slurm/r17_plot_cluster_metadata.slurm`
- Inputs: `outputs/cluster_metadata/cluster_metadata_join.tsv`
- Outputs: `outputs/cluster_metadata/plots/*.pdf` (multiple PDF plots)

13) Calculate alpha diversity metrics
- Script: `scripts/18_calculate_alpha_diversity.R`
- SLURM: `scripts/slurm/r18_calculate_alpha_diversity.slurm`
- Inputs: `data/intermediate/asv_table_soil_only.csv`
- Outputs:
  - `data/processed/alpha_diversity_soil.csv` (Observed OTUs, Chao1, Shannon)
  - `outputs/figures/alpha_diversity_distributions.png`

14) Plot alpha diversity by cluster
- Script: `scripts/19_plot_alpha_diversity_by_cluster.R`
- SLURM: `scripts/slurm/r19_plot_alpha_diversity_by_cluster.slurm`
- Inputs:
  - `data/processed/alpha_diversity_soil.csv`
  - `outputs/cluster_metadata/cluster_assignments.tsv`
- Outputs:
  - `outputs/cluster_metadata/cluster_alpha_diversity.tsv`
  - `outputs/cluster_metadata/cluster_alpha_diversity_summary.tsv`
  - `outputs/cluster_metadata/plots/alpha_div_*.pdf` (7 PDF plots: 3 boxplots, 3 scatter, 1 heatmap)

## Final outputs of record
- `data/intermediate/genus_matrix_soil_only.csv`
- `data/processed/lineage_clr_matrix_soil_filtered.csv`
- `data/processed/genus_prevalence_summary.csv`
- `data/processed/lineage_prevalence_thresholds.csv`
- `data/processed/alpha_diversity_soil.csv` (alpha diversity metrics)
- `models/best_model_emp_k20.RData` (if modeling ran)
- `outputs/cluster_metadata/cluster_metadata_tests.tsv` (if cluster comparison ran)
- `outputs/cluster_metadata/cluster_alpha_diversity.tsv` (alpha diversity by cluster)

## Legacy
- QIIME scripts and artifacts are retained under `legacy/qiime`.
- OS artifacts are under `legacy/osx_artifacts`.

20) Enhanced Cluster CLR Statistics
- Script: `scripts/20_cluster_clr_statistics.R`
- SLURM: `scripts/slurm/r20_cluster_clr_statistics.slurm`
- Inputs:
  - `data/processed/lineage_clr_matrix_soil_filtered.csv`
  - `outputs/cluster_metadata/cluster_assignments.tsv`
- Outputs:
  - `outputs/cluster_metadata/cluster_clr_statistics.csv` (mean, SD, SE, CI per lineage per cluster)
  - `outputs/cluster_metadata/cluster_clr_statistics_summary.tsv`
  - `outputs/ordination/cluster_errorbar_top20_taxa.pdf`
  - `outputs/ordination/cluster_errorbar_selected_phyla.pdf`
  - `outputs/ordination/cluster_variance_heatmap.pdf`

21) PCA/PCoA Ordination Analysis
- Script: `scripts/21_ordination_analysis.R`
- SLURM: `scripts/slurm/r21_ordination_analysis.slurm`
- Inputs:
  - `data/processed/lineage_clr_matrix_soil_filtered.csv`
  - `outputs/cluster_metadata/cluster_assignments.tsv`
- Outputs:
  - Data files: `outputs/ordination/pca_results.rds`, `pca_sample_scores.csv`, `pca_variance_explained.csv`, `pca_loadings.csv`, `pca_loadings_top50_pc{1,2}.csv`, `pcoa_results.rds`, `pcoa_sample_scores.csv`
  - Plots: `outputs/ordination/pca_scree_plot.pdf`, `pca_biplot_pc1_pc2.pdf`, `pca_biplot_pc1_pc3.pdf`, `pca_biplot_pc2_pc3.pdf`, `pca_loadings_pc1_top30.pdf`, `pca_loadings_pc2_top30.pdf`, `pca_loadings_biplot.pdf`, `pcoa_biplot_axis1_axis2.pdf`
