# Soil DADA2 + modeling pipeline (R)

This document describes what the R pipeline does and what outputs are produced.
QIIME artifacts are legacy and live under `legacy/qiime`.

## Step-by-step
All R steps are executed through SLURM wrappers in `scripts/slurm`.

0) Filter classifier (optional)
- Script: `scripts/00_filter_classifier.R`
- SLURM: `scripts/slurm/r00_filter_classifier.slurm`
- Inputs: `data/reference/silva_nr99_v138.2_toGenus_trainset.fa`
- Outputs: `data/reference/silva_nr99_v138.2_toGenus_trainset_filtered.fa`

1) Filter and trim reads
- Script: `scripts/01_filter_and_trim.R`
- SLURM: `scripts/slurm/r01_filter_and_trim.slurm`
- Inputs: `data/raw/fastq_files/*_1.fastq.gz`, `data/raw/fastq_files/*_2.fastq.gz`
- Outputs: `data/intermediate/filtered_fastq/*_F_filt.fastq.gz`,
  `data/intermediate/filtered_fastq/*_R_filt.fastq.gz`,
  `data/intermediate/filtering_summary.rds`

2) Learn error rates
- Script: `scripts/02_learn_errors.R`
- SLURM: `scripts/slurm/r02_learn_errors.slurm`
- Inputs: `data/intermediate/filtered_fastq/*`
- Outputs: `data/intermediate/errF.rds`, `data/intermediate/errR.rds`

3) Denoise reads
- Script: `scripts/03_denoise.R`
- SLURM: `scripts/slurm/r03_denoise.slurm`
- Inputs: `data/intermediate/filtered_fastq/*`, `data/intermediate/errF.rds`,
  `data/intermediate/errR.rds`
- Outputs: `data/intermediate/dadaFs.rds`, `data/intermediate/dadaRs.rds`

4) Merge paired reads
- Script: `scripts/04_merge.R`
- SLURM: `scripts/slurm/r04_merge.slurm`
- Inputs: `data/intermediate/dadaFs.rds`, `data/intermediate/dadaRs.rds`,
  `data/intermediate/filtered_fastq/*`
- Outputs: `data/intermediate/mergers.rds`

5) Build ASV table and remove chimeras
- Script: `scripts/05_asv_table_and_chimera.R`
- SLURM: `scripts/slurm/r05_asv_table_and_chimera.slurm`
- Inputs: `data/intermediate/mergers.rds`
- Outputs: `data/intermediate/seqtab.rds`, `data/intermediate/seqtab_nochim.rds`

6) Remove Eukaryota
- Script: `scripts/06_remove_eukaryota.R`
- SLURM: `scripts/slurm/r06_remove_eukaryota.slurm`
- Inputs: `data/intermediate/seqtab_nochim.rds`, `data/intermediate/taxa.rds`
- Outputs: `data/intermediate/seqtab_nochim_bacteria.rds`,
  `data/intermediate/taxa_bacteria.rds`,
  `data/intermediate/taxonomy_table_bacteria.csv`

7) Read tracking summary
- Script: `scripts/07_read_tracking.R`
- SLURM: `scripts/slurm/r07_read_tracking.slurm`
- Inputs: `data/intermediate/*` from steps 1–5
- Outputs: `outputs/qc_tables/dada2_read_tracking.csv`

8) Assign taxonomy with filtered classifier
- Script: `scripts/08_classify.R`
- SLURM: `scripts/slurm/r08_classify.slurm`
- Inputs: `data/intermediate/seqtab_nochim.rds`,
  `data/reference/silva_nr99_v138.2_toGenus_trainset_filtered.fa`
- Outputs:
  - `data/intermediate/taxa_filtered_classifier.rds`
  - `data/intermediate/taxonomy_table_filtered_classifier.csv`
  - `data/intermediate/sample_taxa_matrix_L6_filtered_classifier.csv`
  - `outputs/figures/genus_level_barplot_filtered_classifier.png`

9) Assign taxonomy (bacteria only)
- Script: `scripts/09_classify_filtered.R`
- SLURM: `scripts/slurm/r09_classify_filtered.slurm`
- Inputs: `data/intermediate/seqtab_nochim_bacteria.rds`,
  `data/intermediate/taxa_bacteria.rds`
- Outputs:
  - `data/intermediate/taxonomy_table_bacteria.csv`
  - `data/intermediate/sample_taxa_matrix_L6_bacteria.csv`
  - `outputs/figures/genus_level_barplot_bacteria.png`

10) QC stats
- Script: `scripts/10_stats.R`
- SLURM: `scripts/slurm/r10_stats.slurm`
- Inputs: `data/intermediate/sample_taxa_matrix_L6_bacteria.csv`
- Outputs: `outputs/qc_tables/*`, `outputs/figures/*`

11) ASV histogram
- Script: `scripts/11_asv_histogram.R`
- SLURM: `scripts/slurm/r11_asv_histogram.slurm`
- Inputs: `data/intermediate/seqtab_nochim.rds`
- Outputs: `outputs/figures/soil_hist_asv_counts_per_sample.png`

12) Taxa prevalence summary
- Script: `scripts/12_taxa_presence_summary.R`
- SLURM: `scripts/slurm/r12_taxa_presence_summary.slurm`
- Inputs: `data/intermediate/sample_taxa_matrix_L6_bacteria.csv`
- Outputs:
  - `data/processed/taxa_presence_summary_soil.csv`
  - `outputs/figures/taxa_presence_histogram_soil.png`

13) Taxa prevalence thresholds
- Script: `scripts/13_taxa_presence_thresholds.R`
- SLURM: `scripts/slurm/r13_taxa_presence_thresholds.slurm`
- Inputs: `data/intermediate/sample_taxa_matrix_L6_bacteria.csv`
- Outputs: `data/processed/soil_fraction_taxa_remaining_1_to_20.csv`

14) Bacteria-only prevalence thresholds
- Script: `scripts/14_taxa_presence_thresholds_bacteria.R`
- SLURM: `scripts/slurm/r14_taxa_thresholds_bacteria.slurm`
- Inputs: `data/intermediate/taxonomy_table_bacteria.csv`
- Outputs: `data/processed/fraction_taxa_remaining_bacteria.csv`

15) Lineage counts
- Script: `scripts/15_taxonomy_lineage_counts.R`
- SLURM: `scripts/slurm/r15_taxonomy_lineage_counts.slurm`
- Inputs: `data/intermediate/taxonomy_table_bacteria.csv`
- Outputs: `data/processed/taxonomy_full_lineage_counts.csv`

16) Lineage thresholds
- Script: `scripts/16_lineage_thresholds.R`
- SLURM: `scripts/slurm/r16_lineage_thresholds.slurm`
- Inputs: `data/processed/taxonomy_full_lineage_counts.csv`
- Outputs: `data/processed/soil_fraction_taxa_remaining_1_to_20.csv`
  - Note: this output filename overlaps with step 13.

17) Lineage prevalence thresholds
- Script: `scripts/17_lineage_prevalence_thresholds.R`
- SLURM: `scripts/slurm/r17_lineage_prevalence_thresholds.slurm`
- Inputs: `data/intermediate/taxonomy_table_bacteria.csv`,
  `data/intermediate/sample_taxa_matrix_L6_bacteria.csv`
- Outputs:
  - `data/processed/asv_fraction_vs_lineage_threshold.csv`
  - `outputs/figures/fraction_lineages_remaining_plot.png`
  - `outputs/figures/fraction_matrix_remaining_plot.png`

18) Build lineage CLR matrix (analysis-ready)
- Script: `scripts/18_build_lineage_clr_matrix.R`
- SLURM: `scripts/slurm/r18_build_lineage_clr_matrix.slurm`
- Inputs: `data/intermediate/seqtab_nochim_bacteria.rds`,
  `data/intermediate/taxa_bacteria.rds`
- Outputs: `data/processed/lineage_clr_matrix_soil_filtered.csv`

19) Modeling
- Script: `scripts/19_model_fit_mixture.R`
- SLURM: `scripts/slurm/r19_model_fit_mixture.slurm`
- Inputs: `data/processed/lineage_clr_matrix_soil_filtered.csv`
- Outputs: `models/best_model_soil_k20.RData`

20) Cluster metadata comparison
- Script: `scripts/20_cluster_metadata_compare.R`
- SLURM: `scripts/slurm/r20_cluster_metadata_compare.slurm`
- Inputs:
  - `models/best_model_soil_k20.RData`
  - `data/raw/sraRunTable.csv`
  - `data/processed/lineage_clr_matrix_soil_filtered.csv`
- Outputs:
  - `outputs/cluster_metadata/cluster_assignments.tsv`
  - `outputs/cluster_metadata/cluster_metadata_join.tsv`
  - `outputs/cluster_metadata/cluster_sizes.tsv`
  - `data/intermediate/soil_metadata.tsv`

21) Cluster CLR averages
- Script: `scripts/21_cluster_clr_averages.R`
- SLURM: `scripts/slurm/r21_cluster_clr_averages.slurm`
- Inputs:
  - `data/processed/lineage_clr_matrix_soil_filtered.csv`
  - `outputs/cluster_metadata/cluster_assignments.tsv`
- Outputs:
  - `outputs/cluster_metadata/cluster_avg_clr_matrix.csv`
  - `outputs/cluster_metadata/cluster_avg_clr_long.tsv`
  - `outputs/cluster_metadata/cluster_avg_clr_summary.tsv`

22) Plot cluster metadata
- Script: `scripts/22_plot_cluster_metadata.py`
- SLURM: `scripts/slurm/r22_plot_cluster_metadata.slurm`
- Inputs: `outputs/cluster_metadata/cluster_metadata_join.tsv`
- Outputs: `outputs/cluster_metadata/plots/*.pdf` (multiple PDF plots)

23) Calculate alpha diversity metrics
- Script: `scripts/23_calculate_alpha_diversity.R`
- SLURM: `scripts/slurm/r23_calculate_alpha_diversity.slurm`
- Inputs: `data/intermediate/seqtab_nochim_bacteria.rds`
- Outputs:
  - `data/processed/alpha_diversity_soil.csv` (Observed OTUs, Chao1, Shannon)
  - `outputs/figures/alpha_diversity_distributions.png`

24) Plot alpha diversity by cluster
- Script: `scripts/24_plot_alpha_diversity_by_cluster.R`
- SLURM: `scripts/slurm/r24_plot_alpha_diversity_by_cluster.slurm`
- Inputs:
  - `data/processed/alpha_diversity_soil.csv`
  - `outputs/cluster_metadata/cluster_assignments.tsv`
- Outputs:
  - `outputs/cluster_metadata/cluster_alpha_diversity.tsv`
  - `outputs/cluster_metadata/cluster_alpha_diversity_summary.tsv`
  - `outputs/cluster_metadata/plots/alpha_div_*.pdf` (7 PDF plots: 3 boxplots, 3 scatter, 1 heatmap)

## Final outputs of record
- `data/processed/lineage_clr_matrix_soil_filtered.csv`
- `data/processed/taxa_presence_summary_soil.csv`
- `data/processed/asv_fraction_vs_lineage_threshold.csv`
- `data/processed/alpha_diversity_soil.csv` (alpha diversity metrics)
- `models/best_model_soil_k20.RData`
- `outputs/cluster_metadata/cluster_assignments.tsv`
- `outputs/cluster_metadata/cluster_alpha_diversity.tsv` (alpha diversity by cluster)

## Legacy
- QIIME scripts and artifacts are retained under `legacy/qiime`.
- OS artifacts are under `legacy/osx_artifacts`.

25) Enhanced Cluster CLR Statistics
- Script: `scripts/25_cluster_clr_statistics.R`
- SLURM: `scripts/slurm/r25_cluster_clr_statistics.slurm`
- Inputs:
  - `data/processed/lineage_clr_matrix_soil_filtered.csv`
  - `outputs/cluster_metadata/cluster_assignments.tsv`
- Outputs:
  - `outputs/cluster_metadata/cluster_clr_statistics.csv` (mean, SD, SE, CI per lineage per cluster)
  - `outputs/cluster_metadata/cluster_clr_statistics_summary.tsv`
  - `outputs/ordination/cluster_errorbar_top20_taxa.pdf`
  - `outputs/ordination/cluster_errorbar_selected_phyla.pdf`
  - `outputs/ordination/cluster_variance_heatmap.pdf`

26) PCA/PCoA Ordination Analysis
- Script: `scripts/26_ordination_analysis.R`
- SLURM: `scripts/slurm/r26_ordination_analysis.slurm`
- Inputs:
  - `data/processed/lineage_clr_matrix_soil_filtered.csv`
  - `outputs/cluster_metadata/cluster_assignments.tsv`
- Outputs:
  - Data files: `outputs/ordination/pca_results.rds`, `pca_sample_scores.csv`, `pca_variance_explained.csv`, `pca_loadings.csv`, `pca_loadings_top50_pc{1,2}.csv`, `pcoa_results.rds`, `pcoa_sample_scores.csv`
  - Plots: `outputs/ordination/pca_scree_plot.pdf`, `pca_biplot_pc1_pc2.pdf`, `pca_biplot_pc1_pc3.pdf`, `pca_biplot_pc2_pc3.pdf`, `pca_loadings_pc1_top30.pdf`, `pca_loadings_pc2_top30.pdf`, `pca_loadings_biplot.pdf`, `pcoa_biplot_axis1_axis2.pdf`
