# global_topsoil_matthew (R-only pipeline)

This folder contains a soil 16S processing and modeling pipeline in R (DADA2 +
post-processing). Legacy QIIME artifacts are preserved under `legacy/qiime` and
are not part of the R pipeline.

## Layout
- `data/raw`: raw inputs (FASTQ, manifests, metadata)
- `data/reference`: classifier/reference files
- `data/intermediate`: DADA2 intermediates (seq tables, taxonomy tables, etc.)
- `data/processed`: analysis-ready outputs (CLR matrices, prevalence summaries)
- `outputs/qc_tables`: QC summaries and tables
- `outputs/figures`: plots and histograms
- `scripts`: R pipeline scripts
- `scripts/slurm`: SLURM runners for the R pipeline
- `models`: model outputs
- `logs`: SLURM logs
- `legacy`: QIIME outputs and OS artifacts
- `docs/PIPELINE.md`: step-by-step pipeline with inputs/outputs
- `docs/FILE_MAP.tsv`: mapping from original names to cleaned names

## R pipeline (high level)
Run R steps through SLURM wrappers in `scripts/slurm`:
1. `scripts/slurm/r00_filter_classifier.slurm` (optional)
2. `scripts/slurm/r01_filter_and_trim.slurm`
3. `scripts/slurm/r02_learn_errors.slurm`
4. `scripts/slurm/r03_denoise.slurm`
5. `scripts/slurm/r04_merge.slurm`
6. `scripts/slurm/r05_asv_table_and_chimera.slurm`
7. `scripts/slurm/r06_remove_eukaryota.slurm`
8. `scripts/slurm/r07_read_tracking.slurm`
9. `scripts/slurm/r08_classify.slurm`
10. `scripts/slurm/r09_classify_filtered.slurm`
11. `scripts/slurm/r10_stats.slurm`
12. `scripts/slurm/r11_asv_histogram.slurm`
13. `scripts/slurm/r12_taxa_presence_summary.slurm`
14. `scripts/slurm/r13_taxa_presence_thresholds.slurm`
15. `scripts/slurm/r14_taxa_thresholds_bacteria.slurm`
16. `scripts/slurm/r15_taxonomy_lineage_counts.slurm`
17. `scripts/slurm/r16_lineage_thresholds.slurm`
18. `scripts/slurm/r17_lineage_prevalence_thresholds.slurm`
19. `scripts/slurm/r18_build_lineage_clr_matrix.slurm`
20. `scripts/slurm/r19_model_fit_mixture.slurm`
21. `scripts/slurm/r20_cluster_metadata_compare.slurm`
22. `scripts/slurm/r21_cluster_clr_averages.slurm`
23. `scripts/slurm/r22_plot_cluster_metadata.slurm`
24. Alpha diversity:
    - `scripts/slurm/r23_calculate_alpha_diversity.slurm`
    - `scripts/slurm/r24_plot_alpha_diversity_by_cluster.slurm`
25. Cluster statistics and ordination analysis:
    - `scripts/slurm/r25_cluster_clr_statistics.slurm`
    - `scripts/slurm/r26_ordination_analysis.slurm`

## Final outputs (typical)
- `data/processed/lineage_clr_matrix_soil_filtered.csv`
- `data/processed/taxa_presence_summary_soil.csv`
- `data/processed/asv_fraction_vs_lineage_threshold.csv`
- `data/processed/alpha_diversity_soil.csv` (alpha diversity metrics)
- `models/best_model_soil_k20.RData`
- `outputs/cluster_metadata/cluster_metadata_tests.tsv`
- `outputs/cluster_metadata/cluster_alpha_diversity.tsv`
- `outputs/cluster_metadata/cluster_clr_statistics.csv` (cluster mean/SD/SE)
- `outputs/figures/alpha_diversity_distributions.png`
- `outputs/cluster_metadata/plots/alpha_div_*.pdf` (7 plots)
- `outputs/ordination/pca_sample_scores.csv` (PCA results)
- `outputs/ordination/pca_variance_explained.csv`
- `outputs/ordination/*.pdf` (12 ordination plots: PCA/PCoA biplots, scree, loadings)

## Notes
- QIIME artifacts and QIIME SLURM scripts are retained under `legacy/`.
- Some QC outputs are placeholders (see `scripts/10_stats.R`).
