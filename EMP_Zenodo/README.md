# EMP_Zenodo_Matthew (R-only pipeline)

This folder contains the EMP soil-only processing pipeline in R. Legacy QIIME
artifacts are preserved under `legacy/qiime` but are not part of the final R
pipeline.

## Layout
- `data/raw`: original EMP inputs
- `data/reference`: reference DBs/trees
- `data/intermediate`: derived tables used between steps
- `data/processed`: analysis-ready outputs (final matrices, prevalence summaries)
- `outputs/qc_tables`: QC summaries
- `outputs/figures`: plots
- `scripts`: R pipeline scripts
- `scripts/slurm`: SLURM runners for the R pipeline (use these to run R scripts)
- `models`: model outputs
- `logs`: SLURM logs
- `legacy`: QIIME outputs and OS artifacts
- `docs/PIPELINE.md`: step-by-step pipeline with inputs/outputs
- `docs/FILE_MAP.tsv`: mapping from original names to cleaned names

## R pipeline (high level)
Run all R steps through SLURM wrappers in `scripts/slurm`:
1. `scripts/slurm/r01_extract_soil_sample_ids.slurm`
2. `scripts/slurm/r01b_extract_soil_metadata.slurm` (inline R; optional helper)
3. `scripts/slurm/r02_assign_taxonomy_and_soil_asv_table.slurm`
4. `scripts/slurm/r03_collapse_genus_all_samples.slurm`
5. `scripts/slurm/r07_collapse_genus_soil_only.slurm`
6. QC and stats:
   - `scripts/slurm/r04_qc_read_counts_vs_genus.slurm`
   - `scripts/slurm/r05_qc_asv_id_match.slurm`
   - `scripts/slurm/r06_stats_all_samples.slurm`
   - `scripts/slurm/r08_stats_soil_only.slurm`
7. Prevalence/thresholds:
   - `scripts/slurm/r10_genus_prevalence_summary.slurm`
   - `scripts/slurm/r11_genus_prevalence_thresholds.slurm`
   - `scripts/slurm/r12_lineage_counts_thresholds.slurm`
   - `scripts/slurm/r13_lineage_prevalence_thresholds.slurm`
8. `scripts/slurm/r09_build_lineage_clr_matrix.slurm`
9. `scripts/slurm/r14_model_fit_mixture.slurm`
10. `scripts/slurm/r15_cluster_metadata_compare.slurm`
11. `scripts/slurm/r16_cluster_clr_averages.slurm`
12. `scripts/slurm/r17_plot_cluster_metadata.slurm`
13. Alpha diversity:
   - `scripts/slurm/r18_calculate_alpha_diversity.slurm`
   - `scripts/slurm/r19_plot_alpha_diversity_by_cluster.slurm`
14. Cluster statistics and ordination analysis:
   - `scripts/slurm/r20_cluster_clr_statistics.slurm`
   - `scripts/slurm/r21_ordination_analysis.slurm`

## Final outputs (typical)
- `data/intermediate/genus_matrix_soil_only.csv` (sample x genus)
- `data/processed/lineage_clr_matrix_soil_filtered.csv`
- `data/processed/genus_prevalence_summary.csv`
- `data/processed/lineage_prevalence_thresholds.csv`
- `data/processed/alpha_diversity_soil.csv` (alpha diversity metrics)
- `models/best_model_emp_k20.RData`
- `outputs/cluster_metadata/cluster_metadata_tests.tsv`
- `outputs/cluster_metadata/cluster_alpha_diversity.tsv`
- `outputs/cluster_metadata/cluster_clr_statistics.csv` (cluster mean/SD/SE)
- `outputs/figures/alpha_diversity_distributions.png`
- `outputs/cluster_metadata/plots/alpha_div_*.pdf` (7 plots)
- `outputs/ordination/pca_sample_scores.csv` (PCA results)
- `outputs/ordination/pca_variance_explained.csv`
- `outputs/ordination/*.pdf` (12 ordination plots: PCA/PCoA biplots, scree, loadings)

## Notes
- `data/intermediate/asv_taxonomy_matrix_all_samples.csv` is a large
  intermediate used by QC scripts. It is assumed to exist from a prior run.
- Legacy QIIME outputs are retained for reference only under `legacy/qiime`.
