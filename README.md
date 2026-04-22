# Thesis — Microbial Community Structure Across Global Soil Environments

Analysis code and supplementary data for my undergraduate thesis examining microbial community structure in global soil environments using 16S amplicon sequencing.

## Datasets

### Global Topsoil (Bahram et al. 2018)
DADA2-based re-analysis of the Bahram et al. 2018 *Nature* dataset (189 topsoil sites, 7,560 subsamples across all terrestrial biomes).

- **EBI Accession:** PRJEB19856 (ERP021922)
- **Reference:** Bahram M. et al. (2018). Structure and function of the global topsoil microbiome. *Nature*, 560, 233–237. https://doi.org/10.1038/s41586-018-0386-6
- See [`global_topsoil/`](global_topsoil/) for pipeline and scripts.

### Earth Microbiome Project (EMP)
Re-analysis of EMP Zenodo soil samples using an equivalent mixture model pipeline for cross-dataset comparison.

- **Source:** EMP Release 1 (Zenodo)
- See [`EMP_Zenodo/`](EMP_Zenodo/) for pipeline and scripts.

## Repository Structure

```
├── figures/                 # All manuscript and supplementary figures
├── global_topsoil/          # Global Topsoil dataset
│   ├── docs/                # PIPELINE.md, FILE_MAP.tsv
│   └── scripts/             # R scripts (00–31) + SLURM wrappers
│       └── slurm/
├── EMP_Zenodo/              # Earth Microbiome Project dataset
│   ├── docs/                # PIPELINE.md, FILE_MAP.tsv
│   └── scripts/             # R scripts (01–26) + SLURM wrappers
│       └── slurm/
├── supplementary_tables/    # Supplementary Tables 1–6
├── make_paper_figures.py    # Generates all manuscript figures
├── rf_biome_stats.csv       # Random forest biome classification statistics
└── main.tex                 # Thesis manuscript (LaTeX)
```

## Pipeline Overview

Both datasets follow the same general pipeline:

1. **Pre-processing** — Quality filtering, DADA2 denoising (Global Topsoil) or biom import (EMP)
2. **Taxonomy** — SILVA v138.2 classification
3. **Lineage filtering** — Prevalence threshold sweep → CLR transformation
4. **Mixture modeling** — Dirichlet mixture model (k=20 components)
5. **Downstream analysis** — Cluster metadata, alpha diversity, ordination (PCA/PCoA), t-SNE, random forest biome classification

See each dataset's `docs/PIPELINE.md` for step-by-step details.

## Requirements

**R packages:** `dada2`, `vegan`, `mixtools`, `randomForest`, `ggplot2`, `phyloseq`

**Python packages:** `matplotlib`, `seaborn`, `pandas`, `numpy`, `scikit-learn`

**HPC:** Scripts in `scripts/slurm/` are written for SLURM job schedulers.
