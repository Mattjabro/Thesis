#!/usr/bin/env python3
"""
Generate publication-quality figures for thesis Results section.
Outputs saved to:
  global_topsoil_matthew/outputs/paper_figures/
  EMP_Zenodo_Matthew/outputs/paper_figures/
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
from collections import Counter

# ── paths ──────────────────────────────────────────────────────────────────
BASE_TOPSOIL = "/hopper/groups/yoosephlab/datasets/global_topsoil_matthew"
BASE_EMP     = "/hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew"

OUT_TOPSOIL  = os.path.join(BASE_TOPSOIL, "outputs", "paper_figures")
OUT_EMP      = os.path.join(BASE_EMP,     "outputs", "paper_figures")
os.makedirs(OUT_TOPSOIL, exist_ok=True)
os.makedirs(OUT_EMP,     exist_ok=True)

# ── shared style ────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.size":        10,
    "axes.titlesize":   11,
    "axes.labelsize":   10,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "legend.fontsize":  8,
    "figure.dpi":       150,
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "axes.linewidth":   0.8,
    "xtick.major.width":0.8,
    "ytick.major.width":0.8,
})

# ── colour palettes ─────────────────────────────────────────────────────────
TOPSOIL_BIOME_COLORS = {
    "TROPICAL_FOREST":       "#2ca02c",
    "TEMPERATE_FOREST":      "#1f77b4",
    "BOREAL_FOREST":         "#9467bd",
    "GRASSLAND_SHRUBLAND":   "#e6ac00",
    "MEDITERRANEAN":         "#d62728",
    "ARCTIC_TUNDRA":         "#17becf",
    "-":                     "#aec7e8",
    "UNKNOWN":               "#aec7e8",
}

EMP_BIOME_COLORS = {
    "AGRICULTURAL":          "#8c564b",
    "AQUATIC_FRESHWATER":    "#17becf",
    "AQUATIC_MARINE":        "#1f77b4",
    "COLD_TUNDRA_POLAR":     "#9467bd",
    "DESERT":                "#e6ac00",
    "FOREST":                "#2ca02c",
    "GRASSLAND_SHRUBLAND":   "#bcbd22",
    "TROPICAL_FOREST":       "#d62728",
    "URBAN_DEVELOPED":       "#7f7f7f",
    "OTHER":                 "#aec7e8",
    "UNKNOWN":               "#aec7e8",
}

# 20-colour cluster palette
CLUSTER_COLORS = [
    "#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f",
    "#edc948","#b07aa1","#ff9da7","#9c755f","#bab0ac",
    "#d4a6c8","#86bcb6","#f1ce63","#a0cbe8","#ffbe7d",
    "#8cd17d","#b6992d","#499894","#e15759","#d37295",
]

def cluster_color_map(clusters):
    unique = sorted(set(clusters))
    return {c: CLUSTER_COLORS[i % len(CLUSTER_COLORS)] for i, c in enumerate(unique)}

# ════════════════════════════════════════════════════════════════════════════
# FIGURE 1  –  PCA biplots (Topsoil PC1/PC2, EMP PC1/PC2), coloured by biome
# ════════════════════════════════════════════════════════════════════════════
print("Generating Figure 1: PCA biplots (both datasets) …")

# ── load data ─────────────────────────────────────────────────────────────
ts_pca  = pd.read_csv(os.path.join(BASE_TOPSOIL,"outputs","ordination","pca_sample_scores.csv"))
ts_tsne = pd.read_csv(os.path.join(BASE_TOPSOIL,"outputs","tsne","tsne_coordinates.csv"))
ts_pca  = ts_pca.merge(
    ts_tsne[["SampleID","env_biome","group_biome"]],
    on="SampleID", how="left")
ts_pca["group_biome"] = ts_pca["group_biome"].fillna("-")

emp_pca  = pd.read_csv(os.path.join(BASE_EMP,"outputs","ordination","pca_sample_scores.csv"))
emp_tsne = pd.read_csv(os.path.join(BASE_EMP,"outputs","tsne","tsne_coordinates.csv"))
emp_pca  = emp_pca.merge(
    emp_tsne[["SampleID","env_biome","group_biome"]],
    on="SampleID", how="left")
emp_pca["group_biome"] = emp_pca["group_biome"].fillna("OTHER")

ts_var  = pd.read_csv(os.path.join(BASE_TOPSOIL,"outputs","ordination","pca_variance_explained.csv"))
emp_var = pd.read_csv(os.path.join(BASE_EMP,"outputs","ordination","pca_variance_explained.csv"))

def biome_label(raw):
    return raw.replace("_"," ").title()

fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

for ax, pca_df, var_df, color_map, title in [
    (axes[0], ts_pca,  ts_var,  TOPSOIL_BIOME_COLORS,
     "(a) Global Topsoil"),
    (axes[1], emp_pca, emp_var, EMP_BIOME_COLORS,
     "(b) EMP Soil Subset"),
]:
    pc1_pct = var_df.loc[0,"prop_variance"]*100
    pc2_pct = var_df.loc[1,"prop_variance"]*100

    for biome, grp in pca_df.groupby("group_biome"):
        col = color_map.get(biome, "#aec7e8")
        ax.scatter(grp["PC1"], grp["PC2"],
                   c=col, s=18, alpha=0.75, linewidths=0,
                   label=biome_label(biome), zorder=3)
        if len(grp) >= 5:
            from matplotlib.patches import Ellipse
            import matplotlib.transforms as transforms
            mu    = grp[["PC1","PC2"]].mean()
            cov   = grp[["PC1","PC2"]].cov()
            vals, vecs = np.linalg.eigh(cov.values)
            order = vals.argsort()[::-1]
            vals, vecs = vals[order], vecs[:,order]
            angle = np.degrees(np.arctan2(*vecs[:,0][::-1]))
            w, h  = 2*np.sqrt(vals * 5.991)   # 95 % CI chi2(2) ≈ 5.991
            ell = Ellipse(xy=(mu["PC1"],mu["PC2"]),
                          width=w, height=h, angle=angle,
                          color=col, fill=False, lw=0.8,
                          linestyle="--", zorder=2)
            ax.add_patch(ell)

    ax.axhline(0, color="grey", lw=0.5, ls=":")
    ax.axvline(0, color="grey", lw=0.5, ls=":")
    ax.set_xlabel(f"PC1 ({pc1_pct:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pc2_pct:.1f}% variance)")
    ax.set_title(title, fontweight="bold")

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),
              loc="best", frameon=True, framealpha=0.8,
              edgecolor="0.7", fontsize=7.5,
              markerscale=1.3, handletextpad=0.4, borderpad=0.5)

fig.suptitle(
    "PCA of CLR-transformed microbial composition profiles\n"
    "coloured by grouped biome",
    fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_pca_both_datasets.pdf"),
            bbox_inches="tight")
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_pca_both_datasets.png"),
            bbox_inches="tight", dpi=200)
plt.close(fig)
print("  → fig_pca_both_datasets.pdf/.png")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 2  –  t-SNE coloured by cluster (both datasets side-by-side)
# ════════════════════════════════════════════════════════════════════════════
print("Generating Figure 2: t-SNE biplots …")

ts_tsne2  = pd.read_csv(os.path.join(BASE_TOPSOIL,"outputs","tsne","tsne_coordinates.csv"))
emp_tsne2 = pd.read_csv(os.path.join(BASE_EMP,    "outputs","tsne","tsne_coordinates.csv"))

fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

for ax, tsne_df, title in [
    (axes[0], ts_tsne2,  "(a) Global Topsoil"),
    (axes[1], emp_tsne2, "(b) EMP Soil Subset"),
]:
    tsne_df = tsne_df.copy()
    tsne_df["cluster"] = tsne_df["cluster"].astype(str)
    cmap = cluster_color_map(tsne_df["cluster"].unique())

    for cl, grp in tsne_df.groupby("cluster"):
        ax.scatter(grp["tSNE1"], grp["tSNE2"],
                   c=cmap[cl], s=14, alpha=0.70, linewidths=0,
                   label=f"C{cl}", zorder=3)

    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
    ax.set_title(title, fontweight="bold")

    ncol = 4 if len(cmap) > 10 else 2
    ax.legend(loc="best", frameon=True, framealpha=0.8,
              edgecolor="0.7", fontsize=7, ncol=ncol,
              markerscale=1.2, handletextpad=0.3, borderpad=0.5)

fig.suptitle(
    "t-SNE projection of CLR-transformed profiles coloured by MixGGM cluster",
    fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_tsne_both_datasets.pdf"),
            bbox_inches="tight")
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_tsne_both_datasets.png"),
            bbox_inches="tight", dpi=200)
plt.close(fig)
print("  → fig_tsne_both_datasets.pdf/.png")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 3  –  Alpha diversity (Shannon) by cluster — both datasets
# ════════════════════════════════════════════════════════════════════════════
print("Generating Figure 3: Alpha diversity by cluster …")

ts_alpha  = pd.read_csv(os.path.join(BASE_TOPSOIL,"outputs","cluster_metadata",
                                      "cluster_alpha_diversity.tsv"), sep="\t")
emp_alpha = pd.read_csv(os.path.join(BASE_EMP,"outputs","cluster_metadata",
                                      "cluster_alpha_diversity.tsv"), sep="\t")

METRICS = [
    ("adiv_observed_otus", "Observed OTUs"),
    ("adiv_chao1",         "Chao1 richness"),
    ("adiv_shannon",       "Shannon diversity index"),
]

fig, axes = plt.subplots(3, 2, figsize=(14, 12))

for col_idx, (alpha_df, col_title) in enumerate([
    (ts_alpha,  "Global Topsoil"),
    (emp_alpha, "EMP Soil Subset"),
]):
    alpha_df = alpha_df.copy()
    alpha_df["cluster"] = alpha_df["cluster"].astype(str)
    order = sorted(alpha_df["cluster"].unique(),
                   key=lambda x: int(x) if x.isdigit() else x)
    ns = [len(alpha_df[alpha_df["cluster"] == c]) for c in order]
    pal = plt.cm.tab20(np.linspace(0, 1, len(order)))

    for row_idx, (metric_col, metric_label) in enumerate(METRICS):
        ax = axes[row_idx][col_idx]
        data = [alpha_df.loc[alpha_df["cluster"] == c, metric_col].dropna().values
                for c in order]

        bp = ax.boxplot(data, patch_artist=True, notch=False,
                        medianprops=dict(color="black", lw=1.5),
                        whiskerprops=dict(lw=0.8),
                        capprops=dict(lw=0.8),
                        flierprops=dict(marker="o", markersize=3,
                                        markerfacecolor="0.5", alpha=0.5))
        for patch, col in zip(bp["boxes"], pal):
            patch.set_facecolor(col)
            patch.set_alpha(0.75)

        ax.set_xticks(range(1, len(order) + 1))
        ax.set_xticklabels(
            [f"C{c} (n={n})" for c, n in zip(order, ns)],
            rotation=45, ha="right", fontsize=7)
        ax.set_ylabel(metric_label)
        ax.yaxis.grid(True, color="0.9", zorder=0)
        ax.set_xlabel("Cluster")

        if row_idx == 0:
            ax.set_title(f"({('a','b')[col_idx]}) {col_title}",
                         fontweight="bold")

fig.suptitle(
    "Alpha diversity per MixGGM cluster — Observed OTUs, Chao1, and Shannon",
    fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_alpha_diversity_by_cluster.pdf"),
            bbox_inches="tight")
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_alpha_diversity_by_cluster.png"),
            bbox_inches="tight", dpi=200)
plt.close(fig)
print("  → fig_alpha_diversity_by_cluster.pdf/.png")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 4  –  Cluster size distributions (both datasets)
# ════════════════════════════════════════════════════════════════════════════
print("Generating Figure 4: Cluster size distributions …")

ts_sizes  = pd.read_csv(os.path.join(BASE_TOPSOIL,"outputs","cluster_metadata",
                                      "cluster_sizes.tsv"), sep="\t")
emp_sizes = pd.read_csv(os.path.join(BASE_EMP,"outputs","cluster_metadata",
                                      "cluster_sizes.tsv"), sep="\t")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, sz, title in [
    (axes[0], ts_sizes,  "(a) Global Topsoil  (n = 193, k = 19 non-empty)"),
    (axes[1], emp_sizes, "(b) EMP Soil Subset  (n = 2,209, k = 10 non-empty)"),
]:
    sz = sz.sort_values("n_samples", ascending=False).reset_index(drop=True)
    pal = plt.cm.tab20(np.linspace(0,1,len(sz)))
    bars = ax.bar(range(len(sz)), sz["n_samples"],
                  color=pal, edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(len(sz)))
    ax.set_xticklabels([f"C{int(c)}" for c in sz["cluster"]],
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Number of samples")
    ax.set_xlabel("Cluster")
    ax.set_title(title, fontweight="bold")
    ax.yaxis.grid(True, color="0.9", zorder=0)
    for bar, n in zip(bars, sz["n_samples"]):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5,
                str(n), ha="center", va="bottom", fontsize=7)

fig.suptitle("MixGGM cluster size distributions", fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_cluster_sizes.pdf"),
            bbox_inches="tight")
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_cluster_sizes.png"),
            bbox_inches="tight", dpi=200)
plt.close(fig)
print("  → fig_cluster_sizes.pdf/.png")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 5  –  RF feature importance (top 20), both datasets side-by-side
# ════════════════════════════════════════════════════════════════════════════
print("Generating Figure 5: RF feature importance (top 20) …")

def short_lineage(lin, maxlen=40):
    parts = lin.split(";")
    genus = parts[-1].strip() if parts else lin
    phylum = parts[1].strip() if len(parts)>1 else ""
    label = f"{genus}  [{phylum}]"
    return label if len(label)<=maxlen else label[:maxlen-1]+"…"

ts_imp  = pd.read_csv(os.path.join(BASE_TOPSOIL,"outputs","ml_models",
                                    "rf_grouped_biomes","feature_importance.csv"))
emp_imp = pd.read_csv(os.path.join(BASE_EMP,"outputs","ml_models",
                                    "rf_grouped_biomes","feature_importance.csv"))

fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))

for ax, imp_df, title, color in [
    (axes[0], ts_imp,  "(a) Global Topsoil (grouped biomes)", "#4e79a7"),
    (axes[1], emp_imp, "(b) EMP Soil Subset (grouped biomes)", "#e15759"),
]:
    top = imp_df.head(20).copy()
    top["label"] = top["lineage"].apply(short_lineage)
    top = top.sort_values("MeanDecreaseAccuracy")

    ax.barh(range(len(top)), top["MeanDecreaseAccuracy"],
            color=color, alpha=0.80, edgecolor="white", linewidth=0.4)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["label"], fontsize=8)
    ax.set_xlabel("Mean decrease in accuracy")
    ax.set_title(title, fontweight="bold")
    ax.xaxis.grid(True, color="0.9", zorder=0)

fig.suptitle(
    "Top 20 lineages by random forest permutation importance\n"
    "(grouped biome classifiers)",
    fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_rf_feature_importance.pdf"),
            bbox_inches="tight")
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_rf_feature_importance.png"),
            bbox_inches="tight", dpi=200)
plt.close(fig)
print("  → fig_rf_feature_importance.pdf/.png")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 6  –  Biome composition per cluster (stacked bar, both datasets)
# ════════════════════════════════════════════════════════════════════════════
print("Generating Figure 6: Biome composition per cluster …")

ts_meta  = pd.read_csv(os.path.join(BASE_TOPSOIL,"outputs","tsne",
                                     "tsne_coordinates.csv"))
emp_meta = pd.read_csv(os.path.join(BASE_EMP,"outputs","tsne",
                                     "tsne_coordinates.csv"))

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, meta, color_map, title in [
    (axes[0], ts_meta,  TOPSOIL_BIOME_COLORS, "(a) Global Topsoil"),
    (axes[1], emp_meta, EMP_BIOME_COLORS,     "(b) EMP Soil Subset"),
]:
    meta = meta.copy()
    meta["cluster"] = meta["cluster"].astype(str)
    meta["group_biome"] = meta["group_biome"].fillna("OTHER")

    ct = pd.crosstab(meta["cluster"], meta["group_biome"])
    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

    # order clusters by size
    order = ct.sum(axis=1).sort_values(ascending=False).index.tolist()
    ct_pct = ct_pct.loc[order]

    biomes = ct_pct.columns.tolist()
    colors = [color_map.get(b, "#cccccc") for b in biomes]

    bottoms = np.zeros(len(ct_pct))
    for biome, col in zip(biomes, colors):
        vals = ct_pct[biome].values
        ax.bar(range(len(ct_pct)), vals, bottom=bottoms,
               color=col, label=biome_label(biome),
               edgecolor="white", linewidth=0.3)
        bottoms += vals

    ax.set_xticks(range(len(ct_pct)))
    cluster_sizes = ct.loc[order].sum(axis=1)
    ax.set_xticklabels(
        [f"C{c}" for c in order],
        rotation=45, ha="right", fontsize=8)
    # annotate n= above each bar
    for i, c in enumerate(order):
        n = int(cluster_sizes[c])
        ax.text(i, 102, f"n={n}", ha="center", va="bottom",
                fontsize=5.5, color="0.3", rotation=90)
    ax.set_ylabel("Fraction of samples (%)")
    ax.set_xlabel("Cluster")
    ax.set_title(title, fontweight="bold")
    ax.set_ylim(0, 118)
    ax.yaxis.grid(True, color="0.9", zorder=0)

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),
              loc="upper right", frameon=True, framealpha=0.85,
              edgecolor="0.7", fontsize=7, ncol=1,
              markerscale=1.1)

fig.suptitle(
    "Biome composition within each MixGGM cluster\n"
    "(clusters ordered by size)",
    fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_biome_per_cluster.pdf"),
            bbox_inches="tight")
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_biome_per_cluster.png"),
            bbox_inches="tight", dpi=200)
plt.close(fig)
print("  → fig_biome_per_cluster.pdf/.png")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 7  –  Scree plots (both datasets, side-by-side)
# ════════════════════════════════════════════════════════════════════════════
print("Generating Figure 7: Scree plots …")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

for ax, var_df, title in [
    (axes[0], ts_var,  "(a) Global Topsoil"),
    (axes[1], emp_var, "(b) EMP Soil Subset"),
]:
    n = 10
    pcs   = var_df["PC"].values[:n]
    props = var_df["prop_variance"].values[:n] * 100
    cum   = var_df["cumulative_variance"].values[:n] * 100

    # single axis — both series in the same % units, no twin-axis distortion
    bars = ax.bar(range(n), props, color="#4e79a7", alpha=0.75,
                  edgecolor="white", linewidth=0.4, label="Per-PC variance")
    ax.plot(range(n), cum, "o-", color="#e15759",
            lw=1.5, ms=4, label="Cumulative variance")

    # % label above each bar, rotated to prevent overlap
    for i, (bar, pct) in enumerate(zip(bars, props)):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.8,
                f"{pct:.1f}%",
                ha="center", va="bottom",
                fontsize=7, rotation=90, color="#4e79a7")

    ax.set_ylim(0, 80)
    ax.set_xticks(range(n))
    ax.set_xticklabels(pcs, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Variance explained (%)")
    ax.set_xlabel("Principal component")
    ax.set_title(title, fontweight="bold")
    ax.yaxis.grid(True, color="0.9", zorder=0)
    ax.legend(loc="upper right", fontsize=7.5, frameon=True,
              framealpha=0.85, edgecolor="0.7")

fig.suptitle("PCA scree plots — variance explained by principal component",
             fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_pca_scree.pdf"),
            bbox_inches="tight")
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_pca_scree.png"),
            bbox_inches="tight", dpi=200)
plt.close(fig)
print("  → fig_pca_scree.pdf/.png")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 8  –  DADA2 read retention (Topsoil only)
# ════════════════════════════════════════════════════════════════════════════
print("Generating Figure 8: DADA2 read retention …")

track = pd.read_csv(os.path.join(BASE_TOPSOIL,"outputs","qc_tables",
                                  "dada2stats_allseqs.csv"))
track = track.sort_values("input").reset_index(drop=True)

stages = ["input","filtered","merged","non_chimeric"]
labels = ["Raw input","Quality-filtered","Merged pairs","Non-chimeric"]
colors = ["#aec7e8","#4e79a7","#f28e2b","#59a14f"]

fig, ax = plt.subplots(figsize=(8, 5))
for i, (col, lab, col_c) in enumerate(zip(stages, labels, colors)):
    if col in track.columns:
        ax.plot(range(len(track)), track[col]/1e3,
                lw=0.8, alpha=0.55, color=col_c)

for i, (col, lab, col_c) in enumerate(zip(stages, labels, colors)):
    if col in track.columns:
        med = track[col].median()/1e3
        ax.axhline(med, color=col_c, lw=2, ls="--",
                   label=f"{lab} (median {med*1e3:,.0f})")

ax.set_xlabel("Samples (sorted by input read count)")
ax.set_ylabel("Read count (×10³)")
ax.set_title("DADA2 read tracking — Global Topsoil dataset\n"
             "Dashed lines show per-stage medians", fontweight="bold")
ax.legend(fontsize=8, frameon=True, framealpha=0.85, edgecolor="0.7")
ax.yaxis.grid(True, color="0.9", zorder=0)
fig.tight_layout()
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_dada2_read_tracking.pdf"),
            bbox_inches="tight")
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_dada2_read_tracking.png"),
            bbox_inches="tight", dpi=200)
plt.close(fig)
print("  → fig_dada2_read_tracking.pdf/.png")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 9  –  Prevalence filter sweep curves (both datasets)
# ════════════════════════════════════════════════════════════════════════════
print("Generating Figure 9: Prevalence sweep curves …")

ts_prev  = pd.read_csv(os.path.join(BASE_TOPSOIL,"data","processed",
                                     "fraction_lineages_and_counts_remaining.csv"))
emp_prev = pd.read_csv(os.path.join(BASE_EMP,"data","processed",
                                     "lineage_prevalence_thresholds.csv"))

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, prev_df, sel_thresh, abund_col, title in [
    (axes[0], ts_prev,  20, "FractionCountsRemaining",
     "(a) Global Topsoil (threshold = 20 samples)"),
    (axes[1], emp_prev, 50, "FractionMatrixRemaining",
     "(b) EMP Soil Subset (threshold = 50 samples)"),
]:
    x     = prev_df["Threshold"]
    flin  = prev_df["FractionLineagesRemaining"] * 100
    fcnt  = prev_df[abund_col] * 100

    ax2 = ax.twinx()
    ax.plot(x, flin,  "o-", color="#4e79a7", lw=1.5, ms=4,
            label="Fraction of lineages retained")
    ax2.plot(x, fcnt, "s--", color="#e15759", lw=1.5, ms=4,
             label="Fraction of read counts retained")
    ax.axvline(sel_thresh, color="black", lw=1.2, ls=":",
               label=f"Selected threshold ({sel_thresh})")

    ax2.set_ylabel("Read counts retained (%)", color="#e15759", fontsize=9)
    ax2.tick_params(axis="y", labelcolor="#e15759", labelsize=8)
    ax2.set_ylim(90, 101)
    ax2.spines["right"].set_visible(True)

    ax.set_xlabel("Minimum prevalence (no. samples)")
    ax.set_ylabel("Lineages retained (%)")
    ax.set_title(title, fontweight="bold")
    ax.yaxis.grid(True, color="0.9", zorder=0)

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1+lines2, labels1+labels2, fontsize=7.5,
              frameon=True, framealpha=0.85, edgecolor="0.7")

fig.suptitle("Prevalence filter threshold sweep", fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_prevalence_sweep.pdf"),
            bbox_inches="tight")
fig.savefig(os.path.join(OUT_TOPSOIL, "fig_prevalence_sweep.png"),
            bbox_inches="tight", dpi=200)
plt.close(fig)
print("  → fig_prevalence_sweep.pdf/.png")


print("\n✓  All figures written to:")
print(f"   {OUT_TOPSOIL}/")
