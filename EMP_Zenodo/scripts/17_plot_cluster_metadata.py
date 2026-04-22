#!/usr/bin/env python3
import argparse
import os
import re
import glob

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


DEFAULT_CATEGORICAL_MAX_UNIQUE = 20


def sanitize_filename(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return value.strip("_") or "variable"


def is_categorical(series: pd.Series, max_unique: int) -> bool:
    if pd.api.types.is_numeric_dtype(series):
        return series.nunique(dropna=True) <= max_unique
    return True


def _cluster_order(series: pd.Series) -> list:
    try:
        numeric = pd.to_numeric(series, errors="coerce")
        if numeric.notna().any():
            return [str(v) for v in sorted(numeric.dropna().unique())]
    except Exception:
        pass
    return sorted(series.astype(str).unique())


def plot_categorical(df: pd.DataFrame, variable: str, output_path: str, col_wrap: int) -> None:
    plot_df = df[["cluster", variable]].dropna()
    plot_df["cluster"] = plot_df["cluster"].astype(str)

    g = sns.catplot(
        data=plot_df,
        y=variable,
        kind="count",
        col="cluster",
        col_wrap=col_wrap,
        height=3.2,
        aspect=1.1,
        sharex=False,
        sharey=False,
    )
    g.set_titles("Cluster {col_name}")
    for ax in g.axes.flat:
        ax.tick_params(axis="y", labelsize=8)
        ax.tick_params(axis="x", labelsize=8)
        ax.set_xlabel("Count")
        ax.set_ylabel(variable)
    g.figure.tight_layout()
    g.figure.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(g.figure)


def plot_numeric(df: pd.DataFrame, variable: str, output_path: str) -> None:
    plot_df = df[["cluster", variable]].dropna()
    plot_df["cluster"] = plot_df["cluster"].astype(str)

    order = _cluster_order(plot_df["cluster"])
    fig, ax = plt.subplots(figsize=(12, 5))
    sns.boxplot(
        data=plot_df,
        x="cluster",
        y=variable,
        ax=ax,
        order=order,
    )
    ax.set_xlabel("Cluster")
    ax.set_ylabel(variable)
    ax.tick_params(axis="x", rotation=45, labelsize=9)
    fig.tight_layout()
    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot metadata distributions by cluster.")
    parser.add_argument("--input", required=True, help="Path to cluster_metadata_join.tsv")
    parser.add_argument("--output-dir", required=True, help="Directory for PDF outputs")
    parser.add_argument("--dataset-name", required=True, help="Label used in output filenames")
    parser.add_argument(
        "--categorical-max-unique",
        type=int,
        default=DEFAULT_CATEGORICAL_MAX_UNIQUE,
        help="Treat numeric columns with <= this many unique values as categorical",
    )
    parser.add_argument("--col-wrap", type=int, default=5, help="Facet column wrap for categorical plots")
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove existing dataset PDF outputs before regenerating",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.clean:
        pattern = os.path.join(
            args.output_dir,
            f"{sanitize_filename(args.dataset_name)}_*_by_cluster.pdf",
        )
        for path in glob.glob(pattern):
            try:
                os.remove(path)
            except OSError:
                pass

    df = pd.read_csv(args.input, sep="\t")
    if "cluster" not in df.columns:
        raise ValueError("Expected a 'cluster' column in the input file.")

    skip_columns = {
        "X.SampleID",
        "cluster",
        "BarcodeSequence",
        "LinkerPrimerSequence",
        "Description",
        "host_subject_id",
        "study_id",
        "title",
        "principal_investigator",
        "doi",
        "ebi_accession",
        "target_gene",
        "target_subfragment",
        "pcr_primers",
        "illumina_technology",
        "extraction_center",
        "run_center",
        "run_date",
        "read_length_bp",
        "sequences_split_libraries",
        "observations_closed_ref_greengenes",
        "observations_closed_ref_silva",
        "observations_open_ref_greengenes",
        "observations_deblur_90bp",
        "observations_deblur_100bp",
        "observations_deblur_150bp",
        "emp_release1",
        "qc_filtered",
        "subset_10k",
        "subset_5k",
        "subset_2k",
    }
    variables = [col for col in df.columns if col not in skip_columns]

    min_non_na = 10  # Minimum non-NA values required to plot
    min_unique = 2   # Minimum unique values required (skip constant columns)

    for variable in variables:
        series = df[variable]
        non_na_series = series.dropna()
        
        # Skip if too few non-NA values
        if len(non_na_series) < min_non_na:
            print(f"Skipping '{variable}': only {len(non_na_series)} non-NA values (need >= {min_non_na})")
            continue
        
        # Skip if only one unique value (uninformative)
        n_unique = non_na_series.nunique()
        if n_unique < min_unique:
            print(f"Skipping '{variable}': only {n_unique} unique value(s) (need >= {min_unique})")
            continue

        file_label = sanitize_filename(variable)
        output_path = os.path.join(
            args.output_dir,
            f"{sanitize_filename(args.dataset_name)}_{file_label}_by_cluster.pdf",
        )

        if is_categorical(series, args.categorical_max_unique):
            plot_categorical(df, variable, output_path, args.col_wrap)
        else:
            plot_numeric(df, variable, output_path)


if __name__ == "__main__":
    sns.set_theme(style="whitegrid")
    main()
