#!/bin/bash
#SBATCH --job-name=soil_r22_plot_cluster
#SBATCH --output=global_topsoil_matthew/logs/soil_r22_plot_cluster_metadata.%j.out
#SBATCH --error=global_topsoil_matthew/logs/soil_r22_plot_cluster_metadata.%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail
mkdir -p global_topsoil_matthew/logs

INPUT="global_topsoil_matthew/outputs/cluster_metadata/cluster_metadata_join.tsv"
OUTDIR="global_topsoil_matthew/outputs/cluster_metadata/plots"
DATASET="topsoil"
VENV_PY="/hopper/groups/yoosephlab/datasets/.venv_cluster_plots/bin/python"

"$VENV_PY" global_topsoil_matthew/scripts/22_plot_cluster_metadata.py \
  --input "$INPUT" \
  --output-dir "$OUTDIR" \
  --dataset-name "$DATASET" \
  --clean
