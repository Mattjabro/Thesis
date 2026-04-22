#!/bin/bash
#SBATCH --job-name=emp_r17_plot_cluster
#SBATCH --output=EMP_Zenodo_Matthew/logs/emp_r17_plot_cluster_metadata.%j.out
#SBATCH --error=EMP_Zenodo_Matthew/logs/emp_r17_plot_cluster_metadata.%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=1028G
#SBATCH --cpus-per-task=100

set -euo pipefail
mkdir -p EMP_Zenodo_Matthew/logs

INPUT="EMP_Zenodo_Matthew/outputs/cluster_metadata/cluster_metadata_join.tsv"
OUTDIR="EMP_Zenodo_Matthew/outputs/cluster_metadata/plots"
DATASET="EMP"
VENV_PY="/hopper/groups/yoosephlab/datasets/.venv_cluster_plots/bin/python"

"$VENV_PY" EMP_Zenodo_Matthew/scripts/17_plot_cluster_metadata.py \
  --input "$INPUT" \
  --output-dir "$OUTDIR" \
  --dataset-name "$DATASET" \
  --clean
