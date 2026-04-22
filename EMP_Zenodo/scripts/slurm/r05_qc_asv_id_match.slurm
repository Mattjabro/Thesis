#!/bin/bash
#SBATCH --job-name=emp_r05_qc_asv_id_match
#SBATCH --output=logs/emp_r05_qc_asv_id_match.out
#SBATCH --error=logs/emp_r05_qc_asv_id_match.err
#SBATCH --time=11:00:00
#SBATCH --mem=1028G
#SBATCH --cpus-per-task=128

module load R
Rscript /hopper/groups/yoosephlab/datasets/EMP_Zenodo_Matthew/scripts/05_qc_asv_id_match.R
