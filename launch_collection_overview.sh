#!/usr/bin/env bash
#SBATCH --job-name=collection_overview
#SBATCH --output=/scratch/rheinnec/collection_overview/collection_overview_%j.out
#SBATCH --error=/scratch/rheinnec/collection_overview/collection_overview_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail

# Edit this file to define a complete collection-table overview run.
# It can be run directly with bash or submitted with sbatch.

repo_dir="/home/rheinnec/repos/trec_biobank"
outdir="/home/rheinnec/collection_table_overview"
container="/home/rheinnec/container/tidyverse_latest.sif"
prefix="collection_table_overview"

mkdir -p "${outdir}"

bash "${repo_dir}/analysis/overview/run_collection_table_overview.sh" \
  --outdir "${outdir}" \
  --container "${container}" \
  --prefix "${prefix}"



## on the cluster environment

repo_dir="/g/schwab/marco/repos/trec_biobank"
outdir="/g/schwab/tem_screen/overview"
container="/g/schwab/marco/container_legacy/tidyverse_latest.sif"
prefix="collection_table_overview"

mkdir -p "${outdir}"

bash "${repo_dir}/analysis/overview/run_collection_table_overview.sh" \
  --outdir "${outdir}" \
  --container "${container}" \
  --prefix "${prefix}"

