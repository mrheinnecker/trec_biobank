#!/usr/bin/env bash
## Runs the TEM annotation analysis scripts in order inside an R container.
##
## usage: run_analysis.sh <annotations.tsv> <outdir> [step ...]
##   steps: 01 02 03 04 05 06 (default: all); 05 needs the output of 04
## env:   R_SIF   container with R + tidyverse, cowplot, patchwork, lme4
##        BINDS   extra apptainer bind paths (default: /g,/scratch,/home)
set -euo pipefail

if [[ $# -lt 2 ]]; then
  sed -n '2,8p' "$0"; exit 1
fi

input=$(realpath "$1"); outdir=$2; shift 2
steps=("$@"); [[ ${#steps[@]} -eq 0 ]] && steps=(01 02 03 04 05 06)
mkdir -p "$outdir"; outdir=$(realpath "$outdir")

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
R_SIF=${R_SIF:-/g/schwab/marco/container_legacy/tidyverse_latest.sif}
BINDS=${BINDS:-/g,/scratch,/home}

for step in "${steps[@]}"; do
  script=$(ls "$script_dir"/"${step}"_*.R)
  echo "== $(basename "$script")"
  apptainer exec -B "$BINDS" "$R_SIF" \
    Rscript "$script" --input "$input" --outdir "$outdir"
done
