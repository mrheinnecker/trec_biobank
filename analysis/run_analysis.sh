#!/usr/bin/env bash
## Runs the TEM annotation analysis scripts in order inside an R container.
##
## usage: run_analysis.sh <outdir> [step ...]
##   steps: 01 02 03 04 05 06 07 08 09 10 11 12 13 (default: all); 05 and 09 need the output of 04,
##          08 that of 07, 10 that of 06 / 09, 12 / 13 that of 11 (13 also the metab table)
##   Every run writes into its own folder <outdir>/<YYYYmmdd_HHMMSS>_round<ROUNDS>/
##   (plots/ and tables/). First the main annotation sheet is downloaded into
##   tables/annotations_main.tsv of that folder (00_download_annotations.R); all steps use it.
## env:   INPUT   use this tsv snapshot instead of downloading
##        ROUNDS  annotation rounds to analyse (default 1 = first annotation only; e.g. 1,2)
##        TEM_GOOGLE_KEY  Google service-account key (default: key in the repo on /g)
##        R_SIF   container with R + tidyverse, googlesheets4, cowplot, patchwork, lme4
##        BINDS   extra apptainer bind paths (default: /g,/scratch,/home)
##        TEM_PNG_DIR, TEM_THUMB_DIR, TEM_SCALE_FACTOR   images for the montages of step 10
set -euo pipefail

if [[ $# -lt 1 ]]; then
  sed -n '2,15p' "$0"; exit 1
fi

base_outdir=$1; shift 1
steps=("$@"); [[ ${#steps[@]} -eq 0 ]] && steps=(01 02 03 04 05 06 07 08 09 10 11 12 13)

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
R_SIF=${R_SIF:-/g/schwab/marco/container_legacy/tidyverse_latest.sif}
BINDS=${BINDS:-/g,/scratch,/home}
ROUNDS=${ROUNDS:-1}

## one folder per run
outdir="$base_outdir/$(date +%Y%m%d_%H%M%S)_round${ROUNDS//,/-}"
mkdir -p "$outdir"; outdir=$(realpath "$outdir")
echo "Output: $outdir"

if [[ -n "${INPUT:-}" ]]; then
  input=$(realpath "$INPUT")
else
  echo "== 00_download_annotations.R"
  apptainer exec -B "$BINDS" "$R_SIF" \
    Rscript "$script_dir/00_download_annotations.R" --outdir "$outdir"
  input="$outdir/tables/annotations_main.tsv"
fi

for step in "${steps[@]}"; do
  script=""
  for candidate in "$script_dir"/"${step}"_*.R "$script_dir"/"${step}"_*.py; do
    [[ -f "$candidate" ]] && { script=$candidate; break; }
  done
  [[ -n "$script" ]] || { echo "No script for step $step"; exit 1; }
  echo "== $(basename "$script")"
  case "$script" in
    *.py) runner=python3 ;;
    *)    runner=Rscript ;;
  esac
  apptainer exec -B "$BINDS" "$R_SIF" \
    "$runner" "$script" --input "$input" --outdir "$outdir" --rounds "$ROUNDS"
done
