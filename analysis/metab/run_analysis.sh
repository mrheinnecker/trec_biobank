#!/usr/bin/env bash
## Runs the metabarcoding analysis scripts in order inside the R container.
##
## usage: run_analysis.sh <outdir> <taxonomy_dir> [step ...]
##   taxonomy_dir: output folder of the taxonomy workflow (asv_table_taxonomy.tsv, asv_taxonomy.tsv)
##   steps: 01 02 03 04 05 06 07 (default: all); every step needs the data/ of 01 from the same run,
##          06 also the table of 05
##   Every run writes into its own folder <outdir>/<YYYYmmdd_HHMMSS>_<ampliseq run>/
##   (data/, plots/, tables/, run_info.txt). The ampliseq run name is the folder above taxonomy/
##   (e.g. ampliseq_2026-11-19).
## env:   METAB_MAPPING  sample_mapping.tsv (default in common.R)
##        METAB_W4_DB    Tara W4 database for 05 (default in common.R)
##        R_SIF          container with R + tidyverse, vegan, sf, rnaturalearth, patchwork, cowplot
##        BINDS          apptainer bind paths (default: /mnt/c,/tmp)
set -euo pipefail

if [[ $# -lt 2 ]]; then
  sed -n '2,14p' "$0"; exit 1
fi

base_outdir=$1; taxonomy_dir=$2; shift 2
steps=("$@"); [[ ${#steps[@]} -eq 0 ]] && steps=(01 02 03 04 05 06 07)

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
R_SIF=${R_SIF:-/mnt/c/container/r_analysis.sif}
BINDS=${BINDS:-/mnt/c,/tmp}

[[ -f "$taxonomy_dir/asv_table_taxonomy.tsv" ]] || { echo "Not a taxonomy folder: $taxonomy_dir"; exit 1; }
taxonomy_dir=$(realpath "$taxonomy_dir")

## run name: folder above taxonomy/, skipping a "results" level (ampliseq_X/results/taxonomy)
run_name=$(basename "$(dirname "$taxonomy_dir")")
[[ "$run_name" == "results" ]] && run_name=$(basename "$(dirname "$(dirname "$taxonomy_dir")")")

## one folder per run
outdir="$base_outdir/$(date +%Y%m%d_%H%M%S)_${run_name}"
mkdir -p "$outdir"; outdir=$(realpath "$outdir")
echo "Output: $outdir"

export METAB_TAXONOMY_DIR=$taxonomy_dir
export METAB_OUTDIR=$outdir

{
  echo "date          $(date -Iseconds)"
  echo "taxonomy_dir  $taxonomy_dir"
  echo "mapping       ${METAB_MAPPING:-default in common.R}"
  echo "w4_db         ${METAB_W4_DB:-default in common.R}"
  echo "steps         ${steps[*]}"
  echo "container     $R_SIF"
  echo "git commit    $(git -C "$script_dir" rev-parse --short HEAD 2>/dev/null || echo unknown)$(git -C "$script_dir" diff --quiet -- . 2>/dev/null || echo ' (+ uncommitted changes in analysis/metab)')"
} > "$outdir/run_info.txt"

for step in "${steps[@]}"; do
  script=""
  for candidate in "$script_dir"/"${step}"_*.R; do
    [[ -f "$candidate" ]] && { script=$candidate; break; }
  done
  [[ -n "$script" ]] || { echo "No script for step $step"; exit 1; }
  echo "== $(basename "$script")"
  apptainer exec -B "$BINDS" "$R_SIF" Rscript "$script" 2>&1 | tee -a "$outdir/run.log"
done

echo "Done: $outdir"
