#!/usr/bin/env bash
## Runs the biobank overview scripts inside the R container (called by analysis/launcher.sh).
##
## usage: run_overview.sh [step ...]
##   data    biobank_overview.R          -> new snapshot <data_root>/analysis/overview/data/<date>_<time>/
##   figure  biobank_overview_figure.R   -> new folder   <data_root>/analysis/overview/runs/<date>_<time>[_<note>]/
##   sheet   export_mapping_sheet.R      -> mapping tables to OVERVIEW_MAPPING_SHEET and <snapshot>/mapping_sheet/
##   (default: figure)
## env:   OVERVIEW_REFRESH, OVERVIEW_DATA_SNAPSHOT, OVERVIEW_RUN_NOTE, OVERVIEW_MAPPING_SHEET,
##        TREC_GOOGLE_KEY, TREC_DATA_ROOT, R_SIF, BINDS (see analysis/launcher.sh)
set -euo pipefail

steps=("$@"); [[ ${#steps[@]} -eq 0 ]] && steps=(figure)

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
R_SIF=${R_SIF:-/mnt/c/container/r_analysis.sif}
BINDS=${BINDS:-/mnt/c,/tmp}

run_r() { (cd "$script_dir" && apptainer exec -B "$BINDS" "$R_SIF" Rscript "$@"); }

for step in "${steps[@]}"; do
  case "$step" in
    data)
      if [[ "${OVERVIEW_REFRESH:-TRUE}" == "TRUE" ]]; then
        [[ -f "${TREC_GOOGLE_KEY:-}" ]] || { echo "Google key not found: ${TREC_GOOGLE_KEY:-}"; exit 1; }
      fi
      echo "== biobank_overview.R"
      run_r biobank_overview.R
      ;;
    figure)
      echo "== biobank_overview_figure.R"
      ## circlize prints a note for every label outside its track; hide those
      run_r biobank_overview_figure.R 2>&1 | { grep -v -E "^Note: |^'[0-9]+'\.$|track$|^track '[0-9]+'\.$" || true; }
      ;;
    sheet)
      [[ -z "${OVERVIEW_MAPPING_SHEET:-}" || -f "${TREC_GOOGLE_KEY:-}" ]] || { echo "Google key not found: ${TREC_GOOGLE_KEY:-}"; exit 1; }
      echo "== export_mapping_sheet.R"
      run_r export_mapping_sheet.R
      ;;
    *)
      echo "Unknown overview step: $step (data, figure, sheet)"; exit 1 ;;
  esac
done
