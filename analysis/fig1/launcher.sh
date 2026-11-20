#!/usr/bin/env bash
## Launcher for the biobank overview figure (fig1; local, WSL).
## Edit the settings below, then run:
##   bash /mnt/c/repos/trec_biobank/analysis/fig1/launcher.sh
##
## Steps:
##   data    biobank_overview.R          -> new snapshot <FIG1_ROOT>/data/<date>_<time>/
##   figure  biobank_overview_figure.R   -> new folder   <FIG1_ROOT>/runs/<date>_<time>[_<RUN_NOTE>]/
##           (plots, copy of the figure script, run_info.txt)
## Colours and legend labels: COLOURS AND LEGEND LABELS section at the top of biobank_overview_figure.R
set -euo pipefail

## ---- settings ----
REPO=/mnt/c/repos/trec_biobank

## steps to run: "figure" (default; uses an existing data snapshot), "data figure" (new snapshot first)
STEPS="figure"

## data step: TRUE = download the sheets / collection tables again, FALSE = reuse the previous snapshot
## (only recompute the event maps). Downloading needs a Google service-account key with read access to
## the sheets (the container cannot do the interactive login).
export FIG1_REFRESH=TRUE
export FIG1_GOOGLE_KEY=$REPO/trec-tem-screen-e98a2e03f58b.json

## figure step: data snapshot (folder name in <FIG1_ROOT>/data; empty = the newest) and an optional
## note appended to the run folder name (e.g. "new_colours")
export FIG1_DATA_SNAPSHOT=""
export FIG1_RUN_NOTE=""

## R container (rebuild with: bash $REPO/container/build_r_analysis.sh). Packages missing in it
## (circlize, ggalluvial, ggnewscale in older builds) are installed once into R_LIB.
R_SIF=/mnt/c/container/r_analysis.sif
R_LIB=/mnt/c/container/r_analysis_libs
BINDS=/mnt/c,/tmp

## time zone for the run folder names (the container runs in UTC otherwise)
export TZ=Europe/Berlin


## ---- run ----
[[ -f "$R_SIF" ]] || { echo "Container not found: $R_SIF"; exit 1; }
mkdir -p "$R_LIB"
export R_LIBS_USER=$R_LIB
fig1_dir=$REPO/analysis/fig1

run_r() {
  apptainer exec -B "$BINDS" --env "R_LIBS_USER=$R_LIB,TZ=$TZ" "$R_SIF" Rscript "$@"
}

## missing packages -> R_LIB (once)
run_r -e '
  pkgs <- c("circlize", "ggalluvial", "ggnewscale")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing) > 0) {
    message("Installing into ", Sys.getenv("R_LIBS_USER"), ": ", paste(missing, collapse=", "))
    install.packages(missing, lib=Sys.getenv("R_LIBS_USER"),
                     repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", quiet=TRUE)
  }'

for step in $STEPS; do
  case "$step" in
    data)
      if [[ "$FIG1_REFRESH" == "TRUE" ]]; then
        [[ -f "$FIG1_GOOGLE_KEY" ]] || { echo "Google key not found: $FIG1_GOOGLE_KEY"; exit 1; }
      fi
      echo "== biobank_overview.R"
      (cd "$fig1_dir" && run_r biobank_overview.R)
      ;;
    figure)
      echo "== biobank_overview_figure.R"
      ## circlize prints a note for every label outside its track; hide those
      (cd "$fig1_dir" && run_r biobank_overview_figure.R 2>&1 | grep -v -E "^Note: |^'[0-9]+'\.$|track$|^track '[0-9]+'\.$")
      ;;
    *)
      echo "Unknown step: $step (data, figure)"; exit 1 ;;
  esac
done

latest=$(ls -d /mnt/c/projects/trec_biobank/fig1/runs/*/ 2>/dev/null | tail -1)
[[ -n "$latest" ]] && echo "Latest figure run: $latest"
