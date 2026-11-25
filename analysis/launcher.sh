#!/usr/bin/env bash
## Launcher for all downstream analyses (local, WSL). Edit the settings below, then run:
##   bash /mnt/c/repos/trec_biobank/analysis/launcher.sh
##
## Modules (run in this order, each in the R container; every run gets its own output folder):
##   overview  biobank overview figure (samples, events, imaging per modality)  analysis/overview/
##             -> <DATA_ROOT>/analysis/overview/{data,runs}/<date>_<time>/
##   metab     metabarcoding analysis (ampliseq taxonomy -> composition, diversity, trophic modes)
##             analysis/metab/  -> <DATA_ROOT>/analysis/metab/<date>_<time>_<ampliseq run>/
##   tem       TEM annotation analysis (QC, organelles, clustering, maps)    analysis/tem/
##             -> <DATA_ROOT>/analysis/tem/<date>_<time>_round<ROUNDS>/
## tem step 13 combines TEM with the metab trophic modes: it uses the metab run of the same launch,
## otherwise the newest metab run that has the table (or TEM_METAB_TROPHIC_TABLE below).
##
## Data locations are defined once in analysis/paths.R (defaults); the values below override them.
## Upstream processing (Nextflow imaging workflows, ampliseq, taxonomy) is not run from here:
## see launch_wf*.sh in the repo root and metab/.
set -euo pipefail

## ======================================================================================
## ---- what to run ----
## modules, in this order; any subset, e.g. "overview" or "metab tem"
MODULES="overview metab tem"

## steps per module; empty = all steps of that module
##   overview: data (new snapshot from the sheets), figure, sheet (mapping sheet); default "figure"
##   metab:    01 prepare data   02 sample QC   03 composition   04 diversity   05 trophic modes
##             06 trophic modes map   07 map of the TEM major groups (all need 01 in the same run,
##             06 needs 05)
##   tem:      01 annotator QC   02 site overview   03 organelles per site   04/05 organelle
##             associations (tests/plots)   06 clustering   07/08 organelles x metadata
##             09 collection tables of the negative pairs -> Google sheet   10 montages of 06 / 09
##             11 chloroplast / food vacuole per site   12 map of 11   13 TEM + metab trophic map
##             14 TEM feature maps (05 and 09 need 04, 08 needs 07, 10 needs 06 / 09, 12-13 need 11)
OVERVIEW_STEPS="figure"
METAB_STEPS="01 02 03 04 05 06 07"
TEM_STEPS="01 02 03 04 05 06 07 08 09 10 11 12 13 14"

## ======================================================================================
## ---- shared ----
REPO=/mnt/c/repos/trec_biobank
export TREC_DATA_ROOT=/mnt/c/projects/trec_biobank

## Google service-account key (read access to all sheets, edit access to the output sheets);
## the service account is nextflow-gsheets@trec-tem-screen.iam.gserviceaccount.com
export TREC_GOOGLE_KEY=$REPO/trec-tem-screen-e98a2e03f58b.json

## R container (rebuild with: bash $REPO/container/build_r_analysis.sh). Packages missing in an
## older build are installed once into R_LIB.
export R_SIF=/mnt/c/container/r_analysis.sif
R_LIB=/mnt/c/container/r_analysis_libs
export BINDS=/mnt/c,/tmp

## time zone for the run folder names (the container runs in UTC otherwise)
export TZ=Europe/Berlin

## ======================================================================================
## ---- inputs (outputs of the upstream pipelines) ----
## metab: taxonomy workflow output (asv_table_taxonomy.tsv, asv_taxonomy.tsv) and sample mapping
## (metab/make_sample_mapping.R)
METAB_TAXONOMY_DIR=$TREC_DATA_ROOT/metab/ampliseq_2026-11-19/taxonomy
export METAB_MAPPING=$TREC_DATA_ROOT/metab/sample_mapping.tsv
## Tara Oceans W4 barcode database with trophic modes (metab 05)
export METAB_W4_DB=$TREC_DATA_ROOT/analysis/databases/Database_W4_barcode_occurences.tsv
## TEM: overview PNGs of wfTEM (one per cell, <TEM_PNG_DIR>/<site>/<source_name>.png) and the cache of
## their reduced copies (tem 10)
export TEM_PNG_DIR=$TREC_DATA_ROOT/pngs
export TEM_THUMB_DIR=$TREC_DATA_ROOT/png_thumbs
## TEM: rerun on an earlier download of the annotation sheet instead of fetching it again:
## <tem run>/tables/annotations_main.tsv; empty = download fresh
TEM_INPUT=""

## ======================================================================================
## ---- overview settings ----
## data step: TRUE = download the sheets / collection tables again, FALSE = reuse the previous
## snapshot (only recompute the event maps)
export OVERVIEW_REFRESH=TRUE
## figure / sheet step: data snapshot (folder name in <DATA_ROOT>/analysis/overview/data; empty = newest)
export OVERVIEW_DATA_SNAPSHOT=""
## note appended to the figure run folder name (e.g. "new_colours")
export OVERVIEW_RUN_NOTE=""
## Google sheet for the mapping tables (sheet step); empty = only tsv files
export OVERVIEW_MAPPING_SHEET="https://docs.google.com/spreadsheets/d/1DabzZ1275Z4Agh_ertKYsTSQXndRIsf-qnUz4VQvj2E/edit"
## colours and legend labels: COLOURS AND LEGEND LABELS section of overview/biobank_overview_figure.R

## ---- metab settings ----
## min reads per sample, top groups in plots: analysis settings section of metab/common.R

## ---- tem settings ----
## annotation rounds: 1 = first annotation only, 1,2 = also the validations
export ROUNDS=1
## spreadsheet for the collection tables of step 09; empty = only local tsv
export TEM_COLLECTION_SHEET_URL="https://docs.google.com/spreadsheets/d/18b2AVslU7sOQOTRRIwxGHKPtrh00ezT-byOrQpiVo8Y/edit?gid=0#gid=0"
## step 06: number of clusters; empty = automatic (k_method in 06). Features: use_features in tem/assoc_common.R
export TEM_CLUSTER_K=5
## step 10: reduction of the PNGs (4 -> 28 nm/px, ~22k px per figure; 2 -> ~40k px)
export TEM_SCALE_FACTOR=4
## step 10: size of the plot part of the cluster tree figure and upper limit of the cell image size
export TEM_TREE_PLOT_SCALE=1.5
export TEM_TREE_IMAGE_SCALE=0.5
## step 13: metab table 06_trophic_modes_site_fraction.tsv; empty = from the metab run of this
## launch, otherwise the newest metab run that has it
TEM_METAB_TROPHIC_TABLE=""


## ======================================================================================
## ---- run ----
[[ -f "$R_SIF" ]] || { echo "Container not found: $R_SIF"; exit 1; }
mkdir -p "$R_LIB"
export R_LIBS_USER=$R_LIB
## passed into the container (apptainer does not forward TZ)
export APPTAINERENV_TZ=$TZ APPTAINERENV_R_LIBS_USER=$R_LIB
export TREC_REPO=$REPO
export TEM_GOOGLE_KEY=$TREC_GOOGLE_KEY
OVERVIEW_STEPS=${OVERVIEW_STEPS:-figure}
analysis_dir=$REPO/analysis
runs_root=$TREC_DATA_ROOT/analysis

newest() { ls -d "$1"/*/ 2>/dev/null | grep -v "/rstudio/$" | sort | tail -1 | sed 's#/$##'; }

## missing packages -> R_LIB (once)
apptainer exec -B "$BINDS" "$R_SIF" Rscript -e '
  pkgs <- c("circlize", "ggalluvial", "patchwork")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing) > 0) {
    message("Installing into ", Sys.getenv("R_LIBS_USER"), ": ", paste(missing, collapse=", "))
    install.packages(missing, lib=Sys.getenv("R_LIBS_USER"),
                     repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", quiet=TRUE)
  }'

summary=()
for module in $MODULES; do
  echo
  echo "######## $module"
  case "$module" in
    overview)
      # shellcheck disable=SC2086
      bash "$analysis_dir/overview/run_overview.sh" $OVERVIEW_STEPS
      if [[ " $OVERVIEW_STEPS " == *" data "* ]]; then summary+=("overview data:   $(newest "$runs_root/overview/data")"); fi
      if [[ " $OVERVIEW_STEPS " == *" figure "* ]]; then summary+=("overview figure: $(newest "$runs_root/overview/runs")"); fi
      ;;
    metab)
      [[ -f "$METAB_MAPPING" ]] || { echo "Mapping not found: $METAB_MAPPING"; exit 1; }
      # shellcheck disable=SC2086
      bash "$analysis_dir/metab/run_analysis.sh" "$runs_root/metab" "$METAB_TAXONOMY_DIR" $METAB_STEPS
      metab_run=$(newest "$runs_root/metab")
      summary+=("metab:           $metab_run")
      ## tem 13 uses this run when it has the trophic table
      if [[ -z "$TEM_METAB_TROPHIC_TABLE" && -f "$metab_run/tables/06_trophic_modes_site_fraction.tsv" ]]; then
        TEM_METAB_TROPHIC_TABLE=$metab_run/tables/06_trophic_modes_site_fraction.tsv
      fi
      ;;
    tem)
      if [[ -n "$TEM_INPUT" ]]; then
        [[ -f "$TEM_INPUT" ]] || { echo "Input not found: $TEM_INPUT"; exit 1; }
        export INPUT=$TEM_INPUT
      else
        [[ -f "$TREC_GOOGLE_KEY" ]] || { echo "Google key not found: $TREC_GOOGLE_KEY"; exit 1; }
        unset INPUT
      fi
      if [[ -n "$TEM_METAB_TROPHIC_TABLE" ]]; then
        export TEM_METAB_TROPHIC_TABLE
        echo "metab table for step 13: $TEM_METAB_TROPHIC_TABLE"
      fi
      # shellcheck disable=SC2086
      bash "$analysis_dir/tem/run_analysis.sh" "$runs_root/tem" $TEM_STEPS
      summary+=("tem:             $(newest "$runs_root/tem")")
      ;;
    *)
      echo "Unknown module: $module (overview, metab, tem)"; exit 1 ;;
  esac
done

echo
echo "######## done"
printf '%s\n' "${summary[@]}"
