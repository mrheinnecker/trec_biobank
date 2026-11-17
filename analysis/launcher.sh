#!/usr/bin/env bash
## Launcher for the TEM annotation analysis (local, WSL).
## Edit the settings below, then run:
##   bash /mnt/c/repos/tem_classification/analysis/launcher.sh
##
## Every run creates a new folder <OUTDIR>/<date>_<time>_round<ROUNDS>/ with plots/ and tables/.
## Which organelles are tested / clustered: use_features in analysis/assoc_common.R
## Clustering settings: top of analysis/06_organelle_clustering.R
set -euo pipefail

## ---- settings ----
REPO=/mnt/c/repos/tem_classification
OUTDIR=/mnt/c/projects/tem_classification/analysis

## R container (rebuild with: bash $REPO/container/build_r_analysis.sh)
export R_SIF=/mnt/c/container/r_analysis.sif
export BINDS=/mnt/c,/tmp

## Google service-account key (read access to the annotation spreadsheet)
export TEM_GOOGLE_KEY=$REPO/trec-tem-screen-e98a2e03f58b.json

## annotation rounds: 1 = first annotation only, 1,2 = also the validations
export ROUNDS=1

## steps to run; empty = all
##   01 annotator QC + uncertain entries   02 site overview   03 organelles per site
##   04/05 organelles x organelles (tests/plots)   06 organelle clustering
##   07/08 organelles x metadata (tests/plots)
##   09 collection tables of the negative organelle pairs -> Google sheet (TEM_COLLECTION_SHEET_URL)
##   10 overview figures (JPEG montages) of the collection tables of 06 (clusters) and 09 (pairs)
##   (05 and 09 need 04, 08 needs 07, 10 needs 09 in the same run)
STEPS="01 02 03 04 05 06 07 08 09 10 11 12"

## spreadsheet for the collection tables of step 09 (the service account
## nextflow-gsheets@trec-tem-screen.iam.gserviceaccount.com needs edit access); empty = only local tsv
export TEM_COLLECTION_SHEET_URL="https://docs.google.com/spreadsheets/d/18b2AVslU7sOQOTRRIwxGHKPtrh00ezT-byOrQpiVo8Y/edit?gid=0#gid=0"

## step 06: number of clusters; empty = the k with the best mean silhouette width
export TEM_CLUSTER_K=6

## step 10: one PNG per cell in <TEM_PNG_DIR>/<site>/<source_name>.png; all PNGs have the same pixel
## size, so the cells are packed side by side true to scale. All are reduced by TEM_SCALE_FACTOR
## (full resolution exceeds the JPEG limit; 4 -> 28 nm/px, ~22k px per figure; 2 -> ~40k px).
## The reduced images are cached in TEM_THUMB_DIR, so the PNGs are read only once.
export TEM_PNG_DIR=/mnt/c/projects/tem_classification/pngs
export TEM_THUMB_DIR=/mnt/c/projects/tem_classification/png_thumbs
export TEM_SCALE_FACTOR=4
## combined tree figure of the clusters: size of the plot part (tree, metadata, organelles, labels;
## also sets the width per cell and with it the whole figure width) and the upper limit of the cell
## image size (images are reduced further until the widest cell of each cluster fits its column)
export TEM_TREE_PLOT_SCALE=1.5
export TEM_TREE_IMAGE_SCALE=0.5

## rerun on an earlier download instead of fetching the main sheet again:
## path to <run folder>/tables/annotations_main.tsv; empty = download fresh
INPUT=""


## ---- run ----
[[ -f "$R_SIF" ]] || { echo "Container not found: $R_SIF"; exit 1; }
if [[ -n "$INPUT" ]]; then
  [[ -f "$INPUT" ]] || { echo "Input not found: $INPUT"; exit 1; }
  export INPUT
else
  [[ -f "$TEM_GOOGLE_KEY" ]] || { echo "Google key not found: $TEM_GOOGLE_KEY"; exit 1; }
  unset INPUT
fi

# shellcheck disable=SC2086
bash "$REPO/analysis/run_analysis.sh" "$OUTDIR" $STEPS
