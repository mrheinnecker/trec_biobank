#!/usr/bin/env bash
## Launcher for the metabarcoding analysis (local, WSL).
## Edit the settings below, then run:
##   bash /mnt/c/repos/trec_biobank/metab/analysis/launcher.sh
##
## Every run creates a new folder <OUTDIR>/<date>_<time>_<ampliseq run>/ with data/, plots/,
## tables/, run_info.txt and run.log. Analysis settings (min reads per sample, top groups
## in plots): analysis settings section of metab/analysis/common.R
set -euo pipefail

## ---- settings ----
REPO=/mnt/c/repos/trec_biobank
BASE=/mnt/c/projects/trec_biobank/metab
OUTDIR=$BASE/analysis

## ampliseq run to analyse: the taxonomy workflow output (asv_table_taxonomy.tsv, asv_taxonomy.tsv)
TAXONOMY_DIR=$BASE/ampliseq_2026-11-19/taxonomy

## sample metadata, one row per run (from metab/make_sample_mapping.R)
export METAB_MAPPING=$BASE/sample_mapping.tsv
## Tara Oceans W4 barcode database with trophic modes (step 05)
export METAB_W4_DB=/mnt/c/projects/trec_biobank/analysis/databases/Database_W4_barcode_occurences.tsv

## R container (rebuild with: bash $REPO/container/build_r_analysis.sh)
export R_SIF=/mnt/c/container/r_analysis.sif
export BINDS=/mnt/c,/tmp

## steps to run; empty = all
##   01 prepare data (taxonomy + mapping -> data/*.rds)   02 sample QC (reads, categories, controls)
##   03 composition per site / size fraction   04 alpha + beta diversity (rarefied / Bray-Curtis)
##   05 trophic modes (PR2 -> Tara W4)   06 trophic modes map per site and size fraction
##   07 map of the 18S reads in the major groups of the TEM screen
##   (every step needs 01 in the same run, 06 needs 05)
##   The TEM figure 13 (analysis/launcher.sh) reads tables/06_trophic_modes_site_fraction.tsv
##   of a run folder: point TEM_METAB_TROPHIC_TABLE there.
STEPS="01 02 03 04 05 06 07"


## ---- run ----
[[ -f "$R_SIF" ]] || { echo "Container not found: $R_SIF"; exit 1; }
[[ -f "$METAB_MAPPING" ]] || { echo "Mapping not found: $METAB_MAPPING"; exit 1; }

# shellcheck disable=SC2086
bash "$REPO/metab/analysis/run_analysis.sh" "$OUTDIR" "$TAXONOMY_DIR" $STEPS
