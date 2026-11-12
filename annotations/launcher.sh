## TEM annotation: distribute images to the annotators and collect the results.
## Run inside the tidyverse container on the cluster:

apptainer shell -B /home -B /scratch -B /g /g/schwab/marco/container_legacy/tidyverse_latest.sif

cd /g/schwab/marco/repos/tem_classification

KEY=/g/schwab/marco/repos/tem_classification/trec-tem-screen-e98a2e03f58b.json
ANNOTATION_SHEET="https://docs.google.com/spreadsheets/d/1NDyVERdrl7nXJrQRWBbwHjyHCMNEZhj1RQnBKUObwuU/edit?gid=0#gid=0"
TEM_COLLECTION="https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=2014112369#gid=2014112369"


## ---- one-off: convert the old main sheet to the round format ----
## 1. dry run: backs up the live sheets, rebuilds the rounds from all backups,
##    writes main_rounds.tsv, status.tsv and issues.tsv to --outdir. Check issues.tsv.
Rscript annotations/migrate_annotation_rounds.R \
  --google_key "$KEY" \
  --collection_table_url "$ANNOTATION_SHEET" \
  --source_collection_table_url "$TEM_COLLECTION" \
  --source_collection_sheet tem_collection_table \
  --log_dir /g/schwab/tem_screen/annotations/log \
  --outdir /g/schwab/tem_screen/annotations/migration_rounds \
  --write FALSE

## 2. same with --write TRUE: overwrites the main and status sheets (personal tabs untouched)


## ---- every round of annotation work ----
Rscript annotations/make_split_collection_tables.R \
  --sheet_mode google \
  --google_key "$KEY" \
  --collection_table_url "$ANNOTATION_SHEET" \
  --annotations_sheet main \
  --status_sheet status \
  --source_collection_table_url "$TEM_COLLECTION" \
  --source_collection_sheet tem_collection_table \
  --annotation_log_dir /g/schwab/tem_screen/annotations/log \
  --assignees "marco,chandni,yannick,karel,ella" \
  --former_assignees "viktoria"

## All open images are split equally over --assignees (first annotations first,
## validations behind them). Optional: --max_rows_per_person "50,20,20,20,20" caps the rows.
## Someone leaves: move the name from --assignees to --former_assignees (their tab is
## collected, backed up and emptied). Someone new: add the name to --assignees.
