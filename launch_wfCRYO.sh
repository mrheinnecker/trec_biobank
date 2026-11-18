#!/usr/bin/env bash
#SBATCH --job-name=wfCRYO
#SBATCH --output=/scratch/rheinnec/cryo_screen/wfCRYO_%j.out
#SBATCH --error=/scratch/rheinnec/cryo_screen/wfCRYO_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
set -euo pipefail

# Edit this file to define a complete CRYO screening run.
# The profile only selects where Nextflow jobs run:
#   interactive = current cluster node
#   cluster     = submit jobs through Slurm

timestamp="$(date +%Y-%m-%d_%H-%M)"

repo_dir="/g/schwab/marco/repos/trec_biobank"
cd /scratch/rheinnec/cryo_screen

bash "${repo_dir}/cryo/cryo_main.sh" interactive \
  --sheet_mode google \
  --sheet_url "https://docs.google.com/spreadsheets/d/1ePRpa56mmMvCeRTLXmwOywOLy5_I3AFrxJepSUYGR1s/edit?gid=1442254503#gid=1442254503" \
  --sheet_name "cryo_input_table" \
  --google_key "${repo_dir}/trec-tem-screen-e98a2e03f58b.json" \
  --collection_table_url "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=199938698#gid=199938698" \
  --collection_table_sheet "cryo_collection_table" \
  --main_dir "/scratch/rheinnec/cryo_screen" \
  --outdir "/scratch/rheinnec/cryo_screen/processed" \
  --logdir "/g/schwab/biobank_logs/cryo/wfCRYO_${timestamp}" \
  --work_dir "/scratch/rheinnec/cryo_screen/work" \
  --persistent_metadata_dir "/g/schwab/marco/central_data_processing/cryo/metadata" \
  --workflow_stage all \
  --dryrun FALSE \
  --dryrun_n 2 \
  --s3_bucket "s3test/imatrec/central_data_processing/cryo" \
  --zarr_format 3 \
  --scale_unit nm \
  --resume FALSE

# If CZI metadata is missing or inconsistent, provide fallback physical scales:
#
#   --default_x_scale 250 \
#   --default_y_scale 250 \
#   --default_z_scale 1000 \
#   --scale_unit nm \
#
# For collection-table regeneration only:
#
# bash "${repo_dir}/cryo/cryo_main.sh" cluster \
#   --sheet_mode google \
#   --sheet_url "https://docs.google.com/spreadsheets/d/REPLACE_WITH_CRYO_INPUT_SHEET_ID/edit?gid=0#gid=0" \
#   --sheet_name "cryo_input_table" \
#   --google_key "${repo_dir}/cryo/trec-tem-screen-e98a2e03f58b.json" \
#   --collection_table_url "https://docs.google.com/spreadsheets/d/REPLACE_WITH_COLLECTION_TABLE_SHEET_ID/edit?gid=0#gid=0" \
#   --collection_table_sheet "cryo_collection_table" \
#   --main_dir "/scratch/rheinnec/cryo_screen" \
#   --logdir "/g/schwab/biobank_logs/cryo/wfCRYO_${timestamp}" \
#   --work_dir "/scratch/rheinnec/cryo_screen/work" \
#   --persistent_metadata_dir "/g/schwab/marco/central_data_processing/cryo/metadata" \
#   --workflow_stage collection \
#   --dryrun FALSE \
#   --s3_bucket "s3test/imatrec/central_data_processing/cryo" \
#   --resume FALSE
