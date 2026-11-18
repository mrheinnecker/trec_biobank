#!/usr/bin/env bash
#SBATCH --job-name=wfEXM
#SBATCH --output=/scratch/rheinnec/exm_screen/wfEXM_%j.out
#SBATCH --error=/scratch/rheinnec/exm_screen/wfEXM_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
set -euo pipefail

# Edit this file to define a complete EXM screening run.
# The profile only selects where Nextflow jobs run:
#   interactive = current cluster node
#   cluster     = submit jobs through Slurm

timestamp="$(date +%Y-%m-%d_%H-%M)"

repo_dir="/g/schwab/marco/repos/trec_biobank"
cd /scratch/rheinnec/exm_screen

bash "${repo_dir}/exm/exm_main.sh" interactive \
  --sheet_mode google \
  --sheet_url "https://docs.google.com/spreadsheets/d/1ePRpa56mmMvCeRTLXmwOywOLy5_I3AFrxJepSUYGR1s/edit?gid=1331814489#gid=1331814489" \
  --sheet_name "exm_input_table" \
  --google_key "${repo_dir}/trec-tem-screen-e98a2e03f58b.json" \
  --collection_table_url "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit" \
  --collection_table_sheet "exm_collection_table" \
  --main_dir "/scratch/rheinnec/exm_screen" \
  --outdir "/scratch/rheinnec/exm_screen/processed" \
  --logdir "/g/schwab/biobank_logs/exm/wfEXM_${timestamp}" \
  --work_dir "/scratch/rheinnec/exm_screen/work" \
  --persistent_metadata_dir "/g/schwab/marco/central_data_processing/exm/metadata" \
  --workflow_stage all \
  --dryrun FALSE \
  --s3_bucket "s3test/imatrec/central_data_processing/exm" \
  --zarr_format 3 \
  --scale_unit nm \
  --resume FALSE
  
  
  
#--dryrun_images "BAR_Img0098,BAR_Img0099" \
# If the .nd2 pixel size should be overridden (e.g. pre-expansion size), provide physical scales:
#
#   --default_x_scale 250 \
#   --default_y_scale 250 \
#   --default_z_scale 1000 \
#   --scale_unit nm \
#
# For collection-table regeneration only:
#
# bash "${repo_dir}/exm/exm_main.sh" cluster \
#   --sheet_mode google \
#   --sheet_url "https://docs.google.com/spreadsheets/d/1ePRpa56mmMvCeRTLXmwOywOLy5_I3AFrxJepSUYGR1s/edit?gid=1331814489#gid=1331814489" \
#   --sheet_name "exm_input_table" \
#   --google_key "${repo_dir}/trec-tem-screen-e98a2e03f58b.json" \
#   --collection_table_url "https://docs.google.com/spreadsheets/d/REPLACE_WITH_COLLECTION_TABLE_SHEET_ID/edit?gid=0#gid=0" \
#   --collection_table_sheet "exm_collection_table" \
#   --main_dir "/scratch/rheinnec/exm_screen" \
#   --logdir "/g/schwab/biobank_logs/exm/wfEXM_${timestamp}" \
#   --work_dir "/scratch/rheinnec/exm_screen/work" \
#   --persistent_metadata_dir "/g/schwab/marco/central_data_processing/exm/metadata" \
#   --workflow_stage collection \
#   --dryrun FALSE \
#   --s3_bucket "s3test/imatrec/central_data_processing/exm" \
#   --resume FALSE
