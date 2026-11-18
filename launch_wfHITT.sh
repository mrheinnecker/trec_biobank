#!/usr/bin/env bash
#SBATCH --job-name=wfHITT
#SBATCH --output=/scratch/rheinnec/hitt_screen/wfHITT_%j.out
#SBATCH --error=/scratch/rheinnec/hitt_screen/wfHITT_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
set -euo pipefail

# Edit this file to define a complete screening run.
# Keep passwords outside the file, for example:
#   export HITT_SSHPASS='PASSWORD'
#
# The profile only selects where Nextflow jobs run:
#   interactive = current cluster node
#   cluster     = submit jobs through Slurm

#repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
timestamp="$(date +%Y-%m-%d_%H-%M)"

#timestamp="2026-06-01_20-54"

repo_dir="/g/schwab/marco/repos/trec_biobank"
cd /scratch/rheinnec/hitt_screen
# HITT workflow
# To use a local table instead, replace the sheet options with:
#   --sheet_mode local \
#   --input_table "/path/to/hitt_images.tsv" \
bash "${repo_dir}/hitt/hitt_main.sh" cluster \
  --sheet_mode google \
  --sheet_url "https://docs.google.com/spreadsheets/d/1ePRpa56mmMvCeRTLXmwOywOLy5_I3AFrxJepSUYGR1s/edit?gid=0#gid=0" \
  --sheet_name "hitt_input_table" \
  --google_key "${repo_dir}/trec-tem-screen-e98a2e03f58b.json" \
  --collection_table_url "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=1062275333#gid=1062275333" \
  --collection_table_sheet "hitt_collection_table" \
  --main_dir "/scratch/rheinnec/hitt_screen" \
  --logdir "/scratch/rheinnec/hitt_screen/logs/wfHITT_${timestamp}" \
  --work_dir "/scratch/rheinnec/hitt_screen/work" \
  --workflow_stage all \
  --dryrun FALSE \
  --s3_bucket "s3embl/imatrec/central_data_processing/hitt" \
  --x_scale 650 \
  --y_scale 650 \
  --z_scale 650 \
  --input_suffix "recon_111_1/tomo" \
  --output_name "omezarr" \
  --overwrite TRUE \
  --convert_uint16 FALSE \
  --uint16_lower_percentile 0.1 \
  --uint16_upper_percentile 99.9 \
  --uint16_sample_values 2000000 \
  --copy_data TRUE \
  --copy_dest_root "/scratch/rheinnec/tmp_hitt" \
  --copy_max_forks 10 \
  --crop_stack TRUE \
  --crop_bright_threshold "auto" \
  --crop_auto_percentile 98.0 \
  --crop_min_bright_fraction 0.0025 \
  --crop_padding_low_slices 30 \
  --crop_padding_high_slices 30 \
  --crop_bridge_gap_slices 3 \
  --crop_min_run_slices 3 \
  --crop_sample_values_per_slice 100000 \
  --remote_user "p3l-yschwab" \
  --remote_host "cerberus.embl-hamburg.de" \
  --remote_port 22443 \
  --password $HITT_SSHPASS \
  --resume FALSE

# Add other workflow launcher calls below when they are ready, for example:
#
# bash "${repo_dir}/tem/main.sh" cluster \
#   --workflow_stage all
#
# bash "${repo_dir}/sem/sem_main.sh" cluster \
#   --workflow_stage all
