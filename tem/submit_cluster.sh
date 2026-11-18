#!/usr/bin/env bash
#SBATCH --job-name=wfTEM
#SBATCH --output=wfTEM_%j.out
#SBATCH --error=wfTEM_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00

set -euo pipefail

PATH="/usr/local/bin:/usr/bin:/bin:${PATH:-}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

bash "${script_dir}/main.sh" \
  --profile cluster \
  --resume "${RESUME:-TRUE}" \
  --dryrun "${DRYRUN:-FALSE}" \
  --dryrun_n "${DRYRUN_N:-10}" \
  --sheet_mode "${SHEET_MODE:-google}" \
  --sheet_url "${SHEET_URL:-https://docs.google.com/spreadsheets/d/143uVeeJ72SQE5eK01lzWYCEiT7pJUF3lX7hJl3R9s9I/edit?gid=258669282#gid=258669282}" \
  --google_key "${GOOGLE_KEY:-/g/schwab/marco/repos/trec_biobank/trec-tem-screen-e98a2e03f58b.json}" \
  --collection_table_url "${COLLECTION_TABLE_URL:-https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=1643802951#gid=1643802951}" \
  --collection_table_sheet "${COLLECTION_TABLE_SHEET:-tem_collection_table}" \
  --workflow_stage "${WORKFLOW_STAGE:-all}" \
  --main_dir "${TEM_SCREEN_DIR:-/g/schwab/tem_screen}" \
  --work_dir "${WORK_DIR:-/scratch/rheinnec/tem_screen/work}" \
  --s3_bucket "${S3_BUCKET:-s3test/imatrec/central_data_processing/tem}" \
  --gradient_chunk_rows "${GRADIENT_CHUNK_ROWS:-512}" \
  --gradient_downsample "${GRADIENT_DOWNSAMPLE:-16}"
