#!/usr/bin/env bash
#SBATCH --job-name=wfSEM
#SBATCH --output=/scratch/rheinnec/sem_screen/wfSEM_%j.out
#SBATCH --error=/scratch/rheinnec/sem_screen/wfSEM_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
set -euo pipefail

# Edit this file to define a complete SEM screening run.
# The profile only selects where Nextflow jobs run:
#   interactive = current cluster node
#   cluster     = submit jobs through Slurm

timestamp="$(date +%Y-%m-%d_%H-%M)"

repo_dir="/g/schwab/marco/repos/trec_biobank"
scratch_dir="/scratch/rheinnec/sem_screen"

mkdir -p "${scratch_dir}"
cd "${scratch_dir}"

bash "${repo_dir}/sem/sem_main.sh" interactive \
  --sheet_mode google \
  --sheet_url "https://docs.google.com/spreadsheets/d/1jcpyMkSR4npSxST3D5cFzkAIi9UPmwbPvzAdr2ws55U/edit?gid=2132397683#gid=2132397683" \
  --google_key "${repo_dir}/trec-tem-screen-e98a2e03f58b.json" \
  --collection_table_url "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=1595862976#gid=1595862976" \
  --main_dir "/g/schwab/sem_screen" \
  --rawdir "/g/schwab/Chandni/SEM/IMATREC SEM" \
  --outdir "/g/schwab/sem_screen/processed" \
  --logdir "/g/schwab/biobank_logs/sem/wfSEM_${timestamp}" \
  --work_dir "${scratch_dir}/work" \
  --local_log "/g/schwab/sem_screen/sem_image_log_local.tsv" \
  --workflow_stage all \
  --dryrun FALSE \
  --s3_bucket "s3test/imatrec/central_data_processing/sem" \
  --resume FALSE
