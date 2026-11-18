#!/usr/bin/env bash

# This file is a run template only.
# Copy one command block into your terminal, edit the values there, and execute it.
echo "Template only: copy one command block from this file and run it manually."
exit 0

# Full cluster production run.
bash /g/schwab/marco/repos/trec_biobank/tem/main.sh \
  --profile cluster \
  --resume TRUE \
  --dryrun FALSE \
  --dryrun_n 10 \
  --sheet_mode google \
  --sheet_url "https://docs.google.com/spreadsheets/d/143uVeeJ72SQE5eK01lzWYCEiT7pJUF3lX7hJl3R9s9I/edit?gid=258669282#gid=258669282" \
  --google_key "/g/schwab/marco/repos/trec_biobank/trec-tem-screen-e98a2e03f58b.json" \
  --collection_table_url "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=1643802951#gid=1643802951" \
  --collection_table_sheet "tem_collection_table" \
  --workflow_stage all \
  --main_dir /g/schwab/tem_screen \
  --work_dir /scratch/rheinnec/tem_screen/work \
  --rawdir /g/schwab/tem_screen/raw \
  --pngdir /g/schwab/tem_screen/pngs \
  --outdir /g/schwab/tem_screen/processed \
  --local_log /g/schwab/tem_screen/image_log_local.tsv \
  --s3_bucket s3test/imatrec/central_data_processing/tem \
  --gradient_chunk_rows 512 \
  --gradient_downsample 16

# Interactive cluster debugging after allocating a node.
bash /g/schwab/marco/repos/trec_biobank/tem/main.sh \
  --profile interactive \
  --resume TRUE \
  --dryrun FALSE \
  --dryrun_n 10 \
  --sheet_mode google \
  --sheet_url "https://docs.google.com/spreadsheets/d/143uVeeJ72SQE5eK01lzWYCEiT7pJUF3lX7hJl3R9s9I/edit?gid=258669282#gid=258669282" \
  --google_key "/g/schwab/marco/repos/trec_biobank/trec-tem-screen-e98a2e03f58b.json" \
  --collection_table_url "https://docs.google.com/spreadsheets/d/15WNNnse7OvlfiJwFOFYbQA4zIp-5nKc0icRZYfJS--o/edit?gid=1643802951#gid=1643802951" \
  --collection_table_sheet "tem_collection_table" \
  --workflow_stage all \
  --main_dir /scratch/rheinnec/tem_screen \
  --work_dir /scratch/rheinnec/tem_screen/work \
  --rawdir /scratch/rheinnec/tem_screen/raw \
  --pngdir /scratch/rheinnec/tem_screen/pngs \
  --outdir /scratch/rheinnec/tem_screen/processed \
  --local_log /scratch/rheinnec/tem_screen/image_log_local.tsv \
  --s3_bucket s3test/imatrec/central_data_processing/tem \
  --gradient_chunk_rows 512 \
  --gradient_downsample 16

# Local Windows discovery run with copied test files.
bash C:/repos/trec_biobank/tem/main.sh \
  --profile local \
  --resume TRUE \
  --dryrun TRUE \
  --dryrun_n 10 \
  --sheet_mode local \
  --workflow_stage discover \
  --main_dir C:/projects/tem_screen \
  --work_dir C:/projects/tem_screen/work \
  --rawdir C:/projects/tem_screen/raw \
  --pngdir C:/projects/tem_screen/pngs \
  --outdir C:/projects/tem_screen/processed \
  --local_log C:/projects/tem_screen/image_log_local.tsv








bash sem_main.sh \
  --profile interactive \
  --main_dir /g/schwab/sem_screen \
  --work_dir /scratch/rheinnec/sem_screen/work \
  --workflow_stage all \
  --dryrun FALSE
  --resume FALSE



