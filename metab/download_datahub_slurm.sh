#!/usr/bin/env bash
#SBATCH --job-name=datahub-download
#SBATCH --array=0-43%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-datahub-%A_%a.out
#SBATCH --error=slurm-datahub-%A_%a.err

# Download the 436 URLs in datahub_links.txt as 44 Slurm array tasks.
# Each task handles at most 10 files; no more than 8 tasks run concurrently.
#
# Submit with the default output directory:
#   sbatch download_datahub_slurm.sh
#
# Or select a different output directory:
#   sbatch --export=ALL,DOWNLOAD_DIR=/path/to/output download_datahub_slurm.sh

set -uo pipefail

readonly FILES_PER_TASK=10
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly LINK_FILE="/g/schwab/marco/repos/trec_biobank/metab/datahub_links.txt"
readonly DOWNLOAD_DIR="${DOWNLOAD_DIR:-/mnt/c/projects/trec_biobank/metab/datahub_downloads}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "Error: submit this script with sbatch; SLURM_ARRAY_TASK_ID is unset." >&2
    exit 2
fi

if [[ ! -r "${LINK_FILE}" ]]; then
    echo "Error: cannot read link file: ${LINK_FILE}" >&2
    exit 2
fi

if ! command -v wget >/dev/null 2>&1; then
    echo "Error: wget is not available on this compute node." >&2
    exit 2
fi

mkdir -p -- "${DOWNLOAD_DIR}"

readonly TOTAL_LINKS="$(awk 'NF { count++ } END { print count + 0 }' "${LINK_FILE}")"
readonly START_INDEX=$((SLURM_ARRAY_TASK_ID * FILES_PER_TASK))
readonly END_INDEX=$((START_INDEX + FILES_PER_TASK))

if (( START_INDEX >= TOTAL_LINKS )); then
    echo "Array task ${SLURM_ARRAY_TASK_ID}: no links assigned; exiting."
    exit 0
fi

echo "Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "Link range:  $((START_INDEX + 1))-$((END_INDEX < TOTAL_LINKS ? END_INDEX : TOTAL_LINKS)) of ${TOTAL_LINKS}"
echo "Destination: ${DOWNLOAD_DIR}"

failed=0
processed=0
line_number=0

while IFS= read -r url || [[ -n "${url}" ]]; do
    url="${url%$'\r'}"
    [[ -z "${url}" ]] && continue

    if (( line_number >= START_INDEX && line_number < END_INDEX )); then
        processed=$((processed + 1))
        echo "Downloading link $((line_number + 1)) of ${TOTAL_LINKS}: ${url}"

        if ! wget \
            --continue \
            --content-disposition \
            --trust-server-names \
            --tries=5 \
            --timeout=60 \
            --waitretry=10 \
            --retry-connrefused \
            --directory-prefix="${DOWNLOAD_DIR}" \
            "${url}"; then
            echo "Error: download failed for link $((line_number + 1)): ${url}" >&2
            failed=$((failed + 1))
        fi
    fi

    line_number=$((line_number + 1))
    (( line_number >= END_INDEX )) && break
done < "${LINK_FILE}"

echo "Array task ${SLURM_ARRAY_TASK_ID} finished: ${processed} processed, ${failed} failed."

if (( failed > 0 )); then
    exit 1
fi
