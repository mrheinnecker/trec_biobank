#!/usr/bin/env bash
#SBATCH --job-name=ampliseq-taxonomy
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-taxonomy-%j.out
#SBATCH --error=slurm-taxonomy-%j.err

# Controller job for wfTAXONOMY.nf: SILVA (16S) and PR2 (18S) classification of
# finished nf-core/ampliseq results. Nextflow submits the actual tasks to Slurm.
#
#   sbatch run_taxonomy_slurm.sh
#
# Every path below can be overridden without editing the file, e.g.
#   sbatch --export=ALL,AMPLISEQ_RESULTS=/scratch/rheinnec/other_run/results run_taxonomy_slurm.sh
# Fresh run without the Nextflow cache:
#   sbatch --export=ALL,RESUME=false run_taxonomy_slurm.sh

set -euo pipefail

# --- adjust these paths ------------------------------------------------------
# sbatch runs a copy of this script from the spool directory, so the repo
# location cannot be derived from the script path.
REPO_TAXONOMY_DIR="${REPO_TAXONOMY_DIR:-/g/schwab/marco/repos/tem_classification/metab/taxonomy}"
AMPLISEQ_RESULTS="${AMPLISEQ_RESULTS:-/scratch/rheinnec/ampliseq_run/results}"
REFDB="${REFDB:-/g/schwab/marco/refdb/metab}"
SILVA_DB="${SILVA_DB:-${REFDB}/silva_nr99_v138.2_toGenus_trainset.fa.gz}"
PR2_DB="${PR2_DB:-${REFDB}/pr2_version_5.1.1_SSU_dada2.fasta.gz}"
DADA2_SIF="${DADA2_SIF:-/g/schwab/marco/container_devel/dada2_1.38.sif}"
RUN_DIR="${RUN_DIR:-/scratch/rheinnec/taxonomy_run}"
APPTAINER_CACHE="${APPTAINER_CACHE:-/scratch/rheinnec/apptainer_cache}"
TAX_CPUS="${TAX_CPUS:-16}"
TAX_MEMORY="${TAX_MEMORY:-32 GB}"
RESUME="${RESUME:-true}"        # false = run every step again, ignore the cache
# -----------------------------------------------------------------------------

for required_file in \
    "${REPO_TAXONOMY_DIR}/wfTAXONOMY.nf" \
    "${REPO_TAXONOMY_DIR}/nextflow.config" \
    "${AMPLISEQ_RESULTS}/dada2/ASV_seqs.fasta" \
    "${AMPLISEQ_RESULTS}/dada2/ASV_table.tsv" \
    "${AMPLISEQ_RESULTS}/barrnap/summary.tsv" \
    "${SILVA_DB}" \
    "${PR2_DB}" \
    "${DADA2_SIF}"; do
    if [[ ! -r "${required_file}" ]]; then
        echo "Error: required file is not readable: ${required_file}" >&2
        exit 2
    fi
done

if ! command -v nextflow >/dev/null 2>&1; then
    echo "Error: Nextflow is not available. Load the appropriate cluster module." >&2
    exit 2
fi

mkdir -p -- "${RUN_DIR}" "${APPTAINER_CACHE}"
export NXF_APPTAINER_CACHEDIR="${APPTAINER_CACHE}"
export NXF_OPTS="${NXF_OPTS:--Xms1g -Xmx2g}"

echo "ampliseq results: ${AMPLISEQ_RESULTS}"
echo "Output:           ${AMPLISEQ_RESULTS}/taxonomy"
echo "SILVA:            ${SILVA_DB}"
echo "PR2:              ${PR2_DB}"
echo "Container:        ${DADA2_SIF}"
echo "Run/work dir:     ${RUN_DIR}"
echo "Resume:           ${RESUME}"

RESUME_ARG=()
[[ "${RESUME}" == "true" ]] && RESUME_ARG=(-resume)

# Run from RUN_DIR so .nextflow/ and the logs land there, not in the repo.
cd -- "${RUN_DIR}"
nextflow run "${REPO_TAXONOMY_DIR}/wfTAXONOMY.nf" \
    -c "${REPO_TAXONOMY_DIR}/nextflow.config" \
    -profile slurm \
    -work-dir "${RUN_DIR}/work" \
    --ampliseq_results "${AMPLISEQ_RESULTS}" \
    --silva_db "${SILVA_DB}" \
    --pr2_db "${PR2_DB}" \
    --container_dada2 "${DADA2_SIF}" \
    --cpus "${TAX_CPUS}" \
    --memory "${TAX_MEMORY}" \
    ${RESUME_ARG[@]+"${RESUME_ARG[@]}"}
