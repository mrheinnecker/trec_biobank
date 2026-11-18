#!/usr/bin/env bash

# Local nf-core/ampliseq pilot (no Slurm). Uses the same samplesheet builder and
# primer-orientation config as the cluster run (main.sh), so a pilot that works
# here checks the real setup: Cutadapt --revcomp, DADA2 per run, controls.
# Edit the variables at the top, then: bash launcher.sh

set -euo pipefail

## -- edit these ------------------------------------------------------------
RAW_DIR="/mnt/c/projects/trec_biobank/metab/datahub_downloads"
ENVIRONMENTAL_SAMPLES="8"          # number or "all"
CONTROL_SAMPLES="all"              # number or "all"
FW_PRIMER="GTGYCAGCMGCCGCGGTAA"
RV_PRIMER="CCGYCAATTYMTTTRAGTTT"
BINNED_QUALITY="2,11,25,37"
MERGEPAIRS_STRATEGY="consensus"   # merge 16S by overlap, concatenate non-overlapping 18S
PIPELINE_VERSION="2.18.0"
MAX_CPUS="16"
MAX_MEMORY="24.GB"
SUBSAMPLE_PAIRS="100000"           # first N read pairs per sample; "all" for full files
RUN_ROOT="/home/rheinnec/ampliseq_pilot"
## ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${RUN_ROOT}/work"
RESULTS_DIR="${RUN_ROOT}/results"
SAMPLESHEET="${RUN_ROOT}/ampliseq_samplesheet.tsv"
PARAMS_FILE="${RUN_ROOT}/ampliseq_params.yaml"
LOCAL_CONFIG="${RUN_ROOT}/local_resources.config"

mkdir -p -- "${WORK_DIR}" "${RESULTS_DIR}"

# Incomplete downloads (empty or truncated .gz) are skipped, not fatal.
SKIP_INCOMPLETE=1 CHECK_GZIP=1 \
    "${SCRIPT_DIR}/make_ampliseq_samplesheet.sh" \
    "${RAW_DIR}" "${SAMPLESHEET}" "${SCRIPT_DIR}/sampleinfo.tsv" \
    "${ENVIRONMENTAL_SAMPLES}" "${CONTROL_SAMPLES}"

# ampliseq counts every FASTQ with Nextflow's countFastq, one file after the
# other and three times per run (input, after Cutadapt, after filtering). That
# takes minutes per full file, so the pilot uses the first N pairs instead.
if [[ "${SUBSAMPLE_PAIRS}" != "all" ]]; then
    SUBSAMPLE_DIR="${RUN_ROOT}/fastq_subsampled"
    mkdir -p -- "${SUBSAMPLE_DIR}"
    lines=$((SUBSAMPLE_PAIRS * 4))
    awk -F'\t' 'NR > 1 { print $2; print $3 }' "${SAMPLESHEET}" | while read -r fastq; do
        out="${SUBSAMPLE_DIR}/$(basename -- "${fastq}")"
        [[ -s "${out}" ]] && continue
        # head closes the pipe early; ignore zcat's resulting SIGPIPE.
        { zcat -- "${fastq}" || true; } | head -n "${lines}" | gzip -1 > "${out}.tmp"
        mv -- "${out}.tmp" "${out}"
    done
    awk -F'\t' -v OFS='\t' -v dir="${SUBSAMPLE_DIR}" '
        NR > 1 { n = split($2, a, "/"); $2 = dir "/" a[n]; n = split($3, a, "/"); $3 = dir "/" a[n] }
        { print }' "${SAMPLESHEET}" > "${SAMPLESHEET}.tmp"
    mv -- "${SAMPLESHEET}.tmp" "${SAMPLESHEET}"
    echo "Subsampled to the first ${SUBSAMPLE_PAIRS} read pairs per sample in ${SUBSAMPLE_DIR}"
fi

{
    printf 'input: "%s"\n' "${SAMPLESHEET}"
    printf 'outdir: "%s"\n' "${RESULTS_DIR}"
    printf 'FW_primer: "%s"\n' "${FW_PRIMER}"
    printf 'RV_primer: "%s"\n' "${RV_PRIMER}"
    printf 'binned_quality: "%s"\n' "${BINNED_QUALITY}"
    printf 'mergepairs_strategy: "%s"\n' "${MERGEPAIRS_STRATEGY}"
    printf 'skip_taxonomy: true\n'
    printf 'min_read_counts: 1\n'
    # Clean blanks can lose all of their few reads; drop them instead of failing.
    printf 'ignore_failed_trimming: true\n'
    printf 'ignore_failed_filtering: true\n'
} > "${PARAMS_FILE}"

# Resource ceilings for this machine; nf-core 2.18 uses resourceLimits.
cat > "${LOCAL_CONFIG}" <<EOF
process {
    resourceLimits = [ cpus: ${MAX_CPUS}, memory: ${MAX_MEMORY}, time: 24.h ]
}
EOF

export NXF_APPTAINER_CACHEDIR="${NXF_APPTAINER_CACHEDIR:-${HOME}/.apptainer/nxf_cache}"
mkdir -p -- "${NXF_APPTAINER_CACHEDIR}"

echo "Samplesheet: ${SAMPLESHEET}"
echo "Params:      ${PARAMS_FILE}"
echo "Results:     ${RESULTS_DIR}"
echo "Work:        ${WORK_DIR}"
echo "Running nf-core/ampliseq ${PIPELINE_VERSION} locally; log in ${RUN_ROOT}/.nextflow.log"

# Run from RUN_ROOT so .nextflow.log and the resume cache live with the run.
cd -- "${RUN_ROOT}"
nextflow run nf-core/ampliseq \
    -r "${PIPELINE_VERSION}" \
    -profile apptainer \
    -params-file "${PARAMS_FILE}" \
    -c "${SCRIPT_DIR}/ampliseq_orientation.config" \
    -c "${SCRIPT_DIR}/nextflow.local.config" \
    -c "${LOCAL_CONFIG}" \
    -work-dir "${WORK_DIR}" \
    -resume
