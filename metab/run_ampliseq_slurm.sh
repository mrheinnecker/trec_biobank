#!/usr/bin/env bash
#SBATCH --job-name=ampliseq-pilot
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm-ampliseq-%j.out
#SBATCH --error=slurm-ampliseq-%j.err

# This long-lived controller job launches nf-core/ampliseq tasks through Slurm.
# Add site-specific #SBATCH --account/--partition directives if required.

set -euo pipefail

# sbatch runs a copy of this script from the Slurm spool directory, so paths
# cannot be resolved relative to it. main.sh passes all of them via --export.
for required_env in PARAMS_FILE CONFIG_FILE RESOURCES_CONFIG WORK_DIR NXF_APPTAINER_CACHEDIR; do
    if [[ -z "${!required_env:-}" ]]; then
        echo "Error: ${required_env} is not set. Submit this job through main.sh." >&2
        exit 2
    fi
done

readonly APPTAINER_CACHE="${NXF_APPTAINER_CACHEDIR}"
readonly PIPELINE_VERSION="${PIPELINE_VERSION:-2.18.0}"
readonly CONTAINER_PROFILE="${CONTAINER_PROFILE:-apptainer}"

for required_file in "${PARAMS_FILE}" "${CONFIG_FILE}" "${RESOURCES_CONFIG}"; do
    if [[ ! -r "${required_file}" ]]; then
        echo "Error: required file is not readable: ${required_file}" >&2
        exit 2
    fi
done

# Load Nextflow from the cluster modules, so the job does not depend on the
# shell it was submitted from. The module brings a Java >= 17, which current
# nf-core plugins need. Pin a version with e.g. NXF_MODULE=Nextflow/25.10.2;
# NXF_MODULE="" skips the module and uses whatever nextflow is on the PATH.
NXF_MODULE="${NXF_MODULE-Nextflow}"
if [[ -n "${NXF_MODULE}" ]]; then
    if ! command -v module >/dev/null 2>&1; then
        for init in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh /usr/share/lmod/lmod/init/bash; do
            if [[ -r "${init}" ]]; then
                # shellcheck disable=SC1090
                set +u; source "${init}"; set -u
                break
            fi
        done
    fi
    set +u; module load "${NXF_MODULE}"; set -u
fi

if ! command -v nextflow >/dev/null 2>&1; then
    echo "Error: Nextflow is not available. Load the appropriate cluster module." >&2
    exit 2
fi
echo "Nextflow:        $(command -v nextflow)"
echo "Java:            $(java -version 2>&1 | head -1)"

if ! command -v apptainer >/dev/null 2>&1 && ! command -v singularity >/dev/null 2>&1; then
    echo "Error: neither Apptainer nor Singularity is available." >&2
    exit 2
fi

mkdir -p -- "${WORK_DIR}" "${APPTAINER_CACHE}"
export NXF_APPTAINER_CACHEDIR="${APPTAINER_CACHE}"
export NXF_OPTS="${NXF_OPTS:--Xms1g -Xmx4g}"

echo "Starting nf-core/ampliseq ${PIPELINE_VERSION}"
echo "Parameters:      ${PARAMS_FILE}"
echo "Custom config:  ${CONFIG_FILE}"
echo "Work directory: ${WORK_DIR}"
echo "Container cache:${NXF_APPTAINER_CACHEDIR}"

nextflow run nf-core/ampliseq \
    -r "${PIPELINE_VERSION}" \
    -profile "${CONTAINER_PROFILE}" \
    -params-file "${PARAMS_FILE}" \
    -c "${CONFIG_FILE}" \
    -c "${RESOURCES_CONFIG}" \
    -work-dir "${WORK_DIR}" \
    -resume
