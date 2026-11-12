#!/usr/bin/env bash
# Download all metab FASTQs from the TREC datahub and verify them by MD5.
#
#   bash download_datahub.sh LINK_FILE SAMPLEINFO OUT_DIR [PARALLEL]
#
# LINK_FILE   one getFile.cgi link per line (the datahub "download all" list;
#             the links carry no file name, the server sends it)
# SAMPLEINFO  datahub file export (TSV with Name and MD5 columns), e.g.
#             raw_input/sampleinfo.tsv
# OUT_DIR     download directory
# PARALLEL    simultaneous downloads (default 4)
#
# Files are written as NAME.part and renamed only after the download finished,
# so an interrupted run never leaves a truncated file under the final name.
# Rerunning skips files that already exist and match their MD5.

set -euo pipefail

if (( $# < 3 )); then
    echo "Usage: download_datahub.sh LINK_FILE SAMPLEINFO OUT_DIR [PARALLEL]" >&2
    exit 1
fi
LINK_FILE="$1"
SAMPLEINFO="$2"
OUT_DIR="$3"
PARALLEL="${4:-4}"

mkdir -p -- "${OUT_DIR}"
EXPECTED="${OUT_DIR}/expected_md5.txt"

# "MD5  Name" list from the export, located by column name.
awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
            { print $col["MD5"] "  " $col["Name"] }' "${SAMPLEINFO}" > "${EXPECTED}"
echo "expected files: $(wc -l < "${EXPECTED}")"

download_one() {
    local url="$1" out_dir="$2" expected="$3" name md5
    # File name from the Content-Disposition header.
    name=$(curl -sS -I -L --retry 3 "${url}" | tr -d '\r' |
           sed -n 's/^[Cc]ontent-[Dd]isposition:.*filename="\{0,1\}\([^";]*\)"\{0,1\}.*/\1/p' | tail -1)
    if [[ -z "${name}" ]]; then
        echo "NO NAME  ${url}" >&2
        return 1
    fi
    md5=$(awk -v n="${name}" '$2 == n { print $1 }' "${expected}")
    if [[ -f "${out_dir}/${name}" && -n "${md5}" ]] &&
       [[ "$(md5sum "${out_dir}/${name}" | cut -d' ' -f1)" == "${md5}" ]]; then
        echo "skip     ${name}"
        return 0
    fi
    if curl -sS -L --fail --retry 3 -o "${out_dir}/${name}.part" "${url}"; then
        mv -- "${out_dir}/${name}.part" "${out_dir}/${name}"
        echo "done     ${name}"
    else
        echo "FAILED   ${name}" >&2
        return 1
    fi
}
export -f download_one

# CRLF line endings are stripped; empty lines ignored.
tr -d '\r' < "${LINK_FILE}" | grep -v '^\s*$' |
    xargs -P "${PARALLEL}" -I{} bash -c 'download_one "$1" "$2" "$3"' _ {} "${OUT_DIR}" "${EXPECTED}" ||
    echo "some downloads failed (see above); rerun to retry them" >&2

echo "checking MD5 ..."
cd -- "${OUT_DIR}"
if md5sum -c --quiet expected_md5.txt; then
    echo "all $(wc -l < expected_md5.txt) files present and correct"
else
    echo "files above are missing or corrupt; rerun this script to fetch them again" >&2
    exit 1
fi
