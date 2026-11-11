#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  bash overview/run_collection_table_overview.sh --outdir /path/to/output [options]

Options:
  --outdir PATH       Output directory for figures and TSVs.
                      Default: ~/collection_table_overview
  --modalities CSV    Optional subset, e.g. TEM,HITT,CRYO.
  --google_key PATH   Optional Google service-account JSON key.
  --authenticated     Use --google_key / GOOGLE_KEY auth instead of anonymous access.
  --container PATH    Singularity image. Default: ~/container/tidyverse_latest.sif
  --prefix VALUE      Output filename prefix. Default: collection_table_overview
  --help              Show this message.
EOF
}

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
outdir="${COLLECTION_OVERVIEW_OUTDIR:-${HOME}/collection_table_overview}"
container="${COLLECTION_OVERVIEW_CONTAINER:-${HOME}/container/tidyverse_latest.sif}"
prefix="collection_table_overview"
modalities=""
google_key="${GOOGLE_KEY:-}"
anonymous="TRUE"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h)
      usage
      exit 0
      ;;
    --outdir)
      outdir="${2:?--outdir requires a path}"
      shift 2
      ;;
    --outdir=*)
      outdir="${1#*=}"
      shift
      ;;
    --modalities)
      modalities="${2:?--modalities requires a comma-separated value}"
      shift 2
      ;;
    --modalities=*)
      modalities="${1#*=}"
      shift
      ;;
    --google_key|--google-key)
      google_key="${2:?--google_key requires a path}"
      anonymous="FALSE"
      shift 2
      ;;
    --google_key=*|--google-key=*)
      google_key="${1#*=}"
      anonymous="FALSE"
      shift
      ;;
    --authenticated)
      anonymous="FALSE"
      shift
      ;;
    --container)
      container="${2:?--container requires a path}"
      shift 2
      ;;
    --container=*)
      container="${1#*=}"
      shift
      ;;
    --prefix)
      prefix="${2:?--prefix requires a value}"
      shift 2
      ;;
    --prefix=*)
      prefix="${1#*=}"
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

mkdir -p "$outdir"

cmd=(
  singularity exec
  -B /g
  -B /home
  "$container"
  Rscript "$repo_dir/overview/collection_table_overview.R"
  --outdir "$outdir"
  --prefix "$prefix"
)

if [[ -n "$modalities" ]]; then
  cmd+=(--modalities "$modalities")
fi

if [[ "$anonymous" == "TRUE" ]]; then
  cmd+=(--anonymous)
elif [[ -n "$google_key" ]]; then
  cmd+=(--google_key "$google_key")
fi

printf 'Writing collection overview to: %s\n' "$outdir"
"${cmd[@]}"
