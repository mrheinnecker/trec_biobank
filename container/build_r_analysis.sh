#!/usr/bin/env bash
## Builds r_analysis.sif from r_analysis.def without running %post inside `apptainer build`.
## Needed on WSL2, where `apptainer build` of any definition with a %post section fails with
##   FATAL: failed to mount binfmt_misc filesystem: invalid argument
## Steps: base image -> writable sandbox (in a Linux filesystem, not /mnt/c) -> the %post section
## of the .def runs there with --fakeroot -> sandbox packed into the .sif.
## On a normal Linux machine `apptainer build --fakeroot r_analysis.sif r_analysis.def` does the same.
##
## usage: bash container/build_r_analysis.sh [output.sif]   (default: /mnt/c/container/r_analysis.sif)
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
def="$script_dir/r_analysis.def"
out=${1:-/mnt/c/container/r_analysis.sif}

base=$(awk '/^From:/ {print $2; exit}' "$def")
work=$(mktemp -d "${TMPDIR:-/tmp}/r_analysis_build.XXXXXX")
trap 'rm -rf "$work"' EXIT

echo "== base image docker://$base -> sandbox $work/sandbox"
apptainer build --sandbox "$work/sandbox" "docker://$base"

echo "== %post of $def"
awk '/^%post/ {f=1; next} /^%/ {f=0} f' "$def" > "$work/post.sh"
## /tmp of the container is the host /tmp, so the script goes to /opt of the sandbox
cp "$work/post.sh" "$work/sandbox/opt/build_post.sh"
apptainer exec --fakeroot --writable --no-home --pwd / "$work/sandbox" bash -e /opt/build_post.sh
rm -f "$work/sandbox/opt/build_post.sh"

echo "== sandbox -> $out"
apptainer build --force "$out" "$work/sandbox"

echo "== test"
apptainer exec "$out" Rscript -e 'library(lme4); library(ComplexHeatmap); library(googlesheets4); cat("ok\n")'
apptainer exec "$out" python3 -c 'import PIL; print("Pillow", PIL.__version__)'
