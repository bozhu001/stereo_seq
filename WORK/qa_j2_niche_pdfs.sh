#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${project_dir}"
output_dir="$(cat .agent/current_j2_unbiased_niche_resolution_out.txt)"
mkdir -p tmp/pdfs/j2_niche_qa

for pdf in "${output_dir}"/figures/0*_J2_*.pdf; do
  stem="$(basename "${pdf}" .pdf)"
  echo "PDF=${pdf}"
  pdfinfo "${pdf}" | grep -E '^(Pages|Page size|File size)'
  pdftoppm \
    -f 1 \
    -singlefile \
    -png \
    -r 110 \
    "${pdf}" \
    "tmp/pdfs/j2_niche_qa/${stem}" \
    >/dev/null 2>&1
done
