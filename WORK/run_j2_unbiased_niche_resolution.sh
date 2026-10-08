#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${project_dir}"

input_dir="results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/J2"
tstamp="$(date +%Y%m%d_%H%M%S)"
output_dir="results/reanalysis/J2_UNBIASED_NICHE_RESOLUTION_${tstamp}"
mkdir -p "${output_dir}"
printf '%s\n' "${output_dir}" > .agent/current_j2_unbiased_niche_resolution_out.txt

{
  echo "[$(date -Is)] START J2 K=3-12 resolution analysis"
  echo "OUTPUT_DIR=${output_dir}"
  bash WORK/run_voltron_r.sh --vanilla --slave \
    -f scripts/1013_j2_unbiased_niche_resolution.R --args \
    "${input_dir}" \
    "${output_dir}"
  echo "[$(date -Is)] COMPLETE J2 K=3-12 resolution analysis"
} 2>&1 | tee "${output_dir}/RUN_LOG.txt"
