#!/usr/bin/env bash
set -euo pipefail

formal_j2_dir="results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/J2"
resolution_dir="results/reanalysis/J2_UNBIASED_NICHE_RESOLUTION_20261008_110522"
timestamp="$(date +%Y%m%d_%H%M%S)"
output_dir="${resolution_dir}/K05_SPATIAL_ARCHITECTURE_${timestamp}"
mkdir -p "${output_dir}"

{
  echo "[$(date -Is)] START"
  echo "SCRIPT=scripts/1018_j2_k5_spatial_architecture.R"
  echo "OUTPUT_DIR=${output_dir}"
  bash WORK/run_voltron_r.sh --vanilla --slave \
    -f scripts/1018_j2_k5_spatial_architecture.R --args \
    "${formal_j2_dir}" \
    "${resolution_dir}" \
    "${output_dir}"
  echo "[$(date -Is)] COMPLETE"
} 2>&1 | tee "${output_dir}/RUN_LOG.txt"

printf '%s\n' "${output_dir}" > .agent/current_j2_k5_spatial_architecture_out.txt
