#!/usr/bin/env bash
set -euo pipefail

existing_l3_dir="results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/L3"
resolution_dir="results/reanalysis/L3_UNBIASED_NICHE_RESOLUTION_20261007_085730"
timestamp="$(date +%Y%m%d_%H%M%S)"
output_dir="${resolution_dir}/K08_SPATIAL_ARCHITECTURE_${timestamp}"
mkdir -p "${output_dir}"

{
  echo "[$(date -Is)] START"
  echo "SCRIPT=scripts/1011_l3_k8_spatial_architecture.R"
  echo "OUTPUT_DIR=${output_dir}"
  bash WORK/run_voltron_r.sh --vanilla --slave \
    -f scripts/1011_l3_k8_spatial_architecture.R --args \
    "${existing_l3_dir}" \
    "${resolution_dir}" \
    "${output_dir}"
  echo "[$(date -Is)] COMPLETE"
} 2>&1 | tee "${output_dir}/RUN_LOG.txt"

