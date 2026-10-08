#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${project_dir}"

input_dir="results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/J2"
output_dir="$(cat .agent/current_j2_unbiased_niche_resolution_out.txt)"

{
  echo "[$(date -Is)] START J2 primary-K composition and UMI sensitivity"
  bash WORK/run_voltron_r.sh --vanilla --slave \
    -f scripts/1014_j2_primaryk_composition_consistency.R --args \
    "${input_dir}" \
    "${output_dir}"
  echo "[$(date -Is)] COMPLETE J2 primary-K composition and UMI sensitivity"
} 2>&1 | tee "${output_dir}/PRIMARYK_RUN_LOG.txt"
