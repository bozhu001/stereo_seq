#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
timestamp="$(date +%y%m%d%H%M%S)"
log_path="${project_dir}/logs/1101_k6_factor_downstream_${timestamp}.log"

mkdir -p "${project_dir}/logs"
cd "${project_dir}"

/home/bozhu/miniconda3/envs/banksy_stereo/bin/python \
  scripts/1101_validate_k6_factor_expression_and_patient_metrics.py \
  --output-root results/reanalysis \
  >"${log_path}" 2>&1

printf '%s\n' "${log_path}"
tail -n 1 "${log_path}"
