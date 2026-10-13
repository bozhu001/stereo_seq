#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
timestamp="$(date +%y%m%d%H%M%S)"
notebook="${project_root}/scripts/1106_b_associated_largest_component_audit_report.ipynb"
output_name="1106_b_associated_largest_component_audit_report_${timestamp}.html"
log_path="${project_root}/logs/1106_b_associated_largest_component_report_${timestamp}.log"

/home/bozhu/miniconda3/envs/stereo/bin/jupyter nbconvert \
  --to html \
  --execute "${notebook}" \
  --ExecutePreprocessor.timeout=600 \
  --output "${output_name}" \
  --output-dir "${project_root}/results" \
  >"${log_path}" 2>&1

printf '%s\n' "${project_root}/results/${output_name}"
printf '%s\n' "${log_path}"
