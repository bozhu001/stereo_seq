#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python_bin="/home/bozhu/miniconda3/envs/stereo/bin/python"
timestamp="$(date +%y%m%d%H%M%S)"
run_stamp="$(date +%Y%m%d_%H%M%S)"
output_dir="${project_root}/results/reanalysis/REFERENCE18_K7_ROBUSTNESS_REPORT_${run_stamp}"
log_path="${project_root}/logs/reference18_k7_robustness_report_${timestamp}.log"
notebook="${project_root}/scripts/1147_reference18_k7_robustness_report.ipynb"
output_name="REFERENCE18_K7_ROBUSTNESS_REPORT_${timestamp}.html"

mkdir -p "${output_dir}" "${project_root}/logs"
cd "${project_root}"

"${python_bin}" -m jupyter nbconvert \
  --to html \
  --execute "${notebook}" \
  --output "${output_name}" \
  --output-dir "${output_dir}" \
  --ExecutePreprocessor.timeout=600 \
  --TemplateExporter.exclude_input=True \
  >"${log_path}" 2>&1

report_path="${output_dir}/${output_name}"
test -s "${report_path}"
sha256sum "${report_path}" >"${output_dir}/REPORT_SHA256.txt"
{
  echo "Script: scripts/1148_render_reference18_k7_robustness_report.sh"
  echo "Notebook: scripts/1147_reference18_k7_robustness_report.ipynb"
  echo "Command: ${python_bin} -m jupyter nbconvert --to html --execute ${notebook} --output ${output_name} --output-dir ${output_dir} --ExecutePreprocessor.timeout=600 --TemplateExporter.exclude_input=True"
  echo "Log: ${log_path}"
  echo "Report: ${report_path}"
  echo "Status: PASS"
} >"${output_dir}/RUN_LOG.txt"
printf '%s\n' "${report_path}"
