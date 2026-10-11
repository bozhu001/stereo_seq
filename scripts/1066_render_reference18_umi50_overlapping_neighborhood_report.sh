#!/usr/bin/env bash
set -euo pipefail

project_win='D:\Lung ST\codex agent\stereo_seq\stereo_seq'
result_win="${project_win}\\results\\reanalysis\\REFERENCE18_UMI50_OVERLAPPING_COMPOSITION_NEIGHBORHOODS_20261011_182142_attempt2"
python_exe="/mnt/d/Lung ST/codex agent/stereo_seq/stereo_seq/envs/stereo_seq_py/Scripts/python.exe"
script_win="${project_win}\\scripts\\1065_render_reference18_umi50_overlapping_neighborhood_report.py"
log_dir='/mnt/d/Lung ST/codex agent/stereo_seq/stereo_seq/logs'
timestamp="$(date +'%y%m%d%H%M%S')"
mkdir -p "${log_dir}"
"${python_exe}" "${script_win}" --result-dir "${result_win}" \
  > "${log_dir}/reference18_umi50_overlapping_neighborhood_report_${timestamp}.log" 2>&1
