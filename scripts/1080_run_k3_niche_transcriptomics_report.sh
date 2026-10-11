#!/usr/bin/env bash
set -euo pipefail

project="/mnt/d/Lung ST/codex agent/stereo_seq/stereo_seq"
windows_project="D:/Lung ST/codex agent/stereo_seq/stereo_seq"
timestamp="$(date +'%y%m%d%H%M%S')"
log_path="${project}/logs/k3_niche_transcriptomics_report_${timestamp}.log"

mkdir -p "${project}/logs"
"${project}/envs/stereo_seq_py/Scripts/python.exe" \
  "${windows_project}/scripts/1080_render_k3_niche_transcriptomics_report.py" \
  --analysis-dir "${windows_project}/results/reanalysis/ALL21_K3_NICHE_PSEUDOBULK_DEG_20261011_201008" \
  --output-dir "${windows_project}/results" \
  >"${log_path}" 2>&1

printf 'RENDER_DONE\n' >>"${log_path}"
