#!/usr/bin/env bash
set -euo pipefail

# Exact render command for scripts/1054_b_lineage_coverage_report.ipynb.
# Run from the project root, preferably through:
# tmux new-session -d -s b_coverage_report 'bash scripts/1054_render_b_lineage_coverage_report.sh'

project_root="$(pwd)"
python_bin="/home/bozhu/miniconda3/envs/stereo/bin/python"
result_dir="$project_root/results/reanalysis/B_LINEAGE_CALIBRATED_COVERAGE_MASK_GE5_20261010_215212"
rule_audit_dir="$project_root/results/reanalysis/B_LINEAGE_BIN_COVERAGE_RULE_AUDIT_20261010_212159"
tstamp="$(date +'%y%m%d%H%M%S')"
output_name="B_LINEAGE_CALIBRATED_COVERAGE_REPORT_${tstamp}.html"
log_path="$project_root/logs/1054_b_lineage_coverage_report_${tstamp}.log"

mkdir -p "$project_root/results" "$project_root/logs"

B_COVERAGE_RESULT_DIR="$result_dir" \
B_COVERAGE_RULE_AUDIT_DIR="$rule_audit_dir" \
B_COVERAGE_PROJECT_ROOT="$project_root" \
"$python_bin" -m nbconvert \
  --to html \
  --execute "$project_root/scripts/1054_b_lineage_coverage_report.ipynb" \
  --ExecutePreprocessor.timeout=600 \
  --HTMLExporter.embed_images=True \
  --output "$output_name" \
  --output-dir "$project_root/results" \
  >"$log_path" 2>&1

printf 'REPORT_PATH=%s\n' "$project_root/results/$output_name"
printf 'LOG_PATH=%s\n' "$log_path"
