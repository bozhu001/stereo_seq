#!/usr/bin/env bash
set -euo pipefail

# Exact render entry point for scripts/1056_b_lineage_subniche_mapping_report.ipynb.
# Run from the project root in tmux:
# tmux new-session -d -s b_subniche_report 'bash scripts/1056_render_b_lineage_subniche_mapping_report.sh'

project_root="$(pwd)"
python_bin="/home/bozhu/miniconda3/envs/stereo/bin/python"
result_dir_file="$project_root/.agent/current_b_lineage_subniche_mapping_out.txt"
result_dir_relative="$(sed '1s/^\xEF\xBB\xBF//' "$result_dir_file" | tr -d '\r\n')"
result_dir="$project_root/$result_dir_relative"
tstamp="$(date +'%y%m%d%H%M%S')"
output_name="B_LINEAGE_SUBNICHE_MAPPING_REPORT_${tstamp}.html"
log_path="$project_root/logs/1056_b_lineage_subniche_mapping_report_${tstamp}.log"

mkdir -p "$project_root/results" "$project_root/logs"

B_SUBNICHE_RESULT_DIR="$result_dir" \
B_SUBNICHE_PROJECT_ROOT="$project_root" \
"$python_bin" -m nbconvert \
  --to html \
  --execute "$project_root/scripts/1056_b_lineage_subniche_mapping_report.ipynb" \
  --ExecutePreprocessor.timeout=600 \
  --HTMLExporter.embed_images=True \
  --output "$output_name" \
  --output-dir "$project_root/results" \
  >"$log_path" 2>&1

printf 'REPORT_PATH=%s\n' "$project_root/results/$output_name"
printf 'LOG_PATH=%s\n' "$log_path"
