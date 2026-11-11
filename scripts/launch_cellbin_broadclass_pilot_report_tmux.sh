#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
tmux kill-session -t cellbin_broadclass_report 2>/dev/null || true
tmux new-session -d -s cellbin_broadclass_report \
  "bash '$script_dir/render_cellbin_broadclass_pilot_report.sh'"
tmux list-sessions
