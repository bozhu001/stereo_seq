#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
windows_script="$(wslpath -w "$script_dir/render_cellbin_broadclass_pilot_report.ps1")"
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "$windows_script"
