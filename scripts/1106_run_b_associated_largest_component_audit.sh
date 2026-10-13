#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
timestamp="$(date +%y%m%d%H%M%S)"
log_path="${project_root}/logs/1106_b_associated_largest_component_audit_${timestamp}.log"

cd "${project_root}"
/home/bozhu/miniconda3/envs/banksy_stereo/bin/python \
  scripts/1106_audit_b_associated_largest_component.py \
  >"${log_path}" 2>&1

tail -n 1 "${log_path}"
printf '%s\n' "${log_path}"
