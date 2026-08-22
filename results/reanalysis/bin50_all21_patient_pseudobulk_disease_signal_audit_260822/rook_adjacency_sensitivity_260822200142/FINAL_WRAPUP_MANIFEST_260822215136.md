# Final wrap-up manifest — 260822215136

## Reproducible scripts

- `scripts/380_rook_adjacency_sensitivity_qc.py`
- `scripts/381_finalize_rook_adjacency_sensitivity.py`
- `scripts/382_rook_adjacency_sensitivity_report.qmd`
- `scripts/383_render_rook_adjacency_report.ps1`

Scripts 380/381 and this rook review were classified by the user as a one-off
analysis audit. They are retained for reproduction; no skill-candidate issue
was created.

## Completed analysis output

- Output directory:
  `results/reanalysis/bin50_all21_patient_pseudobulk_disease_signal_audit_260822/rook_adjacency_sensitivity_260822200142/`
- Run commands: `run_commands.txt`
- Initial completion marker: `SUCCESS`
- Finalization marker: `FINALIZE_SUCCESS`
- Output validation: `output_validation.json`
- Final QC: `06_final_qc_status.json` and `06_final_qc_status.tsv`
- Concise saved analysis report: `20_ROOK_ADJACENCY_SENSITIVITY_REPORT.md`

## Rendered report and QA

- Final standalone HTML:
  `rook_adjacency_sensitivity_report_260822212549.html`
- Size: 20,871,077 bytes
- SHA-256: `BEAA483BA274A1316EEFC81B8BC0AE0EA7EB05886EFBF756B6EE9B6DBCDE984E`
- Render log:
  `../logs/rook_adjacency_sensitivity_report_260822212549.log`
- Visual QA log:
  `visual_qa_260822212549/visual_qa_log.md`
- Primary QA screenshot:
  `visual_qa_260822212549/01_top.png`

Failed render logs and Quarto cleanup residues were retained under
`../logs/failed_render_artifacts*/` for auditability. They are not final
reports.

## ELN and session records

- Approved local ELN:
  `elab/Bin50_frozen_program_rook_adjacency_sensitivity_audit_260822212549.md`
- Project memory: `.agent/memory.md`
- Session summary:
  `.agent/rook_adjacency_sensitivity_wrapup_260822215136.md`

The ELN project tag remains `TODO_PROJECT_TAG`. No input entry IDs were linked,
no remote eLabFTW entry was created, and no Git push was performed.

## Git provenance

- Analysis/report local commit:
  `5ac59c94254242bb4fc667e2fa9e36c7d0e2316c`
- Local wrap-up commit: the commit containing this manifest, the ELN copy,
  project memory update, and session summary; see the final `git log` record.
