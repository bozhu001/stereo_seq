# Bin50 frozen-program rook adjacency sensitivity audit — 260822212549

## Metadata

- **Title:** Bin50 frozen-program rook adjacency sensitivity audit — 260822212549
- **Category:** Bioinformatic
- **Tags:** `TODO_PROJECT_TAG`, `Bin50`, `spatial-omics`,
  `rook-adjacency`, `sensitivity-analysis`
- **Local approval:** Approved by the user on 2026-08-22.
- **Remote status:** Not created because `ELABFTW_API_KEY` is unavailable and
  the user explicitly requested no remote operation.
- **Remote entry ID/URL:** TODO after credentials and project tag become
  available.

## Goal

Assess whether replacing the historical queen eight-neighborhood graph with a
rook four-neighborhood graph changes the completed Bin50 frozen-program spatial
reproducibility calls or patient-level disease conclusions. Preserve the five
frozen gene panels, saved Bin50 program scores, thresholds, permutation count,
patient-level model specification, and 75-test BH family.

## Input

1. Frozen Bin50 program-score table:
   `results/reanalysis/bin50_all21_patient_pseudobulk_disease_signal_audit_260822/minimal_abundance_coordinate_qc_260822193323/04_all21_frozen_program_scores_with_bin_id.tsv.gz`
   (SHA-256 `0a446bff42ba24901276a7389b11424dbe4ec21a19e05fbadd90071ca6c35dee`).
2. Completed queen-adjacency patient/program spatial metrics:
   `results/reanalysis/bin50_all21_patient_pseudobulk_disease_signal_audit_260822/frozen_leading_edge_spatial_mapping_pilot_260822180837/10_patient_program_spatial_metrics.tsv`
   (SHA-256 `0bf84703ee82dc01208bc3db66c4e2bb66ddf87e3ed59e9b945cc9a81609eb9f`).
3. Completed queen-adjacency patient disease tests:
   `results/reanalysis/bin50_all21_patient_pseudobulk_disease_signal_audit_260822/frozen_leading_edge_spatial_mapping_pilot_260822180837/11_chip_adjusted_patient_disease_tests.tsv`
   (SHA-256 `1844c9aa0bcbc4651baa3daedd9b5a114f76be7dd26041c3b3`).
4. Frozen source-component/topology input:
   `results/reanalysis/rctd_bcell_spatial_niche_analysis_260809035034/prepared_niche_inputs.h5`
   (SHA-256 `75815bbd67980001659afe5f242e7aea409c18e09384bdd9c836250742cc756b`).
5. Completed coordinate-QC status and existing K8/J2/L3 validation figures:
   `results/reanalysis/bin50_all21_patient_pseudobulk_disease_signal_audit_260822/minimal_abundance_coordinate_qc_260822193323/11_final_qc_status.tsv`
   (SHA-256 `4027f4047ed1805669f6fab44318979c1a6da5a646b8f4c56001c382c0d88fde`).

**Prior-entry link audit:** no prior ELN markdown copy or entry ID for these
inputs was found in `elab/` or `.agent/memory.md`. The user confirmed that
there is currently no known prior eLabFTW entry ID, so no input association is
recorded. The remote instance was not queried because `ELABFTW_API_KEY` is not
available.

## Script

Analysis commands previously completed:

```text
./envs/cell2location_pilot/bin/python scripts/380_rook_adjacency_sensitivity_qc.py --timestamp 260822200142
./envs/cell2location_pilot/bin/python scripts/381_finalize_rook_adjacency_sensitivity.py --output results/reanalysis/bin50_all21_patient_pseudobulk_disease_signal_audit_260822/rook_adjacency_sensitivity_260822200142
```

Report-only render command used in this wrap-up:

```text
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/383_render_rook_adjacency_report.ps1
```

Key fixed parameters: seed 20260822; 999 permutations; rook offsets `(0, 50)`
and `(50, 0)` native units; high-score threshold `z >= 1`; minimum hotspot size
4 Bin50; Bin50 area 0.000625 mm2. No score, gene panel, GSEA, or differential
expression result was recomputed.

Repository: <https://github.com/bozhu001/stereo_seq>

Committed code/report snapshot:
<https://github.com/bozhu001/stereo_seq/commit/5ac59c94254242bb4fc667e2fa9e36c7d0e2316c>

## Output summary

All three final QC domains passed: bin-score alignment, coordinate orientation,
and spatial adjacency. Rook versus queen adjacency changed 0/5 formal
program-level spatial-QC conclusions and 0/45 adjacency-related patient disease
FDR conclusions. `B_core` remained the only borderline structural program and
remained failed because coherent multigene support was absent under both graph
definitions. J2 retained the identity orientation with no mirror, rotation,
scale, or affine transform.

Result directory:
`results/reanalysis/bin50_all21_patient_pseudobulk_disease_signal_audit_260822/rook_adjacency_sensitivity_260822200142/`

Rendered report:
`rook_adjacency_sensitivity_report_260822212549.html` (20,871,077 bytes;
SHA-256 `BEAA483BA274A1316EEFC81B8BC0AE0EA7EB05886EFBF756B6EE9B6DBCDE984E`).

The report is below 100 MB. It should be attached if this approved local record
is later created remotely after the project tag and credentials are available.

## Deferred remote fields

- Project-specific eLabFTW tag: `TODO_PROJECT_TAG` (unknown; do not infer).
- Remote entry ID and URL: TODO.
- Input-entry associations: none currently confirmed.
- Remote creation and report attachment: not performed in this session.
