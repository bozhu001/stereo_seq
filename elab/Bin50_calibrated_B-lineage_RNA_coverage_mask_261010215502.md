# Bin50 calibrated B-lineage RNA coverage mask — 261010215502

- Category: `Bioinformatic`
- Tags: `Stereo-seq`, `Bin50`, `B-lineage`, `extended-marker`, `RNA-coverage`, `spatial-QC`
- Approval: user approved the ELN draft on 2026-10-10
- Remote eLabFTW ID/URL: unavailable; `ELABFTW_API_KEY` was not provisioned in the agent environment
- Report attachment for manual upload: `results/B_LINEAGE_CALIBRATED_COVERAGE_REPORT_261010215502.html`

## Goal

Construct a Bin50 B-lineage RNA coverage mask that increases sensitivity in
low-UMI spatial data with an extended B-marker panel while limiting calls driven
by a single broad gene or isolated spatial noise. The endpoint is B-lineage RNA
coverage, not a pure B-cell identity, B-cell niche, aggregate, or TLS call.

## Input

1. Frozen 29,645-bin immune candidate pool and K=3/4/5 assignments:
   `results/reanalysis/EXTENDED_IMMUNE_SIGNATURE_VALIDATION_20261010_192044/work/01_FROZEN_K345_ASSIGNMENTS.tsv.gz`.
2. Raw-integer-count-derived center and center-excluded 50 µm neighborhood
   marker scores:
   `results/reanalysis/EXTENDED_IMMUNE_SIGNATURE_VALIDATION_20261010_192044/work/06_CENTER_NEIGHBOR_MARKER_SCORES.tsv.gz`.
3. Frozen extended-signature cluster validation:
   `results/reanalysis/EXTENDED_IMMUNE_SIGNATURE_VALIDATION_20261010_192044/04_CLUSTER_LEVEL_EXTENDED_SIGNATURE_VALIDATION.tsv`.
4. Frozen 50 µm candidate-neighbor edges:
   `results/reanalysis/EXTENDED_IMMUNE_SIGNATURE_VALIDATION_20261010_192044/work/04_EXACT_RADIUS50_CANDIDATE_NEIGHBOR_EDGES.tsv.gz`.
5. Upstream local ELN copy:
   `elab/Extended_Immune_Marker_Library_and_Stereo-seq_RNA_Detectability_Audit_261010185546.md`.
   No remote entry ID exists, so an API provenance link could not be created.

## Script

The documented workflow consists of:

- `scripts/1051_binlevel_b_coverage_pilot.py` and
  `scripts/1051_run_binlevel_b_coverage_pilot.ps1`;
- `scripts/1052_audit_binlevel_b_coverage_rule.py` and
  `scripts/1052_run_binlevel_b_coverage_rule_audit.ps1`;
- `scripts/1053_finalize_calibrated_b_coverage_mask.py` and
  `scripts/1053_run_finalize_calibrated_b_coverage_mask.ps1`;
- `scripts/1054_b_lineage_coverage_report.ipynb` and
  `scripts/1054_render_b_lineage_coverage_report.sh`.

Key parameters and constraints:

- patient/UMI-matched background q90;
- at least two informative extended B genes detected per center-excluded
  50 µm neighborhood;
- primary spatial criterion: at least five connected positive Bin50;
- at least three connected bins retained as a sensitivity result;
- core B genes provide corroboration but do not veto extended coverage;
- disease labels were not used;
- RCTD, graph construction, clustering, and K selection were not rerun.

Analysis/report commit:
`b857c3c223288b693f530dea5fc2eb314f41adb2`
(`https://github.com/bozhu001/stereo_seq/commit/b857c3c223288b693f530dea5fc2eb314f41adb2`).

## Output summary

- Candidate Bin50: 29,645.
- q90 and multi-gene-supported bins before connectivity: 5,376.
- Primary component ≥5 mask: 2,669 Bin50.
- Component ≥3 sensitivity mask: 3,429 Bin50.
- Primary extended gain beyond the core-only comparator: 1,514 Bin50.
- Primary bins with core-marker corroboration: 1,811.
- Primary component ≥5 coverage by chip: J2 545/11,893 (4.58%), L3
  1,386/14,839 (9.34%), and Y40102K8 738/2,913 (25.33%).
- All automated final-mask audit checks passed.

Final output directory:
`results/reanalysis/B_LINEAGE_CALIBRATED_COVERAGE_MASK_GE5_20261010_215212/`.

Final mask:
`results/reanalysis/B_LINEAGE_CALIBRATED_COVERAGE_MASK_GE5_20261010_215212/01_FINAL_CALIBRATED_B_COVERAGE_MASK.tsv.gz`;
SHA256 `F7F9A83A00F74E79AEFB547B0891CA9B0574E9BCBDE86CE9CE2975302CAA2658`.

Rendered HTML:
`results/B_LINEAGE_CALIBRATED_COVERAGE_REPORT_261010215502.html`;
2,158,828 bytes; SHA256
`5BD96C947CA7174627FAE06F4148BC8E26E1E0EF3BF94539DEA1CDD6E2771A9D`.
