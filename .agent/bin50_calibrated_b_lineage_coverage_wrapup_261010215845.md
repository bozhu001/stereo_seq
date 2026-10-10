# Session summary — calibrated Bin50 B-lineage RNA coverage

Date: 2026-10-10

## Objective

Convert the previously validated extended B-marker library into a spatial
Bin50 B-lineage RNA coverage endpoint. The user's intent was coverage rescue in
low-UMI tissues, without interpreting broad single-gene detection as B-cell
identity.

## Starting point

- Frozen Strategy C immune candidate pool: 29,645 Bin50 from 21 patients on
  J2, L3, and Y40102K8.
- Frozen K=3/4/5 immune assignments.
- Raw-integer-count-derived center and center-excluded 50 µm neighborhood
  marker scores.
- Previously validated core and extended-only immune marker panels.
- No disease label was used in calibration.

## Analysis sequence

1. `1051_binlevel_b_coverage_pilot.py` tested a naive rule requiring at least
   two detected extended genes. It called 24,480/29,645 candidate bins (82.58%)
   and was rejected as depth- and broad-gene-sensitive.
2. `1052_audit_binlevel_b_coverage_rule.py` constructed disease-blind positive
   anchors from bins assigned to an independently B-supported cluster at
   least two of K=3/4/5. Zero-resolution bins formed the calibration
   background, and one-resolution bins remained ambiguous. Patient/UMI-matched
   comparison retained ten informative extended genes: ARHGAP24, ITSN2,
   MAP4K4, MGAT5, NCOA3, PRDM2, RBM6, SAMD12, USP34, and ZCCHC7.
3. The calibrated rule required an extended-panel score above the matched
   background q90 and detection of at least two informative genes in the
   center-excluded 50 µm neighborhood. q95 was retained as a score-threshold
   sensitivity.
4. `1053_finalize_calibrated_b_coverage_mask.py` imposed spatial connectivity.
   The user approved a primary minimum component size of five positive Bin50;
   component size three was retained as a sensitivity field.

## Frozen interpretation

- Primary endpoint: q90 + at least two informative extended genes + connected
  component of at least five positive Bin50.
- Core B genes corroborate but do not veto extended support.
- The mask represents B-lineage RNA coverage only, not cell identity, a niche,
  an aggregate, or a TLS.
- RCTD, graph construction, clustering, K selection, and disease testing were
  not rerun.

## Results

- q90 multi-gene-supported before connectivity: 5,376 Bin50.
- Component ≥3 sensitivity: 3,429 Bin50.
- Component ≥5 primary mask: 2,669 Bin50.
- Primary extended gain beyond core-only: 1,514 Bin50.
- Primary core-corroborated: 1,811 Bin50.
- By chip, the primary mask retained J2 545, L3 1,386, and Y40102K8 738 bins.
- All final automated audits passed. Spatial inspection showed connected
  patches rather than uniform whole-chip staining. The large FO22 region was
  stable between q90 and q95 and strongly core-corroborated.

## Outputs and provenance

- Final directory:
  `results/reanalysis/B_LINEAGE_CALIBRATED_COVERAGE_MASK_GE5_20261010_215212/`.
- HTML report:
  `results/B_LINEAGE_CALIBRATED_COVERAGE_REPORT_261010215502.html`.
- Analysis/report commit:
  `b857c3c223288b693f530dea5fc2eb314f41adb2`.
- Approved local ELN copy:
  `elab/Bin50_calibrated_B-lineage_RNA_coverage_mask_261010215502.md`.
- No remote ELN ID/URL is available because `ELABFTW_API_KEY` was not
  provisioned.
- Script-review classification remains pending user choice between recurring
  workflow and one-off analysis; no skill-candidate issue was filed.
