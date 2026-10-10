# Session summary — B-lineage coverage to immune-subniche mapping

Date: 2026-10-10

## Objective

Map the frozen final B-lineage RNA coverage mask onto the frozen K=3/4/5
immune-focused subniche assignments for L3, J2, and Y40102K8. The purpose was
to identify cross-resolution B-enriched branches before expanding the same
coverage framework to other immune lineages.

## Frozen inputs

- Final B-lineage mask:
  `results/reanalysis/B_LINEAGE_CALIBRATED_COVERAGE_MASK_GE5_ONLY_20261010_220617/01_FINAL_CALIBRATED_B_COVERAGE_MASK.tsv.gz`.
- K=3/4/5 assignments:
  `results/reanalysis/EXTENDED_IMMUNE_SIGNATURE_VALIDATION_20261010_192044/work/01_FROZEN_K345_ASSIGNMENTS.tsv.gz`.
- Expected candidate counts were reproduced exactly: L3=14,839, J2=11,893,
  Y40102K8=2,913; total=29,645.
- Expected primary B-mask counts were reproduced exactly: L3=1,386, J2=545,
  Y40102K8=738; total=2,669.

## Method

- Joined every candidate Bin50 to its K=3, K=4, and K=5 assignment.
- Summarized B-mask coverage for each chip × patient × K × cluster.
- Compared each cluster's patient-level B fraction with the remaining clusters
  in the same chip, K, and patient.
- Used patient-direction consistency and median patient enrichment ratios;
  overlapping Bin50 neighborhoods were not treated as independent biological
  replicates.
- Audited core-marker corroboration, extended-only gain, patient contribution,
  connected-component contribution, and K3→K4→K5 transitions.
- No RCTD, graph construction, clustering, K selection, marker calibration, or
  disease comparison was rerun.

## Results

- L3: K3 I1 → K4 I4 → K5 I3. The K3 branch captured 1,032/1,386
  B-mask bins (74.5%); 6/7 patients were directionally enriched and the median
  patient enrichment ratio was 3.08.
- J2: K3 I3 → K4 I1 → K5 I3. The K3 branch captured 249/545 B-mask
  bins (45.7%); 6/7 patients were directionally enriched and the median
  patient enrichment ratio was 1.91.
- Y40102K8: K3 I1 → K4 I4 → K5 I3. The K3 branch captured 734/738
  B-mask bins (99.5%); 5/7 patients were directionally enriched, but 92.2% of
  branch B bins came from the largest contributing patient.
- Eight chip/K clusters met the descriptive recurrent-enrichment rule. All
  retained a patient-imbalance warning; this does not erase cross-patient
  direction consistency but prevents Bin-level pseudo-replication.
- These are mapping labels, not final immune-subniche names. Joint lineage
  annotation should wait for T, NK, Plasma, Myeloid, and Mast/Basophil audits.

## Outputs and provenance

- Final directory:
  `results/reanalysis/B_LINEAGE_SUBNICHE_MAPPING_20261010_225151/`.
- HTML report:
  `results/B_LINEAGE_SUBNICHE_MAPPING_REPORT_261010225750.html`.
- Report SHA256:
  `F94F68D53E8D85095306B32083C694F96B7C291EDE5C24300BC106232D3558D1`.
- Analysis/report commit: `6e94986`.
- User approved the analysis and separately approved the ELN draft.
- Approved local ELN:
  `elab/B-lineage_coverage_mapping_to_immune_subniches_261010225750.md`.
- No remote eLabFTW ID/URL exists because `ELABFTW_API_KEY` was unavailable.
- Script-review recurring-versus-one-off classification remains pending.

## Recommended next analysis

Use the B-lineage workflow as a validated pilot, then build lineage-specific
coverage rules for Plasma, T, NK, inflammatory myeloid, resident macrophage,
and Mast/Basophil. Do not mechanically copy the B q90/≥2 genes/≥5 connected
bins thresholds; calibrate marker specificity, matched background, and spatial
support separately for each lineage. After mapping all lineages to K3/K4/K5,
freeze joint immune-subniche names and only then begin patient-level disease
comparisons.
