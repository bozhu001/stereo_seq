# Session summary: L3 K8 spatial architecture wrap-up

Date: 2026-10-08

## Completed analysis

- Reused the frozen L3 K=8 unbiased-neighborhood labels and existing
  Reference18/coordinate inputs; no clustering or upstream model was rerun.
- Quantified all eight niches symmetrically using patient/ROI-specific Queen
  graphs and assessed all 28 heterotypic niche pairs.
- Generated continuity/component summaries, direct adjacency tables,
  patient/ROI-stratified 999-permutation observed/expected enrichment, and
  cross-patient consistency summaries.
- Generated and visually checked the requested English PDFs and tables in
  `results/reanalysis/L3_UNBIASED_NICHE_RESOLUTION_20261007_085730/K08_SPATIAL_ARCHITECTURE_20261007_121357/`.

## Main findings and cautions

- The analysis covered 39,982 evaluable Bin50s, seven patients, seven ROIs,
  and 149,193 unique Queen-adjacency edges.
- N7, N4, and N5 were the most spatially continuous niches.
- N3 formed spatial components but was relatively fragmented; its median
  largest-component fraction was 0.065.
- N5–N6 showed recurrent positive adjacency enrichment in 6/7 patients.
- Several negative interface patterns were directionally consistent across
  all seven patients, including N4–N7, N2–N7, and N7–N8.
- Because K=8 labels originate from overlapping 50 µm neighborhood
  compositions and the null disrupts spatial autocorrelation, adjacency is
  interpreted only as spatial association/interface preference, not cellular
  communication or independent biological validation.

## Reproducibility and archival state

- Analysis script: `scripts/1011_l3_k8_spatial_architecture.R`.
- Report source: `scripts/1012_l3_k8_spatial_architecture_report.qmd`.
- Rendered report:
  `results/L3_K8_SPATIAL_ARCHITECTURE_REPORT_261008104540.html`.
- Analysis/report commit: `f0d0d1a`.
- User approved the ELN draft on 2026-10-08.
- Approved local ELN copy:
  `elab/L3_K8_spatial_niche_architecture_261008104540.md`.
- No remote eLabFTW ID/URL exists because `ELABFTW_API_KEY` is unavailable in
  the current environment. The approved entry text and attachment path are
  preserved for manual creation in the web UI.

## Scope exclusions

- Frozen lymphoid aggregates were not used.
- Disease labels were not used in spatial definitions or significance tests.
- No DEG, GSEA, pathway, CellChat, or ligand-receptor analysis was run.
