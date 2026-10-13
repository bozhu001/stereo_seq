# Session summary — non-B external microenvironment atlas

Date: 2026-10-13

## Decisions and frozen scope

- The study uses 21 patients on K8, J2 and L3, with seven patients per chip.
- Frozen six-dataset Reference18 RCTD proportions and frozen per-chip K=6
  center-excluded 50 µm Neighborhood NMF results were reused.
- RCTD/NMF were not rerun; the 238 frozen lymphoid aggregates were untouched.
- B-associated analyses remain paused. The external atlas contains only the
  Ciliated, Fibroblast and Lymphatic factor families.
- Disease labels were not used in the external-atlas discovery or testing.

## Preceding factor-metric audits in this session

- The B-associated largest-component audit established that the original
  `largest_high_component_fraction` is the largest high-usage component divided
  by all high-usage Bin50s. The K8 IPF outlier combined 91.5% high-usage
  coverage with 99.6% concentration, so the metric primarily represented
  concentration rather than absolute occupied tissue area. Commit: `5f74a3d`.
- The unified four-family threshold sensitivity analysis compared absolute
  component size, component/valid area and component/high-usage concentration
  under chip q85/q90/q95 and patient-specific q90 thresholds. No association
  passed global FDR across 144 exploratory tests. Patient-specific q90
  concentration was retained as the consistent spatial variable, with absolute
  size and tissue-normalized size as required sensitivities. Commit: `11b6c45`.

## External microenvironment analysis

- Script: `scripts/1108_build_non_b_external_microenvironment_atlas.py`.
- Final directory:
  `results/reanalysis/NON_B_EXTERNAL_MICROENVIRONMENT_ATLAS_20261013_140512/`.
- Main annulus: 100–200 µm, same patient/ROI/tissue component only.
- Primary spatial aggregation: 200 µm blocks; 150 and 250 µm were sensitivity
  scales. Technical adjustment included nUMI, nGene, the reconstructed 0–50 µm
  neighbor count and annulus Bin50 count.
- All 21 patients entered the primary analysis. Patient-level valid-center
  coverage ranged from 0.919 to 1.000 (median 0.994).
- Nine of 54 relationships met the cross-chip reproducible rule and seven were
  partially reproducible. The dominant positive findings were same-state
  continuity for Ciliated, Fibroblast and Lymphatic factors. Ciliated–NK was a
  weaker FDR-significant non-self positive candidate; Ciliated–Secretory and
  Fibroblast–AT2/transitional were partial positive candidates.
- Negative relationships remain relative compositional depletion patterns
  because Reference18 proportions are closed; they are not interpreted as
  biological exclusion.
- The 100–200 µm absolute q90-minus-middle differences were small despite
  reproducible rank correlations. The analysis therefore supports weak,
  repeated spatial gradients rather than a strong multicellular annular niche.

## Reporting and provenance

- Rendered HTML:
  `results/1108_non_b_external_microenvironment_atlas_report_261013141630.html`.
- HTML SHA-256:
  `72251F97E4FABCC87D5E7B9AD2CA9B7297BE6D92DCE7F7413F98EB75726F3FD7`.
- Analysis/report commit: `fbf009a`.
- User approved the scientific analysis and separately approved the ELN draft.
- Approved local ELN:
  `elab/Disease-blind_external_Reference18_microenvironment_atlas_261013141630.md`.
- No remote eLabFTW ID/URL exists because `ELABFTW_API_KEY` is unavailable.
- Script review: this was a bespoke project-specific workflow. Per the user's
  prior instruction, it was not proposed as a reusable skill and no GitHub
  skill-candidate issue was filed.

## Agreed next analysis

The next biological analysis should test disease-associated remodeling inside
the frozen Ciliated, Fibroblast and Lymphatic factor-defined regions. Factor
loadings themselves are fixed per chip and are not patient-specific outcomes.
The primary next layer should therefore be patient × factor weighted
pseudobulk from raw expression, followed by HC/IPF/SSc-ILD comparisons with
chip adjustment and within-chip label permutation. External 100–200 µm
composition-by-disease interactions are secondary. B remains out of scope.
