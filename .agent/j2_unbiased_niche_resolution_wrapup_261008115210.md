# J2 unbiased neighborhood niche resolution wrap-up

Date: 2026-10-08

## Scope

Completed an independent Y40105J2 Bin50-level unbiased neighborhood discovery
using frozen Reference18 continuous weights, UMI >= 50, a 50 µm radius graph,
native VoltRon CLR, and K-means K=3–12 with 20 seeds per K. Frozen lymphoid
aggregate labels were not read, and disease labels did not enter clustering or
K selection.

## Data and methods

- Formal input:
  `results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/J2/01_FROZEN_REFERENCE18_INPUT.tsv.gz`.
- 58,840 eligible Bin50; 58,772 neighborhood-evaluable; 68 empty-neighborhood
  Bin retained as not evaluable.
- Seven patients and seven ROIs: HC=2, IPF=2, SSc-ILD=3.
- Center Bin excluded; 325,776 radius-graph edges; no cross-patient/ROI edges.
- Native CLR retained structural zeros; no pseudocount was added.

## Results

- Primary K=5; lower comparator K=3; higher sensitivity K=7.
- K=5: 20/20 converged, median ARI 0.988, median NMI 0.979, minimum niche
  fraction 15.7%, all niches present in all patients.
- Conservative mixed-niche themes: alveolar-associated, B-enriched,
  airway-associated, lymphatic-endothelial-associated, and
  blood-endothelial-associated.
- Core axes were not single-patient-driven and were directionally retained in
  the prespecified UMI sensitivity checks.

## Outputs and provenance

- Result directory:
  `results/reanalysis/J2_UNBIASED_NICHE_RESOLUTION_20261008_110522/`.
- Standalone report:
  `results/J2_UNBIASED_NICHE_RESOLUTION_REPORT_261008113930.html`.
- Analysis/report commit: `cb6b321`.
- Approved local ELN copy:
  `elab/J2_unbiased_neighborhood_niche_resolution_261008113930.md`.
- No remote eLabFTW entry was created because `ELABFTW_API_KEY` was absent.

## Decisions and next work

- The user approved treating the per-chip L3/J2/upcoming Y40102K8 workflow as
  a recurring workflow. A skill-candidate issue draft is stored in
  `.agent/voltron_unbiased_neighborhood_skill_candidate_261008115210.md`.
- The approved recurring-workflow candidate was filed as GitHub issue #1:
  `https://github.com/bozhu001/stereo_seq/issues/1`.
- Next: J2 primary K=5 simplified spatial architecture audit, then freeze J2,
  then independently run Y40102K8 K=3–12 discovery.
