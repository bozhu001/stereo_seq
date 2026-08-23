# Session summary — B-lineage analysis-line synthesis wrap-up

Date: 2026-08-23 (Europe/Berlin)

## Scope

This session assembled a report-only synthesis of six already completed
B-lineage stages: project-existing marker capture QC, the frozen 49-gene and
six-panel descriptive QC in 16 independently retained patients, FO22
extreme-depth sensitivity, the 15/13-patient disease comparison, the L3 focal
review and the matched spatial-null audit. MS4A4A was excluded.

No new statistic, model, score, randomization, spatial test, image analysis,
marker selection or threshold change was run. The QMD read saved TSV, JSON and
PNG outputs only.

## Integrated evidence

- The project-existing union contains 219 genes, of which 202 were present in
  G15 and 17 absent. Zero and absent rows remained represented. Detected marker
  count was depth-associated descriptively; category CPM-depth correlations
  were weak.
- The literature-driven source froze 49 non-exclusive genes across six
  B-lineage panels and separate `TLS_context`. Independent global QC retained
  16 patients (HC 4, IPF 5, SSc-ILD 7).
- FO22 drove many upper technical extremes. Its removal reduced maxima but did
  not eliminate the residual K8-high capture pattern.
- The cleaned disease comparison was negative: 18/18 primary 95% confidence
  intervals crossed zero and 0/18 reached BH-FDR below 0.05.
- `IPF/FO23-1-06168` and `SSC/15491/14` remain patient-specific
  depth/n_genes/PTPRC-status–matched focal multi-marker B-lineage candidates.
  Strict same-bin is primary; strict rook clustering is construction-dependent
  auxiliary evidence.

## Approved minimal revision

A six-row table now transcribes the exact frozen same-bin rules for
`Mature_B_core`, `ILD_plasma_cell`, `GC_activated_B`,
`IgA_mucosal_plasma`, `Naive_B` and `Memory_tissue_B`. The common rook rule
uses direct orthogonal native-coordinate distance 50 in the same authoritative
source-tissue component; both endpoints contribute evidence and their marker
union must satisfy the corresponding frozen panel rule.

The matched-null interpretation now states that a 100th empirical percentile
means no set among the existing 1,000 matched draws reached the observed
spatial concentration. This is not a formal disease-comparison P value and
cannot support disease-group inference.

The H&E wording is conditional: if corresponding H&E is obtained and reliable
tissue-outline or coarse-region registration is possible, candidate
coordinates may be checked against anatomy and pathology.

Both focal regions remain patient-specific and unvalidated. They cannot be
called B-cell niches, aggregates, TLS or disease-specific structures.

## Output and verification

Final HTML:
`results/reanalysis/bin50_blineage_analysis_line_synthesis_260823215304/report/FINAL_399_bin50_blineage_analysis_line_synthesis_report_260823215304.html`.

Size: 13,794,494 bytes.

SHA-256:
`B97384DE5636CCF29F59D607CD4768AF139FB2A281DB9851A1A754CB72FB52F5`.

Validation confirmed 14 embedded existing PNG images, 12 rendered HTML tables,
all six rule rows, the requested percentile and conditional H&E statements,
one complete closing HTML tag and no escaped kable strings.

Analysis/report commit:
`45f2672f14ec7dacce75930273241236e72f4b7b`.

Local ELN draft:
`elab/DRAFT_Bin50_Blineage_analysis_line_synthesis_260823215304.md`.
Separate ELN approval remains pending; no remote entry was created.

The synthesis and matched spatial-null workflow are classified as one-off.
No reusable skill or GitHub issue was created. No push was performed.

