# Bespoke script review — one-off synthesis

## Classification

The user confirmed that this B-lineage synthesis and the matched spatial-null
workflow are one-off audits. No reusable skill and no GitHub issue will be
created.

## What the QMD does

`399_bin50_blineage_analysis_line_synthesis_report.qmd` reads already saved
TSV, JSON and PNG outputs from six completed B-lineage analysis stages and
assembles a single interpretation-focused HTML report. It does not source or
run any upstream analysis script.

## Inputs

- Project-existing 219-marker capture QC.
- Frozen 49-gene, six-panel QC in the independent 16-sample cohort.
- FO22 extreme-depth 16-versus-15 sensitivity.
- Primary15 and J2/L3-only13 patient-level disease comparison.
- Six-sample L3 focal review.
- L3 depth/n_genes/PTPRC-status matched spatial null.

MS4A4A is explicitly excluded.

## Outputs

- A self-contained DRAFT synthesis HTML.
- Source-analysis and evidence-interpretation manifests.
- HTML validation, render log and local ELN draft.

## Key fixed interpretation parameters

- Strict same-bin is the primary spatial evidence.
- Strict rook clustering is construction-dependent auxiliary evidence.
- Disease comparison is negative: 18/18 confidence intervals cross zero and
  0/18 tests reach BH-FDR below 0.05.
- Two focal results remain patient-specific candidates only.
- Niche, aggregate, TLS and disease-specific language is prohibited.

