# Session summary: Extended Immune Marker Library and RNA Detectability Audit

Date: 2026-10-10

## Scope

Built an atlas-derived immune marker library and audited whole-tissue Bin50 raw-count detectability across the fixed 21-patient Stereo-seq cohort. The session did not rerun RCTD, clustering, neighborhood scoring, disease testing, DEG/GSEA, CellChat, or frozen aggregate analyses.

## Inputs and recovery

- Main atlas: `BI_PF_ILD_atlas_marker_gene_table_260709213053.xlsx`.
- The XLSX contained valid worksheet data but a stale `A1` dimension. The workflow reset reader dimensions in memory only; the source workbook remained unchanged with SHA256 `b6c17a25dd9da4f73b108f3cf3b964397e787b6c08ba97208295f6c38f9c81b8`.
- The optional `260723131417.xlsx` atlas was not found and was not substituted.
- Raw expression input: 178,044 Bin50 by 27,618 genes, sparse `uint32`, covering 21 patients and L3/J2/Y40102K8.
- The prior small panels were traced to `scripts/1049_prepare_immune_subniche_validation.R` and added to `input_manifest.tsv` with SHA256 provenance.

## Main results

- Atlas immune hierarchy: 9 Level 1 and 29 Level 2 categories.
- 1,450 immune marker records; 692 unique genes.
- 27 strict lineage-core genes and 209 eligible extended non-core genes.
- 667 genes were present and detected at least once in the cohort.
- Extended coverage: B 19/22, Plasma 7/8, T 25/27, NK 6/7, Macrophage 68/72.
- Extended panels added 193 detectable lineage-panel gene entries relative to the previous small panels.
- The central interpretation is detection coverage improvement, not demonstrated improvement in immune-cell identification accuracy.

## Outputs and validation

- Final directory: `results/reanalysis/EXTENDED_IMMUNE_MARKER_LIBRARY_20261010/`.
- Ten requested TSV tables, five English PDFs, Chinese findings, README, manifests, reproducibility scripts, and audits were created.
- Twenty-five automated checks passed; all five PDFs were rendered and visually reviewed.
- Final standalone HTML: `reports/EXTENDED_IMMUNE_MARKER_LIBRARY_REPORT_261010185546.html`; SHA256 `367FBA538F808A092475A9A29B5A884AFECFFA542013B876214AAB25474BCA5A`.
- Analysis/report commit: `64a30fee048c59430354f7c7c814832ca0e1389c`.

## ELN and workflow review

- The user approved the scientific analysis and, separately, the eLabFTW draft.
- Approved local ELN copy: `elab/Extended_Immune_Marker_Library_and_Stereo-seq_RNA_Detectability_Audit_261010185546.md`.
- No remote eLabFTW ID/URL was created because `ELABFTW_API_KEY` was unavailable; the standalone HTML is ready for manual attachment.
- Project tag: `Stereo-seq`.
- Script-review classification remains unanswered: no shared-skill candidate or GitHub issue was filed.
