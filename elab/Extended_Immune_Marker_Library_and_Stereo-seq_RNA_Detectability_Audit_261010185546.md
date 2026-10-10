# Extended Immune Marker Library and Stereo-seq RNA Detectability Audit — 261010185546

- Category: `Bioinformatic`
- Tags: `Stereo-seq`, `Bin50`, `immune-marker`, `lung-atlas`, `RNA-detectability`, `QC`
- Approval: user approved the ELN draft on 2026-10-10
- Remote eLabFTW ID/URL: unavailable; `ELABFTW_API_KEY` was not provisioned in the agent environment
- Report attachment for manual upload: `results/reanalysis/EXTENDED_IMMUNE_MARKER_LIBRARY_20261010/reports/EXTENDED_IMMUNE_MARKER_LIBRARY_REPORT_261010185546.html`

## Goal

Construct an extended immune marker library from the BI_PF/ILD lung single-cell atlas and audit marker availability and raw RNA detection across all qualified Bin50 from 21 Stereo-seq lung samples. Compare the previous small marker panels with atlas-derived core and extended panels without rerunning RCTD, clustering, neighborhood scoring, or disease testing.

## Input

1. Primary atlas marker workbook: `results/reanalysis/bin50_input/BI_PF_ILD_atlas_marker_gene_table_260709213053.xlsx`; SHA256 `b6c17a25dd9da4f73b108f3cf3b964397e787b6c08ba97208295f6c38f9c81b8`.
2. Checksum-identical atlas verification copy: `results/visium way/06_hvg_audit/BI_PF_ILD_atlas_marker_gene_table_260709213053.xlsx`.
3. QCed sparse raw-count matrix: `results/reanalysis/bin50_input/BIN50_joint_reanalysis_tissue_raw_counts_common_genes.h5ad`; 178,044 Bin50 by 27,618 genes, sparse `uint32` counts.
4. Ensembl-symbol mapping: `results/reanalysis/bin50_representation_benchmark/ensembl_gene_symbol_mapping.csv`.
5. Previous small-panel source: `scripts/1049_prepare_immune_subniche_validation.R`.

No existing local eLabFTW copy contained these exact inputs, so no prior entry ID was available for API linking. The optional `BI_PF_ILD_atlas_marker_gene_table_260723131417.xlsx` workbook was not found and was not substituted.

## Script

Primary command:

```text
python results/reanalysis/EXTENDED_IMMUNE_MARKER_LIBRARY_20261010/scripts/01_extended_immune_marker_detectability.py --project-root "D:\Lung ST\codex agent\stereo_seq\stereo_seq" --output-dir "D:\Lung ST\codex agent\stereo_seq\stereo_seq\results\reanalysis\EXTENDED_IMMUNE_MARKER_LIBRARY_20261010" --chunk-size 4096
```

Supporting preflight, validation, report-generation, and Quarto-render scripts are stored in `results/reanalysis/EXTENDED_IMMUNE_MARKER_LIBRARY_20261010/scripts/`.

Key parameters and constraints:

- Seed `20261010`; sparse chunk size `4096`.
- All qualified whole-tissue Bin50 were included.
- Detection used raw integer counts greater than zero.
- Disease labels and niche enrichment were not used for marker selection.
- No reclustering, RCTD modification, K selection, neighborhood scoring, disease testing, DEG/GSEA, or CellChat was run.
- Atlas `fold_change` remains definition-unknown and was not assumed to be log2FC.

Repository commit: `64a30fee048c59430354f7c7c814832ca0e1389c` (`https://github.com/bozhu001/stereo_seq/commit/64a30fee048c59430354f7c7c814832ca0e1389c`).

## Output summary

- Atlas immune hierarchy: 9 Level 1 and 29 Level 2 categories.
- Library: 1,450 immune marker records and 692 unique genes.
- Candidate panels: 27 strict lineage-core genes and 209 eligible extended non-core genes.
- Stereo-seq: 667 genes were present and detected in at least one patient.
- Extended-panel coverage: B 19/22, Plasma 7/8, T 25/27, NK 6/7, and Macrophage 68/72.
- Across eight lineage panels, extended panels added 193 detectable lineage-panel gene entries relative to the previous small panels.
- Increased detection coverage does not establish increased immune-cell identification accuracy.
- All 21 patients and L3, J2, and Y40102K8 were included.
- The final scientific audit and all 25 automated output checks passed.

Output directory: `results/reanalysis/EXTENDED_IMMUNE_MARKER_LIBRARY_20261010/`.

Standalone HTML report: `results/reanalysis/EXTENDED_IMMUNE_MARKER_LIBRARY_20261010/reports/EXTENDED_IMMUNE_MARKER_LIBRARY_REPORT_261010185546.html`; 1,787,761 bytes; SHA256 `367FBA538F808A092475A9A29B5A884AFECFFA542013B876214AAB25474BCA5A`.
