# Session summary — L3 Mature_B_core matched spatial null wrap-up

Date: 2026-08-23 (Europe/Berlin)

## Scope and immutable provenance

This analysis reused the persisted Bin-level support flags from
`results/reanalysis/bin50_l3_mature_b_focal_spatial_review_260823202319/`.
It did not read the raw G15 object, reconstruct a score, change a marker or
threshold, or select a sample from B-lineage expression.

The fixed L3 samples were `HC/NL-72`, `IPF/FO23-1-06168`,
`IPF/FO23-1-06170`, `SSC/05957/17B`, `SSC/07998/15A` and
`SSC/15491/14`. The prespecified extremely low-capture sample
`SSC/05957/17B` remained included.

## Matched spatial null

Within each sample, tissue Bins were assigned to fixed joint strata from the
empirical quintile of `log1p(total_counts)`, empirical quintile of `n_genes`,
and PTPRC raw-positive/raw-zero status. For same-bin and rook support
independently, 1,000 without-replacement draws exactly reproduced the observed
stratum counts. Seed `260824` was fixed. All 12 sample/support-mode
combinations passed exact matching validation.

Each observed and random set was summarized by largest rook-connected
component size, largest-component support fraction, component number and
median nearest-neighbor distance. A uniform same-count spatial null was also
generated for like-for-like comparison, while the previous unmatched
depth-overlap table was carried forward unchanged.

## Approved interpretation

Strict same-bin is the primary spatial evidence. In both
`IPF/FO23-1-06168` and `SSC/15491/14`, all four same-bin metrics were at the
clustering-direction 100th empirical percentile of the matched null. Their
signals are not fully explained by the matched technical/immune background.

The two results are retained as
`depth/n_genes/PTPRC-status–matched focal multi-marker B-lineage candidates`:
one patient-specific IPF candidate and one patient-specific SSc-ILD candidate.

Strict rook support uses adjacency in its construction. Its clustering is
construction-dependent auxiliary evidence and not an independent spatial
validation. The candidates are not B-cell niches, aggregates, TLS or
disease-specific structures. No disease inference is made; the technical
interpretation remains limited to “not fully explained by the matched
technical/immune background.”

## Reporting and validation

The final report changed interpretation text only and read existing tables and
figures. It did not rerun the randomization or any B-cell statistic.

Final HTML:

`results/reanalysis/bin50_l3_mature_b_depth_ptprc_matched_spatial_null_260823211809/report/FINAL_398_l3_mature_b_depth_ptprc_matched_spatial_null_report_260823211809.html`

Size: 2,932,894 bytes.

SHA-256:
`D32E7F4E5E854794F38B227D5D393767A7C6AE966FA6854242385401CFFC90B4`.

Validation found one closing HTML tag, four embedded PNG images, 11 rendered
HTML tables, no escaped kable strings, all approved interpretation phrases and
no prohibited depth-independence wording.

Finalization manifest:
`results/reanalysis/bin50_l3_mature_b_depth_ptprc_matched_spatial_null_260823211809/analysis_manifest_260823214027.json`.

Analysis/report commit:
`f4a09b9959ad50f88ca4a097fee37d7d391eaa7a`.

## Documentation state

Local ELN draft:
`elab/L3_Mature_B_core_depth_PTPRC_matched_spatial_null_260823211809.md`.
It remains pending the separate ELN approval gate; no remote entry was created
and no report was uploaded.

The script is bespoke and has a local script-review summary. Recurrence has not
been confirmed, so no new-skill-candidate issue was filed.

No push was performed. No further B-cell statistics or spatial analyses were
run during finalization.
