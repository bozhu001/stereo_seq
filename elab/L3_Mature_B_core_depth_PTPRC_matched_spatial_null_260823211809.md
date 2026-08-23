# DRAFT — L3 Mature_B_core depth/PTPRC-matched spatial null — 260823211809

## Metadata

- **Title:** L3 Mature_B_core depth/PTPRC-matched spatial null — 260823211809
- **Category:** Bioinformatic
- **Tags:** `TODO_PROJECT_TAG`, `Bin50`, `spatial-omics`, `B-lineage`,
  `Mature_B_core`, `L3`, `matched-spatial-null`, `descriptive-QC`
- **Analysis approval:** The calculation and final interpretation were approved
  by the user on 2026-08-23.
- **ELN status:** Local draft only; this entry has not passed the separate ELN
  approval gate.
- **Remote status:** Not created. No API call, attachment upload, push or other
  remote mutation was performed.
- **Remote entry ID/URL:** TODO after separate ELN-draft approval.

## Goal

Determine whether focal frozen `Mature_B_core` strict same-bin and strict rook
support in six fixed L3 samples shows more spatial clustering than random Bin
sets with the same local `total_counts`, `n_genes` and PTPRC raw-positive/raw-
zero background. The analysis was descriptive spatial QC, not disease
inference.

## Input

The immutable Bin-level input was:

`results/reanalysis/bin50_l3_mature_b_focal_spatial_review_260823202319/tables/10_l3_bin_level_frozen_marker_and_support_values.tsv.gz`

MD5: `950a07b5982c359a7be43a822ec9bcfd`.

The fixed sample manifest, frozen aggregate support table and carried-forward
unmatched depth-overlap table came from the same upstream result directory.
The six samples were `HC/NL-72`, `IPF/FO23-1-06168`,
`IPF/FO23-1-06170`, `SSC/05957/17B`, `SSC/07998/15A` and
`SSC/15491/14`. `SSC/05957/17B` was retained as the prespecified extremely
low-capture sample.

**Prior-entry link audit:** no approved local ELN entry ID for the immediate L3
focal-review input was found in `elab/` or `.agent/memory.md`. No remote
instance query was performed, so no input-entry association is claimed.

## Script

Analysis command completed before approval:

```text
Rscript scripts/397_l3_mature_b_depth_ptprc_matched_spatial_null.R --timestamp=260823211809
```

Final report-only render:

```text
$env:L3_MATCHED_NULL_DIR = (Resolve-Path "results/reanalysis/bin50_l3_mature_b_depth_ptprc_matched_spatial_null_260823211809").Path
quarto render scripts/398_l3_mature_b_depth_ptprc_matched_spatial_null_report.qmd --to html
```

The approved final render read existing TSV and PNG outputs only. It did not
rerun spatial randomization, reconstruct a support flag, rescore a marker or
change a threshold.

Key fixed parameters:

- Six prespecified L3 samples; no expression-based sample selection.
- Frozen `strict_same_bin` and `strict_rook` support flags read unchanged.
- Sample-specific empirical quintiles of `log1p(total_counts)` and `n_genes`,
  crossed with PTPRC raw-positive/raw-zero status.
- Exact without-replacement sampling within each observed matching stratum.
- 1,000 draws per sample and support mode; seed `260824`.
- Spatial metrics: largest rook component Bin count and support fraction,
  connected-component count, and median nearest-neighbor distance.
- Empirical percentiles are descriptive ranks, not P values.

Repository: <https://github.com/bozhu001/stereo_seq>

Local analysis/report commit:
`f4a09b9959ad50f88ca4a097fee37d7d391eaa7a`.

Remote status: not pushed as of 2026-08-23.

## Output summary

All 12 sample-by-support-mode combinations completed 1,000 exact matched draws
with no frozen-support count drift. The final interpretation sets strict
same-bin as the primary spatial evidence.

For both `IPF/FO23-1-06168` and `SSC/15491/14`, all four strict same-bin
clustering metrics were at the clustering-direction 100th empirical percentile
of the matched null. Their focal multi-marker signals are therefore **not fully
explained by the matched technical/immune background**. They are retained as
**depth/n_genes/PTPRC-status–matched focal multi-marker B-lineage candidates**:
one patient-specific IPF candidate and one patient-specific SSc-ILD candidate.

Strict rook support itself uses rook adjacency. Its clustering result is
construction-dependent auxiliary evidence and cannot serve as an independent
spatial validation. Neither candidate may be called a B-cell niche, aggregate,
TLS or disease-specific structure. The analysis does not establish a disease
difference; the approved technical interpretation is limited to “not fully
explained by the matched technical/immune background.”

Result directory:

`results/reanalysis/bin50_l3_mature_b_depth_ptprc_matched_spatial_null_260823211809/`

Final standalone HTML:

`report/FINAL_398_l3_mature_b_depth_ptprc_matched_spatial_null_report_260823211809.html`

Size: 2,932,894 bytes.

SHA-256:
`D32E7F4E5E854794F38B227D5D393767A7C6AE966FA6854242385401CFFC90B4`.

The HTML is below 100 MB and should be attached if this draft is separately
approved for remote eLabFTW creation.

Finalization manifest:
`analysis_manifest_260823214027.json`.

## Deferred remote fields

- Separate user approval of this ELN draft: pending.
- Project-specific eLabFTW tag: `TODO_PROJECT_TAG` (unknown; do not infer).
- Remote entry ID and URL: TODO.
- Remote creation and final-HTML attachment: not performed.
- Git push: not performed.
