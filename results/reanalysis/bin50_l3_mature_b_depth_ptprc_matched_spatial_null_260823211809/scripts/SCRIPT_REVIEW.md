# Bespoke script review

## Purpose

`397_l3_mature_b_depth_ptprc_matched_spatial_null.R` performs a prespecified,
descriptive spatial-null QC analysis for frozen Mature_B_core same-bin and rook
support flags in six fixed L3 samples.

## Immutable inputs

- Persisted Bin-level support values from the approved L3 focal review.
- The fixed six-sample manifest and frozen aggregate support table.
- The prior unmatched depth-overlap summary, carried forward unchanged.

The script does not read the raw G15 object and does not reconstruct, rescore,
or alter any support flag.

## Key parameters

- 1,000 draws per sample and support mode.
- Seed: 260824.
- Matching strata: sample-specific empirical quintile of
  `log1p(total_counts)`, empirical quintile of `n_genes`, and PTPRC
  raw-positive/raw-zero status.
- Exact without-replacement sampling within each stratum.
- Rook adjacency: native coordinate offset 50 within topology component.

## Outputs

- Bin-to-stratum audit and exact matching validation.
- Observed spatial metrics.
- Matched and unmatched null draws.
- Empirical-percentile summaries and focus-sample comparison.
- Four report figures and a self-contained DRAFT HTML.

## Scope limits

This is descriptive QC. It performs no disease test, marker/threshold
optimization, Moran's I, hotspot scan, TLS call, clustering, cell2location,
NMF, ELN operation, Git commit, or push.

