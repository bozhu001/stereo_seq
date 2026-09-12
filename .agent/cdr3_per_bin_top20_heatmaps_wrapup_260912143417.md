# Session summary — CDR3 per-Bin Top20 heatmaps

Date: 2026-09-12

## Request and scope

Continued an interrupted Top20 plotting task without restarting or replacing
existing files. The user explicitly prohibited recomputation of upstream
IGH/CDR3 data and requested three RECURRENT main heatmaps, followed by three
FULL supplementary heatmaps, then QC, a MANIFEST, and SHA-256 validation. The
user later approved the analysis and classified the workflow as a one-off.

## Work completed

- Confirmed that `scripts/628_cdr3_per_bin_top20_heatmaps.R` and two prior
  complete, timestamped output directories already existed.
- Added immediate `FILE_COMPLETE` logging after each PNG and PDF is finalized,
  recording the basename, byte size, and completion timestamp.
- Ran a new non-overwriting verification output at
  `results/reanalysis/CDR3_PER_BIN_TOP20_HEATMAP_FINAL_20260912_141500/`.
- Generated three RECURRENT and three FULL figures in the requested order,
  each in PNG and PDF form. The verified run completed in about 25 seconds; no
  FULL heatmap reached the five-minute skip threshold.
- Verified 29/29 QC records, 26/26 MANIFEST rows, file sizes, SHA-256 values,
  PNG decoding, image dimensions, and representative visual appearance.
- Created a Quarto report source and PowerShell timestamped render wrapper.
  Windows/WSL path handling required rendering from `scripts/` with relative
  paths and Quarto `--no-clean`; failed intermediate attempts were preserved.
- Rendered the standalone, six-image HTML report
  `results/CDR3_PER_BIN_TOP20_HEATMAP_REPORT_260912142534.html` (6,327,999
  bytes; SHA-256
  `BBBC42693147166B5B7EC22A060825AF6458EAFD3BAEADFEEB6753E0CB4302EB`).
- Committed the analysis, outputs, logs, and report as
  `7aea87d49493dbd5d95c7ebb7a89c70d466e00e4`.

## Frozen results

- All69: 69 bins, 604 full genes, 93 recurrent genes, 1,183 non-NA entries.
- Shared7: 44 bins, 327 full genes, 55 recurrent genes, 740 non-NA entries.
- Repeated11: 33 bins, 226 full genes, 52 recurrent genes, 587 non-NA entries.
- Frozen upstream counts remained 104 corrected IGH molecules and 79
  patient-specific clonotypes.
- Heatmap values are `log1p(raw_count)` among each physical Bin50's frozen
  Top20 list; blank cells mean not ranked in that Top20, not necessarily
  nondetection. Recurrent means occurrence in at least two physical Bin50s.

## Documentation decisions

- The user approved the analysis before report rendering.
- The workflow is a one-off. No recurring-workflow skill candidate or GitHub
  issue was created; the R script remains for reproduction.
- The user separately approved the ELN draft and selected project tag
  `Stereo-seq`, retaining tags `IGH`, `CDR3`, `Bin50`, `Top20`, `heatmap`, and
  `QC`.
- No direct input had an available submitted eLabFTW entry ID, so the approved
  entry contains no prior-entry links.
- `ELABFTW_API_KEY` was unavailable. No remote entry was created and no remote
  ID/URL exists. The approved local final entry is
  `elab/CDR3_per_Bin_Top20_heatmaps_260912142534.md`; the standalone report is
  ready to attach through the web UI.
