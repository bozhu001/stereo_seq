# CDR3 per-Bin Top20 heatmaps — 260912142534

Status: approved by the user on 2026-09-12. Local final copy; remote eLabFTW
entry not created because `ELABFTW_API_KEY` is unavailable in this environment.

Category: Bioinformatic

Tags: `Stereo-seq`, `IGH`, `CDR3`, `Bin50`, `Top20`, `heatmap`, `QC`

## Goal

Generate publication-oriented recurrent and full heatmaps from the frozen
per-Bin Top20 gene tables for all 69 CDR3-positive physical Bin50s, the 44
cross-patient shared-CDR3 occurrences, and the 33 within-patient repeated-
clonotype Bin50s. Preserve the upstream IGH/CDR3 results unchanged and produce
explicit QC, a file manifest, and SHA-256 provenance.

## Input

1. Frozen All69 per-Bin Top20 table:
   `results/reanalysis/CDR3_PER_BIN_GENE_LIST_20260911_161928/02_ALL69_CDR3_BIN_TOP20_GENES.tsv`.
2. Frozen All69 Bin QC:
   `results/reanalysis/CDR3_PER_BIN_GENE_LIST_20260911_161928/03_ALL69_CDR3_BIN_QC.tsv`.
3. Frozen cross-patient shared-CDR3 Top20 table:
   `results/reanalysis/CDR3_PER_BIN_GENE_LIST_20260911_161928/04_SHARED7_CDR3_OCCURRENCE_TOP20.tsv`.
4. Frozen within-patient repeated-clonotype Top20 table:
   `results/reanalysis/CDR3_PER_BIN_GENE_LIST_20260911_161928/05_WITHIN_PATIENT_REPEATED_TOP20.tsv`.
5. Frozen 104-molecule index, six multi-IGH Bin50 annotations, and Top10 QC
   from `results/reanalysis/CDR3_PER_BIN_TOP10_FINAL_20260911_165003/`.
6. Frozen 69-column annotation table:
   `results/reanalysis/CDR3_PER_BIN_TOP10_HEATMAPS_20260911_172841/01C_ALL69_COLUMN_ANNOTATIONS.tsv`.
7. Top10 heatmap source used to preserve all definitions except TopN:
   `scripts/627_cdr3_per_bin_top10_heatmaps.R`.

The related L3 BCR/CDR3 entry exists only as an unsubmitted local draft, and
no eLabFTW entry IDs or URLs were found for the direct frozen inputs above.
Therefore no prior-entry links are assigned.

## Script

Plotting command:

```text
D:/bb/R/R-4.5.2/bin/x64/Rscript.exe scripts/628_cdr3_per_bin_top20_heatmaps.R results/reanalysis/CDR3_PER_BIN_GENE_LIST_20260911_161928 results/reanalysis/CDR3_PER_BIN_TOP10_FINAL_20260911_165003 results/reanalysis/CDR3_PER_BIN_TOP10_HEATMAPS_20260911_172841 results/reanalysis/CDR3_PER_BIN_TOP20_HEATMAP_FINAL_20260912_141500
```

Report command:

```text
cd scripts
quarto render 629_cdr3_per_bin_top20_heatmap_report.qmd --output CDR3_PER_BIN_TOP20_HEATMAP_REPORT_260912142534.html --output-dir ../results --no-clean
```

Key parameters and definitions:

- TopN is 20; all other definitions are unchanged from the frozen Top10
  heatmap workflow.
- A physical Bin50 is keyed by `sample_id + Bin50_id`.
- Heatmap values are `log1p(raw_count)` only for genes already ranked in that
  Bin50's frozen Top20.
- Blank cells mean “not ranked in this Bin's Top20,” not necessarily
  nondetection.
- Recurrent genes occur in at least two distinct physical Bin50s.
- Rows and columns are not clustered; no gene-wise z-score is applied.
- No raw expression matrix was read and no upstream IGH/CDR3 calculation was
  rerun.

Repository: `https://github.com/bozhu001/stereo_seq`

Analysis/report commit:
`7aea87d49493dbd5d95c7ebb7a89c70d466e00e4`

## Output summary

Output directory:
`results/reanalysis/CDR3_PER_BIN_TOP20_HEATMAP_FINAL_20260912_141500/`

Rendered standalone report:
`results/CDR3_PER_BIN_TOP20_HEATMAP_REPORT_260912142534.html`
(6,327,999 bytes; SHA-256
`BBBC42693147166B5B7EC22A060825AF6458EAFD3BAEADFEEB6753E0CB4302EB`).
The report is under 100 MB and should be attached when the entry is created in
the eLabFTW web UI.

All six requested heatmaps were generated in the required order: three
RECURRENT main figures followed by three FULL supplementary figures, each as
PNG and PDF. No FULL heatmap exceeded five minutes; the complete verified run
finished in approximately 25 seconds.

All69 contained 69 physical Bin50 columns, 604 full genes, 93 recurrent genes,
and 1,183 non-NA Top20 entries. Shared7 contained 44 columns, 327 full genes,
55 recurrent genes, and 740 non-NA entries. Repeated11 contained 33 columns,
226 full genes, 52 recurrent genes, and 587 non-NA entries.

Final QC passed 29/29 checks. Independent verification of all 26 MANIFEST rows
found zero missing paths, size mismatches, or SHA-256 mismatches. All six PNGs
decoded successfully, and visual review found no truncated or blank rendering.
Frozen counts remained 104 corrected IGH molecules and 79 patient-specific
clonotypes.

## Remote-entry status

No remote ID or URL is recorded. The current environment does not contain
`ELABFTW_API_KEY`; create this approved entry through the eLabFTW web UI using
category `Bioinformatic`, the seven approved tags above, and attach the
standalone HTML report.
