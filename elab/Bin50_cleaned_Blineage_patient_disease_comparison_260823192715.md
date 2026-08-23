# DRAFT — Bin50 cleaned B-lineage patient disease comparison — 260823192715

## Metadata

- **Title:** Bin50 cleaned B-lineage patient disease comparison — 260823192715
- **Category:** Bioinformatic
- **Tags:** `TODO_PROJECT_TAG`, `Bin50`, `spatial-omics`, `B-lineage`,
  `patient-level`, `exploratory-disease-comparison`, `HC3`,
  `K8-sensitivity`
- **Analysis approval:** The analysis, statistical interpretation and corrected
  HTML were approved by the user on 2026-08-23.
- **ELN status:** Local draft only; this entry has not passed the separate ELN
  approval gate.
- **Remote status:** Not created. No API call, remote entry, attachment upload,
  commit push or other remote mutation was performed.
- **Remote entry ID/URL:** TODO after separate ELN-draft approval.

## Goal

Compare SSc-ILD, IPF and healthy controls for six frozen B-lineage panels at
the patient level after independent removal of five globally flagged
low-quality samples and one extreme-depth sample. Assess whether adjusted
disease-effect directions are sensitive to the two remaining K8 samples,
without changing samples, markers, thresholds, panel rules or the model in
response to results.

## Input

All inputs were immutable patient-level outputs from:
`results/reanalysis/bin50_hq_frozen_blineage_descriptive_qc_260823180725/`.
No Bin-level score was recalculated.

1. Panel CPM and raw-count summaries:
   `tables/05_sample_panel_detected_marker_counts.tsv`
   (SHA-256 `5402989D4F113D0B40F60029C40892D9746A81565CF733DC9CF2726221C4AAA3`).
2. Frozen strict same-bin support:
   `tables/06_panel_same_bin_strict_support.tsv`
   (SHA-256 `8E0326583DD3E06D00819BE29B2B8C0531CB3BF67774C945B48557FEEB606F69`).
3. Frozen strict rook support:
   `tables/07_panel_rook_strict_support.tsv`
   (SHA-256 `99CD4F28EFB4C9963B05BE8C6C19317AA4847348490026E870634757EA4D7187`).
4. Independently selected high-quality sample metadata:
   `tables/10_high_quality_sample_descriptive_metadata.tsv`
   (SHA-256 `6ED6798104B2FDEC438CA2D0D44940CA75F8D0088BE51AD0BB76086042B0A706`).
5. Frozen six-panel marker manifest and source validation tables from the same
   result directory.

**Prior-entry link audit:** no approved local ELN entry ID for the immediate
frozen B-lineage descriptive-QC input was found in `elab/` or recorded in
`.agent/memory.md`. No remote instance query was performed, so no input-entry
association is claimed.

## Script

Patient-level analysis command completed before approval:

```text
Rscript scripts/393_cleaned_blineage_patient_disease_comparison.R --timestamp=260823192715
```

Final report-only render command:

```text
$env:BLINEAGE_DISEASE_RESULT_DIR = (Resolve-Path "..\results\reanalysis\bin50_cleaned_blineage_patient_disease_comparison_260823192715").Path
quarto render scripts/394_cleaned_blineage_patient_disease_comparison_report.qmd --output FINAL_394_cleaned_blineage_patient_disease_comparison_report_260823192715.html
```

The final render read existing TSV and PNG outputs only. It did not fit or
refit a model.

Key fixed parameters:

- Primary15: HC 4, IPF 4, SSc-ILD 7.
- J2/L3-only13: HC 3, IPF 4, SSc-ILD 6.
- Panels: `Mature_B_core`, `ILD_plasma_cell`, `IgA_mucosal_plasma`,
  `GC_activated_B`, `Naive_B`, `Memory_tissue_B`.
- Endpoint transformation:
  `logit((strict_same_bin_supported_bin_n + 0.5)/(total_bin_n + 1))`.
- Model:
  `transformed_endpoint ~ disease_group + log(sample_total_counts) + chip_id`.
- HC reference; patient as the independent unit; HC3 robust covariance.
- Three prespecified contrasts per panel and one BH correction across all 18
  primary tests.
- Rook support and panel CPM used only for descriptive direction checks.

Repository: <https://github.com/bozhu001/stereo_seq>

Local analysis/report commit:
`446a30841cebf16d09c9094a3f251b6fff916679`.

Remote status: not pushed as of 2026-08-23.

## Output summary

All 18/18 primary 95% confidence intervals crossed zero and 0/18 tests reached
BH-FDR < 0.05. For SSc-ILD versus IPF, all six complete-cohort adjusted effects
were negative in both Primary15 and J2/L3-only13. Leave-one-patient-out refits
retained the direction for five of six panels; `IgA_mucosal_plasma` had three
direction-flipping refits in each cohort. The cross-panel pattern is an
exploratory trend and is not an established disease difference.

Four panels preserved all three complete-cohort effect directions after K8
removal: `GC_activated_B`, `ILD_plasma_cell`, `Memory_tissue_B` and `Naive_B`.
`Mature_B_core` SSc-ILD versus HC and `IgA_mucosal_plasma` IPF versus HC
flipped after K8 removal; both estimates were close to zero. Multiple patients
were flagged across influence diagnostics, so no single patient uniformly
drove every panel, but several near-zero contrasts were leave-one-patient-out
sensitive.

Interpretation limits remain fixed: marker co-capture is not B-cell physical
area; plasma/IgA panels may be influenced by abundant immunoglobulin-related
transcripts; `GC_activated_B` is not a germinal-center or TLS call; and no
nominal P value is treated as a reliable finding. No model, sample, marker or
threshold will be adjusted in response to these results.

Result directory:
`results/reanalysis/bin50_cleaned_blineage_patient_disease_comparison_260823192715/`.

Final standalone HTML:
`report/FINAL_394_cleaned_blineage_patient_disease_comparison_report_260823192715.html`
(3,458,257 bytes; SHA-256
`0CE6E37D4B98CF6358F65772D5129449C951C8F5AB835C6A2116869FEC2F44CC`).

The HTML is below 100 MB and should be attached if this draft is separately
approved for remote eLabFTW creation.

Analysis manifest:
`analysis_manifest_260823200656.json`.

## Deferred remote fields

- Separate user approval of this ELN draft: pending.
- Project-specific eLabFTW tag: `TODO_PROJECT_TAG` (unknown; do not infer).
- Remote entry ID and URL: TODO.
- Remote creation and final-HTML attachment: not performed.
- Git push: not performed.
