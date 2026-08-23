# DRAFT — Bin50 B-lineage analysis-line synthesis — 260823215304

## Metadata

- **Title:** Bin50 B-lineage analysis-line synthesis — 260823215304
- **Category:** Bioinformatic
- **Tags:** `TODO_PROJECT_TAG`, `Bin50`, `spatial-omics`, `B-lineage`,
  `capture-QC`, `patient-level`, `matched-spatial-null`, `synthesis`
- **Analysis status:** Synthesis report and minimal completeness revisions approved on 2026-08-23.
- **ELN status:** Local draft only; this entry has not passed the separate ELN
  approval gate.
- **Remote status:** Not created. No API call, attachment upload or remote
  mutation was performed.
- **Remote entry ID/URL:** TODO after report approval and separate ELN approval.

## Goal

Close the current B-lineage analysis line by integrating existing marker
capture QC, the independently selected 16-patient high-quality cohort, the
FO22 extreme-depth sensitivity, six frozen panels, the 15/13-patient disease
comparison, the L3 focal review and the matched spatial-null audit. The task is
report-only and runs no new statistic, model, scoring, randomization or spatial
analysis. MS4A4A is excluded.

## Input

1. Project-existing marker capture QC:
   `results/reanalysis/bin50_all_existing_bcell_marker_capture_qc_260823170058/`.
2. Frozen 49-gene/six-panel descriptive QC and independent 16-sample cohort:
   `results/reanalysis/bin50_hq_frozen_blineage_descriptive_qc_260823180725/`.
3. FO22 extreme-depth sensitivity:
   `results/reanalysis/bin50_hq_blineage_extreme_depth_loo_sensitivity_260823185322/`.
4. Cleaned 15/13-patient disease comparison:
   `results/reanalysis/bin50_cleaned_blineage_patient_disease_comparison_260823192715/`.
5. Fixed L3 focal review:
   `results/reanalysis/bin50_l3_mature_b_focal_spatial_review_260823202319/`.
6. L3 matched spatial-null review:
   `results/reanalysis/bin50_l3_mature_b_depth_ptprc_matched_spatial_null_260823211809/`.

**Prior-entry link audit:** local ELN drafts exist for the cleaned disease
comparison and matched spatial-null analyses, but neither has an approved
remote eLabFTW entry ID/URL. No approved remote entry IDs were found for the
other immediate inputs. No remote instance query was performed, so no input
entry association is claimed in this draft.

## Script

Report-only render command:

```text
$env:BLINEAGE_SYNTHESIS_DIR = (Resolve-Path "results/reanalysis/bin50_blineage_analysis_line_synthesis_260823215304").Path
quarto render scripts/399_bin50_blineage_analysis_line_synthesis_report.qmd --to html
```

The QMD reads existing TSV, JSON and PNG files only. It does not source or run
the upstream analysis scripts. The workflow was confirmed as a one-off audit;
no reusable skill or GitHub issue was created.

Repository: <https://github.com/bozhu001/stereo_seq>

Local analysis/report commit:
`45f2672f14ec7dacce75930273241236e72f4b7b`.

Existing local commits were preserved. No push was performed.

## Output summary

### Confirmed technical and descriptive results

- The existing project-marker union contains 219 genes: 202 present in G15 and
  17 absent. All zero and absent rows remain represented.
- Detected marker number was strongly depth-associated descriptively across 21
  samples (Spearman rho 0.874 with sample total counts; no P value), while
  category CPM-depth correlations were weak.
- The literature-driven check froze 49 genes across six non-exclusive panels
  plus separate `TLS_context`; `TCL1A` and `FDCSP` were feature-space absent.
- Global technical QC independently retained 16 patients (HC 4, IPF 5,
  SSc-ILD 7). No sample was selected from B-lineage expression.
- FO22 was an extreme technical-depth sample that drove many upper maxima.
  Removing it reduced upper extremes but did not remove the residual K8-high
  capture pattern.

### Negative disease comparison

The primary cohort contained 15 patients and the J2/L3-only sensitivity 13.
All 18/18 primary 95% confidence intervals crossed zero and 0/18 reached
BH-FDR below 0.05. The disease comparison is negative. Point-estimate
directions and leave-one-out patterns do not establish a disease difference.

### Patient-specific focal candidates

`IPF/FO23-1-06168` and `SSC/15491/14` are retained as
**depth/n_genes/PTPRC-status–matched focal multi-marker B-lineage candidates**:
one patient-specific IPF candidate and one patient-specific SSc-ILD candidate.
Strict same-bin is the primary spatial evidence. Their focal signals are not
fully explained by the matched technical/immune background.

Strict rook support uses adjacency in its construction; rook clustering is
construction-dependent auxiliary evidence and not an independent spatial
validation.

The final report now includes a six-row reference table transcribed from the
approved frozen `strict_panel_rules`, `strict_rule` and
`strict_rook_definition`. No marker or threshold was changed. The common rook
auxiliary rule requires a direct native-coordinate orthogonal distance of 50
within the same authoritative source-tissue component, contribution from both
edge endpoints, and a two-endpoint marker union satisfying the corresponding
frozen same-bin rule.

For the two matched focal candidates, the approved report states:
“100th empirical percentile表示在现有1000次匹配随机抽样中，没有随机集合达到观察到的空间集中程度；它不是正式的疾病比较P值，也不能用于疾病组推断。”

### Required next validation and prohibited interpretation

The histology statement is conditional:
“若取得对应H&E且能够完成可靠的组织轮廓或粗区域配准，则核对候选坐标对应的解剖与病理结构。”
Orthogonal RNA/protein confirmation remains required. Current evidence cannot
be called a B-cell niche, B-cell aggregate, TLS or disease-specific B-lineage
structure. Marker co-capture fraction is not actual B-cell area.

Result directory:
`results/reanalysis/bin50_blineage_analysis_line_synthesis_260823215304/`.

FINAL standalone HTML:
`report/FINAL_399_bin50_blineage_analysis_line_synthesis_report_260823215304.html`.

Size: 13,794,494 bytes.

SHA-256:
`B97384DE5636CCF29F59D607CD4768AF139FB2A281DB9851A1A754CB72FB52F5`.

The HTML is below 100 MB and should be attached only after this ELN draft
passes its separate approval gate.

Manifest:
`analysis_manifest_260823215304.json`.

## Deferred fields

- Manual approval of the synthesis report: approved.
- Separate approval of this ELN draft: pending.
- Project-specific eLabFTW tag: `TODO_PROJECT_TAG`.
- Remote entry ID/URL and attachment: not created.
- Analysis/report commit: `45f2672f14ec7dacce75930273241236e72f4b7b`.
- Git push: not performed.
