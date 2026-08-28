# Session summary — L3 7/7 BCR/CDR3 finalization

Date: 2026-08-28

## Objective and scope

Safely finish the L3 seven-patient BCR/CDR3 audit from already generated candidate inputs and completed TRUST4 results, without rerunning any patient, rescanning original FASTQs, generating candidate FASTQs, scanning K8 or producing a PPT.

## Completion checks

- Confirmed exactly one completed TRUST4 attempt for each of seven L3 patients.
- Confirmed nonempty ordinary AIRR, barcode AIRR, assembled-read, assignment and cached-R2 outputs plus `TRUST4_COMPLETE.ok` for every completed attempt.
- Preserved the interrupted `SSC_05957_17B` partial directory; the completed attempt was `SSC_05957_17B_resume_260828_111831`.
- Confirmed the 451 checkpoint launcher ended normally with `TRUST4_REMAINING5_COMPLETE.ok`.
- No 451, TRUST4, 453 or 454 process remained after completion.

## Finalization workflow

Ran in detached tmux sessions:

```bash
python3 scripts/453_finalize_l3_all7_bcr_cdr3.py results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934
python3 scripts/454_summarize_l3_all7_clonal_support.py results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934
```

Script 453 reconnected productive assembled CDR3 QNAMEs to cached R2, candidate R1 and sidecars; checked unique QNAME recovery and R1 structure; applied the frozen J2 quality-aware UMI correction; summarized VDJ/CDR3, fragment, CID/Bin and cross-patient support; and created `L3_ALL7_FINAL_COMPLETE.ok`. Script 454 produced explicit patient-level molecule/clonotype/support metrics and a timestamped HTML report.

## Results

- L3: 65 quality-corrected productive BCR/CDR3 molecules and 55 unique productive clonotypes.
- Classification: 26 singleton, 22 technical candidates, 6 moderate candidates and one formal reliable-rule hit.
- Seven L3 clonotypes had at least two corrected CID-UMI molecules; six had at least two distinct raw UMI sequences; maximum corrected expansion was three molecules.
- Patient corrected molecule / clonotype counts: HC/NL-66 18/17; HC/NL-72 10/8; IPF/FO23-1-06168 4/4; IPF/FO23-1-06170 7/5; SSC/05957/17B 5/5; SSC/07998/15A 6/4; SSC/15491/14 15/12.
- `SSC/05957/17B` had no repeated clonotype. `SSC/07998/15A` had one three-molecule IGK moderate candidate with 2 CID and 2 Bin, but it did not satisfy the reliable IGH rule and carried a cross-patient warning.
- The sole L3 formal hit was `SSC_15491_14_TRUST4_0003`: IGH, 3 corrected molecules, 3 CID, 3 Bin and sufficient fragment support. Because 6/6 supporting reads carried the anomalously frequent `CGCTTGGCCT` motif and the exact CDR3 had a cross-patient warning, the approved wording is **provisional expanded IGH candidate**, never confirmed clonal expansion.
- J2 under the identical frozen classification threshold: 94 corrected molecules, 68 clonotypes, 9 repeated corrected-molecule clonotypes, maximum 11 molecules and zero reliable-rule hits.
- Approved overall conclusion: current unenriched Stereo-seq data provide no reliable evidence of B-cell clonal expansion in either J2 or L3, but this does not prove that biological clonal expansion is absent.

## Key outputs

- `results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934/L3_ALL7_FINAL_COMPLETE.ok`
- `results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934/L3_ALL7_PATIENT_CLONAL_SUPPORT_SUMMARY.tsv`
- `results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934/L3_ALL7_CLONOTYPE_CLASSIFICATION.tsv`
- `results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934/L3_ALL7_CDR3_SUPPORTING_READS.tsv`
- `results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934/L3_ALL7_UMI_CORRECTION_DECISIONS.tsv`
- `results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934/L3_J2_IDENTICAL_THRESHOLD_METRICS.tsv`
- `results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934/L3_ALL7_FINAL_METRICS_REPORT_260828132341.html`

Report SHA-256: `1B4F2942D635C4D15762F5086F742B59F590A5723D1278CB417C65A0BF2F8C9E`.

Analysis/report commit: `c4897b8`.

## Documentation decisions

- The user approved the final analysis and conservative conclusion.
- The ELN draft is saved locally but must not be submitted before the real eLabFTW project tag is confirmed; the tag is omitted and no placeholder is used.
- Scripts 453/454 are retained as a reusable audit workflow within this project. The user explicitly declined promotion to a general bioinfo-wiki skill, so no skill-candidate issue was created.
- K8 remains unscanned and no PPT was created.
