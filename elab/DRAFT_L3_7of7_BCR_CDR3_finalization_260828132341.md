# L3 7/7 BCR/CDR3 finalization — 260828132341

Status: approved local draft only; not submitted to eLabFTW. The project tag is intentionally omitted until its real value is confirmed.

Category: Bioinformatic

## Goal

Complete the L3 seven-patient BCR/CDR3 audit from existing, validated candidate inputs and completed TRUST4 outputs; reconnect productive CDR3-supporting reads to R1 CID/UMI and spatial Bin metadata; apply the frozen quality-aware UMI correction and clonotype classification rules used for J2; and compare J2 and L3 under identical thresholds.

## Input

1. Completed TRUST4 outputs for all seven L3 patients. The five results in `results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934/callers/trust4/` and the two protected earlier pilot results each had exactly one completed attempt with nonempty AIRR, barcode AIRR, assembled-read, assignment and cached-R2 outputs.
2. Existing validated candidate R1 and sidecar files under `results/reanalysis/existing_unenriched_bcr_pilot_260824/streamed_candidates/`. Original FASTQs were not opened or rescanned.
3. Completed two-patient L3 read-level/VDJ/spatial audit inputs reused by script 453.
4. Frozen J2 classification table: `results/reanalysis/l3_j2_final_umi_audit_260827_160202/final_reparsed_all_chunks_260827_234500/CLONOTYPE_RECLASSIFICATION_AFTER_UMI_QC.tsv`.

No prior eLabFTW entry IDs were available for these inputs, so no entry links are assigned in this draft.

## Script

Commands:

```bash
python3 scripts/453_finalize_l3_all7_bcr_cdr3.py results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934
python3 scripts/454_summarize_l3_all7_clonal_support.py results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934
```

Key frozen rules:

- UMI correction is restricted to patient + locus + exact CDR3 nt + CID.
- A quality-corrected molecule is a unique CID plus corrected UMI within that group.
- A formal reliable-rule clonotype requires IGH, at least two corrected molecules, at least two CID, at least two Bin, and at least two independent positional and strict R2 fragment signatures.
- A single UMI, single CID or sequence repetition is never interpreted as clonal expansion.

Repository: `https://github.com/bozhu001/stereo_seq`

Analysis commit: `c4897b8`

## Output summary

Output directory: `results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934/`

Timestamped report: `L3_ALL7_FINAL_METRICS_REPORT_260828132341.html` (7,195 bytes; SHA-256 `1B4F2942D635C4D15762F5086F742B59F590A5723D1278CB417C65A0BF2F8C9E`).

L3 contained 65 quality-corrected productive BCR/CDR3 molecules and 55 unique productive clonotypes: 26 singleton, 22 technical candidates, 6 moderate candidates and one formal reliable-rule hit. The sole formal hit, `SSC_15491_14_TRUST4_0003`, is retained only as a **provisional expanded IGH candidate** because all supporting reads carry the anomalously frequent `CGCTTGGCCT` motif and the exact CDR3 has a cross-patient warning. It is not a confirmed clonal expansion.

Under identical thresholds, J2 had 94 corrected molecules, 68 clonotypes, 9 clonotypes with at least two corrected CID-UMI molecules and no reliable-rule hit; L3 had 65, 55, 7 and one provisional formal hit, respectively.

Approved conservative conclusion: the current unenriched Stereo-seq data provide no reliable evidence of B-cell clonal expansion in either J2 or L3. This absence of reliable evidence does not demonstrate that biological clonal expansion is absent.

K8 was not scanned. No PPT was produced.

## Submission hold

Do not submit this draft or attach the report to eLabFTW until the real project tag is confirmed. No placeholder project tag may be used.
