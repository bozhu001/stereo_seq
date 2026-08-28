# L3 7/7 BCR/CDR3 final report

Generated: 2026-08-28T13:20:55.662916+02:00

## Completion and integrity

L3 is complete for 7/7 manifest patients. Existing candidate FASTQs were used; original FASTQs were not rescanned. The two previously completed TRUST4 patients were read-only and were not rerun. Every productive supporting QNAME was unique in the candidate R1, TRUST4 cached R2 and sidecar manifest; R1 was exactly 35 bp with CID=`R1[0:25]`, UMI=`R1[25:35]`, and complete quality.

Total productive BCR clonotypes: **55**. Classification: **26 singleton, 22 technical, 6 moderate, 1 reliable-rule hit**.

## Moderate/reliable candidates

- SSC/07998/15A / SSC_07998_15A_TRUST4_0002 / IGK / moderate candidate
- SSC/15491/14 / SSC_15491_14_TRUST4_0002 / IGK / moderate candidate
- SSC/15491/14 / SSC_15491_14_TRUST4_0003 / IGH / reliable expanded clonotype
- IPF/FO23-1-06170 / IPF_FO23-1-06170_TRUST4_0001 / IGK / moderate candidate
- HC/NL-66 / HC_NL-66_TRUST4_0002 / IGK / moderate candidate
- HC/NL-72 / HC_NL-72_TRUST4_0004 / IGK / moderate candidate
- HC/NL-72 / HC_NL-72_TRUST4_0005 / IGK / moderate candidate

`SSC_15491_14_TRUST4_0003` remains **spatially supported provisional expanded IGH candidate** after its completed read-level, VDJ and spatial validation. It meets the uniform reliable-rule molecule gate, but is deliberately worded as provisional because all molecules use the anomalously frequent UMI motif and the cross-patient junction warning remains.

## UMI motif

`CGCTTGGCCT` occurs in **430/462 (93.07%)** productive CDR3-supporting reads. Repetition of this motif is never treated as molecule independence; CID, Bin and R2 fragment evidence are primary.

## Cross-patient audit

Exact cross-patient CDR3 warnings: **7**. These are labelled `public-like_or_technical_warning` and are never merged across patients.

## J2 comparison

Both J2 and L3 are complete 7/7. J2: 68 productive clonotypes, 34 singleton, 25 technical, 9 moderate and 0 reliable. L3: 55 productive clonotypes, 26 singleton, 22 technical, 6 moderate and 1 reliable-rule hit (reported as provisional after spatial review).

The comparison is descriptive: unenriched random capture, sparse direct junction support and the chip-level UMI motif abnormality prevent interpreting raw reads as clone size or testing disease prevalence from these seven-patient chips.

## Scope

Conclusions apply only to completed J2 and L3. They do not extend to K8 or all 21 patients. No TRA/TRB/TRD/TRG record is included in BCR counts.
