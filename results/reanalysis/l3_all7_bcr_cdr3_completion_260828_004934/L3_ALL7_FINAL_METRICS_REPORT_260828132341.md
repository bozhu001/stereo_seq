# L3 7/7 final BCR/CDR3 support metrics

Generated: 2026-08-28T13:23:41.664646+02:00

## Definitions

- A quality-corrected molecule is one unique CID plus corrected UMI within patient, chain and exact CDR3 nt.
- Sequence repetition, a single UMI, or a single CID is not clonal expansion.
- The uniform reliable rule requires IGH, at least 2 corrected molecules, at least 2 CID, at least 2 Bin, and at least 2 positional and strict R2 fragment signatures.
- Distinct raw UMI sequence counts are shown separately because identical UMI sequences in distinct CID can still be distinct molecules but are not distinct UMI sequences.

## L3 patient metrics

| patient | quality_corrected_productive_bcr_molecules | unique_productive_bcr_clonotypes | clonotypes_ge2_corrected_cid_umi_molecules | clonotypes_ge2_distinct_raw_umi_sequences | maximum_corrected_molecules_per_clonotype | maximum_clonotype_unique_cid | maximum_clonotype_unique_bin | formal_reliable_rule_clonotypes | reliable_expansion_interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC/NL-66 | 18 | 17 | 1 | 1 | 2 | 1 | 1 | 0 | no_reliable_clonal_expansion |
| HC/NL-72 | 10 | 8 | 2 | 2 | 2 | 2 | 2 | 0 | no_reliable_clonal_expansion |
| IPF/FO23-1-06168 | 4 | 4 | 0 | 0 | 1 | 1 | 1 | 0 | no_reliable_clonal_expansion |
| IPF/FO23-1-06170 | 7 | 5 | 1 | 1 | 3 | 2 | 2 | 0 | no_reliable_clonal_expansion |
| SSC/05957/17B | 5 | 5 | 0 | 0 | 1 | 1 | 1 | 0 | no_reliable_clonal_expansion |
| SSC/07998/15A | 6 | 4 | 1 | 1 | 3 | 2 | 2 | 0 | no_reliable_clonal_expansion |
| SSC/15491/14 | 15 | 12 | 2 | 1 | 3 | 3 | 3 | 1 | formal_reliable_rule_hit_but_provisional_not_unqualified_expansion |

## J2 versus L3 under the frozen identical classification threshold

| chip | patients | quality_corrected_productive_bcr_molecules | unique_productive_bcr_clonotypes | clonotypes_ge2_corrected_cid_umi_molecules | clonotypes_ge2_distinct_raw_umi_sequences | maximum_corrected_molecules_per_clonotype | repeated_clonotypes_with_ge2_cid_and_ge2_bin | singleton | technical_candidate | moderate_candidate | formal_reliable_rule_clonotypes | interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| J2 | 7 | 94 | 68 | 9 | 8 | 11 | 6 | 34 | 25 | 9 | 0 | no_reliable_expanded_clonotype |
| L3 | 7 | 65 | 55 | 7 | 6 | 3 | 4 | 26 | 22 | 6 | 1 | one_formal_hit_reported_as_spatially_supported_provisional_expanded_IGH_candidate |

## Interpretation

L3 has 55 productive BCR clonotypes: 26 singleton, 22 technical, 6 moderate, and 1 formal reliable-rule hit. The sole formal hit is SSC_15491_14_TRUST4_0003 (IGH; 3 corrected molecules; 3 CID; 3 Bin), retained as a spatially supported provisional expanded IGH candidate because all supporting reads carry the anomalously frequent CGCTTGGCCT motif and the exact CDR3 has a cross-patient warning. It is not reported as unqualified reliable clonal expansion.

Python: 3.12.3 (Linux-6.6.87.2-microsoft-standard-WSL2-x86_64-with-glibc2.39)
