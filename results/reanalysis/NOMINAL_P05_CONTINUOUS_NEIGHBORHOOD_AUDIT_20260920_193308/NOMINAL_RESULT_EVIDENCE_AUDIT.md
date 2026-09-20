# Nominal P<0.05 Continuous Neighborhood evidence audit

## Scope and interpretation boundary

This is a read-only re-audit of frozen Primary19 Continuous Neighborhood results. No RCTD, neighborhood construction, permutation, threshold optimization, SPARK-X analysis, or broad new hypothesis testing was run. Rows enter the audit when raw P<0.05 or the already-computed chip-stratified permutation P<0.05.

No result lacking family/global FDR control is described as significant, supported, or disease-specific. Pair metrics are co-enrichment summaries and are not evidence of cell-cell interaction.

## Multiplicity accounting

| Contrast | Tests | Expected chance hits (5%) | Raw P<0.05 | Permutation P<0.05 | Either P<0.05 | Family FDR<0.05 | Global FDR<0.05 |
|---|---:|---:|---:|---:|---:|---:|---:|
| SSc-ILD_vs_HC | 2682 | 134.1 |  54 |  55 |  63 | 0 | 0 |
| IPF_vs_HC | 2682 | 134.1 | 101 | 105 | 117 | 0 | 0 |
| SSc-ILD_vs_IPF | 2682 | 134.1 | 111 | 123 | 129 | 0 | 0 |

Across all contrasts, 309 unique rows had raw or permutation P<0.05. The expected chance-hit count is 5% of the tested family shown above, not a subtraction rule and not an estimate of true discoveries.

## Mutually exclusive classification

- IPF_ASSOCIATED_TREND: 4
- SHARED_ILD_VS_HC_TREND: 2
- SSC_VS_HC_PRIORITY_EXPLORATORY: 0
- SSC_VS_IPF_ONLY_LOW_PRIORITY: 3
- TECHNICALLY_SENSITIVE: 199
- UNSTABLE_OR_DESCRIPTIVE: 101

After the frozen technical/stability filters, 9 rows remain in categories 1-4. This count includes contrast-specific rows; paired SSc-ILD-vs-HC and IPF-vs-HC rows for the same shared-ILD hypothesis are counted separately.

A 2/3 chip result is retained only under the already-frozen K8+J2 direction rule and only when LOPO, LOCO, Primary19-All21, and technical checks pass; its L3 disagreement is kept explicit. No new cutoff was introduced.

## SSc-ILD versus HC nominal candidates

- E04027: 21-state Basal epithelial ↔ NK; fisher_z @50 um; SSc-ILD_vs_HC; effect=0.06433; raw P=0.0262; permutation P=0.0256; family FDR=0.9751; global FDR=0.786.

These are exploratory priorities or shared-ILD trends only. They are not FDR-supported disease effects unless the row-level FDR columns explicitly pass.

## Shared ILD versus HC trends

- E04028: 21-state Basal epithelial ↔ NK; fisher_z @50 um; IPF_vs_HC; effect=0.06926; raw P=0.02521; permutation P=0.023; family FDR=0.9751.
- E04027: 21-state Basal epithelial ↔ NK; fisher_z @50 um; SSc-ILD_vs_HC; effect=0.06433; raw P=0.0262; permutation P=0.0256; family FDR=0.9751; global FDR=0.786.

## IPF-associated trends

- E05396: 21-state Plasma ↔ SPP1 macrophage; fisher_z @50 um; IPF_vs_HC; effect=0.09766; raw P=0.004637; permutation P=0.0021; family FDR=0.07559.
- E01352: 18-state Ciliated epithelial ↔ Dendritic cell; fisher_z @50 um; IPF_vs_HC; effect=0.0949; raw P=0.02103; permutation P=0.0173; family FDR=0.6935.
- E01364: 18-state Blood endothelial ↔ Dendritic cell; fisher_z @50 um; IPF_vs_HC; effect=0.05217; raw P=0.03103; permutation P=0.0449; family FDR=0.8946.
- E03524: 21-state Ciliated epithelial ↔ Dendritic cell; fisher_z @50 um; IPF_vs_HC; effect=0.08871; raw P=0.0427; permutation P=0.0353; family FDR=0.9751.

## SSc-ILD versus IPF only, low priority

- E01365: 18-state Blood endothelial ↔ Dendritic cell; fisher_z @50 um; SSc-ILD_vs_IPF; effect=-0.05677; raw P=0.01087; permutation P=0.009399; family FDR=0.6526.
- E01050: 18-state B ↔ T; fisher_z @50 um; SSc-ILD_vs_IPF; effect=-0.07082; raw P=0.04594; permutation P=0.0249; family FDR=0.4482.
- E03573: 21-state Blood endothelial ↔ Dendritic cell; fisher_z @50 um; SSc-ILD_vs_IPF; effect=-0.04796; raw P=0.0531; permutation P=0.047; family FDR=0.994.

These between-fibrotic-disease rows are deliberately ranked below SSc-ILD-vs-HC evidence and are not SSc-specific findings.

## Technical and stability downgrades

- Technical-sensitive nominal rows: 199. This includes strong association with frozen technical variables, chip-range sensitivity, or HC/NL-66 and SSC/05957/17B low-coverage sensitivity.
- Unstable/descriptive nominal rows: 101. Reasons include secondary/scale-sensitivity tier, insufficient chip direction, LOPO/LOCO reversal, Primary19-All21 disagreement, 18/21-state direction discordance, or a low-reference-state caution.

All row-level technical correlations are retained in the main TSV: mean state weight, valid Bin count, tissue area, coverage fraction, median UMI, and median detected genes.

## Potential links to the pending SPARK-X gene programs

- Basal epithelial ↔ NK: test whether the pending SPARK-X program score spatially covaries with the frozen state/pair metric in a narrowly prespecified follow-up. This is a testable bridge, not current support.
- Plasma ↔ SPP1 macrophage: test whether the pending SPARK-X program score spatially covaries with the frozen state/pair metric in a narrowly prespecified follow-up. This is a testable bridge, not current support.
- Ciliated epithelial ↔ Dendritic cell: test whether the pending SPARK-X program score spatially covaries with the frozen state/pair metric in a narrowly prespecified follow-up. This is a testable bridge, not current support.
- Blood endothelial ↔ Dendritic cell: test whether the pending SPARK-X program score spatially covaries with the frozen state/pair metric in a narrowly prespecified follow-up. This is a testable bridge, not current support.

The SPARK-X run was not read, modified, or used to reorder these results. Any future link must be tested after the program definitions are frozen, with a narrowly prespecified state/program pairing and explicit multiplicity control.

## Files

The five requested TSVs and three PNG/PDF figure pairs are in this directory. `NOMINAL_P05_COUNT_AUDIT.tsv` is an additional machine-readable count summary.
