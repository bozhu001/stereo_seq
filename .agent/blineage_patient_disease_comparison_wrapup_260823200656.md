# Session summary — cleaned B-lineage patient disease comparison wrap-up

Date: 2026-08-23 (Europe/Berlin)

## Scope and frozen provenance

This session recovered and reused immutable patient-level outputs from
`results/reanalysis/bin50_hq_frozen_blineage_descriptive_qc_260823180725/`.
It did not rerun Bin-level scoring. The six frozen panels, same-bin/rook rules,
continuity correction, model, samples and contrasts were not adapted to the
observed results.

Primary15 excluded the five pre-existing global low-quality-watch samples and
the independently specified extreme-depth sample `IPF/FO22-1-09404`, yielding
HC 4, IPF 4 and SSc-ILD 7. J2/L3-only13 additionally excluded the remaining K8
samples `HC/NL-50` and `SSC/24-1-18170A2`, yielding HC 3, IPF 4 and SSc-ILD 6.

## Analysis and results

The primary endpoint was the strict same-bin co-capture fraction with the
prespecified corrected-logit transformation. Each panel used the patient-level
model `transformed_endpoint ~ disease_group + log(sample_total_counts) +
chip_id`, HC reference and HC3 robust covariance. The three contrasts were
SSc-ILD versus HC, IPF versus HC and SSc-ILD versus IPF. All 18 primary P
values were adjusted once as a single BH family.

All 18 confidence intervals crossed zero and 0/18 tests reached BH-FDR < 0.05.
The six SSc-ILD-versus-IPF complete-cohort effects were negative in both 15-
and 13-patient models. Five panels retained that direction across estimable
leave-one-patient-out refits; `IgA_mucosal_plasma` had three flips per cohort.
This is an exploratory trend, not an established disease difference.

`GC_activated_B`, `ILD_plasma_cell`, `Memory_tissue_B` and `Naive_B` retained
all three complete-cohort directions after K8 removal. `Mature_B_core`
SSc-ILD-versus-HC and `IgA_mucosal_plasma` IPF-versus-HC flipped; both were
near-zero effects. Influence flags involved multiple patients rather than a
single uniform driver, but several weak contrasts were leave-one-patient-out
sensitive.

## Reporting and verification

The initial DRAFT table escaping problem was fixed with `knitr::asis_output`.
The approved final report changed only review-status labels and was rendered
from existing tables/figures; it did not rerun a model. Final HTML validation:
complete closing tag, 15 normal HTML tables, zero escaped tables, six embedded
images, no external image source and retained approved summary/interpretation.

Final report:
`results/reanalysis/bin50_cleaned_blineage_patient_disease_comparison_260823192715/report/FINAL_394_cleaned_blineage_patient_disease_comparison_report_260823192715.html`.

Size: 3,458,257 bytes.

SHA-256:
`0CE6E37D4B98CF6358F65772D5129449C951C8F5AB835C6A2116869FEC2F44CC`.

Analysis manifest:
`results/reanalysis/bin50_cleaned_blineage_patient_disease_comparison_260823192715/analysis_manifest_260823200656.json`.

Local analysis/report commit:
`446a30841cebf16d09c9094a3f251b6fff916679`.

## Documentation and boundaries

Local ELN draft:
`elab/Bin50_cleaned_Blineage_patient_disease_comparison_260823192715.md`.
It remains pending the separate ELN approval gate. No remote eLabFTW entry was
created, no report was uploaded, and no push was performed.

This bespoke analysis has a local script-review summary. Recurrence was not
confirmed, so no new-skill-candidate issue was filed.

No further changes to models, samples, markers or thresholds are planned from
these results. Any later spatial review requires a new, separately authorized
task.
