# Bespoke script review — pending recurrence decision

- Script: `393_cleaned_blineage_patient_disease_comparison.R`
- Purpose: exploratory patient-level disease comparison of six frozen B-lineage panels after independently prespecified technical exclusions.
- Immutable inputs: frozen patient-level panel CPM, strict same-bin support, strict rook support, sample metadata, panel manifest and source validation from `bin50_hq_frozen_blineage_descriptive_qc_260823180725`.
- Cohorts: Primary15 (HC 4, IPF 4, SSc-ILD 7); J2/L3-only13 (HC 3, IPF 4, SSc-ILD 6).
- Endpoint: continuity-corrected logit of strict same-bin supported Bin fraction.
- Model: `transformed_endpoint ~ disease_group + log(sample_total_counts) + chip_id`, HC reference, HC3 robust covariance.
- Prespecified inference: three contrasts per panel; one BH correction over all 18 primary tests.
- Sensitivities: unchanged model in J2/L3-only13, descriptive rook/CPM direction checks, conventional influence diagnostics and leave-one-patient-out refits.
- Outputs: patient/group tables, full model results, global FDR table, adjusted predictions, sensitivity/stability tables, diagnostics, six PNG/SVG plot pairs and a self-contained DRAFT HTML.
- Scope exclusions: no Bin-level recomputation, marker/threshold changes, TLS/Moran/hotspot calling, cell2location, NMF, GSEA, differential gene testing or new exclusions.

This is a bespoke analysis rather than an existing skill template. No new-skill-candidate issue was filed because recurrence has not been confirmed by the user.
