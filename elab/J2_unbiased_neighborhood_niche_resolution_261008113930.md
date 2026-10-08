# J2 unbiased neighborhood niche resolution – 261008113930

Status: approved by the user on 2026-10-08. Local final copy; remote eLabFTW
entry was not created because `ELABFTW_API_KEY` is unavailable in this
environment.

Category: Bioinformatic

Tags: `Stereo-seq`, `J2`, `Reference18`, `Bin50`, `VoltRon`, `niche`,
`K-means`, `QC`

## Goal

Perform independent Bin50-level, Reference18-based, 50 µm unbiased
neighborhood niche discovery for Y40105J2. Audit K=3–12 resolution stability
and evaluate the recommended K using Reference18 composition, cross-patient
repeatability, and UMI sensitivity.

## Input

1. Frozen J2 Reference18 Bin50 input:
   `results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/J2/01_FROZEN_REFERENCE18_INPUT.tsv.gz`.
2. Reference18 state order:
   `results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/J2/02_REFERENCE18_STATE_ORDER.tsv`.

No submitted local eLabFTW entry matching these direct inputs was found, so no
prior-entry ID or URL is assigned rather than inventing a provenance link.

## Script and methods

Analysis scripts:

- `scripts/1013_j2_unbiased_niche_resolution.R`
- `scripts/1014_j2_primaryk_composition_consistency.R`
- `scripts/1015_select_j2_k_resolution.R`
- `scripts/1016_finalize_j2_unbiased_niche_resolution.R`

Report source: `scripts/1017_j2_unbiased_niche_resolution_report.qmd`.

Key definitions:

- Bin50 with UMI >= 50 and frozen Reference18 continuous weights.
- Official 50 µm radius neighborhood within patient and ROI; the center Bin
  was excluded from its own neighborhood.
- VoltRon native CLR with structural zeros retained and no added pseudocount.
- K-means K=3–12 with 20 random seeds per K.
- Frozen aggregate labels were not read.
- Disease labels were metadata only and did not participate in neighborhood
  construction, clustering, K selection, or annotation.

Repository: `https://github.com/bozhu001/stereo_seq`

Analysis/report commit: `cb6b321`

## Output summary

Output directory:
`results/reanalysis/J2_UNBIASED_NICHE_RESOLUTION_20261008_110522/`

Rendered standalone report:
`results/J2_UNBIASED_NICHE_RESOLUTION_REPORT_261008113930.html`
(4,462,697 bytes; SHA-256
`FB4336E8A88D01E4A4F05195C531680602054C7EBE0514E1324F4C31C7ED2118`).
The report is below 100 MB and should be attached when the entry is created in
the eLabFTW web UI.

The analysis included 58,840 UMI-eligible Bin50 from seven patients
(HC=2, IPF=2, SSc-ILD=3); 58,772 had a non-empty neighborhood and 68 were
retained as not evaluable. The graph contained 325,776 edges with no
cross-patient or cross-ROI edges.

The audited primary resolution was K=5, with K=3 retained as an
under-resolution comparator and K=7 as a higher-resolution sensitivity. K=5
converged in 20/20 fits, had median pairwise ARI 0.988 and NMI 0.979, covered
all seven patients in every niche, and had a smallest niche fraction of 15.7%.

Conservative K=5 annotations were N1 alveolar-associated mixed, N2 B-enriched
mixed, N3 airway-associated mixed, N4 lymphatic-endothelial-associated mixed,
and N5 blood-endothelial-associated mixed. Core composition axes were not
single-patient-driven and retained direction in the specified UMI sensitivity
checks. These are mixed Reference18 neighborhood-composition patterns, not
confirmed pure cell populations.

No disease-significance test, DEG/GSEA, pathway analysis, CellChat,
ligand-receptor analysis, frozen-aggregate overlap, or cross-chip consensus
was performed.

## Remote-entry status

No remote ID or URL is recorded. The current environment does not contain
`ELABFTW_API_KEY`; create this approved entry through the eLabFTW web UI using
category `Bioinformatic`, the approved tags above, and attach the standalone
HTML report.

