# L3 K8 spatial niche architecture – 261008104540

Status: approved by the user on 2026-10-08. Local final copy; remote eLabFTW
entry was not created because `ELABFTW_API_KEY` is unavailable in this
environment.

Category: Bioinformatic

Tags: `Stereo-seq`, `L3`, `Reference18`, `Bin50`, `niche`, `spatial`, `K8`,
`QC`

## Goal

Characterize the spatial continuity, connected-component structure, direct
Bin50-to-Bin50 adjacency, and patient/ROI-stratified permutation-null
enrichment of the frozen L3 K=8 unbiased neighborhood assignment. Analyze all
eight niches and all 28 heterotypic niche pairs symmetrically before any
targeted interpretation.

## Input

1. Frozen L3 Reference18 Bin50 input:
   `results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/L3/01_FROZEN_REFERENCE18_INPUT.tsv.gz`.
2. Frozen K=8 Bin-to-niche assignment:
   `results/reanalysis/L3_UNBIASED_NICHE_RESOLUTION_20261007_085730/K08/02_BIN_NICHE_ASSIGNMENTS.tsv.gz`.
3. Conservative K=8 annotation table:
   `results/reanalysis/L3_UNBIASED_NICHE_RESOLUTION_20261007_085730/K08_HIGH_RESOLUTION_SENSITIVITY_20261007_104602/08_K08_CONSERVATIVE_NICHE_ANNOTATIONS.tsv`.

No submitted local eLabFTW entry matching these direct inputs was found, so no
prior-entry ID or URL is assigned rather than inventing a provenance link.

## Script

Analysis command:

```text
bash WORK/run_l3_k8_spatial_architecture.sh
```

Analysis script: `scripts/1011_l3_k8_spatial_architecture.R`.

Report source: `scripts/1012_l3_k8_spatial_architecture_report.qmd`.

Key definitions:

- Frozen L3 K=8 labels were reused; clustering was not rerun.
- Physical adjacency was Queen adjacency within patient and ROI, with
  `abs(dx) <= 25 µm` and `abs(dy) <= 25 µm`, excluding self-pairs.
- Each undirected Bin pair was counted once.
- All eight niches and all 28 heterotypic pairs entered the main analysis.
- The null fixed the patient/ROI graph and niche counts and permuted labels 999
  times within patient/ROI.
- Frozen lymphoid aggregates and disease labels were not used in any spatial
  definition or test.

Repository: `https://github.com/bozhu001/stereo_seq`

Analysis/report commit: `f0d0d1a`

## Output summary

Output directory:
`results/reanalysis/L3_UNBIASED_NICHE_RESOLUTION_20261007_085730/K08_SPATIAL_ARCHITECTURE_20261007_121357/`

Rendered standalone report:
`results/L3_K8_SPATIAL_ARCHITECTURE_REPORT_261008104540.html`
(1,266,544 bytes; SHA-256
`0893C945FCA3E2CF48C09F37E73B377FD64320952DDCAD30855A05EE1B7ED8A4`).
The report is under 100 MB and should be attached when the entry is created in
the eLabFTW web UI.

The analysis included 39,982 evaluable Bin50s from seven patients and seven
ROIs and 149,193 unique Queen-adjacency edges. N7, N4, and N5 showed the
highest spatial continuity. N3 was comparatively fragmented, with a median
largest-component fraction of 0.065. N5–N6 had recurrent positive adjacency
enrichment in 6/7 patients (median log2 observed/expected +0.284). N4–N7,
N2–N7, and N7–N8 had negative adjacency enrichment in 7/7 patients.

The permutation null disrupts endogenous spatial autocorrelation, and the
niche labels themselves derive from overlapping 50 µm neighborhoods.
Consequently, these results support spatial association/interface preference,
not independent validation, cell-cell communication, or causal interaction.
No frozen aggregate analysis, disease significance test, DEG/GSEA, pathway
analysis, or ligand-receptor/CellChat analysis was performed.

## Remote-entry status

No remote ID or URL is recorded. The current environment does not contain
`ELABFTW_API_KEY`; create this approved entry through the eLabFTW web UI using
category `Bioinformatic`, the approved tags above, and attach the standalone
HTML report.
