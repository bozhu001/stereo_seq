# J2 K=5 spatial niche architecture audit – 261008121443

Status: approved by the user on 2026-10-08. Local final copy; remote eLabFTW
entry was not created because `ELABFTW_API_KEY` is unavailable in this
environment.

Category: Bioinformatic

Tags: `Stereo-seq`, `Y40105J2`, `Reference18`, `Bin50`, `niche`, `spatial architecture`, `K5`, `QC`

## Objective

Audit the spatial continuity, fragmentation, direct Bin-to-Bin adjacency, and
cross-patient spatial consistency of all five niches in the frozen Y40105J2
K=5 unbiased niche assignment.

## Inputs and fixed definitions

- Frozen Reference18 input:
  `results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/J2/01_FROZEN_REFERENCE18_INPUT.tsv.gz`.
- Frozen K=5 assignment:
  `results/reanalysis/J2_UNBIASED_NICHE_RESOLUTION_20261008_110522/K05/02_BIN_NICHE_ASSIGNMENTS.tsv.gz`.
- Conservative annotation table:
  `results/reanalysis/J2_UNBIASED_NICHE_RESOLUTION_20261008_110522/11_J2_PRIMARYK_CONSERVATIVE_ANNOTATIONS.tsv`.
- 58,772 neighborhood-evaluable Bin50 from seven patients and seven ROIs.
- Queen adjacency was defined within patient and ROI as `|dx| <= 25 µm` and
  `|dy| <= 25 µm`, excluding self-edges and counting each physical edge once.
- Frozen lymphoid aggregate labels, members, and boundaries were not read.

The preceding approved local J2 discovery record is
`elab/J2_unbiased_neighborhood_niche_resolution_261008113930.md`. No remote
entry ID or URL is assigned rather than inventing unavailable provenance.

## Script and methods

- Analysis script: `scripts/1018_j2_k5_spatial_architecture.R`.
- Report source: `scripts/1019_j2_k5_spatial_architecture_report.qmd`.
- Run entry point: `bash WORK/run_j2_k5_spatial_architecture.sh`.
- A completed seven-patient permutation checkpoint was reused after a
  downstream table-ordering error was corrected; permutations were not rerun.
- All five niches and all ten heterotypic pairs were analyzed symmetrically.
- Observed/expected adjacency used 999 label permutations within patient and
  ROI, preserving niche counts, with seed `2026100811`.

Repository: `https://github.com/bozhu001/stereo_seq`

Analysis/report commit: `793085cee696e2432ce055212a286229d4429232`

## Results

The graph contained 218,614 unique direct adjacency edges. By the
cross-patient median largest-component fraction, the three most spatially
continuous labels were N2 (0.169), N4 (0.109), and N1 (0.100), although the
overall K=5 label field remained fragmented.

Nine of ten heterotypic pairs had the same direction in at least six of seven
patients; all were depleted relative to the exchangeable-label null. The
strongest seven-of-seven negative relations included N1–N3, N4–N5, N2–N4,
N2–N5, and N1–N2. No recurrent positive heterotypic adjacency enrichment was
identified.

These results describe spatial continuity, separation, and interface
preference. They are not independent cell-identity validation, cell-cell
communication, or causal interaction. K=5 labels derive from overlapping
50 µm Reference18 neighborhoods, and label permutation disrupts endogenous
spatial autocorrelation.

No disease-significance test, DEG/GSEA, pathway analysis, CellChat,
ligand-receptor analysis, frozen-aggregate overlap, or cross-chip consensus
was performed.

## Outputs

Result directory:
`results/reanalysis/J2_UNBIASED_NICHE_RESOLUTION_20261008_110522/K05_SPATIAL_ARCHITECTURE_20261008_120148/`

Rendered standalone report:
`results/J2_K5_SPATIAL_ARCHITECTURE_REPORT_261008121443.html`
(1,899,368 bytes; SHA-256
`65E128189DE0AD36042E8A00E73C8846EBCB8E9ABE6AD0EC136AFE118E86CFD4`).

## Remote-entry status

No remote ID or URL is recorded because `ELABFTW_API_KEY` is unavailable.
This approved local copy is the authoritative ELN record for the run.
