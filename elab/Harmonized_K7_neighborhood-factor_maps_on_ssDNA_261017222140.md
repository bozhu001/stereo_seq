# Harmonized K7 neighborhood-factor maps on ssDNA — 261017222140

Status: approved by the user on 2026-10-17. Local final copy; remote eLabFTW
entry was not created because `ELABFTW_API_KEY` is unavailable in this
environment.

Category: Bioinformatic

Tags: `Stereo-seq`, `ssDNA`, `Reference18`, `Neighborhood NMF`, `K7`,
`Bin50`, `spatial mapping`, `immune microenvironment`

## Goal

Map the frozen harmonized K=7 neighborhood factors onto the K8, J2 and L3
ssDNA images. Display the dominant K=7 factor across all 21 patients and the
continuous B-enriched/mixed-immune factor usage on the same tissue backgrounds.

## Input

1. Frozen K=7 per-chip model and usage archives:
   `results/reanalysis/REFERENCE18_NEIGHBORHOOD_NMF_K7_K12_EXTENSION_20261015_attempt5/models/neighborhood50_{K8,J2,L3}_K07_MODEL_AND_USAGE.npz`.
2. Frozen Reference18 Bin50 coordinates and patient metadata:
   `results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/{K8,J2,L3}/01_FROZEN_REFERENCE18_INPUT.tsv.gz`.
3. Frozen cross-chip K=7 family mapping:
   `results/reanalysis/REFERENCE18_NEIGHBORHOOD_NMF_K7_STATE_SENSITIVITY_20261016_210302_attempt5/07_CROSS_CHIP_FAMILY_MAP.tsv`.
4. Validated chip-specific ssDNA affine transforms:
   `results/ppt_qc/ssDNA_coordinate_audit.json`.
5. Original K8, J2 and L3 ssDNA TIFF images in the chip QC directories.

The frozen K=7 model has an approved local eLabFTW copy in
`elab/Harmonized_K7_Reference18_Neighborhood_NMF_robustness_261017165550.md`,
but no remote entry ID or URL exists. No remote provenance link was invented.
No remote entry was identified for the coordinate audit or raw ssDNA images.

## Script

Spatial mapping and figure generation:

```text
envs/stereo_seq_py/Scripts/python.exe scripts/1157_render_reference18_k7_ssdna_maps.py
```

Standalone HTML report:

```text
envs/stereo_seq_py/Scripts/python.exe scripts/1158_render_reference18_k7_ssdna_report.py --analysis-dir results/reanalysis/REFERENCE18_K7_SSDNA_MAPS_261017220004
```

Key safeguards:

- The K=7 NMF models, cross-chip family mapping and Reference18 inputs were
  reused without refitting.
- Disease labels were used only to arrange figure rows; no disease hypothesis
  test was run.
- Each Bin50 was transformed into raw ssDNA pixel space using the previously
  validated chip-specific inverse affine transform.
- The categorical map used the maximum relative-usage factor for each Bin50.
- The B-enriched/mixed-immune map used a fixed 0–1 continuous scale across all
  patients and chips.
- B-enriched/mixed-immune usage represents neighborhood exposure rather than
  absolute B-cell abundance.

Repository: `https://github.com/bozhu001/stereo_seq`

Analysis/report commit:
`2f7cf0b5241af4f927658052e5cde034ce08d858`

## Output summary

The workflow mapped 126,199 modelled Bin50 from 21 patients: six HC, six IPF
and nine SSc-ILD. All transformed Bin50 centers fell inside the corresponding
chip-specific ssDNA image bounds; the inside-image fraction was 1.0 for K8,
J2 and L3.

The final output includes:

1. A dominant-factor map containing the seven harmonized K=7 factor families
   across all patients on ssDNA, in PNG and PDF formats.
2. A continuous B-enriched/mixed-immune usage map on a fixed 0–1 scale, in PNG
   and PDF formats.

Final output directory:
`results/reanalysis/REFERENCE18_K7_SSDNA_MAPS_261017220004/`.

Rendered standalone HTML:
`results/reanalysis/REFERENCE18_K7_SSDNA_MAPS_261017220004/reports/REFERENCE18_K7_SSDNA_MAPS_REPORT_261017222140.html`.

Report SHA-256:
`4fd3ccd1117ebe490dcd9b7c1b86b0937c3a7d64f98875c6b542d86933a68a4e`.

The report is approximately 28.5 MB and should be attached when the entry is
created in the eLabFTW web UI.

## Remote-entry status

No remote ID or URL is recorded because `ELABFTW_API_KEY` is unavailable in
this environment. Create this approved entry through the eLabFTW web UI using
category `Bioinformatic`, the approved tags above, and attach the standalone
HTML report.
