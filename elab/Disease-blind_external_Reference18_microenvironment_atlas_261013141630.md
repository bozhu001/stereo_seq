# Disease-blind external Reference18 microenvironment atlas around non-B Neighborhood NMF factors — 261013141630

Status: approved by the user on 2026-10-13. Local final copy; remote eLabFTW
entry was not created because `ELABFTW_API_KEY` is unavailable in this
environment.

Category: Bioinformatic

Tags: `Stereo-seq`, `Reference18`, `Neighborhood NMF`, `Bin50`, `spatial`,
`microenvironment`, `Ciliated`, `Fibroblast`, `Lymphatic`

## Goal

Characterize the 100–200 µm external Reference18 environment surrounding the
frozen Ciliated, Fibroblast and Lymphatic per-chip K=6 Neighborhood NMF factor
families across 21 Stereo-seq lung samples, without using disease labels for
discovery or testing.

## Input

1. Frozen per-chip Reference18 Bin50 proportions, metadata and coordinates:
   `results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/{K8,J2,L3}/01_FROZEN_REFERENCE18_INPUT.tsv.gz`.
2. Frozen per-chip K=6 Neighborhood NMF models and continuous relative usages:
   `results/reanalysis/REFERENCE18_PER_CHIP_DUAL_NMF_20261012_204500_attempt4_alpha0/models/neighborhood50_{K8,J2,L3}_K06_MODEL_AND_USAGE.npz`.
3. Frozen factor-family mapping: Ciliated K8 F03/J2 F04/L3 F02;
   Fibroblast K8 F04/J2 F06/L3 F05; Lymphatic K8 F01/J2 F02/L3 F04.

Existing local eLabFTW copies mention some of the frozen per-chip Reference18
inputs, but no remote entry ID or URL is available for these exact combined
NMF inputs. No provenance link was invented.

## Script

Analysis command:

```text
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/1108_run_non_b_external_microenvironment_atlas.ps1
```

Analysis script:
`scripts/1108_build_non_b_external_microenvironment_atlas.py`.

Report source and render command:

```text
bash scripts/1108_render_non_b_external_microenvironment_atlas.sh
```

Key parameters and safeguards:

- Main annulus: 100–200 µm; 0–50 µm was an internal positive control and
  50–100/200–300 µm were radial sensitivities.
- Annuli were restricted to the same patient, ROI and tissue component.
- The primary annulus required at least 20 Bin50 and at least four of eight
  angular sectors.
- Associations were aggregated into 150, 200 and 250 µm spatial blocks.
- The primary patient-level statistic used 200 µm blocks and adjusted rank
  correlations for nUMI, nGene, local 50 µm neighbor count and annulus size.
- Patients, not Bin50s, were replication units. Disease labels were retained
  only as provenance and did not enter discovery, filtering or testing.
- RCTD and NMF were not rerun. The B-associated family and frozen 238
  lymphoid aggregates were not analyzed or modified.

Repository: `https://github.com/bozhu001/stereo_seq`

Analysis/report commit: `fbf009a`

## Output summary

Output directory:
`results/reanalysis/NON_B_EXTERNAL_MICROENVIRONMENT_ATLAS_20261013_140512/`

Rendered standalone report:
`results/1108_non_b_external_microenvironment_atlas_report_261013141630.html`
(SHA-256
`72251F97E4FABCC87D5E7B9AD2CA9B7297BE6D92DCE7F7413F98EB75726F3FD7`).
The report is under 100 MB and should be attached when the entry is created in
the eLabFTW web UI.

All 21 patients passed the primary analysis coverage requirement; the minimum
patient-level valid-center fraction in the 100–200 µm annulus was 0.919 and
the median was 0.994. Across 54 primary factor-family/state relationships,
nine were classified as cross-chip reproducible and seven as partially
reproducible.

The strongest findings were external continuity of Ciliated epithelial,
Fibroblast and Lymphatic endothelial states around their respective factors.
Each was directionally consistent in all 21 patients and all three chips.
Ciliated–NK was a weaker FDR-significant non-self positive relationship;
Ciliated–Secretory/other airway epithelial and Fibroblast–AT2/alveolar
transitional were partial positive candidates. Negative relationships were
retained only as closure-sensitive relative compositional depletion patterns,
not evidence of cellular exclusion.

These results support spatially continuous tissue domains. They do not by
themselves establish direct cell contact, ligand–receptor communication,
causality or disease effects. No disease comparison was performed.

## Remote-entry status

No remote ID or URL is recorded because `ELABFTW_API_KEY` is unavailable in
this environment. Create this approved entry through the eLabFTW web UI using
category `Bioinformatic`, the approved tags above, and attach the standalone
HTML report.
