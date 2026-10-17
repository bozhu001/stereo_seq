# Harmonized K=7 Reference18 Neighborhood NMF robustness — 261017165550

Status: approved by the user on 2026-10-17. Local final copy; remote eLabFTW
entry was not created because `ELABFTW_API_KEY` is unavailable in this
environment.

Category: Bioinformatic

Tags: `Stereo-seq`, `Reference18`, `Neighborhood NMF`, `Bin50`, `K7`,
`robustness`, `sensitivity`, `immune microenvironment`, `QC`

## Goal

Validate whether a common disease-blind K=7 Neighborhood NMF representation
can be used across K8, J2 and L3, and determine whether the seven factor
families remain compositionally and spatially robust to removal of weak or
potentially ambiguous Reference18 states.

## Input

1. Frozen Reference18 center-excluded 50-µm neighborhood assay for 21 patients:
   `results/reanalysis/ALL21_SHARED_UNBIASED_DISCRETE_NICHE_ATLAS_20261011_191652_attempt4/12_SHARED_NICHE_ASSAY.h5`.
2. Frozen per-chip Reference18 Bin50 inputs:
   `results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/{K8,J2,L3}/01_FROZEN_REFERENCE18_INPUT.tsv.gz`.
3. Persisted per-chip K=7 models from:
   `results/reanalysis/REFERENCE18_NEIGHBORHOOD_NMF_K7_K12_EXTENSION_20261015_attempt5/`.
4. K2–K12 restart, convergence, LOPO and adjacent-rank audit:
   `results/reanalysis/REFERENCE18_NEIGHBORHOOD_NMF_K2_K12_RANK_SELECTION_20261016_110753/`.
5. Reference18 marker-validation conclusions defining B as partially specific
   and Lymphatic endothelial/Basal-like epithelial as sensitivity-qualified.

Existing local eLabFTW copies mention some frozen Reference18 inputs, but no
remote entry ID or URL was identified for the exact combined K=7 inputs. No
provenance link was invented.

## Script

Primary state-removal sensitivity:

```text
/home/bozhu/miniconda3/envs/stereo/bin/python scripts/1143_reference18_neighborhood_nmf_k7_state_sensitivity.py --output-dir results/reanalysis/REFERENCE18_NEIGHBORHOOD_NMF_K7_STATE_SENSITIVITY_20261016_210302_attempt5 --n-jobs 4 --main-restarts 50 --lopo-restarts 3 --max-iter 5000
```

K8 convergence extension:

```text
/home/bozhu/miniconda3/envs/stereo/bin/python scripts/1144_k8_k7_state_sensitivity_convergence_extension.py --output-dir results/reanalysis/REFERENCE18_NEIGHBORHOOD_NMF_K8_K7_CONVERGENCE_EXTENSION_20261016_212100 --n-jobs 4 --restarts 50 --max-iter 20000
```

Independent evidence review and final figure:

```text
envs/stereo_seq_py/Scripts/python.exe scripts/1145_reference18_k7_sensitivity_evidence_review.py --output-dir results/reanalysis/REFERENCE18_NEIGHBORHOOD_NMF_K7_SENSITIVITY_EVIDENCE_REVIEW_20261016_212300
envs/stereo_seq_py/Scripts/python.exe scripts/1146_plot_reference18_k7_robustness_summary.py --output-dir results/reanalysis/REFERENCE18_K7_ROBUSTNESS_SUMMARY_20261017_165039
```

Report render:

```text
bash scripts/1148_render_reference18_k7_robustness_report.sh
```

Key locked parameters and safeguards:

- K=7 for K8, J2 and L3; J2/L3 K=12 was explicitly deferred.
- 50 NMF restarts and three LOPO restarts per state-removal model.
- 495 balanced observations per patient and seven patients per chip.
- Unpenalized Frobenius objective with coordinate descent.
- Sensitivity A removed Basal-like epithelial and Lymphatic endothelial;
  sensitivity B removed B; rows were renormalized before refitting.
- Retained factors were reviewed with independent best matches rather than a
  forced one-to-one assignment after one axis disappeared.
- Disease labels were not used in model fitting, rank review, matching or
  sensitivity analysis.

Repository: `https://github.com/bozhu001/stereo_seq`

Analysis/report commit:
`f2ce760d6b2c0ce85905ce0b88815969cc272864`

## Output summary

All seven primary K=7 factor families were shared across K8, J2 and L3. The
minimum pairwise loading cosine was 0.912. The shared families were
LEC-associated, Blood-endothelial, AT1-associated, Perivascular-stromal,
B-enriched/mixed-immune, Ciliated-associated and Fibroblast-associated.

After removal of Basal plus LEC or removal of B, all six families expected to
remain were compositionally and spatially persistent. Across chips and both
sensitivity models, the minimum independent loading cosine was 0.926, the
minimum median-patient Bin50 score Spearman correlation was 0.852, and the
minimum median-patient Q90 hotspot Jaccard was 0.526.

The K8 20,000-iteration extension converged in all 50 restarts for both
state-removal models. Its minimum matched 5,000-to-20,000-iteration loading
cosine was 0.9986, excluding incomplete optimization as the explanation for
the K8 solution.

After B removal, the original B-enriched factor did not map to one coherent
T/NK successor: the maximum median-patient correlation across all candidate
successors was 0.257 in K8, 0.056 in J2 and 0.125 in L3; the corresponding
maximum Q90 hotspot Jaccard was below 0.076 in every chip.

The frozen interpretation retains `B-enriched/mixed-immune` as a primary
immune-focused neighborhood factor for downstream spatial, patient-level and
molecular characterization. Its score represents an immune neighborhood
composition axis and is not interpreted as pure B-cell abundance. LEC remains
sensitivity-qualified because its upstream RCTD marker support was weak. The
absence of separate T, NK, Mast, Basal or AT2 factors at K=7 is not interpreted
as absence of those biological signals.

Final figure directory:
`results/reanalysis/REFERENCE18_K7_ROBUSTNESS_SUMMARY_20261017_165039/`.

Rendered standalone HTML:
`results/reanalysis/REFERENCE18_K7_ROBUSTNESS_REPORT_20261017_165550/REFERENCE18_K7_ROBUSTNESS_REPORT_261017165550.html`.

Report SHA-256:
`25875d2c41ca164542338b644d8b6a1a069251e8a4d94a8f0feae2f57f931cca`.
The report is under 100 MB and should be attached when the entry is created in
the eLabFTW web UI.

## Remote-entry status

No remote ID or URL is recorded because `ELABFTW_API_KEY` is unavailable in
this environment. Create this approved entry through the eLabFTW web UI using
category `Bioinformatic`, the approved tags above, and attach the standalone
HTML report.
