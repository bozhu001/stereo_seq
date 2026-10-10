# Project memory

## 2026-07-30 — BIN50 top3000_atleast6_1725HVG sensitivity

- Used the existing 21-sample Top3000 HVG consensus; sample-wise HVGs were not
  recomputed.
- Selected `number_of_samples >= 6`: 1,725 unique Ensembl IDs, all matched to
  the fixed G15 raw-count AnnData, with no duplicate or unmatched IDs.
- Fixed workflow: normalize_total(10,000), log1p, sparse-preserving scale,
  PCA50, Harmony PC1-PC50 by `sample_id`, neighbors on Harmony Dim1-Dim12
  (`n_neighbors=30`, Euclidean, seed 0), UMAP, and Leiden resolutions
  0.5/0.8/1.0/1.5/2.0.
- Compared with the existing 3,508-HVG Dim1-Dim12 scheme. The 1,725-HVG graph
  was more fragmented (7 components; 177 bins outside the main component),
  had more quantitative UMAP peripheral bins (2.0768% versus 0.3367%), many
  more <100-bin micro-clusters, and lower mean consecutive-resolution ARI.
  Recommendation: retain the 3,508-HVG scheme.
- Final H5AD passed write-after-read checks for X, PCA, Harmony, UMAP,
  connectivities, distances, shape, and all five Leiden label columns.
- Session summary:
  `.agent/bin50_top3000_atleast6_1725hvg_sensitivity_260730130607.md`.
- eLabFTW entry not created: separate user approval is still required.
- No Git operation was performed, per explicit user instruction.

## 2026-07-31 — BIN50 cNMF + Hotspot methodological pilot

- Audited and reused the fixed G15 Bin50 raw-count object (152,640 bins ×
  27,618 common genes; 21 samples) and the existing metadata/Leiden object.
  No HVG, PCA, Harmony, neighbors, UMAP, Leiden, RCTD, or Git operation was
  run.
- Used official cNMF 1.7.1 in a project-isolated environment. Rank benchmark
  used 28,431 sample-balanced bins, 3,508 recurrent HVGs, K=8/12/16/20 and
  20 initializations per K. K=12 was selected (stability 0.7477); K=16 is the
  adjacent higher-rank sensitivity.
- Final K=12/K=16 consensus used 53,677 sample-balanced bins and 30
  initializations per K. Fixed spectra were projected onto all 152,640 bins.
- Seven of 12 K=12 programs passed the defined 30-seed, multisample,
  nontechnical cNMF reproducibility gate. Provisional post-hoc candidates
  included airway/ciliated/club (P04), AT2/alveolar (P02), ECM/collagen (P01),
  and immunoglobulin/B-plasma-like programs (P07/P11); these are not formal
  cell-type or niche labels.
- Hotspot 1.1.3 completed for all 21 samples at physical-space k=8 and k=12
  (42/42 runs). Strict modules formed in only one sample. No cNMF program had
  strong or moderate cross-sample Hotspot validation; six had weak,
  single-sample support.
- Final decision: `NEED_SPATIAL_MODEL`. Do not proceed to formal program
  co-occurrence or niche analysis on this evidence alone.
- Results:
  `results/reanalysis/bin50_cNMF_Hotspot_pilot_20260731_214558/`.
- Rendered HTML:
  `BIN50_cNMF_Hotspot_pilot_report.html` in the result directory.
- Session summary:
  `.agent/bin50_cnmf_hotspot_pilot_260731224758.md`.
- eLabFTW entry was intentionally not created, per explicit user instruction.
  The bespoke workflow remains project-specific; no shared-skill issue was
  filed. No Git operation was performed.

## 2026-07-31 — BIN50 2274 QC-adjusted r0.3/r0.5 Level1 review

- Reused the completed 2,274-feature object, all 62 saved Level1 program
  scores, saved r=0.3/r=0.5 labels, existing r=0.3 marker/enrichment results,
  and spatial coordinates. No upstream analysis was rerun and all 152,640
  bins were retained.
- `log1p_n_genes` and `log1p_total_counts` were almost collinear (Pearson
  0.9983; VIF 299.15), so the primary residual model used
  `log1p_n_genes`, the counts model was sensitivity-only, and the joint model
  was not used.
- r=0.3 cluster 0 remained unclear; cluster 1 retained airway support;
  cluster 2 retained AT2-rich mixed epithelial-stromal support; cluster 3
  was revised from myofibroblast-rich to mixed immune-stromal because strong
  stromal marker enrichment persisted while adjusted immune programs rose.
- r=0.5 contained 7 major, 2 minor, and 6 micro clusters. Cluster 2 provided
  a reproducible airway split, cluster 4 an AT2-rich mixed split, and cluster
  7 a replicated mixed immune-stromal child of r=0.3 cluster 3. No cluster
  passed the strict stable immune-enriched criteria.
- Recommendation: keep r=0.3 as the main annotation level and use r=0.5 only
  for selected substructure review.
- Results:
  `results/reanalysis/bin50_2274_qc_adjusted_r03_r05_review_260731120038/`.
- Session summary:
  `.agent/bin50_2274_qc_adjusted_r03_r05_review_260731122139.md`.
- The user explicitly confirmed that no eLabFTW entry is needed for now.
- The QC-adjusted program-score/r0.3-r0.5 review workflow remains
  project-specific pending fuller validation; no shared-skill candidate or
  bioinfo-wiki issue was created.
- No Git operation was performed, per explicit user instruction.

## 2026-07-30 — Level1 Rule B additions for 1725 vs 3508 HVGs

- Reused the completed 874-marker Level1 variability assessment and existing
  1,725/3,508 recurrent-HVG gene lists. No H5AD was opened and no HVG, PCA,
  Harmony, neighbors, UMAP, Leiden, or Git operation was run.
- Both schemes used Rule B: detection in at least three samples plus either
  Top3000 recurrence in at least two samples or global normalized-variance
  percentile <=0.15.
- Of 874 deduplicated Level1 markers, 855 were in AnnData and 19 were absent.
- The 1,725-HVG baseline already contained 132 markers; 723 were outside it,
  549 passed Rule B and were added, 174 failed Rule B, and the final feature
  set contained 2,274 Ensembl IDs.
- The 3,508-HVG baseline already contained 223 markers; 632 were outside it,
  458 passed Rule B and were added, 174 failed, and the final feature set
  contained 3,966 IDs.
- New-marker overlap: 458 shared, 91 required only by the 1,725 scheme, and
  zero required only by the 3,508 scheme.
- For two marker symbols with multiple mapped Ensembl IDs, the representative
  ID already selected by the completed variability assessment was used so
  each deduplicated marker contributed exactly one feature.
- Results:
  `results/reanalysis/bin50_level1_ruleB_1725_vs_3508_260730182251/`.
- Session summary:
  `.agent/bin50_level1_ruleb_1725_vs_3508_260730182251.md`.
- eLabFTW entry not created: separate user approval is still required.
- No Git operation was performed.

## 2026-07-30 — BIN50 Level1 marker supplementation review

- Reviewed the newest lung-marker workbook in `bin50_input`:
  `BI_PF_ILD_atlas_marker_gene_table_260709213053.xlsx`, sheet `Sheet 1`.
  Level1 was defined exactly as `level == "cell_type"`.
- Level1 contained 1,550 rows, 31 cell types, and 874 unique marker symbols.
  There were no blank/nonstandard marker rows or repeated marker/celltype rows;
  322 markers were shared across multiple Level1 cell types.
- Reused the fixed 3,508-HVG list and all 21 saved sample-wise Top3000
  checkpoints (63,000 records). No HVG, PCA, Harmony, neighbors, UMAP, or
  Leiden analysis was rerun.
- Of 874 Level1 markers, 223 were already in the 3,508-HVG set, 632 were
  present in AnnData but outside the baseline, and 19 were absent from the
  AnnData gene space. Two symbols mapped to multiple Ensembl IDs; all IDs were
  retained and explicitly flagged.
- Candidate additions: Rule A 556 (final 4,064 features), Rule B 458 (final
  3,966), and Rule C 342 (final 3,850). Rule B is the primary recommendation
  because it raises mean per-cell-type coverage to 75.7% while 91.7% of new
  markers have at least two-sample Top3000 support; Rule A is the sensitivity
  scheme.
- Results:
  `results/reanalysis/bin50_level1_marker_supplementation_review_260730175152/`.
- Main HTML:
  `BIN50_Level1_marker_supplementation_review.html`; timestamped copy:
  `BIN50_Level1_marker_supplementation_review_260730175152.html`.
- Session summary:
  `.agent/bin50_level1_marker_supplementation_review_260730175152.md`.
- eLabFTW entry not created: separate user approval is still required.
- No Git operation was performed, per explicit user instruction.

## 2026-07-30 — BIN50 3508-HVG Leiden r0.3/r0.5 cluster review

- Read the completed 3,508-HVG Dim1-Dim12 AnnData without modifying it and
  reused the existing PCA, Harmony, neighbors, UMAP, and Leiden r0.3/r0.5
  results.
- Input shape was 152,640 bins by 27,618 genes. Marker testing used the full
  normalize_total/log1p `X` matrix (`raw` absent), Wilcoxon cluster-versus-rest,
  with 50 unfiltered and 50 technical-gene-filtered markers per cluster.
- r0.3 contained 3 major, 4 minor, and 6 micro clusters; r0.5 contained
  4 major, 3 minor, and 6 micro clusters. No robust low-QC-driven cluster was
  detected.
- r0.3 clusters 2, 3, 4, 5, 6, 7, 8, 10, and 12 remained stable at r0.5.
  Cluster 1 was clearly split, but no new r0.5 child met the combined size,
  multi-sample, multi-chip, marker-distinctness, spatial-reproducibility, and
  parent-purity criteria.
- Current evidence supports r0.3 as the main resolution. r0.5 adds
  fragmentation without enough replicated, marker-distinct, spatially
  reproducible major structure to justify it.
- Results:
  `results/reanalysis/bin50_top3000_atleast5_3508HVG_r0_3_r0_5_cluster_review_260730160546/`.
- Rendered HTML:
  `BIN50_3508HVG_Leiden_r0.3_r0.5_cluster_review_260730164200.html` within the
  result directory.
- Session summary:
  `.agent/bin50_3508hvg_r03_r05_cluster_review_260730164200.md`.
- eLabFTW entry not created: separate user approval is still required.
- No Git operation was performed, per explicit user instruction.

## 2026-07-31 — BIN50 K=12 cNMF direct spatial validation

- Reused the completed K=12 cNMF usage for all 152,640 fixed G15 Bin50 and
  existing metadata. No cNMF, Hotspot, HVG, PCA, Harmony, neighbors, UMAP,
  Leiden, RCTD, cell2location, mNSF, or Git operation was run.
- Input audit passed: 152,640 exact ordered unique Bin IDs, 21 samples,
  complete physical coordinates, and nonnegative usages.
- Fit sample-specific QC models for all 12 programs: usage ~ log1p_n_genes
  (primary) and usage ~ log1p_total_counts (sensitivity). Matching residual
  correlations were effectively zero (maximum absolute Pearson 4.6e-9).
- Built independent physical-coordinate kNN graphs per sample at k=8 and k=12.
  Moran's I and Geary's C used 999 fixed-seed permutations. All 21 samples
  succeeded; 1,512 Moran rows, 1,512 Geary rows, and 1,512 patch rows were
  produced.
- Stable programs P01, P02, P03, P04, and P07 received moderate direct spatial
  support; P06 was weak and P09 had none. No stable program met the Strong
  definition because compact predefined Q90 patches were not replicated across
  multiple samples.
- P07 and P11 had highly correlated spectra (r=0.924) but uncorrelated usages
  (r=-0.065), low patch overlap, and spatial support only for P07. They were
  classified as related but distinct; P11 remains exploratory.
- Hotspot module absence does not imply zero spatial structure: five stable
  usages retained sample-wise QC-adjusted Moran support. Evidence is more
  consistent with broad spatial gradients plus stringent/mismatched gene-module
  construction than with replicated compact niches.
- Final decision: `NEED_MNSF_PILOT`. Formal program co-occurrence and local
  niche analysis are not yet approved.
- Results:
  `results/reanalysis/bin50_cNMF_all12_direct_spatial_validation_20260731_231529/`.
- Session summary:
  `.agent/bin50_cnmf_all12_direct_spatial_260731233049.md`.
- No eLabFTW, shared-skill filing, or Git operation was performed, following
  explicit user direction.

## 2026-08-01 - BIN50 mNSF technical stop and P07 spatial characterization

- Formally stopped the mNSF route with decision
  `STOP_MNSF_TECHNICAL_ROUTE`: official synthetic smoke test passed, but the
  six-sample CPU pilot had slow initialization and repeated Inf/NaN divergence;
  no K=8 seed completed and K=12 never started. All caches, checkpoints, logs,
  recovery objects, loss plots, and failure records were preserved.
- Reused saved K=12 cNMF P07 usage, per-sample QC residuals, permutation Moran
  results, and patch metrics for all 152,640 fixed G15 Bin50. No upstream model
  or source object was changed.
- P07 phenotypes across 21 samples: 2 focal_high (HC/NL-72 and
  IPF/FO22-1-09404), 0 diffuse_high, 17 scattered, 2 low, and 0 QC_confounded.
- P11 provided no significant positive n_genes-adjusted k=8 Moran support and
  remained a sensitivity-only unstable program.
- Prioritized 10 review ROIs spanning ECM, airway, lymphoid, independent
  P07-high, and P07-low control strata. Seven include Bin20 and eight include
  CellBin review; decision `GO_BOTH_BIN20_AND_CELLBIN`, pending morphology QA.
- Results:
  `results/reanalysis/bin50_P07_Bplasma_spatial_characterization_20260801_194219/`.
- Session summary:
  `.agent/bin50_mnsf_stop_p07_spatial_260801194315.md`.
- No eLabFTW, shared-skill filing, or Git operation was performed, following
  explicit user direction.

## 2026-08-22 — BIN50 frozen-program rook adjacency sensitivity wrap-up

- Continued from the saved rook-adjacency sensitivity output without rerunning
  statistical, GSEA, differential-expression, scoring, or spatial analysis.
- The completed sensitivity changed only the graph from queen eight-neighbor
  adjacency to rook four-neighbor adjacency with a 50-native-unit orthogonal
  step; seed 20260822, 999 permutations, `z >= 1`, minimum 4-bin hotspots, and
  0.000625 mm2 Bin50 area remained fixed.
- All three final QC domains passed. Rook versus queen changed 0/5 formal
  program spatial-QC calls and 0/45 adjacency-related disease FDR conclusions.
  `B_core` remained failed because coherent multigene support was absent under
  both graphs. J2 retained identity orientation.
- Dynamic report source:
  `scripts/382_rook_adjacency_sensitivity_report.qmd`; render wrapper:
  `scripts/383_render_rook_adjacency_report.ps1`.
- Final standalone HTML:
  `results/reanalysis/bin50_all21_patient_pseudobulk_disease_signal_audit_260822/rook_adjacency_sensitivity_260822200142/rook_adjacency_sensitivity_report_260822212549.html`
  (20,871,077 bytes; SHA-256
  `BEAA483BA274A1316EEFC81B8BC0AE0EA7EB05886EFBF756B6EE9B6DBCDE984E`).
- The user approved the scientific report and the separate ELN draft. Local ELN
  copy:
  `elab/Bin50_frozen_program_rook_adjacency_sensitivity_audit_260822212549.md`.
  Remote entry ID/URL was not created because `ELABFTW_API_KEY` is unavailable
  and the user requested no remote operation. Project tag remains
  `TODO_PROJECT_TAG`; no prior input entry IDs are associated.
- Scripts 380/381 and the rook review were classified by the user as a one-off
  audit, so no recurring-workflow skill candidate or GitHub issue was created.
- Analysis/report commit: `5ac59c94254242bb4fc667e2fa9e36c7d0e2316c`.
  No push was performed per user instruction.
- Session summary:
  `.agent/rook_adjacency_sensitivity_wrapup_260822215136.md`.

## 2026-08-23 — Cleaned B-lineage patient-level disease comparison

- Reused frozen patient-level outputs from
  `results/reanalysis/bin50_hq_frozen_blineage_descriptive_qc_260823180725/`;
  no Bin-level scoring, marker selection or threshold optimization was rerun.
- Primary15 contained HC 4, IPF 4 and SSc-ILD 7 after the five independent
  low-quality-watch exclusions plus extreme-depth `IPF/FO22-1-09404`.
  J2/L3-only13 further removed the two remaining K8 samples.
- Six frozen panels used corrected-logit strict same-bin fractions in the
  patient-level model `disease_group + log(sample_total_counts) + chip_id`,
  with HC3 covariance and one BH adjustment across 18 contrasts.
- All 18/18 confidence intervals crossed zero and 0/18 reached BH-FDR < 0.05.
  SSc-ILD versus IPF was negative for all six complete-cohort models in both
  cohorts, but IgA had LOO direction flips; this remains an exploratory trend,
  not an established disease difference.
- Final HTML:
  `results/reanalysis/bin50_cleaned_blineage_patient_disease_comparison_260823192715/report/FINAL_394_cleaned_blineage_patient_disease_comparison_report_260823192715.html`
  (3,458,257 bytes; SHA-256
  `0CE6E37D4B98CF6358F65772D5129449C951C8F5AB835C6A2116869FEC2F44CC`).
- Analysis/report commit:
  `446a30841cebf16d09c9094a3f251b6fff916679`.
- Local ELN draft:
  `elab/Bin50_cleaned_Blineage_patient_disease_comparison_260823192715.md`;
  separate ELN approval remains pending and no remote entry was created.
- Session summary:
  `.agent/blineage_patient_disease_comparison_wrapup_260823200656.md`.
- No push was performed.

## 2026-08-23 — L3 Mature_B_core depth/PTPRC-matched spatial null

- Reused persisted frozen same-bin and rook support flags for all six fixed L3
  samples; no raw-object read, rescoring, marker/threshold change or
  expression-based sample selection was performed.
- Ran 1,000 exact stratum-matched random draws per sample and support mode,
  matching sample-specific `log1p(total_counts)` quintile, `n_genes` quintile
  and PTPRC raw-positive/raw-zero status. All 12 combinations passed exact
  composition validation.
- Final evidence hierarchy: strict same-bin is primary. Strict rook clustering
  is construction-dependent auxiliary evidence because adjacency is part of
  the support definition; it is not an independent spatial validation.
- `IPF/FO23-1-06168` and `SSC/15491/14` are retained as
  depth/n_genes/PTPRC-status–matched focal multi-marker B-lineage candidates:
  one patient-specific IPF candidate and one patient-specific SSc-ILD
  candidate. Their same-bin clustering was not fully explained by the matched
  technical/immune background.
- The candidates are not B-cell niches, aggregates, TLS or disease-specific
  structures; the technical wording is limited to “not fully explained by the
  matched technical/immune background,” with no disease difference claimed.
- Final HTML:
  `results/reanalysis/bin50_l3_mature_b_depth_ptprc_matched_spatial_null_260823211809/report/FINAL_398_l3_mature_b_depth_ptprc_matched_spatial_null_report_260823211809.html`
  (2,932,894 bytes; SHA-256
  `D32E7F4E5E854794F38B227D5D393767A7C6AE966FA6854242385401CFFC90B4`).
- Analysis/report commit:
  `f4a09b9959ad50f88ca4a097fee37d7d391eaa7a`.
- Local ELN draft:
  `elab/L3_Mature_B_core_depth_PTPRC_matched_spatial_null_260823211809.md`;
  separate ELN approval remains pending and no remote entry was created.
- Session summary:
  `.agent/l3_mature_b_depth_ptprc_matched_spatial_null_wrapup_260823214027.md`.
- No push was performed.

## 2026-08-23 — B-lineage analysis-line synthesis

- Integrated six existing B-lineage stages without running any new statistic,
  score, randomization or spatial analysis. MS4A4A was excluded.
- Confirmed technical/descriptive capture feasibility and strong depth/chip
  effects, including FO22 upper-extreme influence and residual K8-high capture.
- The 15/13-patient disease comparison remained negative: all 18 confidence
  intervals crossed zero and 0/18 reached BH-FDR below 0.05.
- `IPF/FO23-1-06168` and `SSC/15491/14` remain patient-specific
  depth/n_genes/PTPRC-status–matched focal multi-marker B-lineage candidates.
  Strict same-bin is primary; rook clustering is construction-dependent
  auxiliary evidence and not independent spatial validation.
- Added a six-row exact frozen-rule table. The 100th empirical percentile is
  limited to the existing 1,000 matched draws and is not a disease P value.
- H&E review is conditional on obtaining corresponding tissue and achieving
  reliable outline or coarse-region registration.
- Niche, aggregate, TLS and disease-specific structure language remains
  prohibited.
- Final HTML:
  `results/reanalysis/bin50_blineage_analysis_line_synthesis_260823215304/report/FINAL_399_bin50_blineage_analysis_line_synthesis_report_260823215304.html`
  (13,794,494 bytes; SHA-256
  `B97384DE5636CCF29F59D607CD4768AF139FB2A281DB9851A1A754CB72FB52F5`).
- Analysis/report commit:
  `45f2672f14ec7dacce75930273241236e72f4b7b`.
- ELN draft:
  `elab/DRAFT_Bin50_Blineage_analysis_line_synthesis_260823215304.md`;
  separate approval remains pending and no remote entry was created.
- Session summary:
  `.agent/blineage_analysis_line_synthesis_wrapup_260823222413.md`.
- One-off classification retained; no skill or GitHub issue. No push.

## 2026-08-28 — L3 7/7 BCR/CDR3 finalization

- Completed L3 for 7/7 patients from existing validated candidate inputs and
  completed TRUST4 outputs. No patient or original FASTQ was rerun, K8 was not
  scanned, and no PPT was created.
- Verified exactly one completed TRUST4 attempt per patient with nonempty AIRR,
  barcode AIRR, assembled-read, assignment and cached-R2 files plus completion
  markers. The interrupted `SSC_05957_17B` directory remains preserved; its
  completed attempt is `SSC_05957_17B_resume_260828_111831`.
- Final L3 counts: 65 quality-corrected productive BCR/CDR3 molecules and 55
  unique productive clonotypes: 26 singleton, 22 technical, 6 moderate and one
  formal reliable-rule hit. Seven clonotypes had at least two corrected CID-UMI
  molecules; maximum corrected expansion was three molecules.
- `SSC/05957/17B` had no repeated clonotype. `SSC/07998/15A` had one
  three-molecule IGK moderate candidate with 2 CID and 2 Bin, not a reliable
  expanded clonotype.
- `SSC_15491_14_TRUST4_0003` is restricted to **provisional expanded IGH
  candidate**, never confirmed clonal expansion: although it had 3 corrected
  molecules, 3 CID, 3 Bin and sufficient fragment support, all 6 supporting
  reads carried the anomalously frequent `CGCTTGGCCT` motif and its exact CDR3
  had a cross-patient warning.
- Identical-threshold comparison: J2 had 94 corrected molecules, 68 clonotypes,
  9 repeated corrected-molecule clonotypes and zero reliable-rule hits; L3 had
  65, 55, 7 and one provisional formal hit. Approved overall conclusion: the
  current unenriched Stereo-seq data provide no reliable evidence of B-cell
  clonal expansion in J2 or L3, but do not prove biological expansion is absent.
- Final report:
  `results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934/L3_ALL7_FINAL_METRICS_REPORT_260828132341.html`
  (SHA-256
  `1B4F2942D635C4D15762F5086F742B59F590A5723D1278CB417C65A0BF2F8C9E`).
- Analysis/report commit: `c4897b8`.
- Approved local ELN draft:
  `elab/DRAFT_L3_7of7_BCR_CDR3_finalization_260828132341.md`. It was not
  submitted because the real eLabFTW project tag is not yet confirmed; the tag
  is omitted and no placeholder is used.
- Scripts 453/454 are retained as a project-reusable audit workflow. Per user
  decision, no general skill-candidate issue was created.
- Session summary:
  `.agent/l3_all7_bcr_cdr3_finalization_260828133100.md`.
## 2026-08-28 — Fixed marker-based B-cell-containing immune aggregate rule

- The project object is a **B-cell-containing immune aggregate**: reliable
  conventional B-lineage marker support is required, but B cells need not be
  the majority. Pure B-rich, B+T/NK, B+Macrophage/myeloid and mixed
  B+T/NK+Macrophage/myeloid regions are eligible.
- T/NK-only or Macrophage/myeloid-only regions are retained as non-B immune
  controls and excluded from primary B-aggregate counts. Candidates must never
  be called TLS or interpreted as confirmed B-cell clonal expansion.
- Raw integer marker counts plus within-sample depth correction define
  candidate bins and regions. Auxiliary HLA/CD74/LYZ/cathepsin/B2M or
  immunoglobulin signals cannot establish a B aggregate alone.
- Primary regions require at least three connected candidate-positive Bin50 on
  the eight-neighbour graph, at least one high-confidence bin, at least three
  distinct immune-core markers, and immune enrichment BH-FDR q<0.05 against at
  least 1,000 same-size, same-sample, UMI-matched random connected regions.
- Primary B eligibility requires either spatially supported conventional B
  evidence in at least two adjacent bins or a focal multi-marker B bin embedded
  in at least two connected T/NK or Macrophage/myeloid-supported bins, plus
  region-level B enrichment BH-FDR q<0.05.
- Report size as Bin50 count, area (`n_bins × 625 µm²`) and equivalent
  diameter; do not report estimated immune-cell count or size categories.
- Coarse14 is composition-only after marker-defined regions are fixed. The 23
  prior abundance-threshold calls are relabelled **depth-confounded preliminary
  abundance-based calls** and used only for overlap comparison.
- Full operational definition:
  `results/reanalysis/marker_based_B_aggregate_pilot_260828_220532/B_AGGREGATE_DEFINITION.md`.

## 2026-08-28 — Marker-based B-aggregate three-sample pilot result

- The valid computational run is
  `results/reanalysis/marker_based_B_aggregate_pilot_260828_222020_attempt2/`;
  it completed with exit code 0 and final status `PASS_WITH_DEPTH_CAUTION`.
  The earlier `..._220532/` attempt is preserved and explicitly marked invalid
  because its signed residual implementation increased depth dependence and
  its random-region matching was incomplete.
- The main method used raw integer counts and within-sample Poisson-offset
  residual enrichment (`>2`) without a forced top-percentile positive call.
  Median absolute marker-depth Spearman changed from 0.175 raw to 0.168 after
  correction. Residual depth dependence therefore remains, especially for
  IPF macrophage/combined signal, but it is far below the prior Coarse14
  abundance-depth confounding.
- Of 101 tested connected regions, 93 passed immune-region BH-FDR q<0.05 and
  25 met the stored non-B immune aggregate control definition. No region passed
  B-region BH-FDR q<0.05 (minimum B q=0.0865), so the main definition and all
  four sensitivity variants produced zero B-cell-containing immune aggregates,
  zero robust calls and zero B-focal mixed calls.
- Five near-candidate regions had at least two conventional B-core markers with
  adjacent B-marker bins, but all failed the B-region FDR gate. Four regions
  were plasma-associated; none had qualifying conventional-B support. The 936
  B-low rows are evidence-insufficient records, not 936 independent regions.
- Ten of the old 23 depth-confounded preliminary abundance-based calls
  overlapped marker-supported immune regions (HC 6, IPF 1, SSC 3); none
  overlapped a primary B aggregate because no primary B aggregate passed.
- Do not expand this exact calling rule to 21 samples yet. The zero-call result
  indicates limited B-marker power and/or an overly stringent region-level
  B-FDR gate at Bin50. Do not relax it silently; any calibration or change of
  inferential unit requires an explicitly approved follow-up analysis.
- No Cell2location/reference training, NMF, K8 BCR/CDR3 scan or 21-sample
  expansion was run. Coarse14 was used only for post-call composition
  description, never for aggregate eligibility or estimated cell counts.
- Recovery/session summary:
  `.agent/marker_based_b_aggregate_pilot_recovery_260828224951.md`.

## 2026-08-28 — Marker-based B-aggregate locked validation

- Reclassified fixed attempt2 immune aggregates without rereading raw counts
  or recomputing discovery/FDR. After retaining the unchanged >=3-Bin immune
  aggregate gate, development calls were two HC/NL-55 B-containing,
  B-suggestive regions and zero in the other two pilot samples.
- Locked definition SHA-256:
  `41D7C3318AF89D385CD8D76DBDDE65E8C2C77B30C22A51554698FB63EC8814D9`.
  Immune q defines immune aggregate status, spatial multi-marker B evidence
  defines B-containing status, and B q only defines enrichment class.
- Independent raw-count validation used IPF/FO23-1-06168, SSC/15491/14 and
  same-chip B-low control IPF/FO23-1-06170. B-containing counts were 3, 9 and
  3; B-enriched counts were 2, 6 and 3, respectively.
- Old strict-B overlap into new B-containing regions was 7/37 for
  IPF/FO23-1-06168, 17/65 for SSC/15491/14 and 0/8 for the control.
- Median absolute marker-depth Spearman was 0.1923 raw and 0.1841 adjusted;
  severe depth recurrence was not detected. Positive-sample counts at the
  >=2/>=3/>=4 Bin gates were 14/12/8.
- Do not expand to 21 samples: the B-low control produced three
  B-containing/B-enriched calls and the size-sensitivity stability criterion
  failed. The locked definition was not modified.
- Results:
  `results/reanalysis/marker_based_B_aggregate_locked_validation_260828_230422/`.
- Session summary:
  `.agent/marker_based_b_aggregate_locked_validation_260828234038.md`.
- No GitHub issue, Git commit/push or eLabFTW operation was performed.

## 2026-08-29 — Locked marker-based B-aggregate all21 analysis

- Completed all 21 Bin50 samples (9 SSc-ILD, 6 IPF, 6 HC; K8/J2/L3) using the
  frozen marker panel and locked definition. Final results are under
  `results/reanalysis/marker_based_B_aggregate_all21_locked_260828_235537/final_all21/`.
- Preserved the first 15 completed upstream results. Corrected all 45 regions
  previously sampled with replacement using unique UMI-matched connected
  backgrounds without replacement; no >=3-Bin primary region had fewer than
  100 unique backgrounds.
- Resumed only samples 16-21. One HC/NL-72 two-Bin sensitivity region retained
  `insufficient_matching_background` with NA P/q; no primary matching failure
  occurred.
- Final counts: 152,640 Bin rows, 1,038 tested regions, 551 immune aggregates,
  and 195 B-containing aggregates (82 enriched, 72 suggestive, 41 present but
  not enriched). Disease totals were HC 56, IPF 59 and SSc-ILD 80.
- Patient-level median B-containing density/mm2 was HC 1.3256, IPF 0.7833 and
  SSc-ILD 0.1785. No disease or descriptive chip-effect comparison survived
  BH correction. K8 nevertheless contributed 151/195 calls and should remain
  a descriptive chip-concentration caution.
- Sensitivity totals were 239/195/145 at >=2/>=3/>=4 Bin; all primary calls
  matched >=2 and 145/195 also met >=4.
- Median absolute depth rho was 0.1731 raw and 0.1609 adjusted; maximum
  adjusted absolute rho was 0.3999 (<0.5 stop threshold).
- Final QA passed, including 21 sample maps, unique region IDs, matching state,
  reports and required tables. No Cell2location, NMF, BCR, K8 scan or CDR3
  overlay ran.
- Session summary:
  `.agent/marker_based_b_aggregate_all21_locked_260829010557.md`.
- No GitHub issue, Git commit/push or eLabFTW operation was performed, per the
  user's explicit instruction.

## 2026-08-29 — B-marker specificity audit of fixed all21 aggregates

- Held all 551 marker-supported immune aggregates and their Bin memberships
  fixed. Recomputed B evidence using only CD79A, CD79B, MS4A1, CD22 and CD19;
  CD37/CD83 were supportive-only and excluded from score, spatial evidence and
  B-region inference.
- Of the original 195 B-containing calls, 73 retained anchor-confirmed status:
  16 Primary B-enriched, 36 B-suggestive and 21 B-present_not_enriched.
  Seventeen original calls were downgraded specifically as CD37/CD83-only.
- Size labels are now 73 `primary_3bin` and 59 `stable_4bin`; the previous
  >=2-based `robust` terminology is not used in this audit.
- K8/J2/L3 anchor-confirmed counts were 56/4/13 and enriched counts 13/0/3.
  K8 anchor density remained 11.86-fold the mean J2/L3 density, but K8 also had
  7.82-fold UMI/mm2 and 7.22-fold median n_genes/Bin; patient-level density
  correlated 0.804 with UMI/mm2 and 0.761 with n_genes. All five anchor rates
  were broadly elevated on K8. Treat K8 concentration as strongly depth/chip
  entangled, not unadjusted disease biology.
- Final valid output:
  `results/reanalysis/B_marker_specificity_audit_260829_013112_attempt2/`.
  QA passed for 551 unchanged regions, 21 samples, 152,640 Bin rows, unique
  no-replacement backgrounds and 30 manual review figures.
- Disease comparison is only conditionally suitable with patient-level chip
  stratification/depth control. CDR3 overlay is suitable for anchor-confirmed,
  especially enriched stable_4bin regions, but was not run.
- Session summary:
  `.agent/b_marker_specificity_audit_all21_260829015234.md`.
- No immune discovery, Cell2location, NMF, BCR, K8 scan, CDR3 overlay, GitHub,
  Git commit/push or eLabFTW operation was performed.

## 2026-08-29 — L3 CDR3 overlap with fixed marker-based B aggregates

- Used the completed Y40105L3 all-seven quality-corrected CDR3 support and the
  immutable `B_marker_specificity_audit_260829_013112_attempt2` boundaries.
  No FASTQ, BCR caller, Cell2location, aggregate discovery or anchor-panel
  modification was performed.
- Reconstructed 65 corrected molecules from 462 supporting reads with the
  frozen `sample+chain+CDR3 nt+CID+corrected UMI` key. All 65 uniquely matched
  the frozen sample mask; 60 mapped exactly to marker-analysis Bin rows. Five
  HC/NL-66 mask-valid bins absent from the marker table were explicitly left
  aggregate-unassigned.
- All three prespecified Level A enriched stable_4bin regions had 0 CDR3
  molecules and 0 clonotypes. All 13 L3 anchor-confirmed Level B regions and
  all 71 Level C anchor-detected downgraded regions also had 0 CDR3 molecules.
- L3 globally contained seven repeated primary clonotypes, four cross-Bin;
  none overlapped any fixed marker-supported immune aggregate. The prior
  SSC/15491/14 three-CID/three-Bin IGH candidate lies outside the fixed regions.
- Final output:
  `results/reanalysis/L3_CDR3_marker_B_aggregate_overlap_260829_022556_attempt6/`.
  QA passed, exit code 0, empty stderr, and three Level A figures were created.
- K8 reads can be CID/spatially restricted during streaming, but no K8 cache
  exists and FASTQ is not spatially random-accessible, so a future run would
  still require one whole-chip sequential scan. The L3 0/13 overlap is not by
  itself a scientific rationale to start it; no K8 access occurred.
- Session summary:
  `.agent/l3_cdr3_marker_b_aggregate_overlap_260829024100.md`.
- No GitHub issue, Git commit/push or eLabFTW operation was performed, following
  the user's explicit scope restrictions.

## 2026-08-29 — L3 repeated CDR3 spatial-position audit

- Kept all 65 quality-corrected L3 molecules, the frozen primary clonotype key,
  551 immune boundaries, 73 anchor-confirmed B boundaries, 13 L3 anchor regions,
  and 71 L3 Level C regions unchanged.
- Used the ALL21 manifest-authoritative raw integer matrix and true member-Bin
  distances. Sixty molecules mapped exactly; five HC/NL-66 records remained
  protected exclusions. No molecule was inside a fixed immune aggregate.
- The seven repeated clonotypes comprised three same-Bin spatially clustered
  marker-weak records and four spatially dispersed molecular-only records.
  Boundary-associated=0, B-anchor-supported=0, and repeats across distinct
  Bin50 positions within three steps=0.
- SSC/15491/14 L3P_e2181118f0348020 had 3 distinct CIDs/3 Bins but one UMI
  sequence, no same/adjacent anchor, nearest immune boundary 127.5 um (5 steps),
  and was classified spatially dispersed molecular-only—not confirmed expansion.
- Final output:
  `results/reanalysis/L3_repeated_CDR3_spatial_audit_260829_025518_attempt7/`;
  QA passed with 13 non-empty figures and empty stderr.
- K8 scanning was not started and is not recommended from this audit alone.
- Session summary:
  `.agent/l3_repeated_cdr3_spatial_audit_260829032200.md`.
- No GitHub, git commit/push, eLabFTW, PPT, K8/J2 scan, Cell2location, NMF, or
  disease comparison was performed.

## 2026-09-07 - IGH/CDR3 continuous cellular-context overlay recovery

- Recovered the interrupted workflow in the existing fixed directory
  `results/reanalysis/IGH_CDR3_CURRENT_PROJECT/05_IGH_CONTINUOUS_CONTEXT_OVERLAY/`;
  no parallel result directory was created and frozen inputs were unchanged.
- The actual task script is `scripts/605_igh_continuous_context_overlay.py`;
  number 604 is occupied by the upstream RCTD continuous-context QC script.
- Final counts passed: 79/68/11/9 clonotypes and 69 deduplicated positions,
  split into 54 direct-RCTD and 15 context-only positions.
- Completed matched controls, 10,000 within-sample permutations,
  sensitivities, 11 evidence cards, figures, report, final validation, and a
  74-row hash-validated manifest.
- No all79 result reached BH-FDR <0.05. Four singleton-stratum FDR hits had
  cluster-bootstrap CIs crossing zero and do not support a robust overall
  claim. Direction-only robustness retained lower T/NK, B+Plasma/immune, and
  T/NK/immune.
- Multiple resume invocations and an invalid K8-flag result were isolated under
  `RUN_HISTORY/20260907_172500/INVALID_K8_FLAG_RUN`; only the corrected run is
  formal output.
- Session summary:
  `.agent/igh_continuous_context_overlay_recovery_260907173407.md`.
- No Git, eLabFTW, PPT, program analysis, or remote operation was performed.

## 2026-09-12 — CDR3 per-Bin Top20 heatmaps

- Continued from the existing frozen per-Bin Top20 inputs without reading a
  raw expression matrix or recomputing upstream IGH/CDR3 data.
- Generated three RECURRENT main heatmaps followed by three FULL supplementary
  heatmaps, each as PNG and PDF, in
  `results/reanalysis/CDR3_PER_BIN_TOP20_HEATMAP_FINAL_20260912_141500/`.
  Every completed file was logged immediately with its size and timestamp; the
  complete verified run took about 25 seconds and no FULL heatmap timed out.
- All69 contained 69 bins/604 full genes/93 recurrent genes; Shared7 contained
  44/327/55; Repeated11 contained 33/226/52. Frozen counts remained 104
  corrected IGH molecules and 79 patient-specific clonotypes.
- Final QC passed 29/29 checks. Independent verification found 0/26 MANIFEST
  path, size, or SHA-256 mismatches; all six PNGs decoded successfully and the
  visual audit found no truncated or blank rendering.
- Standalone HTML:
  `results/CDR3_PER_BIN_TOP20_HEATMAP_REPORT_260912142534.html` (6,327,999
  bytes; SHA-256
  `BBBC42693147166B5B7EC22A060825AF6458EAFD3BAEADFEEB6753E0CB4302EB`).
- Analysis/report commit:
  `7aea87d49493dbd5d95c7ebb7a89c70d466e00e4`.
- The user classified the workflow as a one-off; no skill-candidate issue was
  created. The R script is retained for reproduction.
- The user separately approved the ELN content and tags: project tag
  `Stereo-seq`, plus `IGH`, `CDR3`, `Bin50`, `Top20`, `heatmap`, and `QC`.
  Approved local entry:
  `elab/CDR3_per_Bin_Top20_heatmaps_260912142534.md`. No remote entry ID/URL
  exists because `ELABFTW_API_KEY` was unavailable.
- Session summary:
  `.agent/cdr3_per_bin_top20_heatmaps_wrapup_260912143417.md`.

## 2026-09-30 — Reference18 LYMPHOID annotation-consistency audit

- Audited all 238 frozen LYMPHOID regions in all 12 detected patients without changing boundaries or rerunning upstream models.
- Found a display-only but biologically important column-order error in the previous Chinese atlas: it used the state-definition row order rather than HDF5 `/states`. Its T/NK panel was actually B + Basal-like epithelial, and its Fibroblast panel was actually NK + T. Frozen region membership is unaffected because the generating script used `/states` correctly.
- With corrected mapping, 205/238 regions were B/Plasma-dominant (164/194 small; 41/44 large), 30 mixed, and 3 small regions T/NK-dominant. Three regions had prespecified weight-versus-expression discordance.
- No coordinate-registered B-cell stain/observation exists in the project records; historical B/plasma ROI labels are expression-module-derived, so sampling observation versus RCTD conflict is not directly assessable.
- Final directory: `results/reanalysis/REFERENCE18_LYMPHOID_ANNOTATION_CONSISTENCY_AUDIT_20260930_003136/`.
- Corrected English atlas: `pdf/REFERENCE18_LYMPHOID_DETECTED12_SPATIAL_ATLAS_EN_CORRECTED.pdf`; 12 page PNGs passed visual text/layout QA.
- HTML report: `REFERENCE18_LYMPHOID_ANNOTATION_CONSISTENCY_AUDIT_20260930_005517.html`.
- Session summary: `.agent/reference18_lymphoid_annotation_consistency_20260930_005900.md`.
- No eLabFTW submission, commit, push, or GitHub issue was performed; the worktree was already extensively dirty and `.agent/memory.md` had pre-existing edits.

## 2026-09-30 — Reference18 HDF5 state-order root-cause audit

- Traced the previous 12-patient atlas error to `scripts/894_reference18_lymphoid_spatial_atlas.py`: it interpreted `weights_normalized` columns using the state-definition biological order instead of HDF5 `/states`.
- The error affected only the atlas Reference18 weight panels, large-region layer-weight table, and weight-derived narratives. Frozen regions, mainline Level-1 composition, raw-count module/gene expression, disease detection/burden comparisons, and disease conclusions do not consume those erroneous atlas values and are unaffected.
- Revalidated frozen counts: LYMPHOID 238 regions/12 patients (194 small, 44 large); BP 43 regions/6 patients. All 331 BP member Bin50 satisfy B + Plasma >= 0.40, minimum 0.4003641195; all frozen coordinates match tissue H5AD.
- Corrected independent output: `results/reanalysis/REFERENCE18_H5_STATE_ORDER_ROOT_CAUSE_AUDIT_20260930_124304/`; the prior erroneous atlas remains preserved.
- Corrected English PDF contains 12 patient pages and passed 12/12 visual QA. Timestamped HTML root-cause report was rendered successfully.
- Session summary: `.agent/reference18_h5_state_order_root_cause_audit_20260930_130626.md`.
- No eLabFTW submission (per user instruction), commit, push, or GitHub issue was performed.

## 2026-10-08 – L3 K8 unbiased-niche spatial architecture

- Reused the frozen L3 K=8 assignment and analyzed all eight niches and all 28
  heterotypic pairs symmetrically; no clustering, RCTD, aggregate analysis, or
  disease significance test was rerun.
- The patient/ROI-specific Queen graph contained 39,982 evaluable Bin50s from
  seven patients/seven ROIs and 149,193 unique undirected edges. The null used
  999 label permutations within patient/ROI.
- N7, N4 and N5 were the most spatially continuous. N3 was comparatively
  fragmented (median largest-component fraction 0.065).
- N5–N6 showed recurrent positive adjacency enrichment in 6/7 patients;
  N4–N7, N2–N7 and N7–N8 were negative in 7/7 patients. These are spatial
  associations only: K8 labels derive from overlapping 50 µm neighborhoods,
  and the permutation null disrupts endogenous autocorrelation.
- Final directory:
  `results/reanalysis/L3_UNBIASED_NICHE_RESOLUTION_20261007_085730/K08_SPATIAL_ARCHITECTURE_20261007_121357/`.
- Rendered HTML:
  `results/L3_K8_SPATIAL_ARCHITECTURE_REPORT_261008104540.html`.
- Analysis/report commit: `f0d0d1a`.
- The user approved the ELN draft on 2026-10-08. Approved local copy:
  `elab/L3_K8_spatial_niche_architecture_261008104540.md`. No remote entry
  ID/URL exists because `ELABFTW_API_KEY` is unavailable.
- Session summary:
  `.agent/l3_k8_spatial_architecture_wrapup_261008105343.md`.

## 2026-10-08 – J2 unbiased neighborhood niche resolution

- Independently analyzed Y40105J2 using 58,840 UMI>=50 Bin50, frozen
  Reference18 continuous weights, official 50 µm radius neighborhoods, native
  VoltRon CLR, and K-means K=3–12 with 20 seeds per K.
- 58,772 Bin50 had non-empty neighborhoods; 68 were retained as not
  evaluable. The graph had 325,776 edges and no cross-patient/ROI edges.
- The audited primary resolution is K=5, with K=3 as the lower-resolution
  comparator and K=7 as the higher-resolution sensitivity. K=5 converged in
  20/20 fits, median ARI 0.988, median NMI 0.979, and every niche covered all
  seven patients.
- Conservative K=5 themes are alveolar-associated, B-enriched,
  airway-associated, lymphatic-endothelial-associated, and
  blood-endothelial-associated mixed niches. Core axes were not
  single-patient-driven and retained direction in UMI sensitivity checks.
- Frozen aggregates and disease labels were not used in discovery; no disease
  significance, DEG/GSEA, pathway, communication, or cross-chip consensus was
  performed.
- Final directory:
  `results/reanalysis/J2_UNBIASED_NICHE_RESOLUTION_20261008_110522/`.
- Rendered HTML:
  `results/J2_UNBIASED_NICHE_RESOLUTION_REPORT_261008113930.html`.
- Analysis/report commit: `cb6b321`.
- Approved local ELN:
  `elab/J2_unbiased_neighborhood_niche_resolution_261008113930.md`; no remote
  entry ID/URL because `ELABFTW_API_KEY` is unavailable.
- The user approved promoting the per-chip L3/J2/Y40102K8 framework to a
  recurring workflow. The candidate was filed as GitHub issue #1:
  `https://github.com/bozhu001/stereo_seq/issues/1`; the submitted draft is
  preserved at
  `.agent/voltron_unbiased_neighborhood_skill_candidate_261008115210.md`.
- Session summary:
  `.agent/j2_unbiased_niche_resolution_wrapup_261008115210.md`.

## 2026-10-08 – J2 K=5 spatial niche architecture and freeze

- Reused the frozen Y40105J2 K=5 assignment for 58,772 evaluable Bin50 from
  seven patients and seven ROIs; no clustering or upstream analysis was rerun.
- The within-patient/ROI Queen graph contained 218,614 unique edges. All five
  niches and all ten heterotypic pairs were analyzed symmetrically with 999
  label permutations (seed `2026100811`).
- N2, N4, and N1 had the highest median largest-component fractions, although
  the overall niche fields remained fragmented. Nine of ten heterotypic pairs
  were direction-consistent in at least 6/7 patients, all negative relative to
  the exchangeable-label null; no recurrent positive pair was identified.
- Final directory:
  `results/reanalysis/J2_UNBIASED_NICHE_RESOLUTION_20261008_110522/K05_SPATIAL_ARCHITECTURE_20261008_120148/`.
- Rendered HTML:
  `results/J2_K5_SPATIAL_ARCHITECTURE_REPORT_261008121443.html`.
- Analysis/report commit:
  `793085cee696e2432ce055212a286229d4429232`.
- Approved local ELN:
  `elab/J2_K5_spatial_niche_architecture_261008121443.md`; no remote ID/URL
  because `ELABFTW_API_KEY` is unavailable.
- J2 discovery, composition, consistency, annotation, and spatial outputs are
  frozen by `J2_NICHE_RESULTS_FREEZE_20261008.tsv` and `J2_FROZEN.ok`; the user
  authorized no further J2 analysis.
- Session summary:
  `.agent/j2_k5_spatial_architecture_wrapup_261008122000.md`.

## 2026-10-10 — Extended Immune Marker Library and RNA detectability audit

- Audited `BI_PF_ILD_atlas_v1` and recovered the workbook's full data despite
  a stale `A1` worksheet dimension; the source Excel was not modified.
- The immune hierarchy contained 9 Level 1 and 29 Level 2 categories, with
  1,450 marker records and 692 unique genes. The conservative library retained
  27 strict lineage-core and 209 eligible extended non-core genes.
- Used sparse raw `uint32` counts for all qualified Bin50 from 21 patients on
  L3, J2 and Y40102K8. In total, 667 assessed genes were present and detected
  at least once. Extended panels added 193 detectable lineage-panel entries
  relative to the prior small panels.
- The result supports improved detection coverage only, not improved immune-cell
  identification accuracy. No RCTD, clustering, K selection, neighborhood
  scoring, disease testing, DEG/GSEA, CellChat or frozen aggregate analysis ran.
- Final directory:
  `results/reanalysis/EXTENDED_IMMUNE_MARKER_LIBRARY_20261010/`.
- Standalone HTML:
  `reports/EXTENDED_IMMUNE_MARKER_LIBRARY_REPORT_261010185546.html` within the
  final directory; SHA256
  `367FBA538F808A092475A9A29B5A884AFECFFA542013B876214AAB25474BCA5A`.
- Analysis/report commit:
  `64a30fee048c59430354f7c7c814832ca0e1389c`.
- The user separately approved the ELN draft. Approved local copy:
  `elab/Extended_Immune_Marker_Library_and_Stereo-seq_RNA_Detectability_Audit_261010185546.md`.
  No remote ID/URL exists because `ELABFTW_API_KEY` was unavailable.
- Script-review classification remains pending; no skill-candidate issue was
  filed.
- Session summary:
  `.agent/extended_immune_marker_library_wrapup_261010190419.md`.
