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
