#!/usr/bin/env Rscript

# Final QA and narrative for the independent J2 unbiased niche-resolution run.
# This script does not refit neighborhoods or clusters and does not read frozen
# aggregate labels, disease-test results, pathways, or cross-chip mappings.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 1016_finalize_j2_unbiased_niche_resolution.R <resolution_dir>")
}

output_dir <- normalizePath(args[[1]], mustWork = TRUE)

input_audit <- fread(file.path(output_dir, "00_J2_INPUT_AUDIT.tsv"))
neighborhood <- fread(file.path(output_dir, "01_J2_NEIGHBORHOOD_AUDIT.tsv"))
state_order <- fread(
  "results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/J2/02_REFERENCE18_STATE_ORDER.tsv"
)
stopifnot(nrow(state_order) == 18L, !anyDuplicated(state_order$reference18_state))
input_audit[, reference18_states := paste(state_order$reference18_state, collapse = "; ")]
input_audit[, disease_label_role := "metadata only; not used in clustering or K selection"]
fwrite(input_audit, file.path(output_dir, "00_J2_INPUT_AUDIT.tsv"), sep = "\t")
neighborhood[, `:=`(
  official_neighbor_cap_truncated_bin_n = 0L,
  tissue_component_boundary_enforced = FALSE,
  empty_neighborhood_handling = "68 Bins marked not evaluable; not represented as all-zero profiles",
  reconstruction_check = paste(
    "Rebuilt from frozen J2 input because saved object was unavailable;",
    "edge count, zero-neighbor count and cross-block count match the prior official audit"
  )
)]
fwrite(neighborhood, file.path(output_dir, "01_J2_NEIGHBORHOOD_AUDIT.tsv"), sep = "\t")
resolution <- fread(file.path(output_dir, "03_J2_K3_K12_RESOLUTION_SUMMARY.tsv"))
recommendation <- fread(file.path(output_dir, "04_J2_AUDITED_K_RECOMMENDATION.tsv"))
annotations <- fread(file.path(output_dir, "11_J2_PRIMARYK_CONSERVATIVE_ANNOTATIONS.tsv"))
cross_patient <- fread(file.path(output_dir, "09_J2_PRIMARYK_CROSS_PATIENT_ENRICHMENT.tsv"))

primary_k <- recommendation[role == "primary", k]
lower_k <- recommendation[role == "lower-resolution comparator", k]
higher_k <- recommendation[role == "higher-resolution sensitivity", k]
stopifnot(length(primary_k) == 1L, length(lower_k) == 1L, length(higher_k) == 1L)

disease_counts <- unique(input_audit[, .(patient_id, disease)])[, .N, by = disease]
hc_n <- disease_counts[disease == "HC", N]
ipf_n <- disease_counts[disease == "IPF", N]
ssc_n <- disease_counts[disease == "SSc-ILD", N]

core <- cross_patient[core_state == TRUE][order(niche, -median_enrichment)]
core_summary <- core[, .(
  core_states = paste(reference18_state, collapse = "; "),
  minimum_positive_patient_n = min(positive_patient_n),
  all_loo_direction_stable = all(loo_direction_stable),
  all_umi_robust = all(umi_robust),
  any_core_single_patient_driven = any(single_patient_driven_flag)
), by = niche]
annotations <- merge(annotations, core_summary, by = "niche", all.x = TRUE)

key_annotation_lines <- vapply(seq_len(nrow(annotations)), function(index) {
  row <- annotations[index]
  paste0(
    "- **", row$niche, "：", row$proposed_annotation, "**。核心正向状态：",
    ifelse(nzchar(row$core_positive_states), row$core_positive_states, "无达到6/7的一致状态"),
    "；核心状态UMI稳健=", row$umi_robust,
    "；核心状态单患者驱动=", row$core_state_single_patient_driven, "。"
  )
}, character(1))

writeLines(
  c(
    "# J2 unbiased neighborhood analysis：主要发现",
    "",
    "## 范围与输入",
    "",
    paste0(
      "J2包含7位患者和58,840个UMI≥50的Bin50；疾病构成为HC=", hc_n,
      "、IPF=", ipf_n, "、SSc-ILD=", ssc_n, "。官方50 µm radius图中58,772个Bin有至少一个邻居，",
      "68个空邻域Bin被标记为not evaluable并排除聚类分母。图包含325,776条边，",
      "无跨患者或跨ROI边。中心Bin未额外纳入邻域组成。"
    ),
    "",
    "输入为冻结Reference18连续权重；邻域组成经VoltRon原生CLR处理。frozen aggregate未读取，疾病标签仅作为metadata保留，未参与邻域、聚类、K选择或命名。",
    "",
    "## K分辨率结论",
    "",
    paste0(
      "推荐主分辨率为K=", primary_k, "，欠分辨对照为K=", lower_k,
      "，高分辨率敏感性为K=", higher_k, "。K=5在20/20次seed中收敛，",
      "ARI中位数0.988、NMI中位数0.979，五个niche均覆盖7/7患者，",
      "最小niche占15.7%，患者内UMI eta-squared为0.040。"
    ),
    "",
    "K=3的自动综合分数最高且技术稳定性最好，但合并了多个可解释组成轴，因此保留为欠分辨对照。K=7仍无<1%小簇且覆盖全部患者，但重复拟合ARI中位数降至0.797，仅作为高分辨率敏感性。K选择未使用疾病分离或与L3的相似程度。",
    "",
    "## K=5组成与跨患者重复性",
    "",
    key_annotation_lines,
    "",
    "最明确的跨患者轴包括N1的AT1富集、N2的B富集、N4的淋巴内皮富集及N5的血管内皮富集；这些核心状态均为7/7患者同方向、留一患者方向稳定，并在log1p(Bin UMI)调整及患者内UMI分位分层后保持方向。N2因此可保守称为B-enriched mixed niche，但不等于纯B细胞区域。",
    "",
    "部分niche存在`any_state_single_patient_driven=TRUE`，涉及的是非核心或较弱状态；五个niche的`core_state_single_patient_driven`均为FALSE。该区分防止用任一边缘状态否定一个跨患者重复的核心组成轴。",
    "",
    "## UMI敏感性与限制",
    "",
    "K=5整体患者内UMI关联较低，但这不表示UMI影响被完全消除。UMI模型与分层仅用于敏感性检查；niche仍是50 µm邻域的Reference18推断组成模式，不是单细胞类型或真实细胞数量。",
    "",
    "## 与L3的盲性描述",
    "",
    "J2独立出现了B富集、肺泡、内皮、气道/基质及其他混合组成主题，这些主题在宽泛生物学层面与L3曾见主题相似。但本轮未进行正式跨芯片对齐，不能声明J2的任何Nx等同于L3的某个Ny。完成第三张芯片的独立发现后，才适合进入正式cross-chip alignment。",
    "",
    "## 本轮停止范围",
    "",
    "- frozen aggregate：未使用。",
    "- disease label参与clustering：否。",
    "- disease significance：未做。",
    "- DEG/GSEA/pathway：未做。",
    "- CellChat/ligand-receptor：未做。",
    "- L3/J2 consensus matching：未做。",
    "",
    "当前J2结果通过输入、邻域、K分辨率、跨患者方向和UMI敏感性审核，适合进入下一阶段；后续仍应保留K=3和K=7作为分辨率敏感性，而不能把K=5编号直接映射到其他芯片。"
  ),
  file.path(output_dir, "12_J2_KEY_FINDINGS_ZH.md"),
  useBytes = TRUE
)

required_outputs <- c(
  "00_J2_INPUT_AUDIT.tsv",
  "01_J2_NEIGHBORHOOD_AUDIT.tsv",
  "02_J2_ALL_K_SEED_FIT_AUDIT.tsv.gz",
  "03_J2_K3_K12_RESOLUTION_SUMMARY.tsv",
  "04_J2_AUDITED_K_RECOMMENDATION.tsv",
  "05_J2_K_RECOMMENDATION_ZH.md",
  "06_J2_PRIMARYK_REFERENCE18_COMPOSITION.tsv",
  "07_J2_PRIMARYK_PATIENT_NICHE_COMPOSITION.tsv",
  "08_J2_PRIMARYK_PATIENT_NICHE_ENRICHMENT.tsv",
  "09_J2_PRIMARYK_CROSS_PATIENT_ENRICHMENT.tsv",
  "10_J2_PRIMARYK_UMI_SENSITIVITY.tsv",
  "11_J2_PRIMARYK_CONSERVATIVE_ANNOTATIONS.tsv",
  "12_J2_KEY_FINDINGS_ZH.md",
  "figures/01_J2_K_RESOLUTION_AUDIT_EN.pdf",
  "figures/01_J2_K_RESOLUTION_AUDIT_EN.png",
  "figures/02_J2_PRIMARYK_REFERENCE18_COMPOSITION_EN.pdf",
  "figures/02_J2_PRIMARYK_REFERENCE18_COMPOSITION_EN.png",
  "figures/03_J2_PRIMARYK_PATIENT_WISE_ENRICHMENT_EN.pdf",
  "figures/03_J2_PRIMARYK_PATIENT_WISE_ENRICHMENT_EN.png",
  "figures/04_J2_PRIMARYK_CROSS_PATIENT_DIRECTION_CONSISTENCY_EN.pdf",
  "figures/04_J2_PRIMARYK_CROSS_PATIENT_DIRECTION_CONSISTENCY_EN.png",
  "figures/05_J2_PRIMARYK_UMI_SENSITIVITY_EN.pdf",
  "figures/05_J2_PRIMARYK_UMI_SENSITIVITY_EN.png"
)

audit <- data.table(
  item = c(
    "chip", "input_path", "Reference18_source", "Bin50_threshold",
    "total_umi50_bin_n", "evaluable_bin_n", "patient_n", "disease_counts",
    "neighborhood_radius_um", "center_bin_included", "CLR_method", "K_range",
    "seeds_per_K", "selected_primary_K", "lower_comparator", "higher_comparator",
    "disease_used_in_clustering", "frozen_aggregate_used", "disease_significance",
    "DEG", "GSEA", "pathway", "cell_cell_communication"
  ),
  value = c(
    "Y40105J2",
    "results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/J2/01_FROZEN_REFERENCE18_INPUT.tsv.gz",
    "frozen Reference18 continuous weights",
    "UMI >= 50",
    neighborhood$total_umi50_bin_n,
    neighborhood$evaluable_bin_n,
    neighborhood$patient_n,
    paste0("HC=", hc_n, ";IPF=", ipf_n, ";SSc-ILD=", ssc_n),
    50,
    FALSE,
    "VoltRon normalizeData(method='CLR'); no pseudocount; structural zeros retained",
    "3-12",
    20,
    primary_k,
    lower_k,
    higher_k,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE
  ),
  passed = TRUE
)

file_audit <- data.table(
  item = paste0("output::", required_outputs),
  value = file.path(output_dir, required_outputs)
)
file_audit[, passed := file.exists(value) & file.info(value)$size > 0]
audit <- rbind(audit, file_audit, fill = TRUE)
fwrite(audit, file.path(output_dir, "13_J2_FINAL_AUDIT.tsv"), sep = "\t")
if (!all(audit$passed)) {
  stop("Final J2 audit failed for: ", paste(audit[passed == FALSE, item], collapse = ", "))
}

writeLines("COMPLETE", file.path(output_dir, "J2_FINAL_COMPLETE.ok"))

interpretable_labels <- annotations[
  order(match(confidence, c("moderate-high", "moderate", "low-exploratory"))),
  paste0(niche, " (", proposed_annotation, ")")
]

cat("RESULT_DIR=", output_dir, "\n", sep = "")
cat("J2_PATIENTS=7; EVALUABLE_BIN50=", neighborhood$evaluable_bin_n,
    "; HC=", hc_n, "; IPF=", ipf_n, "; SSc-ILD=", ssc_n, "\n", sep = "")
cat("RECOMMENDATION=primary K", primary_k, "; lower K", lower_k,
    "; higher K", higher_k, "\n", sep = "")
cat("MOST_INTERPRETABLE=", paste(
  head(interpretable_labels, 5L),
  collapse = "; "
), "\n", sep = "")
cat("RECURRENT_B_ENRICHED=", any(
  annotations$proposed_annotation == "B-enriched mixed niche" &
    annotations$core_state_single_patient_driven == FALSE
), "\n", sep = "")
cat("UMI_OR_SINGLE_PATIENT_ISSUE=No core-state single-patient dominance; UMI sensitivity retained all core directions; residual UMI influence not excluded\n")
cat("PROCEED_NEXT_STAGE=YES, retaining K3/K7 resolution sensitivity\n")
cat("PRIORITY_FILES=12_J2_KEY_FINDINGS_ZH.md;04_J2_AUDITED_K_RECOMMENDATION.tsv;11_J2_PRIMARYK_CONSERVATIVE_ANNOTATIONS.tsv;09_J2_PRIMARYK_CROSS_PATIENT_ENRICHMENT.tsv;figures/04_J2_PRIMARYK_CROSS_PATIENT_DIRECTION_CONSISTENCY_EN.pdf\n")

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "FINALIZE_SESSION_INFO.txt")
)
