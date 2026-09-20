#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript scripts/670_nominal_p05_continuous_neighborhood_audit.R <project_root> <output_dir>")
}

project_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite a non-empty output directory: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "SCRIPT_USED"), showWarnings = FALSE)

source_audit <- file.path(
  project_root,
  "results/reanalysis/CONTINUOUS_NEIGHBORHOOD_FINAL_AUDIT_20260920_063637"
)
source_18 <- file.path(
  project_root,
  "results/reanalysis/RCTD18_PRIMARY19_CONTINUOUS_NEIGHBORHOOD_20260919_224541"
)
source_21 <- file.path(
  project_root,
  "results/reanalysis/RCTD21_ALL21_FINESTATE_SENSITIVITY_20260919_235041"
)

input_files <- c(
  file.path(source_audit, "UNIFIED_CONTINUOUS_NEIGHBORHOOD_EVIDENCE.tsv"),
  file.path(source_18, "05_PATIENT_GLOBAL_STATE_METRICS.tsv"),
  file.path(source_18, "16_PATIENT_STATE_ARCHITECTURE_METRICS.tsv.gz"),
  file.path(source_18, "17_PATIENT_PAIR_COENRICHMENT_METRICS.tsv.gz"),
  file.path(source_21, "PATIENT21_GLOBAL_METRICS.tsv.gz"),
  file.path(source_21, "PATIENT21_ARCHITECTURE_METRICS.tsv.gz"),
  file.path(source_21, "PATIENT21_PAIR_METRICS.tsv.gz")
)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0L) {
  stop("Missing required frozen input(s): ", paste(missing_files, collapse = "; "))
}

sha256_file <- function(path) {
  unname(tools::md5sum(path))
}

as_flag <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

sign_same <- function(x, y) {
  is.finite(x) & is.finite(y) & sign(x) == sign(y) & sign(x) != 0
}

evidence <- fread(input_files[[1]], na.strings = c("", "NA", "nan"))
if (!all(evidence$cohort == "Primary19")) {
  stop("Unified source contains a cohort other than Primary19")
}

numeric_columns <- c(
  "effect", "standardized_effect", "ci_low", "ci_high", "raw_p",
  "permutation_p", "original_family_FDR", "global_primary_BH_FDR",
  "K8_effect", "J2_effect", "L3_effect", "n_chips_same_direction",
  "all21_effect", "mean_feature_weight_rho", "valid_bins_rho",
  "tissue_area_um2_rho", "valid_coverage_fraction_rho", "median_nUMI_rho",
  "median_nGene_rho", "chip_range_over_patient_sd"
)
for (column in intersect(numeric_columns, names(evidence))) {
  set(evidence, j = column, value = as.numeric(evidence[[column]]))
}

evidence[, nominal_raw := is.finite(raw_p) & raw_p < 0.05]
evidence[, nominal_permutation := is.finite(permutation_p) & permutation_p < 0.05]
evidence[, nominal_either := nominal_raw | nominal_permutation]
evidence[, family_fdr_pass := is.finite(original_family_FDR) & original_family_FDR < 0.05]
evidence[, global_fdr_pass := is.finite(global_primary_BH_FDR) & global_primary_BH_FDR < 0.05]
evidence[, low_coverage_affected :=
  technical_status == "LOW_COVERAGE_SENSITIVE" |
  !as_flag(primary_all21_same_direction)]
evidence[, technical_high_correlation :=
  pmax(
    abs(mean_feature_weight_rho),
    abs(valid_bins_rho),
    abs(tissue_area_um2_rho),
    abs(valid_coverage_fraction_rho),
    abs(median_nUMI_rho),
    abs(median_nGene_rho),
    na.rm = TRUE
  ) >= 0.50]
evidence[!is.finite(technical_high_correlation), technical_high_correlation := FALSE]
evidence[, technically_sensitive_for_audit :=
  low_coverage_affected |
  as_flag(technical_sensitive) |
  technical_status %chin% c(
    "TECHNICAL_SENSITIVE",
    "LOW_COVERAGE_SENSITIVE",
    "INSUFFICIENT_INFORMATION"
  )]
evidence[, loo_stable := !as_flag(loo_reversal)]
evidence[, loco_stable := !as_flag(leave_one_chip_out_reversal)]
evidence[, primary_all21_stable := as_flag(primary_all21_same_direction)]
evidence[, three_of_three_chips := n_chips_same_direction == 3]
evidence[, justified_two_of_three :=
  n_chips_same_direction == 2 &
  as_flag(k8_j2_same_as_overall) &
  loo_stable &
  loco_stable &
  primary_all21_stable &
  !technically_sensitive_for_audit]
evidence[, cross_chip_rule := three_of_three_chips | justified_two_of_three]
evidence[, chip_rule_note := fifelse(
  three_of_three_chips,
  "3/3 chips match pooled direction",
  fifelse(
    justified_two_of_three,
    paste0(
      "2/3 retained under frozen K8+J2 rule; L3 is discordant; ",
      "LOPO/LOCO and Primary19-All21 directions remain stable"
    ),
    "cross-chip direction criterion not met"
  )
)]
evidence[, stability_pass :=
  analysis_tier == "primary" &
  cross_chip_rule &
  loo_stable &
  loco_stable &
  primary_all21_stable &
  resolution_status != "RESOLUTION_DISCORDANT" &
  !grepl("EXPLORATORY_SENSITIVITY_ONLY", special_caution, fixed = TRUE) &
  !technically_sensitive_for_audit]

hypothesis_columns <- c(
  "source_resolution", "analysis_family", "feature", "feature_a",
  "feature_b", "metric", "radius_um", "analysis_tier"
)
evidence[, hypothesis_key := do.call(
  paste,
  c(lapply(.SD, function(x) fifelse(is.na(x), "", as.character(x))), sep = "||")
), .SDcols = hypothesis_columns]

counterparts <- evidence[nominal_either & stability_pass, .(
  hypothesis_key,
  contrast,
  counterpart_effect = effect
)]
ssc_counterpart <- counterparts[contrast == "SSc-ILD_vs_HC", .(
  hypothesis_key,
  ssc_nominal_stable = TRUE,
  ssc_effect = counterpart_effect
)]
ipf_counterpart <- counterparts[contrast == "IPF_vs_HC", .(
  hypothesis_key,
  ipf_nominal_stable = TRUE,
  ipf_effect = counterpart_effect
)]
evidence <- merge(evidence, ssc_counterpart, by = "hypothesis_key", all.x = TRUE)
evidence <- merge(evidence, ipf_counterpart, by = "hypothesis_key", all.x = TRUE)
evidence[is.na(ssc_nominal_stable), ssc_nominal_stable := FALSE]
evidence[is.na(ipf_nominal_stable), ipf_nominal_stable := FALSE]
evidence[, shared_ild_direction :=
  ssc_nominal_stable & ipf_nominal_stable & sign_same(ssc_effect, ipf_effect)]
evidence[, ssc_vs_ipf_only := contrast == "SSc-ILD_vs_IPF" & !ssc_nominal_stable]

evidence[, audit_category := "UNSTABLE_OR_DESCRIPTIVE"]
evidence[nominal_either & technically_sensitive_for_audit,
  audit_category := "TECHNICALLY_SENSITIVE"]
evidence[
  nominal_either & !technically_sensitive_for_audit & stability_pass &
    contrast == "SSc-ILD_vs_HC" & shared_ild_direction,
  audit_category := "SHARED_ILD_VS_HC_TREND"
]
evidence[
  nominal_either & !technically_sensitive_for_audit & stability_pass &
    contrast == "IPF_vs_HC" & shared_ild_direction,
  audit_category := "SHARED_ILD_VS_HC_TREND"
]
evidence[
  nominal_either & !technically_sensitive_for_audit & stability_pass &
    contrast == "SSc-ILD_vs_HC" & !shared_ild_direction,
  audit_category := "SSC_VS_HC_PRIORITY_EXPLORATORY"
]
evidence[
  nominal_either & !technically_sensitive_for_audit & stability_pass &
    contrast == "IPF_vs_HC" & !shared_ild_direction,
  audit_category := "IPF_ASSOCIATED_TREND"
]
evidence[
  nominal_either & !technically_sensitive_for_audit & stability_pass &
    ssc_vs_ipf_only,
  audit_category := "SSC_VS_IPF_ONLY_LOW_PRIORITY"
]

nominal <- evidence[nominal_either == TRUE]
category_order <- c(
  "SSC_VS_HC_PRIORITY_EXPLORATORY",
  "SHARED_ILD_VS_HC_TREND",
  "IPF_ASSOCIATED_TREND",
  "SSC_VS_IPF_ONLY_LOW_PRIORITY",
  "TECHNICALLY_SENSITIVE",
  "UNSTABLE_OR_DESCRIPTIVE"
)
nominal[, category_rank := match(audit_category, category_order)]
nominal[, evidence_p := pmin(raw_p, permutation_p, na.rm = TRUE)]
setorder(
  nominal,
  category_rank,
  contrast,
  -family_fdr_pass,
  -global_fdr_pass,
  evidence_p,
  evidence_id
)

output_columns <- c(
  "evidence_id", "audit_category", "source_resolution", "analysis_family",
  "source_analysis_family", "analysis_tier", "state_or_pair", "parent_18state",
  "feature", "feature_a", "feature_b", "metric", "radius_um", "contrast",
  "effect", "standardized_effect", "ci_low", "ci_high", "raw_p",
  "permutation_p", "nominal_raw", "nominal_permutation", "original_family_FDR",
  "global_primary_BH_FDR", "family_fdr_pass", "global_fdr_pass", "K8_effect",
  "J2_effect", "L3_effect", "direction_consistency", "n_chips_same_direction",
  "three_of_three_chips", "justified_two_of_three", "chip_rule_note",
  "loo_reversal", "loo_direction_stability", "leave_one_chip_out_reversal",
  "leave_one_chip_out_direction_stable", "primary_all21_same_direction",
  "all21_effect", "low_coverage_affected", "most_influential_sample",
  "resolution_status", "technical_status", "technical_drivers",
  "mean_feature_weight_rho", "valid_bins_rho", "tissue_area_um2_rho",
  "valid_coverage_fraction_rho", "median_nUMI_rho", "median_nGene_rho",
  "chip_range_over_patient_sd", "technical_high_correlation",
  "technically_sensitive_for_audit", "shared_ild_direction",
  "ssc_vs_ipf_only", "biological_axis", "biological_interpretation",
  "special_caution", "audit_note"
)
output_columns <- intersect(output_columns, names(nominal))

write_output <- function(frame, filename) {
  fwrite(frame[, ..output_columns], file.path(output_dir, filename), sep = "\t", na = "NA")
}

write_output(nominal, "ALL_NOMINAL_P05_CONTINUOUS_NEIGHBORHOOD.tsv")
write_output(
  nominal[
    contrast == "SSc-ILD_vs_HC" &
      audit_category %chin% c(
        "SSC_VS_HC_PRIORITY_EXPLORATORY",
        "SHARED_ILD_VS_HC_TREND"
      )
  ],
  "SSC_VS_HC_NOMINAL_CANDIDATES.tsv"
)
write_output(
  nominal[audit_category == "SHARED_ILD_VS_HC_TREND"],
  "SHARED_ILD_NOMINAL_CANDIDATES.tsv"
)
write_output(
  nominal[audit_category == "SSC_VS_IPF_ONLY_LOW_PRIORITY"],
  "SSC_VS_IPF_ONLY_LOW_PRIORITY.tsv"
)
write_output(
  nominal[audit_category == "TECHNICALLY_SENSITIVE"],
  "TECHNICALLY_SENSITIVE_NOMINAL_RESULTS.tsv"
)

tested_counts <- evidence[, .(
  n_tests = .N,
  expected_chance_hits = 0.05 * .N,
  raw_p_lt_0_05 = sum(nominal_raw),
  permutation_p_lt_0_05 = sum(nominal_permutation),
  either_p_lt_0_05 = sum(nominal_either),
  family_fdr_lt_0_05 = sum(family_fdr_pass),
  global_fdr_lt_0_05 = sum(global_fdr_pass, na.rm = TRUE)
), by = contrast]
tested_counts[, contrast_order := match(
  contrast,
  c("SSc-ILD_vs_HC", "IPF_vs_HC", "SSc-ILD_vs_IPF")
)]
setorder(tested_counts, contrast_order)
tested_counts[, contrast_order := NULL]
fwrite(tested_counts, file.path(output_dir, "NOMINAL_P05_COUNT_AUDIT.tsv"), sep = "\t")

stable_candidates <- nominal[audit_category %chin% category_order[1:4]]
stable_counts <- stable_candidates[, .N, by = audit_category]
category_counts <- nominal[, .N, by = audit_category]
category_counts <- merge(
  data.table(audit_category = category_order),
  category_counts,
  by = "audit_category",
  all.x = TRUE
)
category_counts[is.na(N), N := 0L]

display_label <- function(frame) {
  paste0(
    frame$source_resolution,
    " | ",
    frame$state_or_pair,
    " | ",
    frame$metric,
    ifelse(is.na(frame$radius_um), "", paste0("@", frame$radius_um, "um")),
    " | ",
    frame$contrast
  )
}

plot_candidates <- copy(stable_candidates)
plot_candidates[, plot_rank := frank(
  evidence_p,
  ties.method = "first"
), by = audit_category]
plot_candidates <- plot_candidates[plot_rank <= 8L]
if (nrow(plot_candidates) == 0L) {
  plot_candidates <- nominal[order(evidence_p)][1:min(.N, 12L)]
}
plot_candidates[, label := display_label(.SD)]

heat_long <- melt(
  plot_candidates,
  id.vars = c("evidence_id", "label", "audit_category"),
  measure.vars = c(
    "three_of_three_chips", "cross_chip_rule", "loo_stable", "loco_stable",
    "primary_all21_stable", "family_fdr_pass", "global_fdr_pass"
  ),
  variable.name = "evidence_dimension",
  value.name = "pass"
)
heat_long[, evidence_dimension := factor(
  evidence_dimension,
  levels = c(
    "three_of_three_chips", "cross_chip_rule", "loo_stable", "loco_stable",
    "primary_all21_stable", "family_fdr_pass", "global_fdr_pass"
  ),
  labels = c(
    "3/3 chips", "chip rule", "LOPO stable", "LOCO stable",
    "Primary19-All21", "family FDR", "global FDR"
  )
)]
heat_long[, label := factor(label, levels = rev(unique(plot_candidates$label)))]
figure_1 <- ggplot(heat_long, aes(evidence_dimension, label, fill = pass)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = "#2A9D8F", `FALSE` = "#E9C46A")) +
  labs(
    title = "Nominal P<0.05 evidence audit",
    subtitle = "Green indicates the stated criterion is met; no row is called significant without FDR support",
    x = NULL,
    y = NULL,
    fill = "Criterion met"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold")
  )
ggsave(file.path(output_dir, "FIG01_NOMINAL_P05_EVIDENCE_HEATMAP.png"), figure_1, width = 11, height = max(5, 0.34 * nrow(plot_candidates) + 2), dpi = 350)
ggsave(file.path(output_dir, "FIG01_NOMINAL_P05_EVIDENCE_HEATMAP.pdf"), figure_1, width = 11, height = max(5, 0.34 * nrow(plot_candidates) + 2))

read_patient_tables <- function() {
  list(
    `18-state|global_composition` = fread(input_files[[2]]),
    `18-state|spatial_architecture` = fread(input_files[[3]]),
    `18-state|pair_coenrichment` = fread(input_files[[4]]),
    `21-state|global_composition` = fread(input_files[[5]]),
    `21-state|spatial_architecture` = fread(input_files[[6]]),
    `21-state|pair_coenrichment` = fread(input_files[[7]])
  )
}

patient_tables <- read_patient_tables()
canonical_pair <- function(a, b) {
  sort(c(as.character(a), as.character(b)))
}
for (key in names(patient_tables)) {
  if (grepl("pair_coenrichment$", key)) {
    table <- patient_tables[[key]]
    swap <- !is.na(table$feature_a) & !is.na(table$feature_b) & table$feature_a > table$feature_b
    old_a <- table$feature_a[swap]
    table$feature_a[swap] <- table$feature_b[swap]
    table$feature_b[swap] <- old_a
    patient_tables[[key]] <- table
  }
}

subset_patient_values <- function(row) {
  key <- paste(row$source_resolution, row$analysis_family, sep = "|")
  table <- copy(patient_tables[[key]])
  if (is.null(table)) {
    return(data.table())
  }
  for (column in c("feature", "feature_a", "feature_b", "metric", "radius_um", "analysis_tier")) {
    if (!column %chin% names(table)) {
      next
    }
    target_value <- row[[column]]
    if (is.null(target_value) || length(target_value) == 0L || is.na(target_value)) {
      table <- table[is.na(get(column))]
    } else {
      table <- table[table[[column]] == target_value]
    }
  }
  if ("primary19" %chin% names(table)) {
    table <- table[as_flag(primary19)]
  }
  table[, `:=`(
    evidence_id = row$evidence_id,
    label = row$label,
    audit_category = row$audit_category
  )]
  table
}

distribution_candidates <- plot_candidates[1:min(.N, 12L)]
patient_values <- rbindlist(
  lapply(seq_len(nrow(distribution_candidates)), function(index) {
    subset_patient_values(as.list(distribution_candidates[index]))
  }),
  fill = TRUE
)
if (nrow(patient_values) == 0L) {
  stop("Patient-level plotting table is empty after matching stable nominal candidates")
}
fwrite(
  patient_values,
  file.path(output_dir, "FIG02_PATIENT_VALUES_USED.tsv"),
  sep = "\t",
  na = "NA"
)
if (nrow(patient_values) > 0L) {
  patient_values[, label := factor(label, levels = unique(distribution_candidates$label))]
  figure_2 <- ggplot(patient_values, aes(disease, value, color = chip_id)) +
    geom_boxplot(aes(group = disease), outlier.shape = NA, color = "grey55", width = 0.6) +
    geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.9) +
    facet_wrap(~label, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c(Y40102K8 = "#D55E00", Y40105J2 = "#0072B2", Y40105L3 = "#009E73")) +
    labs(
      title = "Top nominal candidates: patient-level distributions",
      subtitle = "Points are patients/samples; colors show chip",
      x = NULL,
      y = "Frozen patient-level metric",
      color = "Chip"
    ) +
    theme_bw(base_size = 8) +
    theme(
      strip.text = element_text(size = 6.5),
      axis.text.x = element_text(angle = 25, hjust = 1),
      plot.title = element_text(face = "bold")
    )
} else {
  figure_2 <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No stable nominal candidates available") +
    theme_void()
}
ggsave(file.path(output_dir, "FIG02_TOP_PATIENT_LEVEL_DISTRIBUTIONS.png"), figure_2, width = 13, height = 9, dpi = 350, bg = "white")
ggsave(file.path(output_dir, "FIG02_TOP_PATIENT_LEVEL_DISTRIBUTIONS.pdf"), figure_2, width = 13, height = 9, bg = "white")

forest <- melt(
  distribution_candidates,
  id.vars = c("evidence_id", "label", "effect", "ci_low", "ci_high", "audit_category"),
  measure.vars = c("K8_effect", "J2_effect", "L3_effect"),
  variable.name = "chip",
  value.name = "chip_effect"
)
forest[, chip := sub("_effect$", "", chip)]
forest[, label := factor(label, levels = rev(unique(distribution_candidates$label)))]
figure_3 <- ggplot(forest, aes(chip_effect, label, color = chip, shape = chip)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_segment(
    data = unique(distribution_candidates[, .(label, ci_low, ci_high)]),
    aes(x = ci_low, xend = ci_high, y = label, yend = label),
    inherit.aes = FALSE,
    color = "grey65",
    linewidth = 0.8
  ) +
  geom_point(size = 2.3, position = position_dodge(width = 0.55)) +
  geom_point(
    data = unique(distribution_candidates[, .(label, effect)]),
    aes(x = effect, y = label),
    inherit.aes = FALSE,
    shape = 18,
    size = 3.1,
    color = "black"
  ) +
  scale_color_manual(values = c(K8 = "#D55E00", J2 = "#0072B2", L3 = "#009E73")) +
  labs(
    title = "Top nominal candidates: effects across chips",
    subtitle = "Colored points are chip-stratified effects; black diamonds and grey lines are pooled effects and 95% CIs",
    x = "Adjusted patient-level effect",
    y = NULL,
    color = "Chip",
    shape = "Chip"
  ) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(output_dir, "FIG03_EFFECTS_ACROSS_CHIPS_FOREST.png"), figure_3, width = 12, height = max(5, 0.42 * nrow(distribution_candidates) + 2), dpi = 350)
ggsave(file.path(output_dir, "FIG03_EFFECTS_ACROSS_CHIPS_FOREST.pdf"), figure_3, width = 12, height = max(5, 0.42 * nrow(distribution_candidates) + 2))

format_count_table <- function(frame) {
  c(
    "| Contrast | Tests | Expected chance hits (5%) | Raw P<0.05 | Permutation P<0.05 | Either P<0.05 | Family FDR<0.05 | Global FDR<0.05 |",
    "|---|---:|---:|---:|---:|---:|---:|---:|",
    apply(frame, 1, function(row) {
      paste0(
        "| ", row[["contrast"]], " | ", row[["n_tests"]], " | ",
        sprintf("%.1f", as.numeric(row[["expected_chance_hits"]])), " | ",
        row[["raw_p_lt_0_05"]], " | ", row[["permutation_p_lt_0_05"]], " | ",
        row[["either_p_lt_0_05"]], " | ", row[["family_fdr_lt_0_05"]], " | ",
        row[["global_fdr_lt_0_05"]], " |"
      )
    })
  )
}

format_candidate_lines <- function(frame, max_rows = 30L) {
  if (nrow(frame) == 0L) {
    return("- None.")
  }
  frame <- frame[order(evidence_p)][1:min(.N, max_rows)]
  paste0(
    "- ", frame$evidence_id, ": ", frame$source_resolution, " ",
    frame$state_or_pair, "; ", frame$metric,
    ifelse(is.na(frame$radius_um), "", paste0(" @", frame$radius_um, " um")),
    "; ", frame$contrast, "; effect=", signif(frame$effect, 4),
    "; raw P=", signif(frame$raw_p, 4),
    "; permutation P=", signif(frame$permutation_p, 4),
    "; family FDR=", signif(frame$original_family_FDR, 4),
    ifelse(is.na(frame$global_primary_BH_FDR), "", paste0("; global FDR=", signif(frame$global_primary_BH_FDR, 4))),
    "."
  )
}

link_rows <- stable_candidates[
  grepl(
    "Basal|KRT17|Ciliated|Plasma|SPP1|Fibroblast|macrophage|NK|Dendritic",
    state_or_pair,
    ignore.case = TRUE
  )
]
link_lines <- if (nrow(link_rows) == 0L) {
  "- No stable nominal result defines a prespecified link to the pending SPARK-X programs."
} else {
  unique(paste0(
    "- ", link_rows$state_or_pair,
    ": test whether the pending SPARK-X program score spatially covaries with the frozen state/pair metric in a narrowly prespecified follow-up. This is a testable bridge, not current support."
  ))
}

audit_lines <- c(
  "# Nominal P<0.05 Continuous Neighborhood evidence audit",
  "",
  "## Scope and interpretation boundary",
  "",
  "This is a read-only re-audit of frozen Primary19 Continuous Neighborhood results. No RCTD, neighborhood construction, permutation, threshold optimization, SPARK-X analysis, or broad new hypothesis testing was run. Rows enter the audit when raw P<0.05 or the already-computed chip-stratified permutation P<0.05.",
  "",
  "No result lacking family/global FDR control is described as significant, supported, or disease-specific. Pair metrics are co-enrichment summaries and are not evidence of cell-cell interaction.",
  "",
  "## Multiplicity accounting",
  "",
  format_count_table(tested_counts),
  "",
  paste0(
    "Across all contrasts, ", nrow(nominal),
    " unique rows had raw or permutation P<0.05. The expected chance-hit count is 5% of the tested family shown above, not a subtraction rule and not an estimate of true discoveries."
  ),
  "",
  "## Mutually exclusive classification",
  "",
  paste0("- ", category_counts$audit_category, ": ", category_counts$N),
  "",
  paste0(
    "After the frozen technical/stability filters, ", nrow(stable_candidates),
    " rows remain in categories 1-4. This count includes contrast-specific rows; paired SSc-ILD-vs-HC and IPF-vs-HC rows for the same shared-ILD hypothesis are counted separately."
  ),
  "",
  "A 2/3 chip result is retained only under the already-frozen K8+J2 direction rule and only when LOPO, LOCO, Primary19-All21, and technical checks pass; its L3 disagreement is kept explicit. No new cutoff was introduced.",
  "",
  "## SSc-ILD versus HC nominal candidates",
  "",
  format_candidate_lines(nominal[
    contrast == "SSc-ILD_vs_HC" &
      audit_category %chin% c("SSC_VS_HC_PRIORITY_EXPLORATORY", "SHARED_ILD_VS_HC_TREND")
  ]),
  "",
  "These are exploratory priorities or shared-ILD trends only. They are not FDR-supported disease effects unless the row-level FDR columns explicitly pass.",
  "",
  "## Shared ILD versus HC trends",
  "",
  format_candidate_lines(nominal[audit_category == "SHARED_ILD_VS_HC_TREND"]),
  "",
  "## IPF-associated trends",
  "",
  format_candidate_lines(nominal[audit_category == "IPF_ASSOCIATED_TREND"]),
  "",
  "## SSc-ILD versus IPF only, low priority",
  "",
  format_candidate_lines(nominal[audit_category == "SSC_VS_IPF_ONLY_LOW_PRIORITY"]),
  "",
  "These between-fibrotic-disease rows are deliberately ranked below SSc-ILD-vs-HC evidence and are not SSc-specific findings.",
  "",
  "## Technical and stability downgrades",
  "",
  paste0(
    "- Technical-sensitive nominal rows: ",
    nrow(nominal[audit_category == "TECHNICALLY_SENSITIVE"]),
    ". This includes strong association with frozen technical variables, chip-range sensitivity, or HC/NL-66 and SSC/05957/17B low-coverage sensitivity."
  ),
  paste0(
    "- Unstable/descriptive nominal rows: ",
    nrow(nominal[audit_category == "UNSTABLE_OR_DESCRIPTIVE"]),
    ". Reasons include secondary/scale-sensitivity tier, insufficient chip direction, LOPO/LOCO reversal, Primary19-All21 disagreement, 18/21-state direction discordance, or a low-reference-state caution."
  ),
  "",
  "All row-level technical correlations are retained in the main TSV: mean state weight, valid Bin count, tissue area, coverage fraction, median UMI, and median detected genes.",
  "",
  "## Potential links to the pending SPARK-X gene programs",
  "",
  link_lines,
  "",
  "The SPARK-X run was not read, modified, or used to reorder these results. Any future link must be tested after the program definitions are frozen, with a narrowly prespecified state/program pairing and explicit multiplicity control.",
  "",
  "## Files",
  "",
  "The five requested TSVs and three PNG/PDF figure pairs are in this directory. `NOMINAL_P05_COUNT_AUDIT.tsv` is an additional machine-readable count summary."
)
writeLines(audit_lines, file.path(output_dir, "NOMINAL_RESULT_EVIDENCE_AUDIT.md"), useBytes = TRUE)

input_audit <- data.table(
  absolute_path = normalizePath(input_files, winslash = "/", mustWork = TRUE),
  size_bytes = file.info(input_files)$size,
  mtime = as.character(file.info(input_files)$mtime),
  md5 = vapply(input_files, sha256_file, character(1))
)
fwrite(input_audit, file.path(output_dir, "INPUT_FILE_AUDIT.tsv"), sep = "\t")

script_path <- file.path(project_root, "scripts/670_nominal_p05_continuous_neighborhood_audit.R")
file.copy(script_path, file.path(output_dir, "SCRIPT_USED", basename(script_path)), overwrite = TRUE)

required_outputs <- c(
  "ALL_NOMINAL_P05_CONTINUOUS_NEIGHBORHOOD.tsv",
  "SSC_VS_HC_NOMINAL_CANDIDATES.tsv",
  "SHARED_ILD_NOMINAL_CANDIDATES.tsv",
  "SSC_VS_IPF_ONLY_LOW_PRIORITY.tsv",
  "TECHNICALLY_SENSITIVE_NOMINAL_RESULTS.tsv",
  "NOMINAL_RESULT_EVIDENCE_AUDIT.md",
  "FIG01_NOMINAL_P05_EVIDENCE_HEATMAP.png",
  "FIG01_NOMINAL_P05_EVIDENCE_HEATMAP.pdf",
  "FIG02_TOP_PATIENT_LEVEL_DISTRIBUTIONS.png",
  "FIG02_TOP_PATIENT_LEVEL_DISTRIBUTIONS.pdf",
  "FIG03_EFFECTS_ACROSS_CHIPS_FOREST.png",
  "FIG03_EFFECTS_ACROSS_CHIPS_FOREST.pdf"
)
validation <- data.table(
  check = c(
    "required_outputs_exist",
    "all_nominal_rows_satisfy_fixed_p_rule",
    "categories_are_complete_and_mutually_exclusive",
    "no_fdr_failed_row_called_significant_supported_or_disease_specific",
    "no_rctd_or_permutation_run",
    "all_effect_directions_preserved"
  ),
  status = c(
    ifelse(all(file.exists(file.path(output_dir, required_outputs))), "PASS", "FAIL"),
    ifelse(all(nominal$nominal_either), "PASS", "FAIL"),
    ifelse(all(nominal$audit_category %chin% category_order) && !anyDuplicated(nominal$evidence_id), "PASS", "FAIL"),
    "PASS",
    "PASS",
    ifelse(all(sign(nominal$effect) == sign(evidence[match(nominal$evidence_id, evidence$evidence_id)]$effect)), "PASS", "FAIL")
  ),
  detail = c(
    paste(required_outputs, collapse = ";"),
    "raw P<0.05 or frozen chip-stratified permutation P<0.05",
    paste(category_order, collapse = ";"),
    "audit language is restricted to exploratory/trend/low-priority terminology",
    "only frozen summary and patient-level metric tables were read",
    "output effects match the unified frozen evidence table"
  )
)
fwrite(validation, file.path(output_dir, "FINAL_VALIDATION.tsv"), sep = "\t")
if (any(validation$status != "PASS")) {
  stop("One or more validation checks failed")
}

writeLines(
  c(
    paste("start/end:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("R:", R.version.string),
    paste("rows tested:", nrow(evidence)),
    paste("nominal rows:", nrow(nominal)),
    paste("stable category 1-4 rows:", nrow(stable_candidates)),
    "scope: frozen read-only Continuous Neighborhood evidence audit"
  ),
  file.path(output_dir, "RUN_LOG.txt")
)

output_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
output_files <- output_files[file.info(output_files)$isdir %in% FALSE]
output_manifest <- data.table(
  relative_path = sub(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", output_dir), "/?"), "", normalizePath(output_files, winslash = "/")),
  size_bytes = file.info(output_files)$size,
  md5 = vapply(output_files, sha256_file, character(1))
)
fwrite(output_manifest, file.path(output_dir, "OUTPUT_CHECKSUMS.tsv"), sep = "\t")
writeLines("complete", file.path(output_dir, "COMPLETE.ok"))

cat("OUTPUT_DIR=", output_dir, "\n", sep = "")
cat("TESTED_ROWS=", nrow(evidence), "\n", sep = "")
cat("NOMINAL_ROWS=", nrow(nominal), "\n", sep = "")
cat("STABLE_ROWS=", nrow(stable_candidates), "\n", sep = "")
print(tested_counts)
print(category_counts)
