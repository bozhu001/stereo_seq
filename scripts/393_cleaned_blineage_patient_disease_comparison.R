#!/usr/bin/env Rscript

# Exploratory patient-level disease comparison using immutable frozen outputs.
# No Bin-level scoring, marker selection, threshold optimization, or spatial
# calling is performed here. The patient is the sole statistical unit.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(sandwich)
  library(jsonlite)
  library(sessioninfo)
})

options(stringsAsFactors = FALSE)
set.seed(260823L)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix, default = NULL) {
  hit <- args[startsWith(args, paste0("--", prefix, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", prefix, "="), "", hit[[1]])
}

script_path <- normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]),
  winslash = "/", mustWork = TRUE
)
repo <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
tstamp <- arg_value("timestamp", format(Sys.time(), "%y%m%d%H%M%S"))
stopifnot(grepl("^[0-9]{12}$", tstamp))

source_dir <- file.path(repo, "results/reanalysis/bin50_hq_frozen_blineage_descriptive_qc_260823180725")
source_tables <- file.path(source_dir, "tables")
out_dir <- file.path(repo, paste0("results/reanalysis/bin50_cleaned_blineage_patient_disease_comparison_", tstamp))
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
log_dir <- file.path(out_dir, "logs")
script_dir <- file.path(out_dir, "scripts")
report_dir <- file.path(out_dir, "report")
for (path in c(table_dir, figure_dir, log_dir, script_dir, report_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(log_dir, "analysis.log")
log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_path, append = TRUE)
}
assert <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)
write_tsv <- function(x, name) fwrite(x, file.path(table_dir, name), sep = "\t", na = "NA")
save_plot <- function(plot, name, width, height) {
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 180, bg = "white")
  ggsave(file.path(figure_dir, paste0(name, ".svg")), plot, width = width, height = height, bg = "white")
}
inv_logit <- function(x) 1 / (1 + exp(-x))

panel_names <- c(
  "Mature_B_core", "ILD_plasma_cell", "IgA_mucosal_plasma",
  "GC_activated_B", "Naive_B", "Memory_tissue_B"
)
group_levels <- c("HC", "IPF", "SSc-ILD")
contrast_specs <- list(
  "SSc-ILD vs HC" = c("SSc-ILD", "HC"),
  "IPF vs HC" = c("IPF", "HC"),
  "SSc-ILD vs IPF" = c("SSc-ILD", "IPF")
)
primary_samples <- c(
  "HC/NL-50", "HC/NL-55", "HC/NL-63", "HC/NL-72",
  "IPF/FO23-1-03100", "IPF/FO23-1-09473", "IPF/FO23-1-06168", "IPF/FO23-1-06170",
  "SSC/24-1-18170A2", "SSC/23-105334B1", "SSC/FO20-I-08293B1",
  "SSC/FO21-1-01756B1", "SSC/05957/17B", "SSC/07998/15A", "SSC/15491/14"
)
k8_samples <- c("HC/NL-50", "SSC/24-1-18170A2")
j2_l3_samples <- setdiff(primary_samples, k8_samples)
watch_excluded <- c("HC/NL-53", "HC/NL-66", "IPF/FO23-1-06474", "SSC/15275/13", "SSC/3342/13")
depth_excluded <- "IPF/FO22-1-09404"

input_files <- c(
  panel_manifest = file.path(source_tables, "00_frozen_panel_manifest.tsv"),
  sample_eligibility = file.path(source_tables, "01_sample_qc_eligibility.tsv"),
  panel_metrics = file.path(source_tables, "05_sample_panel_detected_marker_counts.tsv"),
  same_bin = file.path(source_tables, "06_panel_same_bin_strict_support.tsv"),
  rook = file.path(source_tables, "07_panel_rook_strict_support.tsv"),
  metadata = file.path(source_tables, "10_high_quality_sample_descriptive_metadata.tsv"),
  validation = file.path(source_tables, "12_output_validation.tsv")
)
assert(all(file.exists(input_files)), paste("Missing frozen input:", paste(input_files[!file.exists(input_files)], collapse = ", ")))
log_message("Starting immutable frozen-output audit; source=", source_dir)

panel_manifest_all <- fread(input_files[["panel_manifest"]])
panel_manifest <- panel_manifest_all[panel %chin% panel_names]
same_bin <- fread(input_files[["same_bin"]])[panel %chin% panel_names]
rook <- fread(input_files[["rook"]])[panel %chin% panel_names]
panel_metrics <- fread(input_files[["panel_metrics"]])[panel %chin% panel_names]
metadata <- fread(input_files[["metadata"]])
source_validation <- fread(input_files[["validation"]])

assert(nrow(metadata) == 16L && uniqueN(metadata$sample_id) == 16L, "Frozen metadata is not 16 unique samples")
assert(nrow(same_bin) == 96L && nrow(rook) == 96L, "Frozen strict support grids are not 16 x 6")
assert(nrow(panel_metrics) == 96L, "Frozen panel metric grid is not 16 x 6")
assert(all(source_validation$passed), "Frozen source validation contains a failure")
assert(setequal(metadata$sample_id, c(primary_samples, depth_excluded)), "Frozen 16 sample membership drift")
assert(setequal(primary_samples, setdiff(metadata$sample_id, depth_excluded)), "Primary15 sample membership drift")
assert(setequal(j2_l3_samples, primary_samples[!primary_samples %chin% k8_samples]), "J2/L3-only13 membership drift")
assert(metadata[sample_id == depth_excluded, chip_id] == "Y40102K8", "FO22 chip mismatch")
assert(all(metadata[sample_id %chin% k8_samples, chip_id] == "Y40102K8"), "Remaining K8 identity mismatch")

depth_cols <- c("sample_total_counts", "median_total_counts", "median_n_genes")
depth_rank <- rbindlist(lapply(depth_cols, function(metric) {
  metadata[, .(
    metric = metric,
    sample_id,
    value = get(metric),
    rank_desc = frank(-get(metric), ties.method = "min"),
    excluded_as_extreme_depth = sample_id == depth_excluded
  )]
}))
assert(all(depth_rank[sample_id == depth_excluded, rank_desc] == 1L), "FO22 is not highest on all three frozen depth metrics")

patient <- merge(
  same_bin[, .(
    sample_id, chip_id, group, panel, total_bin_n,
    strict_same_bin_supported_bin_n,
    strict_same_bin_fraction = strict_same_bin_supported_fraction_all_tissue_bins
  )],
  rook[, .(
    sample_id, panel,
    strict_rook_supported_bin_n,
    strict_rook_fraction = strict_rook_supported_fraction_all_tissue_bins
  )],
  by = c("sample_id", "panel"), all = FALSE
)
patient <- merge(
  patient,
  panel_metrics[, .(sample_id, panel, panel_raw_count_total, panel_cpm_total)],
  by = c("sample_id", "panel"), all = FALSE
)
patient <- merge(
  patient,
  metadata[, .(sample_id, sample_total_counts, median_total_counts, median_n_genes)],
  by = "sample_id", all = FALSE
)
assert(nrow(patient) == 96L && uniqueN(patient, by = c("sample_id", "panel")) == 96L, "Merged frozen patient-panel grid drift")
patient[, disease_group := factor(group, levels = group_levels)]
patient[, panel := factor(panel, levels = panel_names)]
patient[, p_corrected := (strict_same_bin_supported_bin_n + 0.5) / (total_bin_n + 1)]
patient[, transformed_endpoint := log(p_corrected / (1 - p_corrected))]
patient[, `:=`(
  same_bin_zero = strict_same_bin_supported_bin_n == 0,
  same_bin_one_bin_only = strict_same_bin_supported_bin_n == 1,
  strict_rook_zero = strict_rook_supported_bin_n == 0,
  panel_cpm_zero = panel_cpm_total == 0,
  extremely_sparse_predefined = strict_same_bin_supported_bin_n <= 1
)]
patient[, cohort_membership := fifelse(
  sample_id %chin% j2_l3_samples, "primary15_and_j2_l3_only13",
  fifelse(sample_id %chin% k8_samples, "primary15_only_K8", "excluded_extreme_depth")
)]

primary <- patient[sample_id %chin% primary_samples]
sensitivity <- patient[sample_id %chin% j2_l3_samples]
assert(nrow(primary) == 90L && uniqueN(primary$sample_id) == 15L, "Primary15 grid drift")
assert(nrow(sensitivity) == 78L && uniqueN(sensitivity$sample_id) == 13L, "Sensitivity13 grid drift")
assert(identical(as.integer(primary[, uniqueN(sample_id), by = group][match(group_levels, group), V1]), c(4L, 4L, 7L)), "Primary group counts drift")
assert(identical(as.integer(sensitivity[, uniqueN(sample_id), by = group][match(group_levels, group), V1]), c(3L, 4L, 6L)), "Sensitivity group counts drift")

source_manifest <- data.table(
  item = c(
    "analysis_label", "source_directory", "statistical_unit", "primary_endpoint",
    "transformation", "primary_model", "robust_se", "primary_multiplicity",
    "sensitivity_model", "panel_scope", "extreme_sparse_flag", "result_adaptation"
  ),
  frozen_value = c(
    "exploratory cleaned B-lineage patient-level disease comparison",
    source_dir, "patient", "strict same-bin co-capture fraction",
    "logit((supported_bin_n + 0.5)/(total_bin_n + 1))",
    "transformed_endpoint ~ disease_group + log(sample_total_counts) + chip_id",
    "HC3", "single BH adjustment over 6 panels x 3 contrasts = 18 tests",
    "same formula in J2/L3-only13", paste(panel_names, collapse = ";"),
    "strict_same_bin_supported_bin_n <= 1; descriptive flag only",
    "none: no sample, marker, threshold, panel, or model changes based on results"
  )
)
write_tsv(source_manifest, "00_data_source_and_frozen_rules_manifest.tsv")
write_tsv(panel_manifest, "01_frozen_six_panel_marker_manifest.tsv")

sample_list <- unique(patient[, .(sample_id, disease_group = group, chip_id, sample_total_counts, median_total_counts, median_n_genes, cohort_membership)])
sample_list[, `:=`(
  primary15_included = sample_id %chin% primary_samples,
  j2_l3_only13_included = sample_id %chin% j2_l3_samples,
  exclusion_reason_here = fifelse(sample_id == depth_excluded, "extreme technical depth", fifelse(sample_id %chin% k8_samples, "excluded only from J2/L3 sensitivity because K8", "included"))
)]
write_tsv(sample_list[primary15_included == TRUE], "02_primary15_sample_list.tsv")
write_tsv(sample_list[j2_l3_only13_included == TRUE], "03_j2_l3_only13_sample_list.tsv")
write_tsv(depth_rank, "04_extreme_depth_exclusion_audit.tsv")

patient_out <- rbind(
  copy(primary)[, cohort := "primary15"],
  copy(sensitivity)[, cohort := "j2_l3_only13"]
)
setcolorder(patient_out, c(
  "cohort", "sample_id", "group", "chip_id", "sample_total_counts", "median_total_counts",
  "median_n_genes", "panel", "total_bin_n", "strict_same_bin_supported_bin_n",
  "strict_same_bin_fraction", "p_corrected", "transformed_endpoint",
  "strict_rook_supported_bin_n", "strict_rook_fraction", "panel_raw_count_total", "panel_cpm_total",
  "same_bin_zero", "same_bin_one_bin_only", "strict_rook_zero", "panel_cpm_zero",
  "extremely_sparse_predefined"
))
write_tsv(patient_out, "05_patient_level_raw_metrics.tsv")

metric_cols <- c("strict_same_bin_fraction", "strict_rook_fraction", "panel_cpm_total")
group_descriptive <- melt(
  patient_out,
  id.vars = c("cohort", "sample_id", "group", "chip_id", "panel"),
  measure.vars = metric_cols,
  variable.name = "metric", value.name = "value"
)[, .(
  patient_n = .N,
  median = median(value),
  q1 = quantile(value, 0.25, names = FALSE),
  q3 = quantile(value, 0.75, names = FALSE),
  iqr = IQR(value),
  minimum = min(value),
  maximum = max(value),
  maximum_patient = paste(sort(sample_id[value == max(value)]), collapse = ";"),
  zero_patient_n = sum(value == 0)
), by = .(cohort, panel, disease_group = group, metric)]
write_tsv(group_descriptive, "06_group_descriptive_statistics.tsv")

fit_panel <- function(dat, cohort_name, panel_name) {
  d <- copy(dat[as.character(panel) == panel_name])
  d[, disease_group := factor(as.character(disease_group), levels = group_levels)]
  d[, chip_id := factor(chip_id)]
  fit <- lm(transformed_endpoint ~ disease_group + log(sample_total_counts) + chip_id, data = d)
  beta <- coef(fit)
  # The primary15 model has six parameters because chip has three levels;
  # the J2/L3-only13 model has five because chip has two levels.
  model_stable <- all(is.finite(beta)) && qr(fit)$rank == length(beta)
  robust_vcov <- tryCatch(vcovHC(fit, type = "HC3"), error = function(e) NULL)
  model_stable <- model_stable && !is.null(robust_vcov) && all(is.finite(robust_vcov))
  p <- length(beta)
  n <- nrow(d)
  df <- df.residual(fit)
  tcrit <- qt(0.975, df = df)
  condition_number <- kappa(model.matrix(fit), exact = TRUE)

  prediction_rows <- list()
  L_by_group <- list()
  if (model_stable) {
    for (g in group_levels) {
      nd <- copy(d)
      nd[, disease_group := factor(g, levels = group_levels)]
      X <- model.matrix(delete.response(terms(fit)), nd, contrasts.arg = fit$contrasts)
      X <- X[, names(beta), drop = FALSE]
      L <- colMeans(X)
      L_by_group[[g]] <- L
      eta <- as.numeric(crossprod(L, beta))
      eta_se <- sqrt(as.numeric(t(L) %*% robust_vcov %*% L))
      prediction_rows[[g]] <- data.table(
        cohort = cohort_name, panel = panel_name, disease_group = g,
        adjusted_logit = eta, adjusted_logit_se_hc3 = eta_se,
        adjusted_logit_ci_low = eta - tcrit * eta_se,
        adjusted_logit_ci_high = eta + tcrit * eta_se,
        adjusted_fraction = inv_logit(eta),
        adjusted_fraction_ci_low = inv_logit(eta - tcrit * eta_se),
        adjusted_fraction_ci_high = inv_logit(eta + tcrit * eta_se),
        standardization = "average over observed cohort technical covariates"
      )
    }
  }

  contrast_rows <- list()
  for (contrast_name in names(contrast_specs)) {
    spec <- contrast_specs[[contrast_name]]
    if (!model_stable) {
      contrast_rows[[contrast_name]] <- data.table(
        cohort = cohort_name, panel = panel_name, endpoint_role = ifelse(panel_name == "Mature_B_core", "primary biological endpoint", "secondary exploratory endpoint"),
        contrast = contrast_name, numerator_group = spec[1], denominator_group = spec[2],
        effect_logit = NA_real_, se_hc3 = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
        adjusted_fraction_numerator = NA_real_, adjusted_fraction_denominator = NA_real_, adjusted_fraction_difference = NA_real_,
        fraction_difference_ci_low = NA_real_, fraction_difference_ci_high = NA_real_, model_stable = FALSE
      )
      next
    }
    L1 <- L_by_group[[spec[1]]]
    L0 <- L_by_group[[spec[2]]]
    C <- L1 - L0
    effect <- as.numeric(crossprod(C, beta))
    se <- sqrt(as.numeric(t(C) %*% robust_vcov %*% C))
    t_value <- effect / se
    p_value <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
    eta1 <- as.numeric(crossprod(L1, beta))
    eta0 <- as.numeric(crossprod(L0, beta))
    frac1 <- inv_logit(eta1)
    frac0 <- inv_logit(eta0)
    frac_diff <- frac1 - frac0
    grad <- frac1 * (1 - frac1) * L1 - frac0 * (1 - frac0) * L0
    frac_se <- sqrt(as.numeric(t(grad) %*% robust_vcov %*% grad))
    contrast_rows[[contrast_name]] <- data.table(
      cohort = cohort_name, panel = panel_name,
      endpoint_role = ifelse(panel_name == "Mature_B_core", "primary biological endpoint", "secondary exploratory endpoint"),
      contrast = contrast_name, numerator_group = spec[1], denominator_group = spec[2],
      effect_logit = effect, se_hc3 = se,
      ci_low = effect - tcrit * se, ci_high = effect + tcrit * se,
      p_value = p_value,
      adjusted_fraction_numerator = frac1, adjusted_fraction_denominator = frac0,
      adjusted_fraction_difference = frac_diff,
      fraction_difference_ci_low = max(-1, frac_diff - tcrit * frac_se),
      fraction_difference_ci_high = min(1, frac_diff + tcrit * frac_se),
      model_stable = TRUE
    )
  }

  infl <- data.table(
    cohort = cohort_name, panel = panel_name, sample_id = d$sample_id,
    disease_group = as.character(d$disease_group), chip_id = as.character(d$chip_id),
    leverage = hatvalues(fit), cooks_distance = cooks.distance(fit),
    studentized_residual = rstudent(fit)
  )
  infl[, `:=`(
    high_leverage = leverage > 2 * p / n,
    high_cook = cooks_distance > 4 / n,
    large_studentized_residual = abs(studentized_residual) > 2
  )]
  infl[, flagged_any := high_leverage | high_cook | large_studentized_residual]
  diag <- data.table(
    cohort = cohort_name, panel = panel_name, n = n, parameter_n = p,
    residual_df = df, model_rank = qr(fit)$rank, full_rank = qr(fit)$rank == p,
    condition_number = condition_number, hc3_finite = !is.null(robust_vcov) && all(is.finite(robust_vcov)),
    model_stable = model_stable, high_leverage_n = sum(infl$high_leverage),
    high_cook_n = sum(infl$high_cook), large_studentized_residual_n = sum(infl$large_studentized_residual)
  )
  list(fit = fit, contrasts = rbindlist(contrast_rows), predictions = rbindlist(prediction_rows), diagnostics = diag, influence = infl)
}

fits_primary <- lapply(panel_names, function(x) fit_panel(primary, "primary15", x))
fits_sensitivity <- lapply(panel_names, function(x) fit_panel(sensitivity, "j2_l3_only13", x))
primary_results <- rbindlist(lapply(fits_primary, `[[`, "contrasts"))
sensitivity_results <- rbindlist(lapply(fits_sensitivity, `[[`, "contrasts"))
assert(nrow(primary_results) == 18L && nrow(sensitivity_results) == 18L, "Expected 18 contrasts per cohort")
primary_results[, bh_fdr_18 := p.adjust(p_value, method = "BH")]
primary_results[, bh_fdr_family := "all 6 panels x 3 contrasts; n=18"]
write_tsv(primary_results, "07_primary15_model_results_with_global_bh_fdr.tsv")
write_tsv(primary_results[, .(panel, endpoint_role, contrast, effect_logit, ci_low, ci_high, p_value, bh_fdr_18)], "08_primary15_18_test_bh_fdr_results.tsv")
write_tsv(sensitivity_results, "09_j2_l3_only13_sensitivity_model_results.tsv")

predictions <- rbindlist(c(
  lapply(fits_primary, `[[`, "predictions"),
  lapply(fits_sensitivity, `[[`, "predictions")
))
write_tsv(predictions, "10_adjusted_group_predictions.tsv")

stability <- merge(
  primary_results[, .(panel, contrast, effect_primary15 = effect_logit, ci_low_primary15 = ci_low, ci_high_primary15 = ci_high, p_primary15 = p_value, bh_fdr_primary15 = bh_fdr_18, stable_primary15 = model_stable)],
  sensitivity_results[, .(panel, contrast, effect_j2_l3_only13 = effect_logit, ci_low_j2_l3_only13 = ci_low, ci_high_j2_l3_only13 = ci_high, p_j2_l3_only13 = p_value, stable_j2_l3_only13 = model_stable)],
  by = c("panel", "contrast")
)
stability[, `:=`(
  direction_primary15 = fifelse(effect_primary15 > 0, "positive", fifelse(effect_primary15 < 0, "negative", "zero")),
  direction_j2_l3_only13 = fifelse(effect_j2_l3_only13 > 0, "positive", fifelse(effect_j2_l3_only13 < 0, "negative", "zero")),
  absolute_effect_ratio_13_to_15 = abs(effect_j2_l3_only13) / abs(effect_primary15)
)]
stability[, stability_class := fifelse(
  !stable_primary15 | !stable_j2_l3_only13 | !is.finite(effect_primary15) | !is.finite(effect_j2_l3_only13), "model_unstable",
  fifelse(sign(effect_primary15) != sign(effect_j2_l3_only13), "direction_flipped",
    fifelse(absolute_effect_ratio_13_to_15 < 0.8, "direction_consistent_noticeably_weakened",
      fifelse(absolute_effect_ratio_13_to_15 > 1.25, "direction_consistent_noticeably_strengthened", "direction_consistent_similar")
    )
  )
)]
write_tsv(stability, "11_primary15_vs_j2_l3_only13_direction_stability.tsv")

diagnostics <- rbindlist(c(lapply(fits_primary, `[[`, "diagnostics"), lapply(fits_sensitivity, `[[`, "diagnostics")))
influence <- rbindlist(c(lapply(fits_primary, `[[`, "influence"), lapply(fits_sensitivity, `[[`, "influence")))
write_tsv(diagnostics, "12_model_diagnostics.tsv")
write_tsv(influence, "13_patient_influence_diagnostics.tsv")

loo_one <- function(dat, panel_name, cohort_name) {
  full <- fit_panel(dat, cohort_name, panel_name)$contrasts[, .(contrast, full_effect = effect_logit)]
  rbindlist(lapply(unique(dat$sample_id), function(drop_id) {
    reduced <- tryCatch(fit_panel(dat[sample_id != drop_id], paste0(cohort_name, "_loo"), panel_name)$contrasts, error = function(e) NULL)
    if (is.null(reduced)) {
      return(data.table(panel = panel_name, cohort = cohort_name, omitted_sample = drop_id, contrast = names(contrast_specs), loo_effect = NA_real_, loo_model_stable = FALSE))
    }
    reduced[, .(panel = panel_name, cohort = cohort_name, omitted_sample = drop_id, contrast, loo_effect = effect_logit, loo_model_stable = model_stable)]
  }))[, full_effect := full$full_effect[match(contrast, full$contrast)]]
}
loo_results <- rbindlist(c(
  lapply(panel_names, function(x) loo_one(primary, x, "primary15")),
  lapply(panel_names, function(x) loo_one(sensitivity, x, "j2_l3_only13"))
))
loo_results[, `:=`(
  direction_flip_vs_full = is.finite(loo_effect) & is.finite(full_effect) & sign(loo_effect) != sign(full_effect),
  absolute_change = abs(loo_effect - full_effect),
  relative_absolute_change = abs(loo_effect - full_effect) / pmax(abs(full_effect), 1e-8)
)]
write_tsv(loo_results, "14_leave_one_patient_out_effects.tsv")
loo_summary <- loo_results[, .(
  loo_refit_n = .N,
  unstable_refit_n = sum(!loo_model_stable | !is.finite(loo_effect)),
  direction_flip_n = sum(direction_flip_vs_full, na.rm = TRUE),
  max_absolute_change = max(absolute_change, na.rm = TRUE),
  max_relative_absolute_change = max(relative_absolute_change, na.rm = TRUE),
  largest_change_patient = paste(sort(omitted_sample[absolute_change == max(absolute_change, na.rm = TRUE)]), collapse = ";")
), by = .(cohort, panel, contrast)]
write_tsv(loo_summary, "15_leave_one_patient_out_summary.tsv")

direction_data <- melt(
  patient_out,
  id.vars = c("cohort", "sample_id", "group", "panel"),
  measure.vars = c("strict_same_bin_fraction", "strict_rook_fraction", "panel_cpm_total"),
  variable.name = "metric", value.name = "value"
)[, .(group_median = median(value)), by = .(cohort, panel, metric, group)]
direction_rows <- rbindlist(lapply(names(contrast_specs), function(contrast_name) {
  spec <- contrast_specs[[contrast_name]]
  wide <- dcast(direction_data, cohort + panel + metric ~ group, value.var = "group_median")
  wide[, .(
    cohort, panel, metric, contrast = contrast_name,
    numerator_median = get(spec[1]), denominator_median = get(spec[2]),
    median_difference = get(spec[1]) - get(spec[2]),
    descriptive_direction = fifelse(get(spec[1]) - get(spec[2]) > 0, "positive", fifelse(get(spec[1]) - get(spec[2]) < 0, "negative", "zero"))
  )]
}))
direction_rows[, same_bin_direction := descriptive_direction[metric == "strict_same_bin_fraction"][match(paste(cohort, panel, contrast), paste(cohort[metric == "strict_same_bin_fraction"], panel[metric == "strict_same_bin_fraction"], contrast[metric == "strict_same_bin_fraction"]))]]
direction_rows[, direction_matches_same_bin := descriptive_direction == same_bin_direction]
write_tsv(direction_rows, "16_same_bin_rook_cpm_descriptive_direction_consistency.tsv")

theme_set(theme_bw(base_size = 10))
group_colors <- c(HC = "#3B82F6", IPF = "#E45756", `SSc-ILD` = "#7A5195")
chip_shapes <- c(Y40102K8 = 17, Y40105J2 = 16, Y40105L3 = 15)

scatter_plot <- function(dat, title) {
  ggplot(dat, aes(x = disease_group, y = strict_same_bin_fraction, color = disease_group, shape = chip_id)) +
    geom_jitter(width = 0.13, height = 0, size = 2.5, alpha = 0.9) +
    facet_wrap(~panel, scales = "free_y", ncol = 3) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    scale_shape_manual(values = chip_shapes, drop = FALSE) +
    labs(title = title, x = NULL, y = "Strict same-bin co-capture fraction", color = "Disease", shape = "Chip") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
}
save_plot(scatter_plot(primary, "Primary cohort (n=15): patient-level strict same-bin fractions"), "01_primary15_patient_scatter", 12, 7.5)
save_plot(scatter_plot(sensitivity, "J2/L3-only sensitivity cohort (n=13)"), "02_j2_l3_only13_patient_scatter", 12, 7.5)

forest_data <- copy(primary_results)
forest_data[, contrast := factor(contrast, levels = rev(names(contrast_specs)))]
forest_data[, panel := factor(panel, levels = rev(panel_names))]
p_forest <- ggplot(forest_data, aes(x = effect_logit, y = contrast, color = bh_fdr_18 < 0.05)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = 2) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.16) +
  geom_point(size = 2.2) +
  facet_wrap(~panel, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c(`TRUE` = "#B22222", `FALSE` = "#333333"), labels = c(`TRUE` = "BH-FDR < 0.05", `FALSE` = "BH-FDR >= 0.05")) +
  labs(title = "Primary15 adjusted disease contrasts", subtitle = "HC3 95% CI; one BH correction across all 18 tests", x = "Adjusted effect on corrected-logit scale", y = NULL, color = NULL) +
  theme(legend.position = "bottom")
save_plot(p_forest, "03_primary15_forest_18_contrasts", 11, 9)

effect_compare <- melt(
  stability,
  id.vars = c("panel", "contrast", "stability_class"),
  measure.vars = c("effect_primary15", "effect_j2_l3_only13"),
  variable.name = "cohort", value.name = "effect"
)
effect_compare[, cohort := factor(cohort, levels = c("effect_primary15", "effect_j2_l3_only13"), labels = c("Primary15", "J2/L3-only13"))]
p_compare <- ggplot(effect_compare, aes(x = cohort, y = effect, group = interaction(panel, contrast), color = stability_class)) +
  geom_hline(yintercept = 0, color = "grey65", linetype = 2) +
  geom_line(alpha = 0.75) + geom_point(size = 2) +
  facet_grid(panel ~ contrast, scales = "free_y") +
  labs(title = "Effect direction comparison: Primary15 versus J2/L3-only13", x = NULL, y = "Adjusted logit-scale effect", color = "Classification") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
save_plot(p_compare, "04_primary15_vs_j2_l3_only13_effect_direction", 14, 12)

direction_plot_data <- copy(direction_rows)
direction_plot_data[, sign_value := sign(median_difference)]
p_direction <- ggplot(direction_plot_data, aes(x = metric, y = panel, fill = sign_value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2g", median_difference)), size = 2.5) +
  facet_grid(cohort ~ contrast) +
  scale_fill_gradient2(low = "#3B82F6", mid = "white", high = "#E45756", midpoint = 0) +
  labs(title = "Descriptive direction consistency across same-bin, rook and CPM", subtitle = "Cells show disease-group median differences; no additional formal tests", x = NULL, y = NULL, fill = "Direction") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p_direction, "05_same_bin_rook_cpm_direction_consistency", 15, 9)

p_influence <- ggplot(influence, aes(x = leverage, y = studentized_residual, size = cooks_distance, color = flagged_any, label = ifelse(flagged_any, sample_id, ""))) +
  geom_hline(yintercept = c(-2, 2), linetype = 2, color = "grey65") +
  geom_point(alpha = 0.75) +
  geom_text(check_overlap = TRUE, size = 2.2, vjust = -0.6) +
  facet_grid(cohort ~ panel, scales = "free") +
  scale_color_manual(values = c(`TRUE` = "#B22222", `FALSE` = "#555555")) +
  labs(title = "Patient influence diagnostics", subtitle = "Flags are descriptive; no patient was removed", x = "Leverage", y = "Studentized residual", size = "Cook distance", color = "Flagged") +
  theme(legend.position = "bottom")
save_plot(p_influence, "06_model_influence_diagnostics", 16, 8.5)

panel_consistency <- stability[, .(
  contrast_n = .N,
  direction_consistent_n = sum(grepl("^direction_consistent", stability_class)),
  direction_flipped_n = sum(stability_class == "direction_flipped"),
  model_unstable_n = sum(stability_class == "model_unstable"),
  all_three_direction_consistent = all(grepl("^direction_consistent", stability_class))
), by = panel]
write_tsv(panel_consistency, "17_panel_level_direction_consistency_summary.tsv")

k8_drive_summary <- stability[, .(
  contrast_n = .N,
  direction_flip_after_k8_removal_n = sum(stability_class == "direction_flipped"),
  weakened_after_k8_removal_n = sum(stability_class == "direction_consistent_noticeably_weakened"),
  strengthened_after_k8_removal_n = sum(stability_class == "direction_consistent_noticeably_strengthened"),
  similar_after_k8_removal_n = sum(stability_class == "direction_consistent_similar")
), by = panel]
write_tsv(k8_drive_summary, "18_k8_sensitivity_summary.tsv")

fdr_count <- primary_results[bh_fdr_18 < 0.05, .N]
single_patient_flags <- loo_summary[cohort == "primary15" & (direction_flip_n > 0 | max_relative_absolute_change > 1)]
summary_json <- list(
  analysis_label = "extreme-depth-cleaned patient-level exploratory disease comparison",
  timestamp = tstamp,
  primary_n = 15,
  sensitivity_n = 13,
  primary_group_counts = as.list(primary[, uniqueN(sample_id), by = group][, setNames(V1, group)]),
  sensitivity_group_counts = as.list(sensitivity[, uniqueN(sample_id), by = group][, setNames(V1, group)]),
  bh_fdr_lt_0_05_n = fdr_count,
  all_three_contrasts_direction_consistent_panels = panel_consistency[all_three_direction_consistent == TRUE, as.character(panel)],
  any_direction_flip_panels = panel_consistency[direction_flipped_n > 0, as.character(panel)],
  primary_loo_sensitive_panel_contrast_n = nrow(single_patient_flags),
  no_bin_level_recalculation = TRUE,
  eln_created = FALSE,
  git_commit_created = FALSE,
  git_push_performed = FALSE
)
write_json(summary_json, file.path(out_dir, "analysis_summary.json"), pretty = TRUE, auto_unbox = TRUE)

validation <- data.table(
  check = c(
    "frozen_source_validation_passed", "primary15_exact_membership", "primary15_group_counts_4_4_7",
    "j2_l3_only13_exact_membership", "j2_l3_group_counts_3_4_6", "six_frozen_panels_only",
    "primary18_contrasts", "sensitivity18_contrasts", "global_bh_over_18",
    "all_models_full_rank", "all_hc3_covariances_finite", "png_svg_pairs_present",
    "no_bin_level_recalculation"
  ),
  passed = c(
    all(source_validation$passed), setequal(unique(primary$sample_id), primary_samples),
    identical(as.integer(primary[, uniqueN(sample_id), by = group][match(group_levels, group), V1]), c(4L, 4L, 7L)),
    setequal(unique(sensitivity$sample_id), j2_l3_samples),
    identical(as.integer(sensitivity[, uniqueN(sample_id), by = group][match(group_levels, group), V1]), c(3L, 4L, 6L)),
    setequal(unique(as.character(patient$panel)), panel_names), nrow(primary_results) == 18L,
    nrow(sensitivity_results) == 18L, length(primary_results$bh_fdr_18) == 18L,
    all(diagnostics$full_rank), all(diagnostics$hc3_finite),
    length(list.files(figure_dir, pattern = "\\.png$")) == 6L && length(list.files(figure_dir, pattern = "\\.svg$")) == 6L,
    TRUE
  )
)
write_tsv(validation, "19_output_validation.tsv")
assert(all(validation$passed), paste("Validation failed:", paste(validation[passed == FALSE, check], collapse = ", ")))

writeLines(capture.output(session_info()), file.path(log_dir, "session_info.txt"))
file.copy(script_path, file.path(script_dir, basename(script_path)), overwrite = TRUE)
writeLines(c(
  "SUCCESS",
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "scope=patient-level only; frozen Bin-level outputs reused",
  "primary_n=15",
  "j2_l3_only_n=13",
  paste0("bh_fdr_lt_0.05_n=", fdr_count),
  "eln_created=false",
  "git_commit=false",
  "git_push=false"
), file.path(out_dir, "SUCCESS"))
log_message("Completed successfully; BH-FDR<0.05 count=", fdr_count)
cat(out_dir, "\n")
