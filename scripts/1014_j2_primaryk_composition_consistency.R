#!/usr/bin/env Rscript

# J2 primary-K patient-specific Reference18 neighborhood composition analysis.
# Reuses the independently selected J2 primary-K assignment and saved radius50
# VoltRon graph/Niche assay from script 1013.
# Frozen lymphoid aggregate labels are intentionally not read.

suppressPackageStartupMessages({
  library(VoltRon)
  library(data.table)
  library(ggplot2)
  library(igraph)
})

options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: 1014_j2_primaryk_composition_consistency.R ",
    "<existing_j2_dir> <resolution_dir>"
  )
}

existing_dir <- normalizePath(args[[1]], mustWork = TRUE)
resolution_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- resolution_dir
fig_dir <- file.path(output_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

status_path <- file.path(output_dir, "RUN_STATUS.tsv")
write_status <- function(stage, detail = "") {
  fwrite(
    data.table(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stage = stage,
      detail = detail
    ),
    status_path,
    sep = "\t",
    append = file.exists(status_path)
  )
  message("[J2 primary-K composition] ", stage, if (nzchar(detail)) paste0(": ", detail))
}

read_gz <- function(path) {
  fread(cmd = paste("gzip -dc", shQuote(path)))
}

save_plot <- function(plot, stem, width, height, dpi = 320) {
  ggsave(
    file.path(fig_dir, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf,
    limitsize = FALSE
  )
  ggsave(
    file.path(fig_dir, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = FALSE
  )
}

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

same_sign <- function(x, y, tolerance = 1e-12) {
  if (!is.finite(x) || !is.finite(y) || abs(x) <= tolerance || abs(y) <= tolerance) {
    return(NA)
  }
  sign(x) == sign(y)
}

write_status("read_fixed_inputs")
input_path <- file.path(existing_dir, "01_FROZEN_REFERENCE18_INPUT.tsv.gz")
object_path <- file.path(resolution_dir, "14_VOLTRON_OBJECT_RADIUS50_CLR.rds")
recommendation <- fread(file.path(resolution_dir, "04_J2_AUDITED_K_RECOMMENDATION.tsv"))
primary_k <- recommendation[role == "primary", k]
stopifnot(length(primary_k) == 1L)
assignment_path <- file.path(
  resolution_dir,
  sprintf("K%02d/02_BIN_NICHE_ASSIGNMENTS.tsv.gz", primary_k)
)

input <- read_gz(input_path)
object <- readRDS(object_path)
assignment <- read_gz(assignment_path)

weight_columns <- grep("^weight__", names(input), value = TRUE)
states <- sub("^weight__", "", weight_columns)
niche_levels <- paste0("N", seq_len(primary_k))

stopifnot(
  nrow(input) == 58840L,
  length(states) == 18L,
  uniqueN(input$patient_id) == 7L,
  nrow(assignment) == 58772L,
  !anyDuplicated(assignment$internal_id),
  setequal(unique(assignment$niche), niche_levels),
  all(assignment$k == primary_k)
)

object_metadata_df <- as.data.frame(Metadata(object))
internal_id <- if ("id" %in% names(object_metadata_df)) {
  as.character(object_metadata_df$id)
} else {
  rownames(object_metadata_df)
}
object_metadata <- as.data.table(object_metadata_df)
object_metadata[, internal_id := internal_id]
entity_to_internal <- setNames(object_metadata$internal_id, object_metadata$entity_key)
input[, internal_id := unname(entity_to_internal[entity_key])]
stopifnot(!anyNA(input$internal_id), !anyDuplicated(input$internal_id))

vrMainFeatureType(object) <- "Niche"
niche_raw <- as.matrix(vrData(object, norm = FALSE))
graph <- vrGraph(object, graph.type = "radius50")
degree_by_id <- degree(graph, v = colnames(niche_raw), mode = "all", loops = FALSE)
evaluable_ids <- names(degree_by_id)[degree_by_id > 0]
stopifnot(
  length(evaluable_ids) == 58772L,
  setequal(evaluable_ids, assignment$internal_id),
  max(abs(colSums(niche_raw) - degree_by_id)) < 1e-6
)

neighbor_average <- sweep(
  niche_raw[, evaluable_ids, drop = FALSE],
  2,
  degree_by_id[evaluable_ids],
  "/"
)
stopifnot(
  identical(rownames(neighbor_average), states),
  max(abs(colSums(neighbor_average) - 1)) < 1e-6
)

assignment <- assignment[match(evaluable_ids, internal_id)]
stopifnot(all(assignment$internal_id == evaluable_ids))
assignment[, `:=`(
  niche = factor(niche, levels = niche_levels),
  log1p_umi = log1p(nUMI),
  umi_quintile = pmin(
    5L,
    floor((frank(nUMI, ties.method = "average") - 1) * 5 / .N) + 1L
  )
), by = patient_id]

write_status("build_bin_level_neighborhood_weights", "58772 Bins x 18 states")
bin_state <- as.data.table(t(neighbor_average), keep.rownames = "internal_id")
bin_state <- melt(
  bin_state,
  id.vars = "internal_id",
  variable.name = "reference18_state",
  value.name = "neighbor_mean_weight"
)
bin_state <- merge(
  bin_state,
  assignment[, .(
    internal_id,
    entity_key,
    patient_id,
    disease,
    roi_id,
    niche,
    nUMI,
    log1p_umi,
    umi_quintile,
    neighbor_count
  )],
  by = "internal_id",
  all.x = TRUE,
  sort = FALSE
)
stopifnot(!anyNA(bin_state$patient_id), nrow(bin_state) == 58772L * 18L)

expected_comparison_n <- 7L * primary_k * 18L
write_status("summarize_patient_niche_composition", paste0("expected ", expected_comparison_n, " rows"))
composition <- bin_state[, .(
  bin_n = .N,
  mean_neighbor_weight = mean(neighbor_mean_weight),
  sd_neighbor_weight = sd(neighbor_mean_weight),
  median_neighbor_weight = median(neighbor_mean_weight),
  q25_neighbor_weight = as.numeric(quantile(neighbor_mean_weight, 0.25)),
  q75_neighbor_weight = as.numeric(quantile(neighbor_mean_weight, 0.75)),
  median_bin_umi = as.numeric(median(nUMI)),
  mean_bin_umi = mean(nUMI),
  median_neighbor_count = as.numeric(median(neighbor_count))
), by = .(patient_id, disease, niche, reference18_state)]
stopifnot(nrow(composition) == expected_comparison_n)

patient_state_totals <- bin_state[, .(
  patient_state_weight_sum = sum(neighbor_mean_weight),
  patient_bin_n = .N
), by = .(patient_id, disease, reference18_state)]
niche_state_totals <- bin_state[, .(
  niche_state_weight_sum = sum(neighbor_mean_weight),
  niche_bin_n = .N,
  niche_mean_log1p_umi = mean(log1p_umi),
  niche_median_umi = as.numeric(median(nUMI)),
  state_umi_spearman_rho_within_niche = safe_cor(neighbor_mean_weight, log1p_umi)
), by = .(patient_id, disease, niche, reference18_state)]

enrichment <- merge(
  niche_state_totals,
  patient_state_totals,
  by = c("patient_id", "disease", "reference18_state")
)
enrichment[, `:=`(
  other_bin_n = patient_bin_n - niche_bin_n,
  niche_mean_weight = niche_state_weight_sum / niche_bin_n,
  other_niches_mean_weight = (
    patient_state_weight_sum - niche_state_weight_sum
  ) / (patient_bin_n - niche_bin_n)
)]
enrichment[, enrichment := niche_mean_weight - other_niches_mean_weight]

patient_umi_totals <- assignment[, .(
  patient_log_umi_sum = sum(log1p_umi),
  patient_bin_n = .N
), by = patient_id]
niche_umi <- assignment[, .(
  niche_log_umi_sum = sum(log1p_umi),
  niche_bin_n = .N
), by = .(patient_id, niche = as.character(niche))]
niche_umi <- merge(niche_umi, patient_umi_totals, by = "patient_id")
niche_umi[, other_niches_mean_log1p_umi := (
  patient_log_umi_sum - niche_log_umi_sum
) / (patient_bin_n - niche_bin_n)]
niche_umi[, niche_mean_log1p_umi := niche_log_umi_sum / niche_bin_n]
niche_umi[, niche_vs_other_log1p_umi_difference := (
  niche_mean_log1p_umi - other_niches_mean_log1p_umi
)]

enrichment[, niche := as.character(niche)]
enrichment <- merge(
  enrichment,
  niche_umi[, .(
    patient_id,
    niche,
    other_niches_mean_log1p_umi,
    niche_vs_other_log1p_umi_difference
  )],
  by = c("patient_id", "niche"),
  all.x = TRUE
)

write_status("compute_umi_sensitivity", paste0(expected_comparison_n, " one-versus-rest comparisons"))
model_rows <- vector("list", expected_comparison_n)
row_index <- 1L
patient_order <- unique(assignment$patient_id)
for (patient in patient_order) {
  patient_assignment <- assignment[patient_id == patient]
  patient_ids <- patient_assignment$internal_id
  patient_weights <- neighbor_average[, patient_ids, drop = FALSE]
  patient_niche <- as.character(patient_assignment$niche)
  patient_log_umi <- patient_assignment$log1p_umi
  patient_quintile <- patient_assignment$umi_quintile

  for (state in states) {
    y <- as.numeric(patient_weights[state, ])
    for (current_niche in niche_levels) {
      target <- as.integer(patient_niche == current_niche)
      design <- cbind(1, target, as.numeric(scale(patient_log_umi)))
      fit <- lm.fit(design, y)
      adjusted_enrichment <- unname(fit$coefficients[[2]])

      strata_rows <- vector("list", 5L)
      for (quintile in 1:5) {
        in_stratum <- patient_quintile == quintile
        target_keep <- in_stratum & target == 1L
        other_keep <- in_stratum & target == 0L
        target_n <- sum(target_keep)
        other_n <- sum(other_keep)
        strata_rows[[quintile]] <- data.table(
          umi_quintile = quintile,
          target_n = target_n,
          other_n = other_n,
          valid = target_n >= 10L && other_n >= 10L,
          difference = if (target_n >= 10L && other_n >= 10L) {
            mean(y[target_keep]) - mean(y[other_keep])
          } else {
            NA_real_
          }
        )
      }
      strata <- rbindlist(strata_rows)
      valid_strata <- strata[valid == TRUE]
      stratified_enrichment <- if (nrow(valid_strata) > 0L) {
        mean(valid_strata$difference)
      } else {
        NA_real_
      }

      model_rows[[row_index]] <- data.table(
        patient_id = patient,
        niche = current_niche,
        reference18_state = state,
        umi_adjusted_enrichment = adjusted_enrichment,
        umi_coefficient = unname(fit$coefficients[[3]]),
        umi_stratified_enrichment = stratified_enrichment,
        valid_umi_quintile_n = nrow(valid_strata),
        minimum_target_bin_n_in_valid_quintile = if (nrow(valid_strata)) {
          min(valid_strata$target_n)
        } else {
          NA_integer_
        },
        minimum_other_bin_n_in_valid_quintile = if (nrow(valid_strata)) {
          min(valid_strata$other_n)
        } else {
          NA_integer_
        }
      )
      row_index <- row_index + 1L
    }
  }
}
umi_sensitivity <- rbindlist(model_rows)
stopifnot(nrow(umi_sensitivity) == expected_comparison_n)

enrichment <- merge(
  enrichment,
  umi_sensitivity,
  by = c("patient_id", "niche", "reference18_state"),
  all.x = TRUE
)
enrichment[, `:=`(
  raw_vs_adjusted_same_direction = mapply(same_sign, enrichment, umi_adjusted_enrichment),
  raw_vs_stratified_same_direction = mapply(same_sign, enrichment, umi_stratified_enrichment),
  adjusted_minus_raw = umi_adjusted_enrichment - enrichment,
  stratified_minus_raw = umi_stratified_enrichment - enrichment
)]

write_status("cross_patient_summary")
cross_patient <- enrichment[, {
  values <- enrichment
  adjusted_values <- umi_adjusted_enrichment
  stratified_values <- umi_stratified_enrichment
  full_median <- median(values)
  loo_medians <- vapply(seq_along(values), function(i) median(values[-i]), numeric(1))
  abs_total <- sum(abs(values))
  driver_index <- which.max(abs(values))
  list(
    patient_n = .N,
    median_enrichment = full_median,
    mean_enrichment = mean(values),
    sd_enrichment = sd(values),
    positive_patient_n = sum(values > 0),
    zero_patient_n = sum(values == 0),
    negative_patient_n = sum(values < 0),
    positive_direction_consistency = mean(values > 0),
    median_direction_consistency = if (full_median > 0) {
      mean(values > 0)
    } else if (full_median < 0) {
      mean(values < 0)
    } else {
      NA_real_
    },
    max_absolute_driver_patient = patient_id[[driver_index]],
    max_absolute_enrichment = abs(values[[driver_index]]),
    max_absolute_share = if (abs_total > 0) abs(values[[driver_index]]) / abs_total else NA_real_,
    single_patient_driven = if (abs_total > 0) {
      abs(values[[driver_index]]) / abs_total >= 0.5
    } else {
      NA
    },
    loo_median_min = min(loo_medians),
    loo_median_max = max(loo_medians),
    loo_direction_stable = if (abs(full_median) > 1e-12) {
      all(sign(loo_medians[abs(loo_medians) > 1e-12]) == sign(full_median))
    } else {
      NA
    },
    median_umi_adjusted_enrichment = median(adjusted_values, na.rm = TRUE),
    median_umi_stratified_enrichment = median(stratified_values, na.rm = TRUE),
    adjusted_direction_consistency = if (median(adjusted_values, na.rm = TRUE) > 0) {
      mean(adjusted_values > 0, na.rm = TRUE)
    } else {
      mean(adjusted_values < 0, na.rm = TRUE)
    },
    stratified_direction_consistency = if (median(stratified_values, na.rm = TRUE) > 0) {
      mean(stratified_values > 0, na.rm = TRUE)
    } else {
      mean(stratified_values < 0, na.rm = TRUE)
    }
  )
}, by = .(niche, reference18_state)]
cross_patient[, direction_class := fifelse(
  median_enrichment > 0 & positive_patient_n >= 6L,
  "Consistent positive (>=6/7)",
  fifelse(
    median_enrichment < 0 & negative_patient_n >= 6L,
    "Consistent negative (>=6/7)",
    "Patient-variable"
  )
)]

write_status("write_tables")
composition[, state_order := match(reference18_state, states)]
setorder(composition, patient_id, niche, state_order)
composition[, state_order := NULL]
enrichment[, state_order := match(reference18_state, states)]
setorder(enrichment, niche, state_order, patient_id)
enrichment[, state_order := NULL]
cross_patient[, state_order := match(reference18_state, states)]
setorder(cross_patient, niche, state_order)
cross_patient[, state_order := NULL]

overall_composition <- bin_state[, .(
  bin_n = .N,
  mean_neighbor_weight = mean(neighbor_mean_weight),
  sd_neighbor_weight = sd(neighbor_mean_weight),
  median_neighbor_weight = median(neighbor_mean_weight)
), by = .(niche, reference18_state)]
overall_composition[, state_order := match(reference18_state, states)]
setorder(overall_composition, niche, state_order)
overall_composition[, state_order := NULL]

umi_summary <- enrichment[, .(
  umi_evaluable_patient_n = .N,
  direction_retained_fraction_adjusted = mean(raw_vs_adjusted_same_direction, na.rm = TRUE),
  direction_retained_fraction_stratified = mean(raw_vs_stratified_same_direction, na.rm = TRUE),
  sign_reversal_fraction_adjusted = mean(raw_vs_adjusted_same_direction == FALSE, na.rm = TRUE),
  sign_reversal_fraction_stratified = mean(raw_vs_stratified_same_direction == FALSE, na.rm = TRUE),
  median_attenuation_adjusted = median(
    1 - abs(umi_adjusted_enrichment[abs(enrichment) > 1e-12]) /
      abs(enrichment[abs(enrichment) > 1e-12]),
    na.rm = TRUE
  ),
  median_attenuation_stratified = median(
    1 - abs(umi_stratified_enrichment[abs(enrichment) > 1e-12]) /
      abs(enrichment[abs(enrichment) > 1e-12]),
    na.rm = TRUE
  )
), by = .(niche, reference18_state)]

cross_patient <- merge(
  cross_patient,
  umi_summary,
  by = c("niche", "reference18_state"),
  all.x = TRUE,
  sort = FALSE
)
cross_patient[, `:=`(
  evaluable_patient_count = patient_n,
  positive_fraction = positive_patient_n / patient_n,
  negative_fraction = negative_patient_n / patient_n,
  maximum_single_patient_contribution = max_absolute_share,
  single_patient_driven_flag = single_patient_driven,
  core_state = median_enrichment > 0 & positive_patient_n >= 6L &
    loo_direction_stable == TRUE,
  umi_robust = median_enrichment * median_umi_adjusted_enrichment > 0 &
    median_enrichment * median_umi_stratified_enrichment > 0 &
    direction_retained_fraction_adjusted >= 5 / 7 &
    direction_retained_fraction_stratified >= 5 / 7
)]

program_map <- list(
  B = "B",
  Plasma = "Plasma",
  `T/NK` = c("T", "NK"),
  Myeloid = c(
    "Alveolar/resident macrophage",
    "Inflammatory monocyte/macrophage",
    "Dendritic cell"
  ),
  Endothelial = c("Blood endothelial", "Lymphatic endothelial"),
  Alveolar = c("AT1", "AT2/alveolar transitional"),
  Airway = c(
    "Ciliated epithelial",
    "Secretory/other airway epithelial",
    "Basal-like epithelial"
  ),
  Stromal = c(
    "Perivascular/mesothelial stromal",
    "Fibroblast",
    "Myofibroblast"
  )
)

annotation_rows <- lapply(niche_levels, function(current_niche) {
  current <- cross_patient[niche == current_niche]
  positive <- current[median_enrichment > 0][order(-median_enrichment)]
  core <- positive[core_state == TRUE]
  program_scores <- vapply(
    program_map,
    function(program_states) sum(
      pmax(current[reference18_state %chin% program_states, median_enrichment], 0),
      na.rm = TRUE
    ),
    numeric(1)
  )
  dominant_program <- names(which.max(program_scores))
  proposed_annotation <- switch(
    dominant_program,
    B = "B-enriched mixed niche",
    Plasma = "plasma-associated mixed niche",
    `T/NK` = "T/NK-associated mixed niche",
    Myeloid = "myeloid-associated mixed niche",
    Endothelial = "endothelial-associated mixed niche",
    Alveolar = "alveolar-associated mixed niche",
    Airway = "airway-associated mixed niche",
    Stromal = "stromal-associated mixed niche",
    "compositionally mixed niche"
  )
  data.table(
    niche = current_niche,
    proposed_annotation = proposed_annotation,
    dominant_positive_states = paste(head(positive$reference18_state, 5), collapse = "; "),
    core_positive_states = paste(core$reference18_state, collapse = "; "),
    dominant_negative_states = paste(
      head(current[median_enrichment < 0][order(median_enrichment)]$reference18_state, 5),
      collapse = "; "
    ),
    positive_patient_consistency = if (nrow(core)) {
      paste0(min(core$positive_patient_n), "/7 or greater")
    } else {
      "No positive state reaches 6/7"
    },
    umi_robust = nrow(core) > 0L && all(core$umi_robust),
    any_state_single_patient_driven = any(current$single_patient_driven == TRUE, na.rm = TRUE),
    core_state_single_patient_driven = any(core$single_patient_driven == TRUE, na.rm = TRUE),
    single_patient_driven_states = paste(
      current[single_patient_driven == TRUE, reference18_state],
      collapse = "; "
    ),
    confidence = if (nrow(core) > 0L && all(core$umi_robust) &&
      !any(core$single_patient_driven == TRUE, na.rm = TRUE)) {
      "moderate-high"
    } else if (nrow(core) > 0L) {
      "moderate"
    } else {
      "low-exploratory"
    },
    interpretation_limit = paste(
      "Relative Reference18 neighborhood composition; mixed Bin50 neighborhoods;",
      "not a pure cell type, independent validation, or disease result"
    )
  )
})
annotations <- rbindlist(annotation_rows)

fwrite(
  overall_composition,
  file.path(output_dir, "06_J2_PRIMARYK_REFERENCE18_COMPOSITION.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  composition,
  file.path(output_dir, "07_J2_PRIMARYK_PATIENT_NICHE_COMPOSITION.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  enrichment,
  file.path(output_dir, "08_J2_PRIMARYK_PATIENT_NICHE_ENRICHMENT.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  cross_patient,
  file.path(output_dir, "09_J2_PRIMARYK_CROSS_PATIENT_ENRICHMENT.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  enrichment,
  file.path(output_dir, "10_J2_PRIMARYK_UMI_SENSITIVITY.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  annotations,
  file.path(output_dir, "11_J2_PRIMARYK_CONSERVATIVE_ANNOTATIONS.tsv"),
  sep = "\t",
  na = "NA"
)

fwrite(
  composition,
  file.path(output_dir, "01_PATIENT_NICHE_REFERENCE18_COMPOSITION.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  enrichment,
  file.path(output_dir, "02_PATIENT_NICHE_STATE_ENRICHMENT.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  cross_patient,
  file.path(output_dir, "03_CROSS_PATIENT_ENRICHMENT_SUMMARY.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  enrichment[, .(
    patient_id,
    disease,
    niche,
    reference18_state,
    enrichment,
    state_umi_spearman_rho_within_niche,
    niche_vs_other_log1p_umi_difference,
    umi_adjusted_enrichment,
    umi_coefficient,
    umi_stratified_enrichment,
    valid_umi_quintile_n,
    minimum_target_bin_n_in_valid_quintile,
    minimum_other_bin_n_in_valid_quintile,
    raw_vs_adjusted_same_direction,
    raw_vs_stratified_same_direction,
    adjusted_minus_raw,
    stratified_minus_raw
  )],
  file.path(output_dir, "04_UMI_SENSITIVITY.tsv"),
  sep = "\t",
  na = "NA"
)

input_audit <- data.table(
  item = c(
    "chip", "K", "eligible_input_bin_n", "evaluable_assignment_bin_n",
    "patient_n", "Reference18_state_n", "patient_niche_state_row_n",
    "enrichment_row_n", "radius_um", "center_bin_included",
    "composition_quantity", "enrichment_baseline", "aggregate_labels_read",
    "disease_testing"
  ),
  value = c(
    "J2", primary_k, nrow(input), nrow(assignment), uniqueN(assignment$patient_id),
    length(states), nrow(composition), nrow(enrichment), 50, FALSE,
    "mean raw Reference18 weight across actual radius50 neighboring Bins",
    "all other evaluable Bin neighborhoods in the same patient, Bin-weighted",
    FALSE, "none"
  )
)
fwrite(input_audit, file.path(output_dir, "00_INPUT_AND_DEFINITION_AUDIT.tsv"), sep = "\t")

write_status("draw_figures")
patient_levels <- unique(assignment$patient_id)
composition_plot_data <- copy(composition)
composition_plot_data[, `:=`(
  patient_id = factor(patient_id, levels = patient_levels),
  niche = factor(niche, levels = niche_levels),
  reference18_state = factor(reference18_state, levels = rev(states))
)]
composition_plot <- ggplot(
  composition_plot_data,
  aes(niche, reference18_state, fill = mean_neighbor_weight)
) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_wrap(~patient_id, ncol = 4) +
  scale_fill_viridis_c(option = "C", limits = c(0, max(composition$mean_neighbor_weight))) +
  labs(
    title = paste0("Patient-specific Reference18 neighborhood composition: J2 K=", primary_k),
    subtitle = "Mean raw weight across actual 50 µm neighboring Bins; center Bin excluded",
    x = "Niche",
    y = "Reference18 state",
    fill = "Mean neighbor weight"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
save_plot(composition_plot, "02_J2_PRIMARYK_REFERENCE18_COMPOSITION_EN", 15, 11)

enrichment_plot_data <- copy(enrichment)
enrichment_plot_data[, `:=`(
  niche = factor(niche, levels = niche_levels),
  reference18_state = factor(reference18_state, levels = rev(states)),
  patient_id = factor(patient_id, levels = patient_levels)
)]
median_plot_data <- cross_patient[, .(
  niche = factor(niche, levels = niche_levels),
  reference18_state = factor(reference18_state, levels = rev(states)),
  median_enrichment
)]
patient_colors <- setNames(hcl.colors(length(patient_levels), "Dark 3"), patient_levels)
enrichment_plot <- ggplot(
  enrichment_plot_data,
  aes(enrichment, reference18_state, color = patient_id)
) +
  geom_vline(xintercept = 0, color = "grey45", linewidth = 0.35) +
  geom_point(size = 1.7, alpha = 0.85, position = position_jitter(height = 0.10, width = 0)) +
  geom_point(
    data = median_plot_data,
    aes(median_enrichment, reference18_state),
    inherit.aes = FALSE,
    shape = 18,
    size = 3,
    color = "black"
  ) +
  facet_wrap(~niche, ncol = 2, scales = "free_x") +
  scale_color_manual(values = patient_colors) +
  labs(
    title = paste0("Within-patient niche enrichment: J2 K=", primary_k),
    subtitle = "Niche mean minus the Bin-weighted mean of all other niches in the same patient; black diamond = median",
    x = "Reference18 weight enrichment",
    y = "Reference18 state",
    color = "Patient"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
save_plot(enrichment_plot, "03_J2_PRIMARYK_PATIENT_WISE_ENRICHMENT_EN", 13, 12)

consistency_plot_data <- copy(cross_patient)
consistency_plot_data[, `:=`(
  niche = factor(niche, levels = niche_levels),
  reference18_state = factor(reference18_state, levels = rev(states)),
  label = paste0(positive_patient_n, "/", patient_n)
)]
consistency_limit <- max(abs(consistency_plot_data$median_enrichment))
consistency_plot <- ggplot(
  consistency_plot_data,
  aes(niche, reference18_state)
) +
  geom_point(
    aes(size = positive_direction_consistency, fill = median_enrichment),
    shape = 21,
    color = "black",
    stroke = 0.25
  ) +
  geom_text(aes(label = label), size = 2.35) +
  scale_size_continuous(range = c(4.5, 9), limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-consistency_limit, consistency_limit)
  ) +
  labs(
    title = paste0("Cross-patient direction consistency: J2 K=", primary_k),
    subtitle = "Label = patients with positive enrichment / evaluable patients",
    x = "Niche",
    y = "Reference18 state",
    size = "Positive fraction",
    fill = "Median enrichment"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )
save_plot(consistency_plot, "04_J2_PRIMARYK_CROSS_PATIENT_DIRECTION_CONSISTENCY_EN", 10, 9)

umi_plot <- ggplot(
  enrichment_plot_data,
  aes(enrichment, umi_adjusted_enrichment, color = patient_id)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(alpha = 0.75, size = 1.5) +
  facet_wrap(~niche, ncol = 2, scales = "free") +
  scale_color_manual(values = patient_colors) +
  labs(
    title = "UMI sensitivity of patient-specific niche enrichment",
    subtitle = "Raw one-versus-rest enrichment versus the coefficient adjusted for log1p(Bin UMI)",
    x = "Raw enrichment",
    y = "UMI-adjusted enrichment",
    color = "Patient"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
save_plot(umi_plot, "05_J2_PRIMARYK_UMI_SENSITIVITY_EN", 11, 9)

write_status("validation")
validation <- data.table(
  check = c(
    "eligible input Bin count", "primary-K assignment Bin count", "patient count",
    "Reference18 state count", "patient-niche-state row count",
    "enrichment row count", "neighbor weights sum to one",
    "all patients have all four niches", "aggregate labels not read"
  ),
  observed = c(
    nrow(input), nrow(assignment), uniqueN(assignment$patient_id), length(states),
    nrow(composition), nrow(enrichment),
    max(abs(colSums(neighbor_average) - 1)),
    all(assignment[, uniqueN(niche), by = patient_id]$V1 == primary_k),
    TRUE
  ),
  expected = c(58840, 58772, 7, 18, expected_comparison_n, expected_comparison_n, 0, TRUE, TRUE)
)
validation[, passed := ifelse(
  check == "neighbor weights sum to one",
  as.numeric(observed) < 1e-6,
  as.character(observed) == as.character(expected)
)]
fwrite(validation, file.path(output_dir, "05_OUTPUT_VALIDATION.tsv"), sep = "\t")
if (!all(validation$passed)) stop("Output validation failed")

config <- data.table(
  parameter = c(
    "chip", "K", "radius_um", "center_bin_included", "patient_n",
    "state_n", "enrichment_baseline", "UMI_model",
    "UMI_stratification", "single_patient_driver_rule",
    "aggregate_handling", "disease_testing"
  ),
  value = c(
    "J2", primary_k, 50, FALSE, 7, 18,
    "target niche mean minus Bin-weighted mean across all other niches in the same patient",
    "per patient/niche/state linear model: neighbor weight ~ niche indicator + scaled log1p(Bin UMI)",
    "within-patient UMI quintiles; require >=10 target and >=10 other Bins per valid stratum; equal mean across valid strata",
    "maximum absolute patient contribution >=50% of total absolute enrichment",
    "not read or used", "none"
  )
)
fwrite(config, file.path(output_dir, "RUN_CONFIG.tsv"), sep = "\t")

manifest <- data.table(
  role = c(
    "frozen Reference18 UMI>=50 input",
    "saved VoltRon object with radius50 graph and Niche assay",
    "independently selected primary-K assignments",
    "analysis script",
    "aggregate input"
  ),
  project_relative_path = c(
    "results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/J2/01_FROZEN_REFERENCE18_INPUT.tsv.gz",
    file.path(resolution_dir, "14_VOLTRON_OBJECT_RADIUS50_CLR.rds"),
    assignment_path,
    "scripts/1014_j2_primaryk_composition_consistency.R",
    "NOT READ"
  ),
  used = c(TRUE, TRUE, TRUE, TRUE, FALSE)
)
fwrite(manifest, file.path(output_dir, "06_REPRODUCIBILITY_MANIFEST.tsv"), sep = "\t")

write_status("complete", paste0(expected_comparison_n, " patient-niche-state comparisons"))
writeLines("COMPLETE", file.path(output_dir, "COMPLETE.ok"))
if (requireNamespace("sessioninfo", quietly = TRUE)) {
  capture.output(
    sessioninfo::session_info(),
    file = file.path(output_dir, "SESSION_INFO.txt")
  )
} else {
  capture.output(sessionInfo(), file = file.path(output_dir, "SESSION_INFO.txt"))
}

cat("Completed J2 primary-K patient-specific Reference18 composition analysis.\n")
