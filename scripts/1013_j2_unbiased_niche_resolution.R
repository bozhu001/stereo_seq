#!/usr/bin/env Rscript

# J2-only unbiased 50 um neighborhood niche resolution analysis.
# This close variant of the audited L3 workflow uses the frozen Reference18
# UMI>=50 input. It reuses a saved official radius50 VoltRon object when
# available; otherwise it reconstructs that object from the same frozen input.
# Frozen lymphoid aggregate labels are intentionally not read.

suppressPackageStartupMessages({
  library(VoltRon)
  library(data.table)
  library(ggplot2)
  library(igraph)
  library(cluster)
})

options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 1013_j2_unbiased_niche_resolution.R <existing_j2_dir> <output_dir>")
}

existing_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
fig_dir <- file.path(output_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

status_path <- file.path(output_dir, "RUN_STATUS.tsv")
seed_base <- 20261008L
k_values <- 3:12
repeat_n <- 20L
iter_max <- 100L
silhouette_n <- 3000L
radius_um <- 50

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
  message("[J2 unbiased resolution] ", stage, if (nzchar(detail)) paste0(": ", detail))
}

read_gz <- function(path) {
  fread(cmd = paste("gzip -dc", shQuote(path)))
}

save_plot <- function(plot, stem, width, height, directory = fig_dir, dpi = 320) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(directory, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf,
    limitsize = FALSE
  )
  ggsave(
    file.path(directory, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = FALSE
  )
}

adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  choose2 <- function(x) x * (x - 1) / 2
  n <- sum(tab)
  if (n < 2L) return(NA_real_)
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  total <- choose2(n)
  expected <- sum_rows * sum_cols / total
  maximum <- (sum_rows + sum_cols) / 2
  if (maximum == expected) return(1)
  (sum_cells - expected) / (maximum - expected)
}

normalized_mutual_information <- function(a, b) {
  tab <- table(a, b)
  pxy <- tab / sum(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)
  expected <- outer(px, py)
  keep <- pxy > 0 & expected > 0
  mutual_information <- sum(pxy[keep] * log(pxy[keep] / expected[keep]))
  entropy_x <- -sum(px[px > 0] * log(px[px > 0]))
  entropy_y <- -sum(py[py > 0] * log(py[py > 0]))
  if (entropy_x == 0 || entropy_y == 0) return(1)
  mutual_information / sqrt(entropy_x * entropy_y)
}

eta_squared <- function(values, groups) {
  keep <- is.finite(values) & !is.na(groups)
  values <- values[keep]
  groups <- as.factor(groups[keep])
  if (length(values) < 2L || nlevels(groups) < 2L || var(values) == 0) return(NA_real_)
  total <- sum((values - mean(values))^2)
  between <- sum(vapply(
    split(values, groups),
    function(x) length(x) * (mean(x) - mean(values))^2,
    numeric(1)
  ))
  between / total
}

cramers_v <- function(x, y) {
  tab <- table(x, y)
  if (nrow(tab) < 2L || ncol(tab) < 2L) return(NA_real_)
  statistic <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic[[1]])
  sqrt(statistic / (sum(tab) * min(nrow(tab) - 1, ncol(tab) - 1)))
}

rank_score <- function(values, higher_is_better = TRUE) {
  output <- rep(NA_real_, length(values))
  keep <- is.finite(values)
  if (!any(keep)) return(output)
  if (length(unique(values[keep])) == 1L) {
    output[keep] <- 0.5
    return(output)
  }
  ranked <- (rank(values[keep], ties.method = "average") - 1) / (sum(keep) - 1)
  output[keep] <- if (higher_is_better) ranked else 1 - ranked
  output
}

capture_kmeans <- function(data, k, seed) {
  warnings <- character()
  error_message <- NA_character_
  fit <- tryCatch(
    withCallingHandlers(
      {
        set.seed(seed)
        stats::kmeans(
          data,
          centers = k,
          iter.max = iter_max,
          nstart = 1,
          algorithm = "Hartigan-Wong"
        )
      },
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      error_message <<- conditionMessage(condition)
      NULL
    }
  )
  list(
    fit = fit,
    warning = paste(unique(warnings), collapse = " | "),
    error = error_message
  )
}

program_map <- list(
  "Alveolar" = c("AT1", "AT2/alveolar transitional"),
  "Airway" = c(
    "Ciliated epithelial",
    "Secretory/other airway epithelial",
    "Basal-like epithelial"
  ),
  "Myeloid" = c(
    "Alveolar/resident macrophage",
    "Inflammatory monocyte/macrophage",
    "Dendritic cell"
  ),
  "Lymphoid" = c("B", "Plasma", "T", "NK", "Mast/Basophil"),
  "Fibrotic/stromal" = c(
    "Perivascular/mesothelial stromal",
    "Fibroblast",
    "Myofibroblast"
  ),
  "Vascular" = c("Blood endothelial", "Lymphatic endothelial")
)

theme_niche <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

mean_sd <- function(x) {
  data.frame(
    y = mean(x),
    ymin = mean(x) - sd(x),
    ymax = mean(x) + sd(x)
  )
}

write_status("read_fixed_inputs")
input <- read_gz(file.path(existing_dir, "01_FROZEN_REFERENCE18_INPUT.tsv.gz"))

weight_columns <- grep("^weight__", names(input), value = TRUE)
states <- sub("^weight__", "", weight_columns)
stopifnot(nrow(input) == 58840L)
stopifnot(length(states) == 18L)
stopifnot(uniqueN(input$patient_id) == 7L)
stopifnot(all(input$nUMI >= 50))
stopifnot(!anyDuplicated(input$entity_key))
stopifnot(max(abs(rowSums(as.matrix(input[, ..weight_columns])) - 1)) < 1e-6)

saved_object_path <- file.path(existing_dir, "14_VOLTRON_OBJECT_WITH_K3_K8.rds")
if (file.exists(saved_object_path)) {
  write_status("reuse_saved_radius50_object", saved_object_path)
  object <- readRDS(saved_object_path)
  object_source <- saved_object_path
} else {
  write_status("rebuild_official_radius50_object", "saved J2 VoltRon object unavailable")
  block_ids <- unique(input$analysis_block_id)
  objects <- vector("list", length(block_ids))
  for (block_index in seq_along(block_ids)) {
    block_id <- block_ids[[block_index]]
    current <- input[analysis_block_id == block_id]
    raw_matrix <- t(as.matrix(current[, ..weight_columns]))
    rownames(raw_matrix) <- states
    colnames(raw_matrix) <- current$entity_key
    coordinates <- data.frame(
      x = current$x_um,
      y = current$y_um,
      row.names = current$entity_key
    )
    metadata <- as.data.frame(current[, .(
      entity_key, bin_id, patient_id, sample_id, roi_id,
      analysis_block_id, disease, chip_id, chip_label, x_raw, y_raw,
      x_um, y_um, nUMI, nGene, tissue_component_id
    )])
    rownames(metadata) <- current$entity_key
    objects[[block_index]] <- formVoltRon(
      data = raw_matrix,
      metadata = metadata,
      coords = coordinates,
      assay.type = "spot",
      sample_name = sprintf("J2_block%02d", block_index),
      feature_name = "Reference18",
      project = "J2_Reference18_UMI50_50um_niche"
    )
  }
  object <- if (length(objects) == 1L) {
    objects[[1]]
  } else {
    merge(objects[[1]], objects[-1], verbose = FALSE)
  }
  rm(objects)
  invisible(gc())
  object <- getSpatialNeighbors(
    object,
    method = "radius",
    radius = radius_um,
    graph.key = "radius50",
    calculate.distances = TRUE,
    verbose = TRUE
  )
  object <- getNicheAssay(
    object,
    graph.type = "radius50",
    new_feature_name = "Niche"
  )
  vrMainFeatureType(object) <- "Niche"
  object <- normalizeData(object, method = "CLR", feat_type = "Niche")
  object_source <- file.path(output_dir, "14_VOLTRON_OBJECT_RADIUS50_CLR.rds")
  saveRDS(object, object_source, compress = "xz")
}

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
niche_clr <- as.matrix(vrData(object, norm = TRUE))
graph <- vrGraph(object, graph.type = "radius50")
degree_by_id <- degree(graph, v = colnames(niche_raw), mode = "all", loops = FALSE)
stopifnot(max(abs(colSums(niche_raw) - degree_by_id)) < 1e-6)
zero_neighbor_ids <- names(degree_by_id)[degree_by_id == 0]
evaluable_ids <- names(degree_by_id)[degree_by_id > 0]
stopifnot(length(zero_neighbor_ids) == 68L, length(evaluable_ids) == 58772L)
stopifnot(all(niche_raw[, zero_neighbor_ids, drop = FALSE] == 0))
stopifnot(all(niche_clr[, zero_neighbor_ids, drop = FALSE] == 0))

cluster_data <- t(niche_clr[, evaluable_ids, drop = FALSE])
neighbor_average <- sweep(
  niche_raw[, evaluable_ids, drop = FALSE],
  2,
  degree_by_id[evaluable_ids],
  "/"
)
metadata <- input[match(evaluable_ids, internal_id)]
stopifnot(all(metadata$internal_id == evaluable_ids))

graph_edges <- as_edgelist(graph, names = TRUE)
edge_keep <- graph_edges[, 1] %chin% evaluable_ids & graph_edges[, 2] %chin% evaluable_ids
graph_edges <- graph_edges[edge_keep, , drop = FALSE]

scope_audit <- data.table(
  item = c(
    "chip",
    "reference_states",
    "all_eligible_bin_n",
    "evaluable_nonempty_neighborhood_bin_n",
    "zero_neighbor_bin_n",
    "patient_n",
    "radius_um",
    "center_bin_included",
    "clustering_input",
    "aggregate_labels_read",
    "aggregate_labels_used_for_clustering",
    "aggregate_labels_used_for_k_selection",
    "aggregate_labels_used_for_naming"
  ),
  value = c(
    "J2",
    length(states),
    nrow(input),
    length(evaluable_ids),
    length(zero_neighbor_ids),
    uniqueN(input$patient_id),
    radius_um,
    FALSE,
    "Existing VoltRon native CLR of radius50 Reference18 Niche assay",
    FALSE,
    FALSE,
    FALSE,
    FALSE
  )
)
fwrite(scope_audit, file.path(output_dir, "00_SCOPE_AND_INPUT_AUDIT.tsv"), sep = "\t")

patient_input_audit <- input[, .(
  total_umi50_bin_n = .N,
  evaluable_radius50_bin_n = sum(internal_id %chin% evaluable_ids),
  roi_n = uniqueN(roi_id),
  missing_coordinate_n = sum(!is.finite(x_um) | !is.finite(y_um)),
  missing_reference18_weight_n = sum(!complete.cases(.SD))
), by = .(patient_id, disease, chip_id), .SDcols = weight_columns]
patient_input_audit[, `:=`(
  total_chip_bin_n = nrow(input),
  reference18_state_n = length(states),
  disease_used_in_clustering = FALSE,
  frozen_aggregate_read = FALSE
)]
fwrite(patient_input_audit, file.path(output_dir, "00_J2_INPUT_AUDIT.tsv"), sep = "\t")

neighborhood_audit <- data.table(
  chip = "Y40105J2",
  radius_um = radius_um,
  center_bin_included = FALSE,
  total_umi50_bin_n = nrow(input),
  evaluable_bin_n = length(evaluable_ids),
  zero_neighbor_bin_n = length(zero_neighbor_ids),
  edge_n = ecount(graph),
  patient_n = uniqueN(input$patient_id),
  roi_n = uniqueN(input$analysis_block_id),
  cross_patient_or_roi_edge_n = sum(
    object_metadata$analysis_block_id[match(graph_edges[, 1], internal_id)] !=
      object_metadata$analysis_block_id[match(graph_edges[, 2], internal_id)],
    na.rm = TRUE
  ),
  niche_aggregation = "sum of neighboring Bin50 Reference18 weights; center excluded",
  clr_method = "VoltRon normalizeData(method='CLR'); no pseudocount; structural zeros retained",
  object_source = object_source
)
fwrite(neighborhood_audit, file.path(output_dir, "01_J2_NEIGHBORHOOD_AUDIT.tsv"), sep = "\t")

zero_neighbor_table <- input[internal_id %chin% zero_neighbor_ids, .(
  internal_id,
  entity_key,
  patient_id,
  disease,
  roi_id,
  x_um,
  y_um,
  nUMI,
  reason = "radius50 graph degree = 0; excluded from clustering and marked not evaluable"
)]
fwrite(
  zero_neighbor_table,
  file.path(output_dir, "00_ZERO_NEIGHBOR_NOT_EVALUABLE.tsv"),
  sep = "\t"
)

write_status("fit_k3_to_k12", paste0(repeat_n, " seeds per K"))
all_assignments <- list()
all_composition <- list()
all_niche_summary <- list()
all_patient_proportions <- list()
all_disease_summary <- list()
all_umi_summary <- list()
all_seed_audit <- list()
all_pairwise_stability <- list()
all_interpretability <- list()
k_summary_rows <- list()

for (k in k_values) {
  write_status("fit_k", paste0("K=", k))
  k_dir <- file.path(output_dir, sprintf("K%02d", k))
  k_fig_dir <- file.path(k_dir, "figures")
  dir.create(k_fig_dir, recursive = TRUE, showWarnings = FALSE)
  seeds <- seed_base + k * 1000L + seq_len(repeat_n)
  fits <- vector("list", repeat_n)
  seed_rows <- vector("list", repeat_n)

  for (seed_index in seq_along(seeds)) {
    captured <- capture_kmeans(cluster_data, k, seeds[[seed_index]])
    fits[[seed_index]] <- captured$fit
    seed_rows[[seed_index]] <- data.table(
      k = k,
      seed_index = seed_index,
      seed = seeds[[seed_index]],
      fit_success = !is.null(captured$fit),
      converged = !is.null(captured$fit) && captured$fit$ifault == 0L,
      iter = if (is.null(captured$fit)) NA_integer_ else captured$fit$iter,
      ifault = if (is.null(captured$fit)) NA_integer_ else captured$fit$ifault,
      tot_withinss = if (is.null(captured$fit)) NA_real_ else captured$fit$tot.withinss,
      betweenss_fraction = if (is.null(captured$fit)) NA_real_ else captured$fit$betweenss / captured$fit$totss,
      warning = captured$warning,
      error = captured$error
    )
  }
  seed_audit <- rbindlist(seed_rows)
  all_seed_audit[[as.character(k)]] <- seed_audit
  converged_indices <- which(seed_audit$converged)
  if (length(converged_indices) == 0L) {
    stop("No converged fit for K=", k)
  }
  selected_index <- converged_indices[[which.min(seed_audit$tot_withinss[converged_indices])]]
  selected_fit <- fits[[selected_index]]
  labels <- selected_fit$cluster
  names(labels) <- rownames(cluster_data)
  niche_levels <- paste0("N", seq_len(k))
  niche <- factor(paste0("N", labels), levels = niche_levels)

  pair_rows <- list()
  pair_counter <- 1L
  if (length(converged_indices) >= 2L) {
    for (left_position in seq_len(length(converged_indices) - 1L)) {
      for (right_position in (left_position + 1L):length(converged_indices)) {
        left <- converged_indices[[left_position]]
        right <- converged_indices[[right_position]]
        pair_rows[[pair_counter]] <- data.table(
          k = k,
          seed_left = seeds[[left]],
          seed_right = seeds[[right]],
          ari = adjusted_rand_index(fits[[left]]$cluster, fits[[right]]$cluster),
          nmi = normalized_mutual_information(fits[[left]]$cluster, fits[[right]]$cluster)
        )
        pair_counter <- pair_counter + 1L
      }
    }
  }
  pairwise_stability <- rbindlist(pair_rows)
  all_pairwise_stability[[as.character(k)]] <- pairwise_stability

  assignment <- metadata[, .(
    internal_id,
    entity_key,
    bin_id,
    patient_id,
    disease,
    roi_id,
    analysis_block_id,
    x_um,
    y_um,
    nUMI,
    nGene
  )]
  assignment[, `:=`(
    k = k,
    niche_number = labels,
    niche = as.character(niche),
    selected_seed = seeds[[selected_index]],
    neighbor_count = as.integer(degree_by_id[internal_id])
  )]
  all_assignments[[as.character(k)]] <- assignment

  composition_rows <- list()
  program_rows <- list()
  for (niche_index in seq_len(k)) {
    keep <- labels == niche_index
    mean_weights <- rowMeans(neighbor_average[, keep, drop = FALSE])
    sd_weights <- apply(neighbor_average[, keep, drop = FALSE], 1, sd)
    composition_rows[[niche_index]] <- data.table(
      k = k,
      niche = paste0("N", niche_index),
      reference18_state = rownames(neighbor_average),
      mean_neighbor_weight = mean_weights,
      sd_neighbor_weight = sd_weights,
      bin_n = sum(keep)
    )
    program_values <- vapply(
      program_map,
      function(program_states) sum(mean_weights[program_states]),
      numeric(1)
    )
    ordered_programs <- sort(program_values, decreasing = TRUE)
    dominant <- names(ordered_programs)[[1]]
    margin <- ordered_programs[[1]] - ordered_programs[[2]]
    descriptor <- if (ordered_programs[[1]] >= 0.35 && margin >= 0.05) {
      dominant
    } else {
      "Mixed"
    }
    top_states <- sort(mean_weights, decreasing = TRUE)[seq_len(3)]
    program_rows[[niche_index]] <- data.table(
      k = k,
      niche = paste0("N", niche_index),
      descriptor = descriptor,
      dominant_program = dominant,
      dominant_program_weight = ordered_programs[[1]],
      second_program = names(ordered_programs)[[2]],
      second_program_weight = ordered_programs[[2]],
      dominance_margin = margin,
      top_3_reference18_states = paste(
        paste0(names(top_states), "=", sprintf("%.3f", top_states)),
        collapse = "; "
      ),
      naming_rule = paste0(
        "descriptor only if dominant broad-program weight >=0.35 and top-vs-second margin >=0.05; ",
        "otherwise Mixed; descriptors are composition summaries, not cell identities"
      )
    )
  }
  composition <- rbindlist(composition_rows)
  interpretability <- rbindlist(program_rows)
  all_composition[[as.character(k)]] <- composition
  all_interpretability[[as.character(k)]] <- interpretability

  patient_counts <- assignment[, .N, by = .(k, patient_id, disease, niche)]
  patient_totals <- assignment[, .(patient_evaluable_bin_n = .N), by = .(k, patient_id, disease)]
  patient_proportions <- merge(
    CJ(
      k = k,
      patient_id = unique(assignment$patient_id),
      niche = niche_levels,
      unique = TRUE
    ),
    unique(assignment[, .(patient_id, disease)]),
    by = "patient_id",
    allow.cartesian = TRUE
  )
  patient_proportions <- merge(
    patient_proportions,
    patient_counts,
    by = c("k", "patient_id", "disease", "niche"),
    all.x = TRUE
  )
  patient_proportions[is.na(N), N := 0L]
  patient_proportions <- merge(
    patient_proportions,
    patient_totals,
    by = c("k", "patient_id", "disease")
  )
  patient_proportions[, proportion := N / patient_evaluable_bin_n]
  all_patient_proportions[[as.character(k)]] <- patient_proportions

  niche_summary <- assignment[, .(
    bin_n = .N,
    full_chip_fraction = .N / nrow(assignment),
    patient_coverage_n = uniqueN(patient_id),
    disease_coverage_n = uniqueN(disease),
    median_umi = as.numeric(median(nUMI)),
    mean_umi = mean(nUMI),
    q25_umi = as.numeric(quantile(nUMI, 0.25)),
    q75_umi = as.numeric(quantile(nUMI, 0.75))
  ), by = .(k, niche)]
  niche_summary[, `:=`(
    below_one_percent = full_chip_fraction < 0.01,
    single_patient_only = patient_coverage_n == 1L
  )]
  all_niche_summary[[as.character(k)]] <- niche_summary

  umi_patient_summary <- assignment[, .(
    bin_n = .N,
    median_umi = as.numeric(median(nUMI)),
    mean_umi = mean(nUMI),
    q25_umi = as.numeric(quantile(nUMI, 0.25)),
    q75_umi = as.numeric(quantile(nUMI, 0.75))
  ), by = .(k, patient_id, disease, niche)]
  all_umi_summary[[as.character(k)]] <- umi_patient_summary

  disease_summary <- patient_proportions[, .(
    patient_n = .N,
    mean_patient_proportion = mean(proportion),
    sd_patient_proportion = sd(proportion),
    median_patient_proportion = median(proportion),
    min_patient_proportion = min(proportion),
    max_patient_proportion = max(proportion)
  ), by = .(k, disease, niche)]
  all_disease_summary[[as.character(k)]] <- disease_summary

  patient_centered_log_umi <- log1p(assignment$nUMI) -
    ave(log1p(assignment$nUMI), assignment$patient_id, FUN = mean)
  umi_eta_global <- eta_squared(log1p(assignment$nUMI), assignment$niche)
  umi_eta_within <- eta_squared(patient_centered_log_umi, assignment$niche)
  umi_association <- if (umi_eta_within >= 0.14) {
    "high"
  } else if (umi_eta_within >= 0.06) {
    "moderate"
  } else {
    "low"
  }

  edge_left <- labels[graph_edges[, 1]]
  edge_right <- labels[graph_edges[, 2]]
  spatial_edge_agreement <- mean(edge_left == edge_right)
  patient_cramers_v <- cramers_v(assignment$niche, assignment$patient_id)

  set.seed(seed_base + k * 100L)
  silhouette_index <- sample(seq_len(nrow(cluster_data)), min(silhouette_n, nrow(cluster_data)))
  silhouette_value <- mean(cluster::silhouette(
    labels[silhouette_index],
    dist(cluster_data[silhouette_index, , drop = FALSE])
  )[, "sil_width"])

  non_mixed_fraction <- mean(interpretability$descriptor != "Mixed")
  interpretability_index <- 0.5 * mean(interpretability$dominance_margin) +
    0.5 * non_mixed_fraction
  k_summary_rows[[as.character(k)]] <- data.table(
    k = k,
    evaluable_bin_n = nrow(assignment),
    zero_neighbor_not_evaluable_n = length(zero_neighbor_ids),
    seed_fit_n = nrow(seed_audit),
    converged_seed_n = sum(seed_audit$converged),
    converged_fraction = mean(seed_audit$converged),
    selected_seed = seeds[[selected_index]],
    selected_iter = selected_fit$iter,
    selected_ifault = selected_fit$ifault,
    betweenss_fraction = selected_fit$betweenss / selected_fit$totss,
    silhouette_mean_subsample = silhouette_value,
    repeat_ari_mean = mean(pairwise_stability$ari),
    repeat_ari_median = median(pairwise_stability$ari),
    repeat_ari_min = min(pairwise_stability$ari),
    repeat_nmi_mean = mean(pairwise_stability$nmi),
    repeat_nmi_median = median(pairwise_stability$nmi),
    repeat_nmi_min = min(pairwise_stability$nmi),
    min_niche_bin_n = min(niche_summary$bin_n),
    min_niche_fraction = min(niche_summary$full_chip_fraction),
    min_patient_coverage_n = min(niche_summary$patient_coverage_n),
    below_one_percent_niche_n = sum(niche_summary$below_one_percent),
    single_patient_niche_n = sum(niche_summary$single_patient_only),
    spatial_edge_agreement = spatial_edge_agreement,
    patient_cramers_v = patient_cramers_v,
    global_umi_eta_squared = umi_eta_global,
    within_patient_umi_eta_squared = umi_eta_within,
    umi_association_level = umi_association,
    non_mixed_descriptor_fraction = non_mixed_fraction,
    mean_dominance_margin = mean(interpretability$dominance_margin),
    interpretability_index = interpretability_index,
    unique_non_mixed_descriptor_n = uniqueN(
      interpretability[descriptor != "Mixed", descriptor]
    )
  )

  fwrite(composition, file.path(k_dir, "01_REFERENCE18_NEIGHBORHOOD_COMPOSITION.tsv"), sep = "\t")
  fwrite(assignment, file.path(k_dir, "02_BIN_NICHE_ASSIGNMENTS.tsv.gz"), sep = "\t")
  fwrite(niche_summary, file.path(k_dir, "03_NICHE_SIZE_AND_PATIENT_COVERAGE.tsv"), sep = "\t")
  fwrite(patient_proportions, file.path(k_dir, "04_PATIENT_NICHE_PROPORTIONS.tsv"), sep = "\t")
  fwrite(umi_patient_summary, file.path(k_dir, "05_PATIENT_NICHE_UMI_SUMMARY.tsv"), sep = "\t")
  fwrite(disease_summary, file.path(k_dir, "06_DISEASE_PATIENT_ABUNDANCE_SUMMARY.tsv"), sep = "\t")
  fwrite(seed_audit, file.path(k_dir, "07_SEED_FIT_AUDIT.tsv"), sep = "\t")
  fwrite(pairwise_stability, file.path(k_dir, "08_PAIRWISE_SEED_STABILITY.tsv"), sep = "\t")
  fwrite(interpretability, file.path(k_dir, "09_COMPOSITION_INTERPRETABILITY.tsv"), sep = "\t")

  niche_colors <- setNames(hcl.colors(k, palette = "Dark 3"), niche_levels)
  color_key <- data.table(niche = niche_levels, color = unname(niche_colors))
  fwrite(color_key, file.path(k_dir, "00_NICHE_COLOR_KEY.tsv"), sep = "\t")

  composition_plot_data <- copy(composition)
  composition_plot_data[, niche := factor(niche, levels = niche_levels)]
  assignment_plot_data <- copy(assignment)
  assignment_plot_data[, niche := factor(niche, levels = niche_levels)]
  patient_plot_data <- copy(patient_proportions)
  patient_plot_data[, niche := factor(niche, levels = niche_levels)]

  composition_plot <- ggplot(
    composition_plot_data,
    aes(niche, factor(reference18_state, levels = rev(states)), fill = mean_neighbor_weight)
  ) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.3f", mean_neighbor_weight)), size = if (k <= 6) 2.4 else 1.7) +
    scale_fill_viridis_c(option = "C", limits = c(0, max(composition$mean_neighbor_weight))) +
    labs(
      title = paste0("J2 unbiased radius50 niches: K=", k),
      subtitle = "Mean Reference18 weight across neighboring Bins; center Bin excluded",
      x = "Niche number (composition-based; not a cell type)",
      y = "Reference18 state",
      fill = "Mean neighbor weight"
    ) +
    theme_niche +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(composition_plot, "01_REFERENCE18_COMPOSITION_EN", 11, 8, k_fig_dir)

  spatial_pdf <- file.path(k_fig_dir, "02_ALL_PATIENT_SPATIAL_NICHE_MAPS_EN.pdf")
  cairo_pdf(spatial_pdf, width = 9, height = 8, onefile = TRUE)
  for (patient in unique(assignment_plot_data$patient_id)) {
    current <- assignment_plot_data[patient_id == patient]
    spatial_plot <- ggplot(current, aes(x_um, y_um, color = niche)) +
      geom_point(size = 0.32, alpha = 0.85) +
      scale_color_manual(values = niche_colors, drop = FALSE) +
      scale_y_reverse() +
      coord_fixed() +
      labs(
        title = paste0(patient, ": J2 unbiased radius50 niches (K=", k, ")"),
        subtitle = paste0("", nrow(current), " evaluable Bins; disease: ", unique(current$disease)),
        x = "x (um)",
        y = "y (um)",
        color = "Niche"
      ) +
      theme_niche +
      theme(panel.grid = element_blank())
    print(spatial_plot)
    ggsave(
      file.path(k_fig_dir, paste0("02_SPATIAL_", gsub("[^A-Za-z0-9]+", "_", patient), "_EN.png")),
      plot = spatial_plot,
      width = 9,
      height = 8,
      dpi = 240,
      limitsize = FALSE
    )
  }
  dev.off()

  patient_plot <- ggplot(
    patient_plot_data,
    aes(patient_id, proportion, fill = niche)
  ) +
    geom_col(width = 0.82) +
    scale_fill_manual(values = niche_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = paste0("Patient niche proportions: K=", k),
      subtitle = "All evaluable whole-tissue J2 Bins; patients are descriptive units",
      x = "Patient",
      y = "Niche proportion",
      fill = "Niche"
    ) +
    theme_niche +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_plot(patient_plot, "03_PATIENT_NICHE_PROPORTIONS_EN", 11, 7, k_fig_dir)

  umi_plot <- ggplot(assignment_plot_data, aes(niche, nUMI, fill = niche)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
    scale_fill_manual(values = niche_colors, guide = "none") +
    scale_y_log10() +
    facet_wrap(~patient_id, scales = "free_y") +
    labs(
      title = paste0("Within-patient UMI distributions by niche: K=", k),
      subtitle = paste0(
        "Within-patient eta-squared = ", sprintf("%.3f", umi_eta_within),
        " (", umi_association, " descriptive association)"
      ),
      x = "Niche number",
      y = "Bin UMI (log10 scale)"
    ) +
    theme_niche +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(umi_plot, "04_WITHIN_PATIENT_UMI_DISTRIBUTION_EN", 14, 9, k_fig_dir)

  disease_plot <- ggplot(
    patient_plot_data,
    aes(disease, proportion, color = disease)
  ) +
    geom_jitter(width = 0.10, height = 0, size = 2.2, alpha = 0.85) +
    stat_summary(fun = mean, geom = "point", color = "black", size = 2.4) +
    stat_summary(
      fun.data = mean_sd,
      geom = "errorbar",
      color = "black",
      width = 0.18,
      linewidth = 0.55
    ) +
    facet_wrap(~niche, scales = "free_y") +
    scale_color_manual(values = c("HC" = "#4DAF4A", "IPF" = "#377EB8", "SSc-ILD" = "#E41A1C")) +
    labs(
      title = paste0("Disease-group abundance overview: K=", k),
      subtitle = "One point per patient; black mean +/- SD; descriptive only, no disease testing",
      x = "Disease group",
      y = "Patient niche proportion",
      color = "Disease"
    ) +
    theme_niche +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(disease_plot, "05_DISEASE_PATIENT_ABUNDANCE_EN", 14, 8, k_fig_dir)
}

write_status("combine_and_score")
assignments_all <- rbindlist(all_assignments)
composition_all <- rbindlist(all_composition)
niche_summary_all <- rbindlist(all_niche_summary)
patient_proportions_all <- rbindlist(all_patient_proportions)
disease_summary_all <- rbindlist(all_disease_summary)
umi_summary_all <- rbindlist(all_umi_summary)
seed_audit_all <- rbindlist(all_seed_audit)
pairwise_stability_all <- rbindlist(all_pairwise_stability)
interpretability_all <- rbindlist(all_interpretability)
k_summary <- rbindlist(k_summary_rows)

k_summary[, `:=`(
  score_stability = rank_score((repeat_ari_mean + repeat_nmi_mean) / 2),
  score_cluster_size = rank_score(min_niche_fraction),
  score_patient_coverage = rank_score(min_patient_coverage_n),
  score_low_umi = rank_score(within_patient_umi_eta_squared, higher_is_better = FALSE),
  score_interpretability = rank_score(interpretability_index),
  score_silhouette = rank_score(silhouette_mean_subsample),
  score_spatial = rank_score(spatial_edge_agreement)
)]
k_summary[, recommendation_score :=
  0.25 * score_stability +
  0.15 * score_cluster_size +
  0.10 * score_patient_coverage +
  0.15 * score_low_umi +
  0.20 * score_interpretability +
  0.10 * score_silhouette +
  0.05 * score_spatial
]
k_summary[, eligibility :=
  converged_fraction >= 0.90 &
  below_one_percent_niche_n == 0L &
  single_patient_niche_n == 0L &
  min_patient_coverage_n >= 3L
]
if (any(k_summary$eligibility)) {
  recommended_k <- k_summary[eligibility == TRUE, k[which.max(recommendation_score)]]
} else {
  recommended_k <- k_summary$k[[which.max(k_summary$recommendation_score)]]
}
k_summary[, recommended := k == recommended_k]
k_summary[, recommendation_basis := paste0(
  "Score weights: stability 0.25, minimum niche size 0.15, patient coverage 0.10, ",
  "low within-patient UMI association 0.15, composition interpretability 0.20, ",
  "silhouette 0.10, spatial edge agreement 0.05. Eligibility requires >=90% converged ",
  "seeds, no niche <1%, no single-patient niche, and minimum patient coverage >=3."
)]

fwrite(k_summary, file.path(output_dir, "01_K3_K12_RECOMMENDATION.tsv"), sep = "\t")
fwrite(seed_audit_all, file.path(output_dir, "02_ALL_K_SEED_FIT_AUDIT.tsv.gz"), sep = "\t")
fwrite(pairwise_stability_all, file.path(output_dir, "03_ALL_K_PAIRWISE_STABILITY.tsv.gz"), sep = "\t")
fwrite(niche_summary_all, file.path(output_dir, "04_ALL_K_NICHE_SIZE_PATIENT_COVERAGE.tsv"), sep = "\t")
fwrite(patient_proportions_all, file.path(output_dir, "05_ALL_K_PATIENT_NICHE_PROPORTIONS.tsv.gz"), sep = "\t")
fwrite(disease_summary_all, file.path(output_dir, "06_ALL_K_DISEASE_PATIENT_ABUNDANCE.tsv"), sep = "\t")
fwrite(umi_summary_all, file.path(output_dir, "07_ALL_K_PATIENT_NICHE_UMI.tsv.gz"), sep = "\t")
fwrite(composition_all, file.path(output_dir, "08_ALL_K_REFERENCE18_COMPOSITION.tsv"), sep = "\t")
fwrite(interpretability_all, file.path(output_dir, "09_ALL_K_COMPOSITION_INTERPRETABILITY.tsv"), sep = "\t")
fwrite(assignments_all, file.path(output_dir, "10_ALL_K_BIN_ASSIGNMENTS.tsv.gz"), sep = "\t")
fwrite(seed_audit_all, file.path(output_dir, "02_J2_ALL_K_SEED_FIT_AUDIT.tsv.gz"), sep = "\t")
fwrite(k_summary, file.path(output_dir, "03_J2_K3_K12_RESOLUTION_SUMMARY.tsv"), sep = "\t")

metric_long <- melt(
  k_summary,
  id.vars = "k",
  measure.vars = c(
    "repeat_ari_mean",
    "repeat_nmi_mean",
    "silhouette_mean_subsample",
    "min_niche_fraction",
    "within_patient_umi_eta_squared",
    "interpretability_index",
    "recommendation_score"
  ),
  variable.name = "metric",
  value.name = "value"
)
recommendation_plot <- ggplot(metric_long, aes(k, value)) +
  geom_line(color = "#3B6FB6", linewidth = 0.7) +
  geom_point(aes(fill = k == recommended_k), shape = 21, size = 2.8, color = "black") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#E69F00"), guide = "none") +
  scale_x_continuous(breaks = k_values) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  labs(
    title = "J2 unbiased radius50 niche-resolution comparison",
    subtitle = paste0(
      "K=", recommended_k,
      " is recommended by a predeclared multi-criterion score; silhouette is only one component"
    ),
    x = "K",
    y = "Metric value"
  ) +
  theme_niche
save_plot(recommendation_plot, "01_K3_K12_RESOLUTION_SUMMARY_EN", 12, 11)

small_plot <- ggplot(niche_summary_all, aes(factor(k), full_chip_fraction, color = factor(k))) +
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "#B2182B") +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 2, alpha = 0.8) +
  scale_y_log10() +
  guides(color = "none") +
  labs(
    title = "Niche-size distribution across candidate K",
    subtitle = "Dashed line marks 1% of all evaluable J2 Bins",
    x = "K",
    y = "Whole-chip niche fraction (log10 scale)"
  ) +
  theme_niche
save_plot(small_plot, "02_K3_K12_NICHE_SIZE_EN", 10, 7)

scope_manifest <- data.table(
  role = c(
    "frozen Reference18 UMI>=50 input",
    "saved or reconstructed VoltRon object with radius50 Niche assay",
    "analysis script",
    "aggregate input"
  ),
  path = c(
    file.path(existing_dir, "01_FROZEN_REFERENCE18_INPUT.tsv.gz"),
    object_source,
    normalizePath("scripts/1013_j2_unbiased_niche_resolution.R", mustWork = TRUE),
    "NOT READ"
  )
)
fwrite(scope_manifest, file.path(output_dir, "11_INPUT_AND_SCRIPT_MANIFEST.tsv"), sep = "\t")

config <- data.table(
  parameter = c(
    "chip", "radius_um", "K_values", "seed_base", "seeds_per_K",
    "kmeans_algorithm", "iter_max", "silhouette_subsample_n",
    "zero_neighbor_handling", "aggregate_handling", "disease_testing"
  ),
  value = c(
    "J2", radius_um, "3:12", seed_base, repeat_n,
    "stats::kmeans Hartigan-Wong on existing VoltRon native CLR", iter_max,
    silhouette_n, "excluded from clustering; retained as not evaluable",
    "not read or used", "none"
  )
)
fwrite(config, file.path(output_dir, "RUN_CONFIG.tsv"), sep = "\t")

validation <- data.table(
  check = c(
    "eligible_bin_n",
    "evaluable_bin_n",
    "zero_neighbor_bin_n",
    "patient_n",
    "reference_state_n",
    "K_candidate_n",
    "seed_fit_row_n",
    "assignment_row_n",
    "unique_assignment_key",
    "no_aggregate_input",
    "all_K_have_recommendation_rows",
    "recommended_K_unique"
  ),
  observed = c(
    nrow(input),
    length(evaluable_ids),
    length(zero_neighbor_ids),
    uniqueN(input$patient_id),
    length(states),
    length(k_values),
    nrow(seed_audit_all),
    nrow(assignments_all),
    anyDuplicated(assignments_all[, .(k, internal_id)]),
    TRUE,
    nrow(k_summary),
    sum(k_summary$recommended)
  ),
  expected = c(
    58840, 58772, 68, 7, 18, 10, length(k_values) * repeat_n,
    length(k_values) * length(evaluable_ids), 0, TRUE, length(k_values), 1
  )
)
validation[, passed := as.character(observed) == as.character(expected)]
fwrite(validation, file.path(output_dir, "12_OUTPUT_VALIDATION.tsv"), sep = "\t")
if (!all(validation$passed)) stop("Output validation failed")

write_status("complete", paste0("recommended K=", recommended_k))
writeLines("COMPLETE", file.path(output_dir, "COMPLETE.ok"))
if (requireNamespace("sessioninfo", quietly = TRUE)) {
  print(sessioninfo::session_info())
} else {
  print(sessionInfo())
}
