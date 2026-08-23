#!/usr/bin/env Rscript

# Prespecified depth/n_genes/PTPRC-matched spatial null analysis for the six
# fixed L3 Mature_B_core samples. This script consumes the persisted frozen
# Bin-level support flags from analysis 395. It never reconstructs or changes
# marker expression, same-bin/rook scores, thresholds, or sample membership.
#
# Exact command:
# Rscript scripts/397_l3_mature_b_depth_ptprc_matched_spatial_null.R --timestamp=YYMMDDHHMMSS

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(RANN)
  library(sessioninfo)
  library(scales)
})

options(stringsAsFactors = FALSE)

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

source_dir <- file.path(repo, "results/reanalysis/bin50_l3_mature_b_focal_spatial_review_260823202319")
source_tables <- file.path(source_dir, "tables")
bin_path <- file.path(source_tables, "10_l3_bin_level_frozen_marker_and_support_values.tsv.gz")
sample_path <- file.path(source_tables, "01_fixed_l3_sample_manifest.tsv")
support_path <- file.path(source_tables, "05_mature_b_strict_support_summary.tsv")
prior_unmatched_path <- file.path(source_tables, "07_supported_bin_depth_overlap_vs_random.tsv")

out_dir <- file.path(repo, paste0("results/reanalysis/bin50_l3_mature_b_depth_ptprc_matched_spatial_null_", tstamp))
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
save_png <- function(plot, name, width, height, dpi = 180) {
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot, width = width, height = height,
         dpi = dpi, bg = "white", limitsize = FALSE)
}

sample_order <- c(
  "HC/NL-72", "IPF/FO23-1-06168", "IPF/FO23-1-06170",
  "SSC/05957/17B", "SSC/07998/15A", "SSC/15491/14"
)
focus_samples <- c("IPF/FO23-1-06168", "SSC/15491/14")
support_modes <- c("same_bin", "rook")
metric_order <- c(
  "max_rook_component_bin_n", "max_component_support_fraction",
  "connected_component_n", "median_nearest_neighbor_distance"
)
metric_labels <- c(
  max_rook_component_bin_n = "Largest rook component (Bins)",
  max_component_support_fraction = "Largest component / support",
  connected_component_n = "Connected components",
  median_nearest_neighbor_distance = "Median nearest-neighbor distance"
)
random_draw_n <- 1000L
random_seed <- 260824L
set.seed(random_seed)

required <- c(bin_path, sample_path, support_path, prior_unmatched_path)
assert(all(file.exists(required)), paste("Missing immutable input:", paste(required[!file.exists(required)], collapse = ", ")))
input_md5 <- unname(tools::md5sum(required))
log_message("Start; immutable Bin input=", bin_path)

bins <- fread(bin_path)
samples <- fread(sample_path)
frozen_support <- fread(support_path)
prior_unmatched <- fread(prior_unmatched_path)

required_columns <- c(
  "bin_id", "sample_id", "group", "chip_id", "x", "y", "component",
  "total_counts", "n_genes", "PTPRC", "strict_same_bin", "strict_rook"
)
assert(all(required_columns %chin% names(bins)), "Persisted Bin-level input columns are incomplete")
assert(identical(as.character(samples$sample_id), sample_order), "Fixed sample manifest order or membership drift")
assert(setequal(unique(bins$sample_id), sample_order), "Bin table sample membership drift")
assert(nrow(bins) == sum(samples$total_bin_n), "Bin table row count differs from fixed manifest")
assert(all(samples$chip_id == "Y40105L3") && all(bins$chip_id == "Y40105L3"), "Non-L3 sample detected")

support_check <- bins[, .(
  total_bin_n = .N,
  strict_same_bin_supported_bin_n = sum(strict_same_bin),
  strict_rook_supported_bin_n = sum(strict_rook)
), by = sample_id]
support_check <- merge(
  support_check,
  frozen_support[, .(sample_id, frozen_total_bin_n = total_bin_n,
                     frozen_same_n = strict_same_bin_supported_bin_n,
                     frozen_rook_n = strict_rook_supported_bin_n)],
  by = "sample_id", all.x = TRUE
)
support_check[, `:=`(
  total_n_exact = total_bin_n == frozen_total_bin_n,
  same_n_exact = strict_same_bin_supported_bin_n == frozen_same_n,
  rook_n_exact = strict_rook_supported_bin_n == frozen_rook_n
)]
assert(all(support_check$total_n_exact & support_check$same_n_exact & support_check$rook_n_exact),
       "Persisted support flags do not match frozen aggregate results")

# Deterministic empirical quintiles. Average ranks keep tied values together;
# empty nominal cells are allowed and reported rather than forcing five bins.
bins[, depth_quintile := pmin(5L, pmax(1L, ceiling(5 * frank(log1p(total_counts), ties.method = "average") / .N))),
     by = sample_id]
bins[, n_genes_quintile := pmin(5L, pmax(1L, ceiling(5 * frank(n_genes, ties.method = "average") / .N))),
     by = sample_id]
bins[, ptprc_status := fifelse(PTPRC > 0, "raw_positive", "raw_zero")]
bins[, matching_stratum := paste0("D", depth_quintile, "_G", n_genes_quintile, "_P", ptprc_status)]
bins[, row_in_sample := seq_len(.N), by = sample_id]

strata_summary <- bins[, .(
  stratum_bin_n = .N,
  median_total_counts = as.numeric(median(total_counts)),
  median_n_genes = as.numeric(median(n_genes)),
  ptprc_positive_fraction = mean(PTPRC > 0),
  observed_same_bin_n = sum(strict_same_bin),
  observed_rook_n = sum(strict_rook)
), by = .(sample_id, group, chip_id, depth_quintile, n_genes_quintile, ptprc_status, matching_stratum)]

build_rook_neighbors <- function(dt) {
  keys <- paste(dt$x, dt$y, sep = "_")
  assert(!anyDuplicated(keys), paste("Duplicate coordinates:", dt$sample_id[[1]]))
  neighbors <- vector("list", nrow(dt))
  for (offset in list(c(0L, 50L), c(50L, 0L))) {
    j <- match(paste(dt$x + offset[[1]], dt$y + offset[[2]], sep = "_"), keys)
    i <- which(!is.na(j))
    j <- j[i]
    keep <- dt$component[i] == dt$component[j]
    i <- i[keep]
    j <- j[keep]
    for (k in seq_along(i)) {
      neighbors[[i[[k]]]] <- c(neighbors[[i[[k]]]], j[[k]])
      neighbors[[j[[k]]]] <- c(neighbors[[j[[k]]]], i[[k]])
    }
  }
  neighbors
}

spatial_metrics <- function(selected, coords, neighbors) {
  selected <- sort(as.integer(selected))
  n <- length(selected)
  if (!n) {
    return(c(
      max_rook_component_bin_n = 0,
      max_component_support_fraction = NA_real_,
      connected_component_n = 0,
      median_nearest_neighbor_distance = NA_real_
    ))
  }
  selected_flag <- logical(nrow(coords))
  selected_flag[selected] <- TRUE
  visited <- logical(nrow(coords))
  component_sizes <- integer(0)
  for (start in selected) {
    if (visited[start]) next
    queue <- start
    visited[start] <- TRUE
    size <- 0L
    while (length(queue)) {
      current <- queue[[1]]
      queue <- queue[-1]
      size <- size + 1L
      next_nodes <- neighbors[[current]]
      if (length(next_nodes)) {
        next_nodes <- next_nodes[selected_flag[next_nodes] & !visited[next_nodes]]
        if (length(next_nodes)) {
          visited[next_nodes] <- TRUE
          queue <- c(queue, next_nodes)
        }
      }
    }
    component_sizes <- c(component_sizes, size)
  }
  median_nnd <- if (n <= 1L) {
    NA_real_
  } else {
    median(RANN::nn2(as.matrix(coords[selected, .(x, y)]), k = 2)$nn.dists[, 2])
  }
  c(
    max_rook_component_bin_n = max(component_sizes),
    max_component_support_fraction = max(component_sizes) / n,
    connected_component_n = length(component_sizes),
    median_nearest_neighbor_distance = median_nnd
  )
}

draw_results <- list()
observed_results <- list()
validation_results <- list()
bin_strata <- list()
result_index <- 0L
observed_index <- 0L
validation_index <- 0L

for (sid in sample_order) {
  dt <- copy(bins[sample_id == sid])
  setorder(dt, row_in_sample)
  coords <- dt[, .(x, y)]
  neighbors <- build_rook_neighbors(dt)
  strata_indices <- split(seq_len(nrow(dt)), dt$matching_stratum)
  log_message("Null draws: ", sid, "; tissue Bins=", nrow(dt))

  bin_strata[[sid]] <- dt[, .(
    bin_id, sample_id, group, chip_id, x, y, component, total_counts, n_genes,
    PTPRC, strict_same_bin, strict_rook, depth_quintile, n_genes_quintile,
    ptprc_status, matching_stratum
  )]

  for (mode in support_modes) {
    flag_col <- if (mode == "same_bin") "strict_same_bin" else "strict_rook"
    observed <- which(dt[[flag_col]])
    observed_stratum_n <- table(factor(dt$matching_stratum[observed], levels = names(strata_indices)))
    obs_metrics <- spatial_metrics(observed, coords, neighbors)
    observed_index <- observed_index + 1L
    observed_results[[observed_index]] <- data.table(
      sample_id = sid,
      group = dt$group[[1]],
      chip_id = dt$chip_id[[1]],
      support_mode = mode,
      supported_bin_n = length(observed),
      metric = names(obs_metrics),
      observed_value = as.numeric(obs_metrics)
    )

    matched_matrix <- matrix(NA_real_, nrow = random_draw_n, ncol = length(metric_order), dimnames = list(NULL, metric_order))
    unmatched_matrix <- matched_matrix
    exact_match <- logical(random_draw_n)
    for (draw in seq_len(random_draw_n)) {
      matched <- unlist(Map(function(pool, n_take) {
        if (!n_take) integer(0) else sample(pool, n_take, replace = FALSE)
      }, strata_indices, as.integer(observed_stratum_n)), use.names = FALSE)
      unmatched <- if (length(observed)) sample.int(nrow(dt), length(observed), replace = FALSE) else integer(0)
      matched_matrix[draw, ] <- spatial_metrics(matched, coords, neighbors)[metric_order]
      unmatched_matrix[draw, ] <- spatial_metrics(unmatched, coords, neighbors)[metric_order]
      drawn_counts <- table(factor(dt$matching_stratum[matched], levels = names(strata_indices)))
      exact_match[[draw]] <- identical(as.integer(drawn_counts), as.integer(observed_stratum_n))
    }
    assert(all(exact_match), paste("Exact stratum matching failed:", sid, mode))

    for (null_type in c("depth_n_genes_PTPRC_matched", "unmatched_uniform_tissue")) {
      values <- if (null_type == "depth_n_genes_PTPRC_matched") matched_matrix else unmatched_matrix
      result_index <- result_index + 1L
      draw_results[[result_index]] <- data.table(
        sample_id = sid,
        group = dt$group[[1]],
        chip_id = dt$chip_id[[1]],
        support_mode = mode,
        null_type = null_type,
        draw = rep(seq_len(random_draw_n), each = length(metric_order)),
        metric = rep(metric_order, times = random_draw_n),
        null_value = as.numeric(t(values))
      )
    }
    validation_index <- validation_index + 1L
    validation_results[[validation_index]] <- data.table(
      sample_id = sid,
      support_mode = mode,
      supported_bin_n = length(observed),
      random_draw_n = random_draw_n,
      exact_stratum_composition_draw_n = sum(exact_match),
      all_draws_exact = all(exact_match)
    )
  }
}

observed_metrics <- rbindlist(observed_results)
null_draws <- rbindlist(draw_results)
matching_validation <- rbindlist(validation_results)
bin_strata <- rbindlist(bin_strata)
assert(nrow(null_draws) == 6L * 2L * 2L * 1000L * 4L, "Unexpected null draw table size")

summary_rows <- null_draws[, .(
  null_draw_n = sum(!is.na(null_value)),
  null_min = suppressWarnings(min(null_value, na.rm = TRUE)),
  null_q1 = quantile(null_value, 0.25, na.rm = TRUE, names = FALSE),
  null_median = median(null_value, na.rm = TRUE),
  null_q3 = quantile(null_value, 0.75, na.rm = TRUE, names = FALSE),
  null_max = suppressWarnings(max(null_value, na.rm = TRUE))
), by = .(sample_id, group, chip_id, support_mode, null_type, metric)]
summary_rows[!is.finite(null_min), `:=`(null_min = NA_real_, null_max = NA_real_)]
summary_rows <- merge(summary_rows, observed_metrics, by = c("sample_id", "group", "chip_id", "support_mode", "metric"), all.x = TRUE)

percentile_rows <- null_draws[observed_metrics, on = .(sample_id, group, chip_id, support_mode, metric), allow.cartesian = TRUE]
percentile_rows <- percentile_rows[, .(
  empirical_cdf_percentile = if (supported_bin_n[[1]] < 2L || all(is.na(null_value)) || is.na(observed_value[[1]])) NA_real_ else 100 * mean(null_value <= observed_value[[1]], na.rm = TRUE),
  clustering_direction_percentile = if (supported_bin_n[[1]] < 2L || all(is.na(null_value)) || is.na(observed_value[[1]])) NA_real_ else {
    if (metric[[1]] %chin% c("max_rook_component_bin_n", "max_component_support_fraction")) {
      100 * mean(null_value <= observed_value[[1]], na.rm = TRUE)
    } else {
      100 * mean(null_value >= observed_value[[1]], na.rm = TRUE)
    }
  }
), by = .(sample_id, group, chip_id, support_mode, null_type, metric, observed_value, supported_bin_n)]
null_summary <- merge(
  summary_rows,
  percentile_rows,
  by = c("sample_id", "group", "chip_id", "support_mode", "null_type", "metric", "observed_value", "supported_bin_n")
)
null_summary[, metric_label := unname(metric_labels[metric])]
setcolorder(null_summary, c(
  "sample_id", "group", "chip_id", "support_mode", "supported_bin_n", "null_type",
  "metric", "metric_label", "observed_value", "null_draw_n", "null_min", "null_q1",
  "null_median", "null_q3", "null_max", "empirical_cdf_percentile",
  "clustering_direction_percentile"
))

comparison <- dcast(
  null_summary,
  sample_id + group + chip_id + support_mode + supported_bin_n + metric + metric_label + observed_value ~ null_type,
  value.var = c("null_median", "clustering_direction_percentile")
)
comparison[, percentile_change_matched_minus_unmatched :=
             clustering_direction_percentile_depth_n_genes_PTPRC_matched -
             clustering_direction_percentile_unmatched_uniform_tissue]

candidate_summary <- comparison[sample_id %chin% focus_samples]
candidate_summary[, matched_null_interpretation := fifelse(
  is.na(clustering_direction_percentile_depth_n_genes_PTPRC_matched),
  "not_estimable",
  fifelse(
    clustering_direction_percentile_depth_n_genes_PTPRC_matched == 100,
    "observed clustering is at or beyond the clustering-direction maximum of all 1000 matched draws",
    "report empirical percentile without a significance cutoff"
  )
)]

manifest <- data.table(
  field = c(
    "analysis_label", "status", "immutable_bin_input", "immutable_input_md5",
    "fixed_sample_n", "fixed_samples", "support_flags", "matching_depth",
    "matching_complexity", "matching_immune_context", "quintile_method",
    "matched_sampling", "unmatched_sampling", "random_draw_n", "random_seed",
    "rook_definition", "nearest_neighbor_definition", "inference_scope",
    "forbidden_interpretations", "forbidden_analyses"
  ),
  value = c(
    "depth/PTPRC-matched spatial null analysis",
    "DRAFT — pending manual review",
    normalizePath(bin_path, winslash = "/", mustWork = TRUE),
    input_md5[[1]],
    "6",
    paste(sample_order, collapse = "; "),
    "read unchanged from persisted strict_same_bin and strict_rook columns",
    "sample-specific empirical quintile of log1p(total_counts)",
    "sample-specific empirical quintile of n_genes",
    "PTPRC raw-positive versus raw-zero",
    "average rank; ties remain together; nominal empty cells allowed",
    "without replacement within every observed stratum; exact observed stratum counts",
    "uniform without replacement from all tissue Bins; same total support count",
    as.character(random_draw_n),
    as.character(random_seed),
    "native x/y offset 50, same topology component",
    "Euclidean distance in native coordinates among selected Bins, irrespective of component",
    "descriptive spatial QC only; empirical percentiles are not P values",
    "B-cell niche; aggregate; TLS; established disease difference",
    "disease testing; threshold optimization; Moran's I; hotspot scan; new marker selection"
  )
)

input_audit <- data.table(
  input_file = required,
  md5 = input_md5,
  size_bytes = file.info(required)$size,
  exists = file.exists(required)
)

write_tsv(manifest, "00_analysis_manifest.tsv")
write_tsv(input_audit, "01_immutable_input_audit.tsv")
write_tsv(samples, "02_fixed_l3_sample_manifest.tsv")
write_tsv(support_check, "03_frozen_support_flag_validation.tsv")
write_tsv(strata_summary, "04_matching_strata_summary.tsv")
write_tsv(observed_metrics, "05_observed_spatial_metrics.tsv")
fwrite(bin_strata, file.path(table_dir, "06_bin_matching_strata.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")
fwrite(null_draws[null_type == "depth_n_genes_PTPRC_matched"], file.path(table_dir, "07_matched_null_draws.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")
fwrite(null_draws[null_type == "unmatched_uniform_tissue"], file.path(table_dir, "08_unmatched_spatial_null_draws.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")
write_tsv(null_summary, "09_empirical_percentile_summary.tsv")
write_tsv(matching_validation, "10_exact_matching_validation.tsv")
write_tsv(comparison, "11_unmatched_vs_matched_comparison.tsv")
write_tsv(candidate_summary, "12_focus_candidate_summary.tsv")
write_tsv(prior_unmatched, "13_prior_unmatched_depth_overlap_carryforward.tsv")

# Plot 1: all samples, clustering-direction empirical percentiles.
plot_percentile <- copy(null_summary)
plot_percentile[, sample_id := factor(sample_id, levels = rev(sample_order))]
plot_percentile[, metric_label := factor(metric_label, levels = unname(metric_labels))]
p1 <- ggplot(plot_percentile, aes(x = clustering_direction_percentile, y = sample_id, color = null_type, shape = support_mode)) +
  geom_vline(xintercept = 95, linetype = 2, color = "grey55") +
  geom_point(size = 2.2, position = position_dodge(width = 0.55), na.rm = TRUE) +
  facet_wrap(~metric_label, scales = "free_x", ncol = 2) +
  scale_color_manual(values = c(depth_n_genes_PTPRC_matched = "#0072B2", unmatched_uniform_tissue = "#D55E00")) +
  labs(x = "Clustering-direction empirical percentile (%)", y = NULL,
       color = "Null", shape = "Support", title = "Observed spatial clustering relative to matched and unmatched nulls") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
save_png(p1, "01_all_samples_empirical_percentiles", 12, 9)

# Plot 2: focused null distributions with observed values.
focus_draws <- null_draws[sample_id %chin% focus_samples]
focus_obs <- observed_metrics[sample_id %chin% focus_samples]
focus_obs <- focus_obs[, .(null_type = c("depth_n_genes_PTPRC_matched", "unmatched_uniform_tissue")),
                       by = .(sample_id, group, chip_id, support_mode, supported_bin_n,
                              metric, observed_value)]
focus_draws[, metric_label := factor(unname(metric_labels[metric]), levels = unname(metric_labels))]
focus_obs[, metric_label := factor(unname(metric_labels[metric]), levels = unname(metric_labels))]
p2 <- ggplot(focus_draws, aes(x = null_type, y = null_value, fill = null_type)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, na.rm = TRUE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.3, na.rm = TRUE) +
  geom_point(data = focus_obs, aes(x = null_type, y = observed_value), inherit.aes = FALSE,
             color = "black", shape = 23, fill = "#F0E442", size = 2.6,
             position = position_nudge(x = c(-0.0)), na.rm = TRUE) +
  facet_grid(sample_id + support_mode ~ metric_label, scales = "free_y") +
  scale_fill_manual(values = c(depth_n_genes_PTPRC_matched = "#56B4E9", unmatched_uniform_tissue = "#E69F00")) +
  scale_x_discrete(labels = c(depth_n_genes_PTPRC_matched = "Matched", unmatched_uniform_tissue = "Unmatched")) +
  labs(x = NULL, y = "Null metric value", fill = "Null", title = "Focused candidate distributions; diamonds are observed") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none",
        strip.text = element_text(size = 8))
save_png(p2, "02_focus_candidate_null_distributions", 16, 10)

# Plot 3: effect of matching on the clustering-direction percentile.
p3_dt <- copy(comparison)
p3_dt[, metric_label := factor(metric_label, levels = unname(metric_labels))]
p3 <- ggplot(p3_dt, aes(
  x = clustering_direction_percentile_unmatched_uniform_tissue,
  y = clustering_direction_percentile_depth_n_genes_PTPRC_matched,
  color = sample_id, shape = support_mode
)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = 2) +
  geom_hline(yintercept = 95, color = "grey75", linetype = 3) +
  geom_vline(xintercept = 95, color = "grey75", linetype = 3) +
  geom_point(size = 2.2, na.rm = TRUE) +
  facet_wrap(~metric_label, ncol = 2) +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(x = "Unmatched clustering percentile (%)", y = "Matched clustering percentile (%)",
       color = "Sample", shape = "Support", title = "Does depth/PTPRC matching alter the spatial-QC conclusion?") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
save_png(p3, "03_unmatched_vs_matched_percentiles", 11, 9)

# Plot 4: fixed matching strata and observed support composition.
strata_plot <- melt(
  strata_summary,
  id.vars = c("sample_id", "matching_stratum", "stratum_bin_n"),
  measure.vars = c("observed_same_bin_n", "observed_rook_n"),
  variable.name = "support_mode", value.name = "observed_supported_n"
)
strata_plot[, support_mode := fifelse(support_mode == "observed_same_bin_n", "same_bin", "rook")]
strata_plot[, observed_composition_fraction := {
  denominator <- sum(observed_supported_n)
  if (denominator > 0) observed_supported_n / denominator else rep(NA_real_, .N)
}, by = .(sample_id, support_mode)]
p4 <- ggplot(strata_plot[observed_supported_n > 0], aes(x = sample_id, y = observed_composition_fraction,
                                                        fill = matching_stratum)) +
  geom_col(width = 0.75) +
  facet_wrap(~support_mode, ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Observed support composition", fill = "Fixed stratum",
       title = "Depth / n_genes / PTPRC composition held exactly in every matched draw") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
save_png(p4, "04_observed_support_matching_strata", 13, 8)

# Reproducibility and output-level validation.
capture.output(sessioninfo::session_info(), file = file.path(log_dir, "session_info.txt"))
expected <- c(
  file.path(table_dir, sprintf("%02d_%s", 0:13, c(
    "analysis_manifest.tsv", "immutable_input_audit.tsv", "fixed_l3_sample_manifest.tsv",
    "frozen_support_flag_validation.tsv", "matching_strata_summary.tsv",
    "observed_spatial_metrics.tsv", "bin_matching_strata.tsv.gz", "matched_null_draws.tsv.gz",
    "unmatched_spatial_null_draws.tsv.gz", "empirical_percentile_summary.tsv",
    "exact_matching_validation.tsv", "unmatched_vs_matched_comparison.tsv",
    "focus_candidate_summary.tsv", "prior_unmatched_depth_overlap_carryforward.tsv"
  ))),
  file.path(figure_dir, sprintf("%02d_%s.png", 1:4, c(
    "all_samples_empirical_percentiles", "focus_candidate_null_distributions",
    "unmatched_vs_matched_percentiles", "observed_support_matching_strata"
  ))),
  file.path(log_dir, c("analysis.log", "session_info.txt"))
)
validation <- data.table(
  check = c(
    "all_expected_outputs_exist", "fixed_six_samples", "frozen_support_exact",
    "all_1000_matched_draws_exact", "null_rows_complete", "focus_samples_retained",
    "extremely_low_capture_sample_retained"
  ),
  passed = c(
    all(file.exists(expected)),
    setequal(unique(bins$sample_id), sample_order),
    all(support_check$total_n_exact & support_check$same_n_exact & support_check$rook_n_exact),
    all(matching_validation$all_draws_exact & matching_validation$exact_stratum_composition_draw_n == random_draw_n),
    nrow(null_draws) == 96000L,
    all(focus_samples %chin% unique(null_summary$sample_id)),
    "SSC/05957/17B" %chin% unique(null_summary$sample_id)
  ),
  detail = c(
    paste(length(expected), "expected pre-report files"),
    paste(sample_order, collapse = "; "),
    "persisted Bin flags equal frozen aggregate counts",
    paste(random_draw_n, "draws for each of 12 sample/support-mode combinations"),
    as.character(nrow(null_draws)),
    paste(focus_samples, collapse = "; "),
    "retained despite zero same-bin support"
  )
)
write_tsv(validation, "14_output_validation.tsv")
assert(all(validation$passed), "One or more output validations failed")
file.create(file.path(out_dir, "SUCCESS"))
log_message("SUCCESS; output=", out_dir)
cat(out_dir, "\n")
