#!/usr/bin/env Rscript

# Symmetric J2 K=5 spatial architecture audit. Reuses frozen K=5 labels and
# physical Bin50 coordinates. It does not read frozen aggregate inputs, rerun
# clustering, or use disease in any spatial definition or significance test.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: 1018_j2_k5_spatial_architecture.R ",
    "<formal_j2_dir> <j2_resolution_dir> <output_dir>"
  )
}

formal_dir <- normalizePath(args[[1]], mustWork = TRUE)
resolution_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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
  message("[J2 K5 spatial architecture] ", stage, if (nzchar(detail)) paste0(": ", detail))
}

save_figure <- function(plot, stem, width, height) {
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf,
    limitsize = FALSE
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    limitsize = FALSE
  )
}

read_table <- function(path) {
  if (grepl("[.]gz$", path)) {
    return(fread(cmd = paste("gzip -cd", shQuote(path))))
  }
  fread(path)
}

connected_membership <- function(node_ids, edge_i, edge_j) {
  node_ids <- unique(node_ids)
  node_n <- length(node_ids)
  if (!node_n) return(integer())
  parent <- seq_len(node_n)
  find_root <- function(index) {
    while (parent[[index]] != index) {
      parent[[index]] <<- parent[[parent[[index]]]]
      index <- parent[[index]]
    }
    index
  }
  if (length(edge_i)) {
    local_i <- match(edge_i, node_ids)
    local_j <- match(edge_j, node_ids)
    keep <- !is.na(local_i) & !is.na(local_j)
    for (edge_index in which(keep)) {
      root_i <- find_root(local_i[[edge_index]])
      root_j <- find_root(local_j[[edge_index]])
      if (root_i != root_j) parent[[root_j]] <- root_i
    }
  }
  roots <- vapply(seq_len(node_n), find_root, integer(1))
  as.integer(factor(roots, levels = unique(roots)))
}

niche_levels <- paste0("N", 1:5)
pair_grid <- CJ(niche_a_num = 1:5, niche_b_num = 1:5)[niche_a_num <= niche_b_num]
pair_grid[, pair_id := (niche_a_num - 1L) * 5L + niche_b_num]
pair_grid[, `:=`(
  niche_a = paste0("N", niche_a_num),
  niche_b = paste0("N", niche_b_num),
  niche_pair = paste0("N", niche_a_num, "-N", niche_b_num),
  pair_type = fifelse(niche_a_num == niche_b_num, "homotypic", "heterotypic")
)]

count_pairs <- function(labels, edge_i, edge_j) {
  left <- pmin(labels[edge_i], labels[edge_j])
  right <- pmax(labels[edge_i], labels[edge_j])
  ids <- (left - 1L) * 5L + right
  counts <- tabulate(ids, nbins = 25L)
  counts[pair_grid$pair_id]
}

write_status("read_and_validate_inputs")
input_path <- file.path(formal_dir, "01_FROZEN_REFERENCE18_INPUT.tsv.gz")
assignment_path <- file.path(resolution_dir, "K05", "02_BIN_NICHE_ASSIGNMENTS.tsv.gz")
color_path <- file.path(resolution_dir, "K05", "00_NICHE_COLOR_KEY.tsv")
zero_path <- file.path(resolution_dir, "00_ZERO_NEIGHBOR_NOT_EVALUABLE.tsv")
annotation_path <- file.path(resolution_dir, "11_J2_PRIMARYK_CONSERVATIVE_ANNOTATIONS.tsv")

input <- read_table(input_path)
assignment <- read_table(assignment_path)
colors_table <- read_table(color_path)
zero_neighbor <- read_table(zero_path)
annotations <- read_table(annotation_path)
niche_colors <- setNames(colors_table$color, colors_table$niche)

stopifnot(
  nrow(input) == 58840L,
  nrow(assignment) == 58772L,
  nrow(zero_neighbor) == 68L,
  uniqueN(input$patient_id) == 7L,
  !anyDuplicated(input$entity_key),
  !anyDuplicated(assignment$entity_key),
  all(assignment$entity_key %chin% input$entity_key),
  setequal(unique(assignment$niche), niche_levels),
  all(input$x_um * 2 == input$x_raw),
  all(input$y_um * 2 == input$y_raw),
  all(annotations$niche %chin% niche_levels)
)

map_data <- merge(
  input[, .(
    entity_key, bin_id, patient_id, disease, chip_id, roi_id,
    analysis_block_id, x_um, y_um, nUMI, nGene
  )],
  assignment[, .(entity_key, internal_id, niche, selected_seed, neighbor_count)],
  by = "entity_key",
  all.x = TRUE,
  sort = FALSE
)
map_data[, evaluable := !is.na(niche)]
map_data[, niche := factor(niche, levels = niche_levels)]
setorder(map_data, patient_id, roi_id, y_um, x_um)

stopifnot(
  map_data[evaluable == TRUE, .N] == 58772L,
  map_data[evaluable == FALSE, .N] == 68L
)

nodes <- map_data[evaluable == TRUE, .(
  node_id = .I,
  entity_key,
  internal_id,
  bin_id,
  patient_id,
  disease,
  roi_id,
  x_um,
  y_um,
  nUMI,
  nGene,
  niche = as.character(niche),
  niche_num = as.integer(niche)
)]
stopifnot(!anyDuplicated(nodes[, .(patient_id, roi_id, x_um, y_um)]))

write_status("build_queen_adjacency", "same patient/ROI; four unique 25 um offsets")
right_nodes <- nodes[, .(
  patient_id,
  roi_id,
  target_x = x_um,
  target_y = y_um,
  node_j = node_id
)]
offsets <- data.table(dx = c(25, 0, 25, 25), dy = c(0, 25, 25, -25))
edge_list <- lapply(seq_len(nrow(offsets)), function(offset_index) {
  candidates <- nodes[, .(
    patient_id,
    roi_id,
    target_x = x_um + offsets$dx[[offset_index]],
    target_y = y_um + offsets$dy[[offset_index]],
    node_i = node_id
  )]
  merge(
    candidates,
    right_nodes,
    by = c("patient_id", "roi_id", "target_x", "target_y"),
    all = FALSE,
    sort = FALSE
  )[, .(patient_id, roi_id, node_i, node_j)]
})
edges <- unique(rbindlist(edge_list))
edges <- merge(
  edges,
  nodes[, .(node_i = node_id, niche_i = niche, niche_i_num = niche_num)],
  by = "node_i",
  all.x = TRUE,
  sort = FALSE
)
edges <- merge(
  edges,
  nodes[, .(node_j = node_id, niche_j = niche, niche_j_num = niche_num)],
  by = "node_j",
  all.x = TRUE,
  sort = FALSE
)
stopifnot(
  !anyDuplicated(edges[, .(pmin(node_i, node_j), pmax(node_i, node_j))]),
  all(edges$node_i != edges$node_j),
  all(edges$patient_id == nodes$patient_id[edges$node_i]),
  all(edges$roi_id == nodes$roi_id[edges$node_i])
)

patient_order <- unique(map_data[
  order(factor(disease, levels = c("HC", "IPF", "SSc-ILD")), patient_id),
  patient_id
])

write_status("all_patient_spatial_maps", paste0(length(patient_order), " patient pages"))
plot_patient_map <- function(patient) {
  current <- copy(map_data[patient_id == patient])
  current[, display_label := as.character(niche)]
  current[evaluable == FALSE, display_label := "Neighborhood not evaluable"]
  palette <- c(niche_colors, "Neighborhood not evaluable" = "#F2F2F2")
  ggplot(current, aes(x_um, y_um, fill = display_label)) +
    geom_tile(width = 25, height = 25, linewidth = 0) +
    scale_fill_manual(values = palette, drop = FALSE, name = "K=5 niche") +
    scale_y_reverse() +
    coord_fixed() +
    facet_wrap(~roi_id) +
    labs(
      title = paste0(patient, ": J2 K=5 spatial niche assignment"),
      subtitle = paste0(
        "Disease group: ", unique(current$disease),
        "; evaluable Bin50: ", current[evaluable == TRUE, .N],
        "; Bin50 size: 25 x 25 um"
      ),
      x = "x coordinate (um)",
      y = "y coordinate (um)"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = grid::unit(0.8, "cm"),
      strip.text = element_text(size = 8),
      plot.title = element_text(face = "bold")
    )
}

map_pdf <- file.path(figure_dir, "01_K05_ALL_PATIENT_SPATIAL_NICHE_MAPS_EN.pdf")
grDevices::cairo_pdf(map_pdf, width = 12, height = 8)
for (patient in patient_order) print(plot_patient_map(patient))
grDevices::dev.off()

write_status("connected_components", "symmetric patient-by-niche continuity")
continuity_rows <- list()
for (patient in patient_order) {
  patient_nodes <- nodes[patient_id == patient]
  patient_edges <- edges[patient_id == patient]
  for (current_niche in niche_levels) {
    current_nodes <- patient_nodes[niche == current_niche]
    current_edges <- patient_edges[niche_i == current_niche & niche_j == current_niche]
    membership <- connected_membership(
      current_nodes$node_id,
      current_edges$node_i,
      current_edges$node_j
    )
    if (!length(membership)) next
    sizes <- data.table(component = membership)[, .N, by = component]$N
    continuity_rows[[length(continuity_rows) + 1L]] <- data.table(
      patient_id = patient,
      disease = unique(current_nodes$disease),
      niche = current_niche,
      patient_evaluable_bin_n = nrow(patient_nodes),
      niche_bin_n = nrow(current_nodes),
      niche_area_um2 = nrow(current_nodes) * 625,
      connected_component_n = length(sizes),
      largest_component_bin_n = max(sizes),
      largest_component_fraction = max(sizes) / sum(sizes),
      median_component_size = median(sizes),
      mean_component_size = mean(sizes),
      singleton_component_n = sum(sizes == 1L),
      singleton_component_fraction = mean(sizes == 1L),
      singleton_bin_fraction = sum(sizes == 1L) / sum(sizes),
      fragmentation_index = 1 - max(sizes) / sum(sizes),
      spatial_coverage = nrow(current_nodes) / nrow(patient_nodes)
    )
  }
}
continuity <- rbindlist(continuity_rows)
setorder(continuity, patient_id, niche)
fwrite(
  continuity,
  file.path(output_dir, "01_K05_NICHE_SPATIAL_CONTINUITY.tsv"),
  sep = "\t",
  na = "NA"
)

p_component <- ggplot(continuity, aes(niche, median_component_size, color = niche)) +
  geom_point(
    aes(shape = disease),
    position = position_jitter(width = 0.1, height = 0, seed = 2026100811L),
    size = 2.5,
    alpha = 0.9
  ) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.45) +
  scale_color_manual(values = niche_colors, guide = "none") +
  labs(
    title = "J2 K=5 niche connected-component size across patients",
    subtitle = "Points are patients; black bars are cross-patient medians",
    x = "K=5 niche",
    y = "Median connected-component size (Bin50)",
    shape = "Disease group"
  ) +
  theme_bw(base_size = 10)
save_figure(p_component, "02_K05_NICHE_COMPONENT_SIZE_EN", 9, 5.8)

p_largest <- ggplot(continuity, aes(niche, largest_component_fraction, color = niche)) +
  geom_point(
    aes(shape = disease),
    position = position_jitter(width = 0.1, height = 0, seed = 2026100811L),
    size = 2.5,
    alpha = 0.9
  ) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.45) +
  scale_color_manual(values = niche_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Largest connected-component fraction for J2 K=5 niches",
    subtitle = "Points are patients; black bars are cross-patient medians",
    x = "K=5 niche",
    y = "Largest component / all niche Bin50",
    shape = "Disease group"
  ) +
  theme_bw(base_size = 10)
save_figure(p_largest, "03_K05_NICHE_LARGEST_COMPONENT_FRACTION_EN", 9, 5.8)

write_status("direct_niche_adjacency")
patient_count_rows <- lapply(patient_order, function(patient) {
  current_edges <- edges[patient_id == patient]
  observed <- count_pairs(nodes$niche_num, current_edges$node_i, current_edges$node_j)
  cbind(
    data.table(
      scope = "patient",
      patient_id = patient,
      disease = unique(nodes[patient_id == patient, disease]),
      total_unique_edge_n = nrow(current_edges)
    ),
    pair_grid[, .(niche_a, niche_b, niche_pair, pair_type)],
    unique_edge_count = observed
  )
})
adjacency_counts <- rbindlist(patient_count_rows)
overall_counts <- adjacency_counts[, .(
  scope = "all_J2",
  patient_id = "ALL",
  disease = "ALL",
  total_unique_edge_n = sum(unique(total_unique_edge_n)),
  unique_edge_count = sum(unique_edge_count)
), by = .(niche_a, niche_b, niche_pair, pair_type)]
adjacency_counts <- rbindlist(list(adjacency_counts, overall_counts), use.names = TRUE)
setcolorder(adjacency_counts, c(
  "scope", "patient_id", "disease", "niche_a", "niche_b", "niche_pair",
  "pair_type", "unique_edge_count", "total_unique_edge_n"
))
fwrite(adjacency_counts, file.path(output_dir, "02_K05_NICHE_ADJACENCY_COUNTS.tsv"), sep = "\t")

make_fraction_rows <- function(count_table) {
  rows <- list()
  for (row_index in seq_len(nrow(count_table))) {
    current <- count_table[row_index]
    if (current$niche_a == current$niche_b) {
      rows[[length(rows) + 1L]] <- data.table(
        scope = current$scope,
        patient_id = current$patient_id,
        disease = current$disease,
        source_niche = current$niche_a,
        neighboring_niche = current$niche_b,
        unique_edge_count = current$unique_edge_count,
        source_incident_edge_end_n = 2 * current$unique_edge_count,
        total_unique_edge_n = current$total_unique_edge_n
      )
    } else {
      rows[[length(rows) + 1L]] <- data.table(
        scope = current$scope,
        patient_id = current$patient_id,
        disease = current$disease,
        source_niche = c(current$niche_a, current$niche_b),
        neighboring_niche = c(current$niche_b, current$niche_a),
        unique_edge_count = current$unique_edge_count,
        source_incident_edge_end_n = current$unique_edge_count,
        total_unique_edge_n = current$total_unique_edge_n
      )
    }
  }
  result <- rbindlist(rows)
  result[, source_total_incident_edge_end_n := sum(source_incident_edge_end_n), by = .(scope, patient_id, source_niche)]
  result[, row_normalized_adjacency_fraction := source_incident_edge_end_n / source_total_incident_edge_end_n]
  result[, total_edge_normalized_fraction := unique_edge_count / total_unique_edge_n]
  result[, same_niche_adjacency_fraction := row_normalized_adjacency_fraction[source_niche == neighboring_niche], by = .(scope, patient_id, source_niche)]
  result[, heterotypic_adjacency_proportion := 1 - same_niche_adjacency_fraction]
  result
}

adjacency_fractions <- make_fraction_rows(adjacency_counts)
setorder(adjacency_fractions, scope, patient_id, source_niche, neighboring_niche)
fwrite(adjacency_fractions, file.path(output_dir, "03_K05_NICHE_ADJACENCY_FRACTIONS.tsv"), sep = "\t")

heatmap_data <- adjacency_fractions[scope == "all_J2"]
heatmap_data[, `:=`(
  source_niche = factor(source_niche, levels = niche_levels),
  neighboring_niche = factor(neighboring_niche, levels = niche_levels)
)]
p_adjacency <- ggplot(
  heatmap_data,
  aes(neighboring_niche, source_niche, fill = row_normalized_adjacency_fraction)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", row_normalized_adjacency_fraction)), size = 3.2) +
  scale_fill_viridis_c(option = "C", name = "Row-normalized\nadjacency fraction") +
  coord_equal() +
  labs(
    title = "Direct J2 K=5 niche adjacency",
    subtitle = "Queen adjacency within patient and ROI; each physical edge counted once",
    x = "Neighboring niche",
    y = "Source niche"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())
save_figure(p_adjacency, "04_K05_NICHE_ADJACENCY_HEATMAP_EN", 8, 7)

write_status("permutation_null", "999 label permutations within patient and ROI")
permutation_n <- 999L
permutation_seed <- 2026100811L
small_offset <- 0.5
permutation_checkpoint <- file.path(
  output_dir,
  "04_K05_PATIENT_NICHE_ADJACENCY_ENRICHMENT.tsv"
)
if (file.exists(permutation_checkpoint) && file.info(permutation_checkpoint)$size > 0) {
  patient_enrichment <- fread(permutation_checkpoint)
  stopifnot(
    nrow(patient_enrichment) == length(patient_order) * nrow(pair_grid),
    uniqueN(patient_enrichment$patient_id) == length(patient_order)
  )
  write_status("permutation_checkpoint_reused", "complete 7-patient table")
} else {
  set.seed(permutation_seed)
  permutation_rows <- list()
  for (patient_index in seq_along(patient_order)) {
    patient <- patient_order[[patient_index]]
    current_nodes <- nodes[patient_id == patient]
    current_edges <- edges[patient_id == patient]
    global_to_local <- setNames(seq_len(nrow(current_nodes)), current_nodes$node_id)
    edge_i <- unname(global_to_local[as.character(current_edges$node_i)])
    edge_j <- unname(global_to_local[as.character(current_edges$node_j)])
    labels <- current_nodes$niche_num
    observed <- count_pairs(labels, edge_i, edge_j)
    null_counts <- matrix(0L, nrow = permutation_n, ncol = nrow(pair_grid))
    roi_indices <- split(seq_len(nrow(current_nodes)), current_nodes$roi_id)
    for (permutation_index in seq_len(permutation_n)) {
      permuted <- labels
      for (indices in roi_indices) {
        permuted[indices] <- sample(labels[indices], replace = FALSE)
      }
      null_counts[permutation_index, ] <- count_pairs(permuted, edge_i, edge_j)
    }
    expected <- colMeans(null_counts)
    null_sd <- apply(null_counts, 2, sd)
    p_upper <- (1 + colSums(sweep(null_counts, 2, observed, FUN = ">="))) / (permutation_n + 1)
    p_lower <- (1 + colSums(sweep(null_counts, 2, observed, FUN = "<="))) / (permutation_n + 1)
    p_two_sided <- pmin(1, 2 * pmin(p_upper, p_lower))
    current <- cbind(
      data.table(
        patient_id = patient,
        disease = unique(current_nodes$disease),
        patient_evaluable_bin_n = nrow(current_nodes),
        patient_unique_edge_n = nrow(current_edges)
      ),
      pair_grid[, .(niche_a, niche_b, niche_pair, pair_type)],
      observed_adjacency_count = observed,
      expected_adjacency_count = expected,
      null_sd = null_sd,
      log2_observed_expected_enrichment = log2((observed + small_offset) / (expected + small_offset)),
      empirical_p_enrichment = p_upper,
      empirical_p_depletion = p_lower,
      empirical_p_two_sided = p_two_sided
    )
    current[, bh_family := paste0("within_patient_", pair_type)]
    current[, bh_fdr := p.adjust(empirical_p_two_sided, method = "BH"), by = pair_type]
    permutation_rows[[patient_index]] <- current
    write_status(
      "permutation_patient_complete",
      paste0(patient_index, "/", length(patient_order), ": ", patient)
    )
  }

  patient_enrichment <- rbindlist(permutation_rows)
  setorder(patient_enrichment, patient_id, niche_a, niche_b)
  fwrite(patient_enrichment, permutation_checkpoint, sep = "\t", na = "NA")
}

write_status("cross_patient_spatial_consistency")
cross_patient <- patient_enrichment[pair_type == "heterotypic", {
  values <- log2_observed_expected_enrichment[is.finite(log2_observed_expected_enrichment)]
  patients <- patient_id[is.finite(log2_observed_expected_enrichment)]
  absolute_total <- sum(abs(values))
  driver_index <- if (length(values)) which.max(abs(values)) else NA_integer_
  positive_n <- sum(values > 0)
  negative_n <- sum(values < 0)
  evaluable_n <- length(values)
  list(
    median_log2_enrichment = median(values),
    mean_log2_enrichment = mean(values),
    sd_log2_enrichment = sd(values),
    positive_patient_n = positive_n,
    negative_patient_n = negative_n,
    zero_patient_n = evaluable_n - positive_n - negative_n,
    evaluable_patient_n = evaluable_n,
    positive_fraction = positive_n / evaluable_n,
    median_permutation_p_two_sided = median(empirical_p_two_sided, na.rm = TRUE),
    direction_consistency = fifelse(
      max(positive_n, negative_n) >= 6L,
      paste0(ifelse(positive_n >= negative_n, "Positive", "Negative"), " in >=6/7 patients"),
      fifelse(
        max(positive_n, negative_n) >= 5L,
        paste0(ifelse(positive_n >= negative_n, "Positive", "Negative"), " in 5/7 patients"),
        "Patient-variable direction"
      )
    ),
    max_absolute_driver_patient = if (length(values)) patients[[driver_index]] else NA_character_,
    max_absolute_enrichment = if (length(values)) abs(values[[driver_index]]) else NA_real_,
    max_absolute_share = if (absolute_total > 0) abs(values[[driver_index]]) / absolute_total else NA_real_,
    single_patient_driven = if (absolute_total > 0) abs(values[[driver_index]]) / absolute_total >= 0.5 else NA
  )
}, by = .(niche_a, niche_b, niche_pair)]
cross_patient[, direction_patient_n := pmax(positive_patient_n, negative_patient_n)]
cross_patient[, absolute_median_log2_enrichment := abs(median_log2_enrichment)]
setorder(cross_patient, -direction_patient_n, -absolute_median_log2_enrichment)
fwrite(
  cross_patient,
  file.path(output_dir, "05_K05_CROSS_PATIENT_ADJACENCY_SUMMARY.tsv"),
  sep = "\t",
  na = "NA"
)

plot_enrichment <- copy(patient_enrichment[pair_type == "heterotypic"])
plot_enrichment[, niche_pair := factor(niche_pair, levels = sort(unique(niche_pair)))]
p_cross <- ggplot(
  plot_enrichment,
  aes(niche_pair, log2_observed_expected_enrichment, color = disease)
) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 2026100811L),
    size = 2.2,
    alpha = 0.85
  ) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.45) +
  coord_flip() +
  labs(
    title = "Cross-patient J2 K=5 heterotypic adjacency enrichment",
    subtitle = "Observed versus 999 within-patient/ROI label permutations; each point is one patient",
    x = "Heterotypic niche pair",
    y = "log2 observed / expected direct adjacency",
    color = "Disease group"
  ) +
  theme_bw(base_size = 10)
save_figure(p_cross, "05_K05_CROSS_PATIENT_ADJACENCY_ENRICHMENT_EN", 9.5, 7.5)

continuity_summary <- continuity[, .(
  median_largest_component_fraction = median(largest_component_fraction),
  median_fragmentation_index = median(fragmentation_index),
  median_component_size = median(median_component_size)
), by = niche][order(-median_largest_component_fraction)]
fwrite(
  continuity_summary,
  file.path(output_dir, "06_K05_CROSS_PATIENT_CONTINUITY_SUMMARY.tsv"),
  sep = "\t"
)

top_continuous <- continuity_summary[1:min(3L, .N), paste0(
  niche, " (median largest-component fraction ",
  sprintf("%.3f", median_largest_component_fraction), ")"
)]
recurrent_pairs <- cross_patient[
  pmax(positive_patient_n, negative_patient_n) >= 6L
][order(-pmax(positive_patient_n, negative_patient_n), -abs(median_log2_enrichment))]

findings <- c(
  "# J2 K=5 spatial architecture audit：主要发现",
  "",
  "## 范围与方法",
  "",
  paste0(
    "本轮复用冻结的J2 K=5 assignment，对", format(nrow(nodes), big.mark = ","),
    "个neighborhood-evaluable Bin50进行对称空间审核。Queen邻接仅在同一患者、同一ROI内建立，",
    "要求|dx|和|dy|均不超过25 µm且每条物理边只计一次。"
  ),
  paste0(
    "共纳入", nrow(edges), "条唯一直接邻接边；五个niche和全部10个heterotypic niche pair均进入分析。",
    "observed/expected参照来自每位患者/ROI内999次标签置换。"
  ),
  "",
  "## 空间连续性",
  "",
  paste0("按跨患者最大连通分量占比中位数排序，连续性最高的三个niche为：", paste(top_continuous, collapse = "；"), "。"),
  "连通分量结果描述K=5标签场在真实坐标中的连续或碎片化程度；它不等于对应细胞类型形成了纯组织区室。",
  "",
  "## niche–niche直接邻接",
  "",
  if (nrow(recurrent_pairs)) {
    paste0(
      "共有", nrow(recurrent_pairs), "个heterotypic pair在至少6/7患者中保持同一方向。方向最一致的关系包括：",
      paste(head(paste0(
        recurrent_pairs$niche_pair, " [", recurrent_pairs$direction_consistency,
        ", median log2(O/E)=", sprintf("%.3f", recurrent_pairs$median_log2_enrichment), "]"
      ), 5), collapse = "；"), "。"
    )
  } else {
    "没有heterotypic pair在至少6/7患者中保持同一方向；现有直接邻接差异主要表现为患者可变。"
  },
  "高同类邻接不是正文的主要生物学证据。K=5标签来自重叠的50 µm邻域组成，置换参照又破坏内生空间自相关，因此这里仅支持空间连续性或界面偏好，不能视为独立验证、细胞互作或因果关系。",
  "",
  "## 解释边界",
  "",
  "- 未读取或使用frozen aggregate标签、成员或边界。",
  "- disease仅用于患者图注，未参与空间定义，未做疾病显著性检验。",
  "- 未做DEG/GSEA、pathway、CellChat或ligand-receptor分析。",
  "- 未进行跨芯片consensus或将J2 niche编号映射到L3。"
)
writeLines(findings, file.path(output_dir, "07_K05_SPATIAL_ARCHITECTURE_FINDINGS_ZH.md"), useBytes = TRUE)

readme <- c(
  "# J2 K=5 spatial architecture audit",
  "",
  "This directory reuses the frozen J2 K=5 assignment and evaluates all five niches symmetrically.",
  "",
  "- Queen adjacency: same patient and ROI; |dx| <= 25 um and |dy| <= 25 um; self-pairs excluded.",
  "- Each physical edge is counted once.",
  "- The patient/ROI null uses 999 label permutations and preserves niche counts.",
  "- Frozen aggregate inputs are not read.",
  "- Disease is metadata only; no disease-significance test is performed.",
  "- Adjacency is spatial association, not cell-cell communication."
)
writeLines(readme, file.path(output_dir, "README_ZH.md"), useBytes = TRUE)

audit <- data.table(
  item = c(
    "chip", "K assignment source", "evaluable Bin50", "not-evaluable Bin50",
    "patients", "ROIs", "adjacency definition", "unique adjacency edges",
    "heterotypic pair count", "permutation number", "permutation seed",
    "frozen aggregate used", "disease used in spatial definition",
    "disease significance", "DEG/GSEA", "pathway", "cell-cell communication"
  ),
  value = c(
    "Y40105J2", assignment_path, nrow(nodes), nrow(map_data[evaluable == FALSE]),
    uniqueN(nodes$patient_id), uniqueN(nodes[, .(patient_id, roi_id)]),
    "Queen: same patient/ROI; |dx|<=25 um; |dy|<=25 um; no self-edge",
    nrow(edges), nrow(pair_grid[pair_type == "heterotypic"]), permutation_n,
    permutation_seed, "FALSE", "FALSE", "FALSE", "FALSE", "FALSE", "FALSE"
  ),
  passed = TRUE
)

required_outputs <- c(
  "01_K05_NICHE_SPATIAL_CONTINUITY.tsv",
  "02_K05_NICHE_ADJACENCY_COUNTS.tsv",
  "03_K05_NICHE_ADJACENCY_FRACTIONS.tsv",
  "04_K05_PATIENT_NICHE_ADJACENCY_ENRICHMENT.tsv",
  "05_K05_CROSS_PATIENT_ADJACENCY_SUMMARY.tsv",
  "06_K05_CROSS_PATIENT_CONTINUITY_SUMMARY.tsv",
  "07_K05_SPATIAL_ARCHITECTURE_FINDINGS_ZH.md",
  "README_ZH.md",
  "figures/01_K05_ALL_PATIENT_SPATIAL_NICHE_MAPS_EN.pdf",
  "figures/02_K05_NICHE_COMPONENT_SIZE_EN.pdf",
  "figures/03_K05_NICHE_LARGEST_COMPONENT_FRACTION_EN.pdf",
  "figures/04_K05_NICHE_ADJACENCY_HEATMAP_EN.pdf",
  "figures/05_K05_CROSS_PATIENT_ADJACENCY_ENRICHMENT_EN.pdf"
)
output_audit <- rbindlist(lapply(required_outputs, function(relative_path) {
  path <- file.path(output_dir, relative_path)
  data.table(
    item = paste0("output::", relative_path),
    value = normalizePath(path, mustWork = FALSE),
    passed = file.exists(path) && file.info(path)$size > 0
  )
}))
audit <- rbindlist(list(audit, output_audit), fill = TRUE)
fwrite(audit, file.path(output_dir, "08_K05_SPATIAL_ARCHITECTURE_FINAL_AUDIT.tsv"), sep = "\t")
stopifnot(all(audit$passed))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "SESSION_INFO.txt"))
writeLines("COMPLETE", file.path(output_dir, "COMPLETE.ok"))
write_status("complete", paste0("output_dir=", output_dir))

cat("RESULT_DIR=", output_dir, "\n", sep = "")
cat("EVALUABLE_BIN50=", nrow(nodes), "\n", sep = "")
cat("UNIQUE_QUEEN_EDGES=", nrow(edges), "\n", sep = "")
cat("TOP_CONTINUOUS=", paste(continuity_summary$niche[1:3], collapse = ","), "\n", sep = "")
cat("RECURRENT_HETEROTYPIC_PAIRS=", nrow(recurrent_pairs), "\n", sep = "")
