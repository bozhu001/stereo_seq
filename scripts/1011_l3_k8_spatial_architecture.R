#!/usr/bin/env Rscript

# L3 K=8 spatial architecture analysis.
# Reuses the frozen K=8 assignment and physical Bin50 coordinates. It does not
# read frozen aggregate inputs, rerun clustering, or use disease in any spatial
# definition. Direct Queen adjacency is defined within patient and ROI only.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: 1011_l3_k8_spatial_architecture.R ",
    "<existing_l3_dir> <resolution_dir> <output_dir>"
  )
}

existing_dir <- normalizePath(args[[1]], mustWork = TRUE)
resolution_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
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
  message("[L3 K8 spatial architecture] ", stage, if (nzchar(detail)) paste0(": ", detail))
}

save_pdf <- function(plot, filename, width, height) {
  ggsave(
    file.path(fig_dir, filename),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf,
    limitsize = FALSE
  )
}

read_tsv <- function(path) {
  if (grepl("[.]gz$", path)) {
    return(fread(cmd = paste("gzip -cd", shQuote(path))))
  }
  fread(path)
}

connected_membership <- function(node_ids, edge_i, edge_j) {
  node_ids <- unique(node_ids)
  n <- length(node_ids)
  if (!n) return(integer())
  parent <- seq_len(n)
  names(parent) <- node_ids
  find_root <- function(x) {
    while (parent[[x]] != x) {
      parent[[x]] <<- parent[[parent[[x]]]]
      x <- parent[[x]]
    }
    x
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
  roots <- vapply(seq_len(n), find_root, integer(1))
  as.integer(factor(roots, levels = unique(roots)))
}

pair_grid <- CJ(niche_a_num = 1:8, niche_b_num = 1:8)[niche_a_num <= niche_b_num]
pair_grid[, pair_id := (niche_a_num - 1L) * 8L + niche_b_num]
pair_grid[, `:=`(
  niche_a = paste0("N", niche_a_num),
  niche_b = paste0("N", niche_b_num),
  niche_pair = paste0("N", niche_a_num, "-N", niche_b_num),
  pair_type = fifelse(niche_a_num == niche_b_num, "homotypic", "heterotypic")
)]

count_pairs <- function(labels, edge_i, edge_j) {
  a <- pmin(labels[edge_i], labels[edge_j])
  b <- pmax(labels[edge_i], labels[edge_j])
  ids <- (a - 1L) * 8L + b
  counts <- tabulate(ids, nbins = 64L)
  counts[pair_grid$pair_id]
}

write_status("read_and_validate_frozen_inputs")
input_path <- file.path(existing_dir, "01_FROZEN_REFERENCE18_INPUT.tsv.gz")
k8_assignment_path <- file.path(resolution_dir, "K08", "02_BIN_NICHE_ASSIGNMENTS.tsv.gz")
color_path <- file.path(resolution_dir, "K08", "00_NICHE_COLOR_KEY.tsv")
zero_neighbor_path <- file.path(resolution_dir, "00_ZERO_NEIGHBOR_NOT_EVALUABLE.tsv")
annotation_candidates <- Sys.glob(file.path(
  resolution_dir,
  "K08_HIGH_RESOLUTION_SENSITIVITY_*",
  "08_K08_CONSERVATIVE_NICHE_ANNOTATIONS.tsv"
))
if (length(annotation_candidates) != 1L) {
  stop("Expected exactly one finalized K8 annotation table; found ", length(annotation_candidates))
}
annotation_path <- annotation_candidates[[1]]

input <- read_tsv(input_path)
k8_assignment <- read_tsv(k8_assignment_path)
k8_colors_dt <- fread(color_path)
zero_neighbor <- fread(zero_neighbor_path)
k8_annotations <- fread(annotation_path)
k8_levels <- paste0("N", 1:8)
k8_colors <- setNames(k8_colors_dt$color, k8_colors_dt$niche)

stopifnot(
  nrow(input) == 40139L,
  nrow(k8_assignment) == 39982L,
  nrow(zero_neighbor) == 157L,
  uniqueN(input$patient_id) == 7L,
  !anyDuplicated(input$entity_key),
  !anyDuplicated(k8_assignment$entity_key),
  all(k8_assignment$entity_key %chin% input$entity_key),
  setequal(unique(k8_assignment$niche), k8_levels),
  all(input$x_um * 2 == input$x_raw),
  all(input$y_um * 2 == input$y_raw)
)

map_data <- merge(
  input[, .(
    entity_key, bin_id, patient_id, disease, chip_id, roi_id,
    analysis_block_id, x_um, y_um, nUMI, nGene
  )],
  k8_assignment[, .(entity_key, internal_id, niche, selected_seed, neighbor_count)],
  by = "entity_key",
  all.x = TRUE,
  sort = FALSE
)
map_data[, evaluable := !is.na(niche)]
map_data[, niche := factor(niche, levels = k8_levels)]
setorder(map_data, patient_id, roi_id, y_um, x_um)

if (map_data[evaluable == TRUE, .N] != 39982L || map_data[evaluable == FALSE, .N] != 157L) {
  stop("Evaluable/non-evaluable Bin counts do not match frozen inputs")
}

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

write_status("build_queen_adjacency", "four unique offsets; no cross-patient or cross-ROI edges")
right_nodes <- nodes[, .(
  patient_id,
  roi_id,
  target_x = x_um,
  target_y = y_um,
  node_j = node_id
)]
offsets <- data.table(dx = c(25, 0, 25, 25), dy = c(0, 25, 25, -25))
edge_list <- lapply(seq_len(nrow(offsets)), function(offset_index) {
  candidate <- nodes[, .(
    patient_id,
    roi_id,
    target_x = x_um + offsets$dx[[offset_index]],
    target_y = y_um + offsets$dy[[offset_index]],
    node_i = node_id
  )]
  merge(
    candidate,
    right_nodes,
    by = c("patient_id", "roi_id", "target_x", "target_y"),
    all = FALSE,
    sort = FALSE
  )[, .(patient_id, roi_id, node_i, node_j)]
})
edges <- unique(rbindlist(edge_list))
edges[, edge_id := .I]
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

write_status("spatial_maps", "three multipage PDFs; one page per patient/ROI")
patient_order <- unique(map_data[order(factor(disease, levels = c("HC", "IPF", "SSc-ILD")), patient_id), patient_id])

plot_patient_map <- function(patient, mode = c("all", "n3", "key")) {
  mode <- match.arg(mode)
  current <- copy(map_data[patient_id == patient])
  current[, display_label := as.character(niche)]
  if (mode == "n3") {
    current[evaluable == TRUE & display_label != "N3", display_label := "Other evaluable niche"]
  } else if (mode == "key") {
    current[evaluable == TRUE & !display_label %chin% c("N2", "N3", "N4", "N7"), display_label := "Other evaluable niche"]
  }
  current[evaluable == FALSE, display_label := "Neighborhood not evaluable"]
  palette <- c(k8_colors, "Other evaluable niche" = "#D3D3D3", "Neighborhood not evaluable" = "#F2F2F2")
  title_suffix <- switch(
    mode,
    all = "K=8 spatial niche assignment",
    n3 = "N3 B-enriched mixed niche",
    key = "Highlighted K=8 niches: N2, N3, N4 and N7"
  )
  ggplot(current, aes(x_um, y_um, fill = display_label)) +
    geom_tile(width = 25, height = 25, linewidth = 0) +
    scale_fill_manual(values = palette, drop = FALSE, name = "K=8 niche") +
    scale_y_reverse() +
    coord_fixed() +
    facet_wrap(~roi_id) +
    labs(
      title = paste0(patient, ": ", title_suffix),
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

map_files <- c(
  all = "01_K08_ALL_PATIENT_SPATIAL_NICHE_MAPS_EN.pdf",
  n3 = "02_K08_N3_SPATIAL_MAPS_EN.pdf",
  key = "03_K08_KEY_NICHES_SPATIAL_MAPS_EN.pdf"
)
for (mode in names(map_files)) {
  grDevices::cairo_pdf(file.path(fig_dir, map_files[[mode]]), width = 12, height = 8)
  for (patient in patient_order) print(plot_patient_map(patient, mode))
  grDevices::dev.off()
}

write_status("connected_component_continuity", "56 patient-by-niche summaries")
continuity_rows <- list()
component_memberships <- list()
row_index <- 0L
for (patient in patient_order) {
  patient_nodes <- nodes[patient_id == patient]
  patient_edges <- edges[patient_id == patient]
  for (current_niche in k8_levels) {
    current_nodes <- patient_nodes[niche == current_niche]
    current_edges <- patient_edges[niche_i == current_niche & niche_j == current_niche]
    membership <- connected_membership(
      current_nodes$node_id,
      current_edges$node_i,
      current_edges$node_j
    )
    if (!length(membership)) next
    membership_dt <- data.table(
      node_id = current_nodes$node_id,
      patient_id = patient,
      roi_id = current_nodes$roi_id,
      niche = current_niche,
      component_number = membership
    )
    component_memberships[[length(component_memberships) + 1L]] <- membership_dt
    sizes <- membership_dt[, .N, by = component_number]$N
    row_index <- row_index + 1L
    continuity_rows[[row_index]] <- data.table(
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
component_membership <- rbindlist(component_memberships)
setorder(continuity, patient_id, niche)
fwrite(continuity, file.path(output_dir, "01_K08_NICHE_SPATIAL_CONTINUITY.tsv"), sep = "\t", na = "NA")

p_component_size <- ggplot(
  continuity,
  aes(niche, median_component_size, color = niche)
) +
  geom_point(aes(shape = disease), position = position_jitter(width = 0.1, height = 0, seed = 2026100711L), size = 2.4, alpha = 0.9) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.45) +
  scale_color_manual(values = k8_colors, guide = "none") +
  labs(
    title = "K=8 niche connected-component size across L3 patients",
    subtitle = "Points are patients; black bars are cross-patient medians",
    x = "K=8 niche",
    y = "Median connected-component size (Bin50)",
    shape = "Disease group"
  ) +
  theme_bw(base_size = 10)
save_pdf(p_component_size, "04_K08_NICHE_COMPONENT_SIZE_EN.pdf", 9, 5.8)

p_largest <- ggplot(
  continuity,
  aes(niche, largest_component_fraction, color = niche)
) +
  geom_point(aes(shape = disease), position = position_jitter(width = 0.1, height = 0, seed = 2026100711L), size = 2.4, alpha = 0.9) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.45) +
  scale_color_manual(values = k8_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Largest connected-component fraction for K=8 niches",
    subtitle = "Points are patients; black bars are cross-patient medians",
    x = "K=8 niche",
    y = "Largest component / all niche Bin50",
    shape = "Disease group"
  ) +
  theme_bw(base_size = 10)
save_pdf(p_largest, "05_K08_NICHE_LARGEST_COMPONENT_FRACTION_EN.pdf", 9, 5.8)

write_status("direct_adjacency_counts_and_fractions")
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
  scope = "all_L3",
  patient_id = "ALL",
  disease = "ALL",
  total_unique_edge_n = sum(
    adjacency_counts[, unique(total_unique_edge_n), by = patient_id]$V1
  ),
  unique_edge_count = sum(unique_edge_count)
), by = .(niche_a, niche_b, niche_pair, pair_type)]
adjacency_counts <- rbindlist(list(adjacency_counts, overall_counts), use.names = TRUE)
setcolorder(adjacency_counts, c(
  "scope", "patient_id", "disease", "niche_a", "niche_b", "niche_pair",
  "pair_type", "unique_edge_count", "total_unique_edge_n"
))
fwrite(adjacency_counts, file.path(output_dir, "02_K08_NICHE_ADJACENCY_COUNTS.tsv"), sep = "\t")

make_fraction_rows <- function(count_table) {
  output <- list()
  for (row_number in seq_len(nrow(count_table))) {
    current <- count_table[row_number]
    a <- current$niche_a
    b <- current$niche_b
    count <- current$unique_edge_count
    if (a == b) {
      output[[length(output) + 1L]] <- data.table(
        scope = current$scope,
        patient_id = current$patient_id,
        disease = current$disease,
        source_niche = a,
        neighboring_niche = b,
        unique_edge_count = count,
        source_incident_edge_end_n = 2 * count,
        total_unique_edge_n = current$total_unique_edge_n
      )
    } else {
      output[[length(output) + 1L]] <- data.table(
        scope = current$scope,
        patient_id = current$patient_id,
        disease = current$disease,
        source_niche = c(a, b),
        neighboring_niche = c(b, a),
        unique_edge_count = count,
        source_incident_edge_end_n = count,
        total_unique_edge_n = current$total_unique_edge_n
      )
    }
  }
  result <- rbindlist(output)
  result[, source_total_incident_edge_end_n := sum(source_incident_edge_end_n), by = .(scope, patient_id, source_niche)]
  result[, row_normalized_adjacency_fraction := source_incident_edge_end_n / source_total_incident_edge_end_n]
  result[, total_edge_normalized_fraction := unique_edge_count / total_unique_edge_n]
  result[, same_niche_adjacency_fraction := row_normalized_adjacency_fraction[source_niche == neighboring_niche], by = .(scope, patient_id, source_niche)]
  result[, heterotypic_adjacency_proportion := 1 - same_niche_adjacency_fraction]
  result
}
adjacency_fractions <- make_fraction_rows(adjacency_counts)
setorder(adjacency_fractions, scope, patient_id, source_niche, neighboring_niche)
fwrite(adjacency_fractions, file.path(output_dir, "03_K08_NICHE_ADJACENCY_FRACTIONS.tsv"), sep = "\t")

heatmap_data <- adjacency_fractions[scope == "all_L3"]
heatmap_data[, `:=`(
  source_niche = factor(source_niche, levels = k8_levels),
  neighboring_niche = factor(neighboring_niche, levels = k8_levels)
)]
p_adjacency <- ggplot(
  heatmap_data,
  aes(neighboring_niche, source_niche, fill = row_normalized_adjacency_fraction)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", row_normalized_adjacency_fraction)), size = 3) +
  scale_fill_viridis_c(option = "C", name = "Row-normalized\nadjacency fraction") +
  coord_equal() +
  labs(
    title = "Direct K=8 niche adjacency across evaluable L3 Bin50",
    subtitle = "Queen adjacency within patient and ROI; each physical edge counted once",
    x = "Neighboring niche",
    y = "Source niche"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())
save_pdf(p_adjacency, "06_K08_NICHE_ADJACENCY_HEATMAP_EN.pdf", 8, 7)

write_status("patient_roi_permutation_null", "999 label permutations within each patient and ROI")
permutation_n <- 999L
permutation_seed <- 2026100711L
small_offset <- 0.5
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
    for (indices in roi_indices) permuted[indices] <- sample(labels[indices], replace = FALSE)
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
  write_status("permutation_patient_complete", paste0(patient_index, "/", length(patient_order), ": ", patient))
}
patient_enrichment <- rbindlist(permutation_rows)
setorder(patient_enrichment, patient_id, niche_a, niche_b)
fwrite(
  patient_enrichment,
  file.path(output_dir, "04_K08_PATIENT_NICHE_ADJACENCY_ENRICHMENT.tsv"),
  sep = "\t",
  na = "NA"
)

write_status("cross_patient_adjacency_consistency")
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
      fifelse(max(positive_n, negative_n) >= 5L, "Same direction in 5/7 patients", "Patient-variable direction")
    ),
    max_absolute_driver_patient = if (length(values)) patients[[driver_index]] else NA_character_,
    max_absolute_enrichment = if (length(values)) abs(values[[driver_index]]) else NA_real_,
    max_absolute_share = if (absolute_total > 0) abs(values[[driver_index]]) / absolute_total else NA_real_,
    single_patient_driven = if (absolute_total > 0) abs(values[[driver_index]]) / absolute_total >= 0.5 else NA
  )
}, by = .(niche_a, niche_b, niche_pair)]
setorder(cross_patient, -positive_patient_n, -median_log2_enrichment)
fwrite(cross_patient, file.path(output_dir, "05_K08_CROSS_PATIENT_ADJACENCY_SUMMARY.tsv"), sep = "\t", na = "NA")

n3_summary <- copy(cross_patient[niche_a == "N3" | niche_b == "N3"])
n3_summary[, neighboring_niche := fifelse(niche_a == "N3", niche_b, niche_a)]
n3_raw <- adjacency_counts[scope == "all_L3" & (niche_a == "N3" | niche_b == "N3"), .(
  niche_pair,
  all_l3_unique_direct_edge_count = unique_edge_count,
  all_l3_total_edge_fraction = unique_edge_count / total_unique_edge_n
)]
n3_summary <- merge(n3_summary, n3_raw, by = "niche_pair", all.x = TRUE)
setorder(n3_summary, -positive_patient_n, -median_log2_enrichment)
fwrite(n3_summary, file.path(output_dir, "06_K08_N3_NEIGHBOR_SUMMARY.tsv"), sep = "\t", na = "NA")

# Unbiased global prioritization is computed before the N3-targeted follow-up.
# All 28 heterotypic pairs remain in the table; ranking is descriptive and does
# not alter the spatial analysis or remove any pair.
continuity_for_priority <- continuity[, .(
  endpoint_median_largest_component_fraction = median(largest_component_fraction),
  endpoint_median_fragmentation_index = median(fragmentation_index)
), by = niche]
annotation_for_priority <- k8_annotations[, .(
  niche,
  proposed_annotation,
  confidence,
  UMI_robust,
  core_state_single_patient_driven
)]
global_priority <- copy(cross_patient)
global_priority <- merge(
  global_priority,
  continuity_for_priority,
  by.x = "niche_a",
  by.y = "niche",
  all.x = TRUE
)
setnames(
  global_priority,
  c("endpoint_median_largest_component_fraction", "endpoint_median_fragmentation_index"),
  c("niche_a_median_largest_component_fraction", "niche_a_median_fragmentation_index")
)
global_priority <- merge(
  global_priority,
  continuity_for_priority,
  by.x = "niche_b",
  by.y = "niche",
  all.x = TRUE
)
setnames(
  global_priority,
  c("endpoint_median_largest_component_fraction", "endpoint_median_fragmentation_index"),
  c("niche_b_median_largest_component_fraction", "niche_b_median_fragmentation_index")
)
global_priority <- merge(global_priority, annotation_for_priority, by.x = "niche_a", by.y = "niche", all.x = TRUE)
setnames(
  global_priority,
  c("proposed_annotation", "confidence", "UMI_robust", "core_state_single_patient_driven"),
  paste0(c("proposed_annotation", "confidence", "UMI_robust", "core_state_single_patient_driven"), "_a")
)
global_priority <- merge(global_priority, annotation_for_priority, by.x = "niche_b", by.y = "niche", all.x = TRUE)
setnames(
  global_priority,
  c("proposed_annotation", "confidence", "UMI_robust", "core_state_single_patient_driven"),
  paste0(c("proposed_annotation", "confidence", "UMI_robust", "core_state_single_patient_driven"), "_b")
)
global_priority[, `:=`(
  mean_endpoint_median_largest_component_fraction = (
    niche_a_median_largest_component_fraction + niche_b_median_largest_component_fraction
  ) / 2,
  absolute_median_log2_enrichment = abs(median_log2_enrichment),
  consistent_direction_patient_n = pmax(positive_patient_n, negative_patient_n),
  both_endpoints_umi_robust = UMI_robust_a %in% TRUE & UMI_robust_b %in% TRUE,
  endpoint_core_single_patient_driven = core_state_single_patient_driven_a %in% TRUE |
    core_state_single_patient_driven_b %in% TRUE,
  pair_or_endpoint_single_patient_driven = single_patient_driven %in% TRUE |
    core_state_single_patient_driven_a %in% TRUE |
    core_state_single_patient_driven_b %in% TRUE,
  both_endpoints_high_confidence = confidence_a == "High" & confidence_b == "High"
)]
setorder(
  global_priority,
  -consistent_direction_patient_n,
  pair_or_endpoint_single_patient_driven,
  -absolute_median_log2_enrichment,
  -mean_endpoint_median_largest_component_fraction,
  -both_endpoints_umi_robust,
  -both_endpoints_high_confidence,
  niche_pair
)
global_priority[, global_priority_rank := .I]
global_priority[, priority_interpretation := fifelse(
  consistent_direction_patient_n >= 6L &
    !pair_or_endpoint_single_patient_driven &
    both_endpoints_umi_robust,
  "High-priority reproducible spatial association",
  fifelse(
    consistent_direction_patient_n >= 5L & !pair_or_endpoint_single_patient_driven,
    "Recurrent but requires cautious biological characterization",
    "Exploratory or patient-variable spatial association"
  )
)]
fwrite(
  global_priority,
  file.path(output_dir, "07A_K08_GLOBAL_SPATIAL_PRIORITY.tsv"),
  sep = "\t",
  na = "NA"
)

patient_enrichment[, pair_label := factor(niche_pair, levels = sort(unique(niche_pair)))]
p_cross <- ggplot(
  patient_enrichment[pair_type == "heterotypic"],
  aes(pair_label, log2_observed_expected_enrichment, color = patient_id)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 2026100711L), size = 1.7, alpha = 0.85) +
  coord_flip() +
  labs(
    title = "Patient-level heterotypic K=8 niche adjacency enrichment",
    subtitle = "Observed versus 999 within-patient/ROI label permutations; each point is one patient",
    x = "Heterotypic niche pair",
    y = "log2(observed / expected adjacency)",
    color = "Patient"
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")
save_pdf(p_cross, "07_K08_CROSS_PATIENT_ADJACENCY_ENRICHMENT_EN.pdf", 10, 11)

n3_patient <- patient_enrichment[pair_type == "heterotypic" & (niche_a == "N3" | niche_b == "N3")]
n3_patient[, neighboring_niche := fifelse(niche_a == "N3", niche_b, niche_a)]
n3_patient[, neighboring_niche := factor(neighboring_niche, levels = setdiff(k8_levels, "N3"))]
p_n3 <- ggplot(
  n3_patient,
  aes(neighboring_niche, log2_observed_expected_enrichment, color = patient_id)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 2026100711L), size = 2.3, alpha = 0.9) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.45) +
  labs(
    title = "N3 direct-neighbor adjacency enrichment across L3 patients",
    subtitle = "N3 is the B-enriched mixed niche; points are patients and black bars are medians",
    x = "N3 neighboring niche",
    y = "log2(observed / expected adjacency)",
    color = "Patient"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")
save_pdf(p_n3, "08_K08_N3_NEIGHBOR_ENRICHMENT_EN.pdf", 9.5, 6)

write_status("n3_component_boundary_context")
n3_members <- component_membership[niche == "N3"]
n3_nodes <- merge(
  n3_members,
  nodes[, .(node_id, disease, x_um, y_um)],
  by = "node_id",
  all.x = TRUE
)
n3_nodes[, component_id := paste0(
  gsub("[^A-Za-z0-9]+", "_", patient_id), "__",
  gsub("[^A-Za-z0-9]+", "_", roi_id), "__N3C",
  sprintf("%03d", component_number)
)]
node_component <- setNames(n3_nodes$component_id, n3_nodes$node_id)

n3_contact_i <- edges[niche_i == "N3" & niche_j != "N3", .(
  patient_id,
  roi_id,
  n3_node_id = node_i,
  neighbor_node_id = node_j,
  neighboring_niche = niche_j
)]
n3_contact_j <- edges[niche_j == "N3" & niche_i != "N3", .(
  patient_id,
  roi_id,
  n3_node_id = node_j,
  neighbor_node_id = node_i,
  neighboring_niche = niche_i
)]
n3_contacts <- rbindlist(list(n3_contact_i, n3_contact_j))
n3_contacts[, component_id := unname(node_component[as.character(n3_node_id)])]

component_base <- n3_nodes[, .(
  patient_id = first(patient_id),
  disease = first(disease),
  roi_id = first(roi_id),
  n3_bin_n = .N
), by = component_id]
component_boundary <- n3_contacts[, .(
  boundary_bin_n = uniqueN(n3_node_id),
  total_direct_heterotypic_contact_edge_n = .N
), by = component_id]
component_base <- merge(component_base, component_boundary, by = "component_id", all.x = TRUE)
component_base[is.na(boundary_bin_n), `:=`(
  boundary_bin_n = 0L,
  total_direct_heterotypic_contact_edge_n = 0L
)]

n3_context <- CJ(component_id = component_base$component_id, neighboring_niche = k8_levels)
n3_context <- merge(n3_context, component_base, by = "component_id", all.x = TRUE)
contact_summary <- n3_contacts[, .(
  direct_contact_edge_n = .N,
  boundary_bin_touch_n = uniqueN(n3_node_id)
), by = .(component_id, neighboring_niche)]
n3_context <- merge(n3_context, contact_summary, by = c("component_id", "neighboring_niche"), all.x = TRUE)
n3_context[is.na(direct_contact_edge_n), `:=`(direct_contact_edge_n = 0L, boundary_bin_touch_n = 0L)]
n3_context[, direct_contact_edge_fraction := fifelse(
  total_direct_heterotypic_contact_edge_n > 0,
  direct_contact_edge_n / total_direct_heterotypic_contact_edge_n,
  NA_real_
)]
n3_context[, boundary_bin_touch_fraction := fifelse(
  boundary_bin_n > 0,
  boundary_bin_touch_n / boundary_bin_n,
  NA_real_
)]
n3_context[, nearest_neighboring_niche := {
  max_count <- max(direct_contact_edge_n)
  if (max_count == 0L) "None" else paste(neighboring_niche[direct_contact_edge_n == max_count], collapse = "; ")
}, by = component_id]
setorder(n3_context, patient_id, roi_id, component_id, neighboring_niche)
fwrite(n3_context, file.path(output_dir, "07_K08_N3_COMPONENT_CONTEXT.tsv"), sep = "\t", na = "NA")

p_boundary <- ggplot(
  n3_context[neighboring_niche != "N3" & !is.na(boundary_bin_touch_fraction)],
  aes(neighboring_niche, boundary_bin_touch_fraction, color = patient_id)
) +
  geom_point(aes(size = n3_bin_n), alpha = 0.55, position = position_jitter(width = 0.12, height = 0, seed = 2026100711L)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.45) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Direct niche context at N3 connected-component boundaries",
    subtitle = "Each point is one N3 component; a boundary Bin may touch more than one niche",
    x = "Neighboring niche",
    y = "Fraction of N3 boundary Bin50 touching niche",
    color = "Patient",
    size = "N3 component\nBin50 count"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")
save_pdf(p_boundary, "09_K08_N3_COMPONENT_BOUNDARY_CONTEXT_EN.pdf", 10, 6.5)

write_status("write_findings_and_audit")
continuity_rank <- continuity[, .(
  median_largest_component_fraction = median(largest_component_fraction),
  median_fragmentation_index = median(fragmentation_index),
  patients_with_component = uniqueN(patient_id)
), by = niche][order(-median_largest_component_fraction, median_fragmentation_index)]
top_continuous <- continuity_rank[1:3]
top_fragmented <- continuity_rank[order(-median_fragmentation_index)][1:3]

n3_continuity <- continuity[niche == "N3"]
n3_raw_neighbor_rank <- n3_summary[order(-all_l3_unique_direct_edge_count)]
n3_consistency_rank <- n3_summary[
  order(-pmax(positive_patient_n, negative_patient_n), -abs(median_log2_enrichment))
]
top_global_pairs <- global_priority[1:3]
stable_pairs <- cross_patient[
  positive_patient_n >= 6L | negative_patient_n >= 6L
][order(-pmax(positive_patient_n, negative_patient_n), -abs(median_log2_enrichment))]
driver_pairs <- cross_patient[single_patient_driven %in% TRUE]

format_niche_rank <- function(dt, value_column) {
  paste0(dt$niche, " (", sprintf("%.3f", dt[[value_column]]), ")", collapse = "; ")
}

findings <- c(
  "# L3 K=8 niche内部空间架构：主要发现",
  "",
  "## 范围与解释边界",
  "",
  paste0("本分析复用冻结K=8标签，覆盖", nrow(nodes), "个neighborhood-evaluable Bin50、", uniqueN(nodes$patient_id), "位患者和", uniqueN(nodes$roi_id), "个ROI。"),
  "Queen邻接仅在同一患者、同一ROI内建立；中心距离在x、y方向均不超过25 um，且每条物理边只计一次。",
  paste0("置换零模型在每个patient/ROI内部固定坐标和niche数量，随机打乱标签", permutation_n, "次。"),
  "K8标签本身来自50 um邻域Reference18组成。因此，标签的局部相似不是独立验证；本轮新增信息是连通分量、真实Bin-to-Bin界面及跨患者接触方向的一致性。",
  "",
  "## 1. 空间连续性",
  "",
  paste0("按患者最大连通分量占比的中位数排序，连续性最高的三个niche为：", format_niche_rank(top_continuous, "median_largest_component_fraction"), "。"),
  paste0("碎片化指数中位数最高的三个niche为：", format_niche_rank(top_fragmented, "median_fragmentation_index"), "。"),
  paste0(
    "N3在", uniqueN(n3_continuity$patient_id), "/7位患者均可评估；其最大连通分量占比中位数为",
    sprintf("%.3f", median(n3_continuity$largest_component_fraction)),
    "，每位患者的连通分量数范围为", min(n3_continuity$connected_component_n), "-", max(n3_continuity$connected_component_n),
    "。这说明N3可形成真实空间patch，同时也可能包含多个分离component，而不是单一连续区块。"
  ),
  "",
  "## 2. 全部28个heterotypic pair的无偏比较",
  "",
  "N1-N8均以完全相同的方法进入空间连续性、直接邻接、patient/ROI内置换及跨患者一致性分析；没有使用N3或既往B-cell问题筛选pair。",
  paste0(
    "综合跨患者方向一致性、非单患者驱动、observed/expected效应量、两端niche连续性、Reference18可解释性及UMI稳健性后，排序最前的三个pair为：",
    paste0(
      top_global_pairs$niche_pair, " (", top_global_pairs$consistent_direction_patient_n, "/7 same direction; median |log2 O/E|=",
      sprintf("%.3f", top_global_pairs$absolute_median_log2_enrichment), ")",
      collapse = "; "
    ), "。完整28个pair及全部判据见07A表。"
  ),
  "该排序仅用于选择后续生物学表征对象，不改变任何空间结果，也不是新的显著性门槛。",
  "",
  "## 3. N3 targeted follow-up",
  "",
  paste0(
    "按未校正的直接接触边数，N3最常接触的三个niche为：",
    paste0(
      n3_raw_neighbor_rank$neighboring_niche[1:3], " (",
      n3_raw_neighbor_rank$all_l3_unique_direct_edge_count[1:3], " unique edges)",
      collapse = "; "
    ), "。"
  ),
  paste0(
    "校正各patient/ROI内niche丰度后，N3方向最一致的三个空间关系为：",
    paste0(
      n3_consistency_rank$neighboring_niche[1:3], " (",
      pmax(n3_consistency_rank$positive_patient_n[1:3], n3_consistency_rank$negative_patient_n[1:3]),
      "/7 same direction; median log2 O/E=",
      sprintf("%.3f", n3_consistency_rank$median_log2_enrichment[1:3]), ")",
      collapse = "; "
    ), "。负值表示低于随机标签期望的直接界面，而不是“没有空间结构”。"
  ),
  paste0(
    "N3-N2、N3-N4、N3-N7的结果分别为：",
    paste0(
      n3_summary[match(c("N2", "N4", "N7"), neighboring_niche), neighboring_niche], " ",
      n3_summary[match(c("N2", "N4", "N7"), neighboring_niche), positive_patient_n], "/7正向、",
      n3_summary[match(c("N2", "N4", "N7"), neighboring_niche), negative_patient_n], "/7负向",
      collapse = "; "
    ), "。应结合06表中的效应量与置换P值逐项解释。"
  ),
  "直接接触次数受niche丰度影响；正文判断优先使用patient/ROI内置换后的observed/expected，而不只看raw edge count。",
  "该置换保持每个ROI内niche数量，但会破坏标签的空间自相关，也不保留相邻中心Bin之间50 um输入邻域重叠带来的相关性。因此，广泛的heterotypic负向富集可反映标签patch化与构建尺度，不能全部解释为特异性生物学排斥。",
  "",
  "## 4. 跨患者重复性与单患者驱动",
  "",
  paste0("共有", nrow(stable_pairs), "/28个heterotypic pair在至少6/7患者中方向一致。"),
  if (nrow(stable_pairs)) paste0("这些pair包括：", paste(stable_pairs$niche_pair, collapse = "; "), "。") else "没有heterotypic pair达到6/7方向一致。",
  if (nrow(driver_pairs)) paste0("按单患者贡献>=50%规则，可能单患者驱动的pair为：", paste(driver_pairs$niche_pair, collapse = "; "), "。") else "没有heterotypic pair达到预设的单患者驱动条件。",
  "方向一致不等同于每位患者置换检验显著；两者在05表中分别保留。",
  "",
  "## 5. 是否存在有组织的spatial architecture",
  "",
  if (nrow(stable_pairs)) {
    "K8标签场相对于可交换标签零模型呈现明显空间组织：多个niche形成连通component，部分heterotypic界面在多个患者中保持同向observed/expected偏移。但K8由重叠的50 um邻域组成构建，置换又破坏了这种内生空间自相关；因此这里只能确认有组织的niche label field，不能把全部偏移视为独立生物学验证、细胞互作或因果关系。"
  } else {
    "K8标签可形成连通component，但heterotypic界面尚未显示至少6/7患者同向的重复模式；因此目前只能支持局部空间连续性，不能主张稳定的跨患者界面偏好。"
  },
  paste0(
    "全局绝对效应排序优先的pair为：", paste(top_global_pairs$niche_pair, collapse = "; "),
    "；当前三者均为负向observed/expected，主要代表重复的空间分离/少接触。正向界面偏好应单独查看05与07A表，不能与负向分离混为同一故事。"
  ),
  "",
  "## 6. 明确未做事项",
  "",
  "- 未读取或使用frozen aggregate标签、成员或边界。",
  "- disease仅作为患者图注，未参与空间定义，也未进行疾病显著性检验。",
  "- 未做DEG、GSEA、pathway、CellChat、ligand-receptor或cell-cell communication。",
  "- adjacency只表示直接Bin50空间接触或界面偏好，不等于communication。",
  "- 同一Bin50可包含多个细胞；niche label不是单细胞类型。"
)
writeLines(findings, file.path(output_dir, "08_K08_SPATIAL_ARCHITECTURE_KEY_FINDINGS_ZH.md"), useBytes = TRUE)

readme <- c(
  "# L3 K=8 spatial architecture",
  "",
  "本目录复用冻结K=8 assignment和x_um/y_um物理坐标。",
  "",
  "- Queen adjacency: same patient and ROI; |dx| <= 25 um and |dy| <= 25 um; not the same Bin.",
  "- Only the four unique offsets (25,0), (0,25), (25,25), (25,-25) are constructed, so each physical edge is counted once.",
  "- Fragmentation index = 1 - largest component fraction.",
  "- Spatial coverage = niche Bin count / all evaluable Bin count in that patient.",
  "- Bin50 area is reported as 625 um2 because the frozen physical conversion defines each Bin50 as 25 x 25 um.",
  "- Row-normalized adjacency uses incident edge ends: a homotypic edge contributes two ends to its source niche; a heterotypic edge contributes one end to each source niche.",
  paste0("- Null: ", permutation_n, " label permutations within each patient/ROI; seed ", permutation_seed, "; log2 offset ", small_offset, "."),
  "- Empirical P is two-sided from upper/lower permutation tails with +1 correction; BH is performed separately within each patient for 28 heterotypic and 8 homotypic pairs.",
  "- The exchangeable-label null does not preserve niche-label spatial autocorrelation or overlap among adjacent 50 um input neighborhoods; enrichment therefore describes departure from label exchangeability, not independent validation of a biological interface.",
  "- Frozen aggregate inputs were not read. Disease was not used in graph construction or testing."
)
writeLines(readme, file.path(output_dir, "README_ZH.md"), useBytes = TRUE)

manifest <- data.table(
  role = c(
    "frozen Reference18 UMI>=50 input", "fixed K8 assignment", "K8 color key",
    "zero-neighbor audit", "K8 conservative annotation", "analysis script", "aggregate input"
  ),
  path = c(
    input_path, k8_assignment_path, color_path, zero_neighbor_path,
    annotation_path, "scripts/1011_l3_k8_spatial_architecture.R", "NOT READ"
  ),
  used = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
)
fwrite(manifest, file.path(output_dir, "INPUT_AND_SCRIPT_MANIFEST.tsv"), sep = "\t")

required_outputs <- c(
  "01_K08_NICHE_SPATIAL_CONTINUITY.tsv",
  "02_K08_NICHE_ADJACENCY_COUNTS.tsv",
  "03_K08_NICHE_ADJACENCY_FRACTIONS.tsv",
  "04_K08_PATIENT_NICHE_ADJACENCY_ENRICHMENT.tsv",
  "05_K08_CROSS_PATIENT_ADJACENCY_SUMMARY.tsv",
  "06_K08_N3_NEIGHBOR_SUMMARY.tsv",
  "07_K08_N3_COMPONENT_CONTEXT.tsv",
  "07A_K08_GLOBAL_SPATIAL_PRIORITY.tsv",
  "08_K08_SPATIAL_ARCHITECTURE_KEY_FINDINGS_ZH.md",
  "README_ZH.md",
  "figures/01_K08_ALL_PATIENT_SPATIAL_NICHE_MAPS_EN.pdf",
  "figures/02_K08_N3_SPATIAL_MAPS_EN.pdf",
  "figures/03_K08_KEY_NICHES_SPATIAL_MAPS_EN.pdf",
  "figures/04_K08_NICHE_COMPONENT_SIZE_EN.pdf",
  "figures/05_K08_NICHE_LARGEST_COMPONENT_FRACTION_EN.pdf",
  "figures/06_K08_NICHE_ADJACENCY_HEATMAP_EN.pdf",
  "figures/07_K08_CROSS_PATIENT_ADJACENCY_ENRICHMENT_EN.pdf",
  "figures/08_K08_N3_NEIGHBOR_ENRICHMENT_EN.pdf",
  "figures/09_K08_N3_COMPONENT_BOUNDARY_CONTEXT_EN.pdf"
)

audit <- rbindlist(list(
  data.table(
    item = c(
      "K8 assignment source", "evaluable Bin50", "patients", "ROIs",
      "adjacency definition", "spatial distance definition", "unique physical edges",
      "permutation number", "permutation seed", "frozen aggregate used",
      "disease used in spatial definition", "Bin-level disease significance test",
      "pathway", "cell-cell communication", "all K8 niches analyzed symmetrically",
      "all heterotypic pairs in main summary"
    ),
    value = c(
      k8_assignment_path, nrow(nodes), uniqueN(nodes$patient_id), uniqueN(nodes$roi_id),
      "Queen adjacency within patient and ROI; each edge counted once",
      "abs(dx)<=25 um and abs(dy)<=25 um; excluding self",
      nrow(edges), permutation_n, permutation_seed, FALSE, FALSE, FALSE, FALSE, FALSE,
      uniqueN(continuity$niche), nrow(cross_patient)
    ),
    expected = c(
      "frozen K08 assignment", 39982L, 7L, 7L, "fixed", "25 um Queen", ">0",
      999L, "recorded", FALSE, FALSE, FALSE, FALSE, FALSE, 8L, 28L
    ),
    passed = c(
      TRUE, nrow(nodes) == 39982L, uniqueN(nodes$patient_id) == 7L,
      uniqueN(nodes$roi_id) == 7L, TRUE, TRUE, nrow(edges) > 0L,
      permutation_n == 999L, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      uniqueN(continuity$niche) == 8L, nrow(cross_patient) == 28L
    ),
    detail = c(
      "No clustering rerun", "Frozen evaluable set", "L3 only", "No cross-ROI edges",
      "Four unique offsets", "Uses x_um/y_um", "Graph validation", "Fixed null size",
      "Reproducible stochastic step", "Aggregate files not read", "Disease annotation only",
      "No disease test", "Not run", "Not run", "N1-N8 use identical definitions",
      "All unordered heterotypic pairs retained"
    )
  ),
  rbindlist(lapply(required_outputs, function(relative_path) {
    full_path <- file.path(output_dir, relative_path)
    data.table(
      item = paste0("output: ", relative_path),
      value = if (file.exists(full_path)) file.info(full_path)$size else 0,
      expected = "exists and non-empty",
      passed = file.exists(full_path) && file.info(full_path)$size > 0,
      detail = "File-level check"
    )
  }))
), use.names = TRUE, fill = TRUE)
fwrite(audit, file.path(output_dir, "09_K08_SPATIAL_ARCHITECTURE_FINAL_AUDIT.tsv"), sep = "\t", na = "NA")
if (!all(audit$passed)) stop("Final audit contains failed checks")

writeLines(capture.output(sessionInfo()), file.path(output_dir, "SESSION_INFO.txt"))
write_status("complete", paste0(length(required_outputs), " required outputs passed"))
writeLines("COMPLETE", file.path(output_dir, "COMPLETE.ok"))

cat("RESULT_DIR=", output_dir, "\n", sep = "")
cat("TOP_CONTINUOUS_NICHES=", paste(top_continuous$niche, collapse = ";"), "\n", sep = "")
cat(
  "N3_COMPONENTS=YES; median_largest_component_fraction=",
  sprintf("%.3f", median(n3_continuity$largest_component_fraction)), "\n", sep = ""
)
cat(
  "N3_STABLE_RELATIONS=",
  paste0(
    head(n3_consistency_rank$neighboring_niche, 3), "(",
    ifelse(head(n3_consistency_rank$median_log2_enrichment, 3) > 0, "enriched", "depleted"), ")",
    collapse = ";"
  ),
  "\n",
  sep = ""
)
cat("PRIORITY_PAIRS=", paste(top_global_pairs$niche_pair, collapse = ";"), "\n", sep = "")
cat(
  "PRIORITY_FILES=08_K08_SPATIAL_ARCHITECTURE_KEY_FINDINGS_ZH.md;",
  "06_K08_N3_NEIGHBOR_SUMMARY.tsv;05_K08_CROSS_PATIENT_ADJACENCY_SUMMARY.tsv;",
  "figures/03_K08_KEY_NICHES_SPATIAL_MAPS_EN.pdf;",
  "figures/08_K08_N3_NEIGHBOR_ENRICHMENT_EN.pdf\n",
  sep = ""
)
