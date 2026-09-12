#!/usr/bin/env Rscript

# Derived from scripts/627_cdr3_per_bin_top10_heatmaps.R.
# Only the frozen per-Bin TopN input and Top10 -> Top20 labels are changed.
# This script does not read a raw expression matrix or recompute Top20 ranks.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(digest)
  library(grid)
  library(jsonlite)
  library(magick)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
top20_dir <- if (length(args) >= 1L) args[[1L]] else
  "results/reanalysis/CDR3_PER_BIN_GENE_LIST_20260911_161928"
top10_dir <- if (length(args) >= 2L) args[[2L]] else
  "results/reanalysis/CDR3_PER_BIN_TOP10_FINAL_20260911_165003"
top10_heatmap_dir <- if (length(args) >= 3L) args[[3L]] else
  "results/reanalysis/CDR3_PER_BIN_TOP10_HEATMAPS_20260911_172841"
output_dir <- if (length(args) >= 4L) args[[4L]] else file.path(
  "results/reanalysis",
  paste0("CDR3_PER_BIN_TOP20_HEATMAP_FINAL_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)

paths <- c(
  all = file.path(top20_dir, "02_ALL69_CDR3_BIN_TOP20_GENES.tsv"),
  bin_qc = file.path(top20_dir, "03_ALL69_CDR3_BIN_QC.tsv"),
  shared = file.path(top20_dir, "04_SHARED7_CDR3_OCCURRENCE_TOP20.tsv"),
  repeated = file.path(top20_dir, "05_WITHIN_PATIENT_REPEATED_TOP20.tsv"),
  molecules = file.path(top10_dir, "03_ALL104_IGH_MOLECULE_TOP10_INDEX.tsv"),
  multi = file.path(top10_dir, "06_MULTI_IGH_BIN50_ANNOTATION.tsv"),
  frozen_qc = file.path(top10_dir, "07_TOP10_QC.tsv"),
  top10_column_annotations = file.path(top10_heatmap_dir, "01C_ALL69_COLUMN_ANNOTATIONS.tsv"),
  source_top10_script = "scripts/627_cdr3_per_bin_top10_heatmaps.R"
)

stop_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

read_tsv <- function(path) {
  stop_if(!file.exists(path), paste("Missing input:", path))
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    na.strings = c("NA", ""),
    fileEncoding = "UTF-8"
  )
}

as_flag <- function(x) {
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
}

collapse_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  paste(sort(unique(x)), collapse = "; ")
}

physical_key <- function(sample_id, bin_id) paste(sample_id, bin_id, sep = "||")

all_long <- read_tsv(paths[["all"]])
bin_qc <- read_tsv(paths[["bin_qc"]])
shared_long <- read_tsv(paths[["shared"]])
repeated_long <- read_tsv(paths[["repeated"]])
molecules <- read_tsv(paths[["molecules"]])
multi_bins <- read_tsv(paths[["multi"]])
frozen_qc <- read_tsv(paths[["frozen_qc"]])
base_annotations <- read_tsv(paths[["top10_column_annotations"]])

for (object_name in c("all_long", "shared_long", "repeated_long")) {
  object <- get(object_name)
  required <- c("sample_id", "Bin50_id", "gene_symbol", "raw_count", "within_bin_rank")
  stop_if(!all(required %in% names(object)), paste(object_name, "is missing required columns"))
  object$physical_bin_key <- physical_key(object$sample_id, object$Bin50_id)
  stop_if(any(as.numeric(object$within_bin_rank) > 20), paste(object_name, "contains rank > 20"))
  stop_if(
    any(duplicated(paste(object$physical_bin_key, object$gene_symbol, sep = "||"))),
    paste(object_name, "contains duplicate physical Bin50-gene keys")
  )
  per_bin_n <- table(object$physical_bin_key)
  stop_if(any(per_bin_n > 20L), paste(object_name, "contains >20 genes in a Bin50"))
  assign(object_name, object)
}

expected <- data.frame(
  analysis = c("All69", "Shared7", "Repeated11"),
  input_rows = c(1183L, 740L, 587L),
  bins = c(69L, 44L, 33L),
  full_genes = c(604L, 327L, 226L),
  recurrent_genes = c(93L, 55L, 52L),
  stringsAsFactors = FALSE
)

audit_input <- function(data, label) {
  gene_frequency <- table(data$gene_symbol)
  data.frame(
    analysis = label,
    input_rows = nrow(data),
    bins = length(unique(data$physical_bin_key)),
    full_genes = length(gene_frequency),
    recurrent_genes = sum(gene_frequency >= 2L),
    duplicate_bin_gene = sum(duplicated(paste(
      data$physical_bin_key,
      data$gene_symbol,
      sep = "||"
    ))),
    maximum_genes_per_bin = max(table(data$physical_bin_key)),
    rank_above_20 = sum(as.numeric(data$within_bin_rank) > 20),
    stringsAsFactors = FALSE
  )
}

input_audit <- rbind(
  audit_input(all_long, "All69"),
  audit_input(shared_long, "Shared7"),
  audit_input(repeated_long, "Repeated11")
)
input_audit <- merge(input_audit, expected, by = "analysis", suffixes = c("", "_expected"), sort = FALSE)
input_audit$passed <- with(
  input_audit,
  input_rows == input_rows_expected &
    bins == bins_expected &
    full_genes == full_genes_expected &
    recurrent_genes == recurrent_genes_expected &
    duplicate_bin_gene == 0L &
    maximum_genes_per_bin <= 20L &
    rank_above_20 == 0L
)
stop_if(!all(input_audit$passed), "Top20 input QC failed")
stop_if(nrow(bin_qc) != 69L, "Bin QC does not contain 69 physical Bin50s")
stop_if(nrow(molecules) != 104L, "Frozen molecule index does not contain 104 molecules")
stop_if(length(unique(molecules$corrected_molecule_id)) != 104L, "Molecule IDs are not unique")
stop_if(nrow(multi_bins) != 6L, "Multi-IGH annotation does not contain 6 physical Bin50s")
stop_if(nrow(base_annotations) != 69L, "Frozen All69 column annotation does not contain 69 bins")
stop_if(any(duplicated(base_annotations$physical_bin_key)), "Frozen column annotations duplicate physical bins")

frozen_lookup <- setNames(as.character(frozen_qc$observed), frozen_qc$metric)
stop_if(frozen_lookup[["corrected_molecules"]] != "104", "Frozen corrected-molecule count is not 104")
stop_if(frozen_lookup[["patient_clonotypes"]] != "79", "Frozen patient-clonotype count is not 79")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stop_if(length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L,
        "Output directory exists and is not empty")

log_path <- file.path(output_dir, "RUN_LOG.txt")
log_line <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\t", paste(..., collapse = " "))
  cat(line, "\n", file = log_path, append = TRUE)
  message(line)
}

report_completed_file <- function(path) {
  info <- file.info(path)
  stop_if(is.na(info$size), paste("Completed file is missing:", path))
  log_line(
    "FILE_COMPLETE",
    basename(path),
    paste0("size_bytes=", info$size),
    paste0("completed_at=", format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"))
  )
}
log_line("START", "Top20 heatmaps from frozen per-Bin Top20 TSVs")
log_line("DEFINITION", "Only TopN changes from 10 to 20; no upstream recomputation")

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else NA_character_
if (!is.na(script_path) && file.exists(script_path)) {
  file.copy(
    script_path,
    file.path(output_dir, "CDR3_PER_BIN_TOP20_HEATMAPS.R"),
    overwrite = FALSE
  )
}

run_config <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  top_n = 20,
  value = "log1p(raw_count) for genes already ranked in each Bin50 Top20",
  blank_cell = "NA: not ranked in this Bin's Top20; not necessarily undetected",
  recurrent_rule = "gene enters Top20 of at least 2 distinct physical Bin50s",
  physical_bin_key = "sample_id + Bin50_id",
  column_clustering = FALSE,
  gene_wise_z_score = FALSE,
  source_top10_script = normalizePath(paths[["source_top10_script"]], winslash = "/", mustWork = TRUE),
  frozen_counts = list(corrected_IGH_molecules = 104, patient_specific_clonotypes = 79),
  output_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE)
)
write_json(run_config, file.path(output_dir, "run_config.json"), pretty = TRUE, auto_unbox = TRUE)

base_annotations$IGH_molecule_count <- as.numeric(base_annotations$IGH_molecule_count)
base_annotations$multi_molecule_bin <- as_flag(base_annotations$multi_molecule_bin)
base_annotations$hotspot <- as_flag(base_annotations$hotspot)
base_annotations$direct_RCTD_match <- as_flag(base_annotations$direct_RCTD_match)
base_annotations$context_only <- as_flag(base_annotations$context_only)

disease_levels <- c("HC", "IPF", "SSc-ILD")

make_metadata <- function(data, mode) {
  keys <- unique(data$physical_bin_key)
  rows <- lapply(keys, function(key) {
    z <- data[data$physical_bin_key == key, , drop = FALSE]
    base <- base_annotations[base_annotations$physical_bin_key == key, , drop = FALSE]
    stop_if(nrow(base) != 1L, paste("No unique frozen annotation for", key))
    identity <- if (mode == "shared") {
      collapse_unique(z$CDR3_aa)
    } else if (mode == "repeated") {
      collapse_unique(z$primary_clonotype_id)
    } else {
      collapse_unique(z$CDR3_aa)
    }
    data.frame(
      physical_bin_key = key,
      sample_id = base$sample_id,
      patient_id = base$patient_id,
      disease = base$disease,
      chip = base$chip,
      Bin50_id = base$Bin50_id,
      Bin50_x = as.numeric(base$Bin50_x),
      Bin50_y = as.numeric(base$Bin50_y),
      identity_label = identity,
      direct_RCTD_match = base$direct_RCTD_match,
      context_only = base$context_only,
      RCTD_context = ifelse(base$direct_RCTD_match, "Direct RCTD", "Context only"),
      dominant_lineage = base$dominant_lineage,
      IGH_molecule_count = base$IGH_molecule_count,
      multi_molecule_bin = base$multi_molecule_bin,
      hotspot = base$hotspot,
      stringsAsFactors = FALSE
    )
  })
  metadata <- do.call(rbind, rows)
  if (mode == "all") {
    ordering <- order(
      match(metadata$disease, disease_levels),
      metadata$sample_id,
      metadata$Bin50_x,
      metadata$Bin50_y
    )
    metadata$group_label <- metadata$disease
  } else if (mode == "shared") {
    ordering <- order(
      metadata$identity_label,
      metadata$Bin50_x,
      metadata$Bin50_y,
      metadata$sample_id
    )
    metadata$group_label <- metadata$identity_label
  } else {
    ordering <- order(
      metadata$sample_id,
      metadata$Bin50_x,
      metadata$Bin50_y,
      metadata$identity_label
    )
    metadata$group_label <- metadata$sample_id
  }
  metadata <- metadata[ordering, , drop = FALSE]
  metadata$heatmap_column_id <- sprintf("BINCOL_%03d", seq_len(nrow(metadata)))
  metadata$coordinate <- paste0(metadata$Bin50_x, ":", metadata$Bin50_y)
  metadata
}

build_bundle <- function(data, label, mode) {
  metadata <- make_metadata(data, mode)
  frequency <- table(data$gene_symbol)
  genes <- names(frequency)
  genes <- genes[order(-as.integer(frequency[genes]), genes)]
  categories <- vapply(genes, function(gene) {
    collapse_unique(data$gene_category_all[data$gene_symbol == gene])
  }, character(1))
  gene_annotation <- data.frame(
    gene_symbol = genes,
    number_of_bins_top20 = as.integer(frequency[genes]),
    gene_category_all = categories,
    recurrent_in_at_least_2_bins = as.integer(frequency[genes]) >= 2L,
    stringsAsFactors = FALSE
  )
  matrix_full <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = nrow(metadata),
    dimnames = list(genes, metadata$heatmap_column_id)
  )
  row_index <- match(data$gene_symbol, genes)
  col_index <- match(data$physical_bin_key, metadata$physical_bin_key)
  matrix_full[cbind(row_index, col_index)] <- log1p(as.numeric(data$raw_count))
  recurrent <- gene_annotation$gene_symbol[gene_annotation$recurrent_in_at_least_2_bins]
  list(
    label = label,
    mode = mode,
    data = data,
    metadata = metadata,
    gene_annotation = gene_annotation,
    full = matrix_full,
    recurrent = matrix_full[recurrent, , drop = FALSE]
  )
}

bundles <- list(
  All69 = build_bundle(all_long, "All69", "all"),
  Shared7 = build_bundle(shared_long, "Shared7", "shared"),
  Repeated11 = build_bundle(repeated_long, "Repeated11", "repeated")
)

category_palette <- c(
  "immunoglobulin" = "#6A3D9A",
  "mitochondrial" = "#666666",
  "ribosomal" = "#BDBDBD",
  "Fibroblast marker" = "#1B9E77",
  "Endothelial marker" = "#E78AC3",
  "Epithelial marker" = "#A6761D",
  "other" = "#FEE08B"
)
disease_palette <- c("HC" = "#4DAF4A", "IPF" = "#E41A1C", "SSc-ILD" = "#377EB8")
sample_values <- sort(unique(all_long$sample_id))
sample_palette <- setNames(grDevices::hcl.colors(length(sample_values), "Dynamic"), sample_values)
lineage_values <- sort(unique(na.omit(base_annotations$dominant_lineage)))
lineage_palette <- setNames(grDevices::hcl.colors(length(lineage_values), "Dark 3"), lineage_values)
rctd_palette <- c("Direct RCTD" = "#2166AC", "Context only" = "#BDBDBD")
binary_palette <- c("FALSE" = "#F2F2F2", "TRUE" = "#7F0000")
hotspot_palette <- c("FALSE" = "#F2F2F2", "TRUE" = "#000000")

make_heatmap <- function(bundle, matrix_value, title, png_path, pdf_path, full = FALSE) {
  metadata <- bundle$metadata
  gene_info <- bundle$gene_annotation[
    match(rownames(matrix_value), bundle$gene_annotation$gene_symbol),
    ,
    drop = FALSE
  ]
  categories <- gene_info$gene_category_all
  categories[is.na(categories) | !categories %in% names(category_palette)] <- "other"
  left_annotation <- rowAnnotation(
    `Gene category` = categories,
    col = list(`Gene category` = category_palette),
    width = unit(4, "mm"),
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(
      title_gp = gpar(fontface = "bold", fontsize = 8),
      labels_gp = gpar(fontsize = 7)
    )
  )
  identity_values <- unique(metadata$identity_label)
  identity_palette <- setNames(grDevices::hcl.colors(length(identity_values), "Set 3"), identity_values)
  top_annotation <- HeatmapAnnotation(
    Disease = metadata$disease,
    Sample = metadata$sample_id,
    Identity = metadata$identity_label,
    `RCTD match` = metadata$RCTD_context,
    `Dominant lineage` = metadata$dominant_lineage,
    `IGH molecules` = anno_barplot(
      metadata$IGH_molecule_count,
      gp = gpar(fill = "#5E3C99", col = NA),
      height = unit(10, "mm"),
      axis_param = list(gp = gpar(fontsize = 6))
    ),
    `Multi-IGH Bin` = as.character(metadata$multi_molecule_bin),
    Hotspot = as.character(metadata$hotspot),
    col = list(
      Disease = disease_palette,
      Sample = sample_palette,
      Identity = identity_palette,
      `RCTD match` = rctd_palette,
      `Dominant lineage` = lineage_palette,
      `Multi-IGH Bin` = binary_palette,
      Hotspot = hotspot_palette
    ),
    show_legend = c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
    simple_anno_size = unit(3, "mm"),
    annotation_name_gp = gpar(fontsize = 8)
  )
  values <- matrix_value[is.finite(matrix_value)]
  color_function <- colorRamp2(
    c(0, unname(quantile(values, 0.65)), max(values)),
    c("#FFF7EC", "#FDBB84", "#B30000")
  )
  split_factor <- factor(metadata$group_label, levels = unique(metadata$group_label))
  heatmap <- Heatmap(
    matrix_value,
    name = "log1p(raw count)",
    col = color_function,
    na_col = "#F2F2F2",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = split_factor,
    cluster_column_slices = FALSE,
    show_column_names = TRUE,
    column_labels = metadata$coordinate,
    column_names_rot = 90,
    column_names_gp = gpar(fontsize = if (full) 3.5 else 5),
    row_names_gp = gpar(fontsize = if (full) 4 else 6.5),
    column_title_gp = gpar(fontsize = if (full) 5 else 7, fontface = "bold"),
    top_annotation = top_annotation,
    left_annotation = left_annotation,
    border = TRUE,
    heatmap_legend_param = list(
      title = "Color: log1p(raw count)\namong per-Bin Top20 genes",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)
    )
  )
  width_inches <- max(14, 6 + ncol(matrix_value) * 0.23)
  height_inches <- if (full) {
    max(12, 6.5 + nrow(matrix_value) * 0.075)
  } else {
    max(9, 6.5 + nrow(matrix_value) * 0.12)
  }
  draw_plot <- function() {
    grid.newpage()
    draw(
      heatmap,
      newpage = FALSE,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = TRUE,
      column_title = title,
      column_title_gp = gpar(fontsize = 13, fontface = "bold"),
      padding = unit(c(8, 2, 14, 2), "mm")
    )
    grid.text(
      "Color: log1p(raw count) among per-Bin Top20 genes",
      x = unit(0.01, "npc"), y = unit(0.025, "npc"),
      just = c("left", "bottom"), gp = gpar(fontsize = 8, fontface = "bold")
    )
    grid.text(
      "Blank: not ranked in this Bin's Top20; not necessarily undetected.",
      x = unit(0.01, "npc"), y = unit(0.007, "npc"),
      just = c("left", "bottom"), gp = gpar(fontsize = 8)
    )
  }
  png(
    png_path,
    width = width_inches,
    height = height_inches,
    units = "in",
    res = 300,
    type = "cairo-png",
    bg = "white"
  )
  draw_plot()
  dev.off()
  png_image <- image_read(png_path)
  image_write(png_image, path = png_path, format = "png", density = "300x300")
  rm(png_image)
  report_completed_file(png_path)
  pdf(pdf_path, width = width_inches, height = height_inches, onefile = TRUE, useDingbats = FALSE)
  draw_plot()
  dev.off()
  report_completed_file(pdf_path)
}

figure_specs <- list(
  list(bundle = "All69", matrix = "recurrent", stem = "FIG1_ALL69_RECURRENT_TOP20_HEATMAP", full = FALSE,
       title = "All CDR3-positive physical Bin50s: recurrent per-Bin Top20 genes"),
  list(bundle = "Shared7", matrix = "recurrent", stem = "FIG2_SHARED7_RECURRENT_TOP20_HEATMAP", full = FALSE,
       title = "Cross-patient shared CDR3 occurrences: recurrent per-Bin Top20 genes"),
  list(bundle = "Repeated11", matrix = "recurrent", stem = "FIG3_REPEATED11_RECURRENT_TOP20_HEATMAP", full = FALSE,
       title = "Within-patient repeated clonotypes: recurrent per-Bin Top20 genes"),
  list(bundle = "All69", matrix = "full", stem = "SUPP_FIG1_ALL69_FULL_TOP20_HEATMAP", full = TRUE,
       title = "All CDR3-positive physical Bin50s: all per-Bin Top20 genes"),
  list(bundle = "Shared7", matrix = "full", stem = "SUPP_FIG2_SHARED7_FULL_TOP20_HEATMAP", full = TRUE,
       title = "Cross-patient shared CDR3 occurrences: all per-Bin Top20 genes"),
  list(bundle = "Repeated11", matrix = "full", stem = "SUPP_FIG3_REPEATED11_FULL_TOP20_HEATMAP", full = TRUE,
       title = "Within-patient repeated clonotypes: all per-Bin Top20 genes")
)

for (spec in figure_specs) {
  matrix_value <- bundles[[spec$bundle]][[spec$matrix]]
  log_line("PLOT", spec$stem, paste(dim(matrix_value), collapse = "x"))
  make_heatmap(
    bundles[[spec$bundle]],
    matrix_value,
    spec$title,
    file.path(output_dir, paste0(spec$stem, ".png")),
    file.path(output_dir, paste0(spec$stem, ".pdf")),
    spec$full
  )
}

validate_png <- function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  signature <- readBin(connection, "raw", n = 8L)
  identical(signature, as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))) && file.info(path)$size > 10000
}

validate_pdf <- function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  header <- readBin(connection, "raw", n = 5L)
  size <- file.info(path)$size
  seek(connection, where = max(0, size - 2048), origin = "start")
  tail_raw <- readBin(connection, "raw", n = 2048L)
  eof_raw <- charToRaw("%%EOF")
  starts <- seq_len(length(tail_raw) - length(eof_raw) + 1L)
  has_eof <- any(vapply(starts, function(index) {
    identical(tail_raw[index:(index + length(eof_raw) - 1L)], eof_raw)
  }, logical(1)))
  identical(header, charToRaw("%PDF-")) && has_eof && size > 10000
}

qc_records <- list()
add_qc <- function(analysis, check, observed, expected_value, passed, notes) {
  qc_records[[length(qc_records) + 1L]] <<- data.frame(
    analysis = analysis,
    check = check,
    observed = as.character(observed),
    expected = as.character(expected_value),
    passed = passed,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

for (label in names(bundles)) {
  bundle <- bundles[[label]]
  exp_row <- expected[expected$analysis == label, ]
  checks <- list(
    input_rows = c(nrow(bundle$data), exp_row$input_rows),
    unique_physical_bins = c(ncol(bundle$full), exp_row$bins),
    duplicate_heatmap_columns = c(sum(duplicated(bundle$metadata$physical_bin_key)), 0),
    maximum_genes_per_bin = c(max(table(bundle$data$physical_bin_key)), 20),
    duplicate_bin_gene_keys = c(sum(duplicated(paste(bundle$data$physical_bin_key, bundle$data$gene_symbol))), 0),
    full_gene_count = c(nrow(bundle$full), exp_row$full_genes),
    recurrent_gene_count = c(nrow(bundle$recurrent), exp_row$recurrent_genes),
    full_matrix_non_NA = c(sum(!is.na(bundle$full)), exp_row$input_rows)
  )
  for (name in names(checks)) {
    observed <- checks[[name]][[1L]]
    expected_value <- checks[[name]][[2L]]
    passed <- if (name == "maximum_genes_per_bin") observed <= expected_value else observed == expected_value
    add_qc(label, name, observed, expected_value, passed, "physical key: sample_id + Bin50_id")
  }
}

png_paths <- file.path(output_dir, paste0(vapply(figure_specs, `[[`, character(1), "stem"), ".png"))
pdf_paths <- file.path(output_dir, paste0(vapply(figure_specs, `[[`, character(1), "stem"), ".pdf"))
add_qc("Frozen", "corrected_IGH_molecules", nrow(molecules), 104, nrow(molecules) == 104L, "unchanged frozen count")
add_qc("Frozen", "patient_specific_clonotypes", frozen_lookup[["patient_clonotypes"]], 79,
       frozen_lookup[["patient_clonotypes"]] == "79", "read from frozen QC; not recomputed")
add_qc("Frozen", "multi_IGH_Bin50", nrow(multi_bins), 6, nrow(multi_bins) == 6L, "annotation only")
add_qc("Output", "PNG_files_valid", sum(vapply(png_paths, validate_png, logical(1))), 6, all(vapply(png_paths, validate_png, logical(1))), "signature and non-empty")
add_qc("Output", "PDF_files_valid", sum(vapply(pdf_paths, validate_pdf, logical(1))), 6, all(vapply(pdf_paths, validate_pdf, logical(1))), "PDF header, EOF and non-empty")
qc <- do.call(rbind, qc_records)
write.table(qc, file.path(output_dir, "TOP20_HEATMAP_QC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
stop_if(!all(qc$passed), "Final Top20 heatmap QC failed")

summary_lines <- c(
  "TOP20 HEATMAP QC SUMMARY",
  paste("Generated:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "Definition: each heatmap column is one unique sample_id + Bin50_id physical Bin50.",
  "Value: log1p(raw_count) only for genes already present in that Bin50's frozen Top20 list.",
  "Blank cells are NA: not ranked in this Bin's Top20; not necessarily undetected.",
  sprintf("All69: %d input rows; %d bins; Full %d genes; Recurrent %d genes.",
          nrow(all_long), ncol(bundles$All69$full), nrow(bundles$All69$full), nrow(bundles$All69$recurrent)),
  sprintf("Shared7: %d input rows; %d bins; Full %d genes; Recurrent %d genes.",
          nrow(shared_long), ncol(bundles$Shared7$full), nrow(bundles$Shared7$full), nrow(bundles$Shared7$recurrent)),
  sprintf("Repeated11: %d input rows; %d bins; Full %d genes; Recurrent %d genes.",
          nrow(repeated_long), ncol(bundles$Repeated11$full), nrow(bundles$Repeated11$full), nrow(bundles$Repeated11$recurrent)),
  "All bins contain at most 20 distinct gene symbols; duplicate physical Bin50-gene keys: 0.",
  "Frozen counts retained: 104 corrected IGH molecules and 79 patient-specific clonotypes.",
  "Compared with the Top10 heatmaps, only TopN changed from 10 to 20; all other definitions are unchanged.",
  "QC status: PASS"
)
writeLines(summary_lines, file.path(output_dir, "TOP20_HEATMAP_QC_SUMMARY.txt"), useBytes = TRUE)

log_line("QC", "All69 Full/Recurrent", nrow(bundles$All69$full), nrow(bundles$All69$recurrent))
log_line("QC", "Shared7 Full/Recurrent", nrow(bundles$Shared7$full), nrow(bundles$Shared7$recurrent))
log_line("QC", "Repeated11 Full/Recurrent", nrow(bundles$Repeated11$full), nrow(bundles$Repeated11$recurrent))
log_line("COMPLETE", normalizePath(output_dir, winslash = "/", mustWork = TRUE))

output_files <- list.files(output_dir, full.names = TRUE, recursive = FALSE)
output_files <- output_files[basename(output_files) != "MANIFEST.tsv"]
manifest <- rbind(
  data.frame(
    role = "input",
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    size_bytes = file.info(paths)$size,
    modified_time = format(file.info(paths)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    sha256 = vapply(paths, digest, character(1), file = TRUE, algo = "sha256"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    role = "output",
    path = normalizePath(output_files, winslash = "/", mustWork = TRUE),
    size_bytes = file.info(output_files)$size,
    modified_time = format(file.info(output_files)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    sha256 = vapply(output_files, digest, character(1), file = TRUE, algo = "sha256"),
    stringsAsFactors = FALSE
  )
)
write.table(manifest, file.path(output_dir, "MANIFEST.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("OUTPUT_DIR\t", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
for (spec in figure_specs) {
  cat("FIGURE\t", normalizePath(file.path(output_dir, paste0(spec$stem, ".png")), winslash = "/", mustWork = TRUE), "\n", sep = "")
  cat("FIGURE\t", normalizePath(file.path(output_dir, paste0(spec$stem, ".pdf")), winslash = "/", mustWork = TRUE), "\n", sep = "")
}
for (label in names(bundles)) {
  cat(
    "QC\t", label,
    "\tbins=", ncol(bundles[[label]]$full),
    "\tfull_genes=", nrow(bundles[[label]]$full),
    "\trecurrent_genes=", nrow(bundles[[label]]$recurrent),
    "\tnon_NA=", sum(!is.na(bundles[[label]]$full)),
    "\n",
    sep = ""
  )
}
cat("QC\t104_molecules=PASS\t79_patient_clonotypes=PASS\tdefinitions_except_TopN_unchanged=TRUE\n")
