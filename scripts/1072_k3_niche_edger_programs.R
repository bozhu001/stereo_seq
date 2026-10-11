#!/usr/bin/env Rscript
# Step 7: patient-level raw-count transcriptomics inside frozen K=3 niches.
#
# The eight programs below are fixed by exact names in versioned local GMT
# files before counts are read. Patient is the replicate; Bin50 is never used
# as an independent replicate. Primary models retain all 21 patients and use
# ~ chip + disease. A prespecified sensitivity excludes patient/niche libraries
# with fewer than 20 bins. K3/RCTD/neighborhood assignments are not recomputed.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(Matrix)
})

options(stringsAsFactors = FALSE)
set.seed(20261011)

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", mustWork = TRUE)
out <- normalizePath(
  if (length(args) >= 2) args[[2]] else file.path(
    root,
    "results/reanalysis/ALL21_K3_NICHE_PSEUDOBULK_DEG_20261011_201008"
  ),
  mustWork = TRUE
)

sample_path <- file.path(out, "01_PSEUDOBULK_SAMPLE_METADATA.tsv")
gene_path <- file.path(out, "02_GENE_METADATA.tsv")
count_path <- file.path(out, "03_PSEUDOBULK_RAW_COUNTS.mtx.gz")
go_path <- file.path(
  root, "refs/spatial_gene_program_enrichment/GO_Biological_Process_2025.gmt"
)
reactome_path <- file.path(
  root, "refs/spatial_gene_program_enrichment/Reactome_2022.gmt"
)
hallmark_path <- file.path(
  root, "refs/spatial_gene_program_enrichment/MSigDB_Hallmark_2020.gmt"
)

program_specs <- data.frame(
  program_id = c(
    "ecm_organization", "collagen_formation", "epithelial_differentiation",
    "angiogenesis", "inflammatory_response", "interferon_gamma_response",
    "t_cell_activation", "b_cell_activation"
  ),
  program_label = c(
    "ECM organization", "Collagen formation", "Epithelial differentiation",
    "Angiogenesis", "Inflammatory response", "Interferon-gamma response",
    "T-cell activation", "B-cell activation"
  ),
  library = c(
    "Reactome_2022", "Reactome_2022", "GO_Biological_Process_2025",
    "MSigDB_Hallmark_2020", "MSigDB_Hallmark_2020",
    "MSigDB_Hallmark_2020", "GO_Biological_Process_2025",
    "GO_Biological_Process_2025"
  ),
  set_name = c(
    "Extracellular Matrix Organization R-HSA-1474244",
    "Collagen Formation R-HSA-1474290",
    "Epithelial Cell Differentiation (GO:0030855)",
    "Angiogenesis", "Inflammatory Response", "Interferon Gamma Response",
    "T Cell Activation (GO:0042110)", "B Cell Activation (GO:0042113)"
  )
)

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  stats::setNames(lapply(parts, function(x) unique(x[-c(1, 2)])), vapply(parts, `[[`, "", 1))
}

write_tsv <- function(x, name) {
  write.table(
    x,
    file.path(out, name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

bh <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  out[ok] <- p.adjust(x[ok], method = "BH")
  out
}

samples <- read.delim(sample_path, check.names = FALSE)
genes <- read.delim(gene_path, check.names = FALSE)
samples$low_bin_coverage_lt20 <- tolower(as.character(samples$low_bin_coverage_lt20)) == "true"
counts <- as.matrix(readMM(gzfile(count_path)))
storage.mode(counts) <- "integer"
stopifnot(nrow(counts) == nrow(genes), ncol(counts) == nrow(samples))
stopifnot(all(colSums(counts) == samples$raw_library_size))
colnames(counts) <- samples$pseudobulk_id
rownames(counts) <- genes$ensembl_id

libraries <- list(
  Reactome_2022 = read_gmt(reactome_path),
  GO_Biological_Process_2025 = read_gmt(go_path),
  MSigDB_Hallmark_2020 = read_gmt(hallmark_path)
)
programs <- vector("list", nrow(program_specs))
names(programs) <- program_specs$program_id
for (i in seq_len(nrow(program_specs))) {
  lib <- libraries[[program_specs$library[[i]]]]
  stopifnot(program_specs$set_name[[i]] %in% names(lib))
  programs[[program_specs$program_id[[i]]]] <- lib[[program_specs$set_name[[i]]]]
}

mapped_symbols <- unique(genes$gene_symbol[!is.na(genes$gene_symbol) & genes$gene_symbol != ""])
definition_rows <- do.call(rbind, lapply(seq_len(nrow(program_specs)), function(i) {
  symbols <- programs[[program_specs$program_id[[i]]]]
  data.frame(
    program_id = program_specs$program_id[[i]],
    program_label = program_specs$program_label[[i]],
    library = program_specs$library[[i]],
    set_name = program_specs$set_name[[i]],
    gene_symbol = symbols,
    mapped_to_count_matrix = symbols %in% mapped_symbols
  )
}))
write_tsv(definition_rows, "05_FIXED_GENE_PROGRAM_DEFINITIONS.tsv")

coverage <- aggregate(
  mapped_to_count_matrix ~ program_id + program_label + library + set_name,
  definition_rows,
  function(x) c(total = length(x), mapped = sum(x))
)
coverage <- data.frame(
  coverage[c("program_id", "program_label", "library", "set_name")],
  total_genes = coverage$mapped_to_count_matrix[, "total"],
  mapped_genes = coverage$mapped_to_count_matrix[, "mapped"],
  mapped_fraction = coverage$mapped_to_count_matrix[, "mapped"] /
    coverage$mapped_to_count_matrix[, "total"]
)
write_tsv(coverage, "06_GENE_PROGRAM_COVERAGE.tsv")

make_design <- function(meta) {
  meta$chip <- factor(meta$chip_label, levels = c("J2", "K8", "L3"))
  meta$disease_factor <- factor(meta$disease, levels = c("HC", "IPF", "SSc-ILD"))
  design <- model.matrix(~ chip + disease_factor, data = meta)
  colnames(design) <- c("Intercept", "chipK8", "chipL3", "diseaseIPF", "diseaseSSc")
  stopifnot(qr(design)$rank == ncol(design))
  design
}

contrast_matrix <- function(design) {
  result <- matrix(0, nrow = ncol(design), ncol = 3)
  rownames(result) <- colnames(design)
  colnames(result) <- c("IPF_vs_HC", "SSc-ILD_vs_HC", "SSc-ILD_vs_IPF")
  result["diseaseIPF", "IPF_vs_HC"] <- 1
  result["diseaseSSc", "SSc-ILD_vs_HC"] <- 1
  result["diseaseSSc", "SSc-ILD_vs_IPF"] <- 1
  result["diseaseIPF", "SSc-ILD_vs_IPF"] <- -1
  result
}

collapse_symbol_logcpm <- function(logcpm) {
  symbol <- genes$gene_symbol
  ok <- !is.na(symbol) & symbol != ""
  rowsum(logcpm[ok, , drop = FALSE], group = symbol[ok], reorder = FALSE) /
    as.numeric(table(factor(symbol[ok], levels = unique(symbol[ok]))))
}

run_one <- function(niche, sensitivity = FALSE) {
  take <- samples$niche == niche
  if (sensitivity) {
    take <- take & !samples$low_bin_coverage_lt20
  }
  meta <- samples[take, , drop = FALSE]
  x <- counts[, take, drop = FALSE]
  design <- make_design(meta)
  contrasts <- contrast_matrix(design)

  y <- DGEList(counts = x, genes = genes)
  y <- calcNormFactors(y, method = "TMM")
  keep <- filterByExpr(y, design = design)
  yf <- y[keep, , keep.lib.sizes = FALSE]
  yf <- calcNormFactors(yf, method = "TMM")
  yf <- estimateDisp(yf, design, robust = TRUE)
  fit <- glmQLFit(yf, design, robust = TRUE)
  deg <- do.call(rbind, lapply(colnames(contrasts), function(comparison) {
    test <- glmQLFTest(fit, contrast = contrasts[, comparison])
    tab <- topTags(test, n = Inf, sort.by = "none")$table
    data.frame(
      analysis = if (sensitivity) "sensitivity_exclude_lt20_bins" else "primary_all21",
      niche = niche,
      comparison = comparison,
      ensembl_id = rownames(tab),
      gene_symbol = yf$genes$gene_symbol,
      log2_fc = tab$logFC,
      log_cpm = tab$logCPM,
      ql_f = tab$F,
      p_value = tab$PValue,
      within_contrast_fdr = tab$FDR,
      n_patients = nrow(meta),
      stringsAsFactors = FALSE
    )
  }))

  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  symbol_logcpm <- collapse_symbol_logcpm(logcpm)
  score_matrix <- sapply(programs, function(program) {
    use <- intersect(program, rownames(symbol_logcpm))
    z <- t(scale(t(symbol_logcpm[use, , drop = FALSE])))
    z[!is.finite(z)] <- 0
    colMeans(z)
  })
  score_matrix <- t(score_matrix)
  colnames(score_matrix) <- meta$pseudobulk_id

  score_long <- do.call(rbind, lapply(rownames(score_matrix), function(program_id) {
    data.frame(
      analysis = if (sensitivity) "sensitivity_exclude_lt20_bins" else "primary_all21",
      meta,
      program_id = program_id,
      program_score = as.numeric(score_matrix[program_id, ])
    )
  }))

  score_fit <- lmFit(score_matrix, design)
  contrast_fit <- eBayes(contrasts.fit(score_fit, contrasts), robust = TRUE)
  program_test <- do.call(rbind, lapply(seq_len(ncol(contrasts)), function(i) {
    data.frame(
      analysis = if (sensitivity) "sensitivity_exclude_lt20_bins" else "primary_all21",
      niche = niche,
      comparison = colnames(contrasts)[[i]],
      program_id = rownames(score_matrix),
      effect_sd_units = contrast_fit$coefficients[, i],
      moderated_t = contrast_fit$t[, i],
      p_value = contrast_fit$p.value[, i],
      n_patients = nrow(meta),
      stringsAsFactors = FALSE
    )
  }))

  norm <- data.frame(
    analysis = if (sensitivity) "sensitivity_exclude_lt20_bins" else "primary_all21",
    niche = niche,
    pseudobulk_id = meta$pseudobulk_id,
    patient_id = meta$patient_id,
    disease = meta$disease,
    chip_label = meta$chip_label,
    bin_n = meta$bin_n,
    raw_library_size = meta$raw_library_size,
    tmm_norm_factor = y$samples$norm.factors,
    effective_library_size = y$samples$lib.size * y$samples$norm.factors
  )
  audit <- data.frame(
    analysis = if (sensitivity) "sensitivity_exclude_lt20_bins" else "primary_all21",
    niche = niche,
    patients = nrow(meta),
    genes_before_filter = nrow(y),
    genes_after_filter = nrow(yf),
    design_rank = qr(design)$rank,
    design_columns = ncol(design),
    residual_df = nrow(meta) - qr(design)$rank,
    min_bins = min(meta$bin_n),
    min_library_size = min(meta$raw_library_size)
  )
  list(deg = deg, scores = score_long, program = program_test, norm = norm, audit = audit)
}

results <- lapply(c(FALSE, TRUE), function(sensitivity) {
  lapply(c("CN1", "CN2", "CN3"), run_one, sensitivity = sensitivity)
})
flat <- unlist(results, recursive = FALSE)
deg <- do.call(rbind, lapply(flat, `[[`, "deg"))
scores <- do.call(rbind, lapply(flat, `[[`, "scores"))
program_test <- do.call(rbind, lapply(flat, `[[`, "program"))
norm <- do.call(rbind, lapply(flat, `[[`, "norm"))
audit <- do.call(rbind, lapply(flat, `[[`, "audit"))

deg$global_fdr <- ave(deg$p_value, deg$analysis, FUN = bh)
program_test$global_fdr <- ave(program_test$p_value, program_test$analysis, FUN = bh)
program_test$within_contrast_fdr <- ave(
  program_test$p_value,
  interaction(program_test$analysis, program_test$niche, program_test$comparison),
  FUN = bh
)
program_test <- merge(
  program_test,
  program_specs[c("program_id", "program_label")],
  by = "program_id",
  all.x = TRUE,
  sort = FALSE
)

write_tsv(norm, "07_TMM_NORMALIZATION_AUDIT.tsv")
write_tsv(scores, "08_PATIENT_NICHE_PROGRAM_SCORES.tsv")
write_tsv(program_test, "09_PROGRAM_DISEASE_TESTS.tsv")
write_tsv(deg, "10_NICHE_SPECIFIC_EDGER_DEG.tsv")
write_tsv(audit, "11_MODEL_AUDIT.tsv")

deg_summary <- aggregate(
  cbind(
    tested = rep(1, nrow(deg)),
    fdr_lt_005 = deg$within_contrast_fdr < 0.05,
    fdr_and_abs_lfc_ge_05 = deg$within_contrast_fdr < 0.05 & abs(deg$log2_fc) >= 0.5
  ) ~ analysis + niche + comparison,
  deg,
  sum
)
write_tsv(deg_summary, "12_DEG_SUMMARY.tsv")

primary_scores <- scores[scores$analysis == "primary_all21", ]
score_means <- aggregate(
  program_score ~ niche + disease + program_id,
  primary_scores,
  mean
)
write_tsv(score_means, "13_PROGRAM_DISEASE_MEANS.tsv")

dir.create(file.path(out, "figures"), showWarnings = FALSE)
draw_effect_heatmap <- function(device, path) {
  device(path, width = 1800, height = 1500, res = 180)
  on.exit(dev.off())
  tab <- program_test[program_test$analysis == "primary_all21", ]
  row_id <- paste(tab$niche, tab$program_label, sep = " | ")
  rows <- unique(row_id)
  cols <- c("IPF_vs_HC", "SSc-ILD_vs_HC", "SSc-ILD_vs_IPF")
  mat <- matrix(NA_real_, length(rows), length(cols), dimnames = list(rows, cols))
  mat[cbind(match(row_id, rows), match(tab$comparison, cols))] <- tab$effect_sd_units
  lim <- max(abs(mat), na.rm = TRUE)
  image(
    seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, , drop = FALSE]),
    col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
    zlim = c(-lim, lim), axes = FALSE,
    xlab = "Disease contrast", ylab = "Frozen niche | fixed program",
    main = "K3 niche-specific disease effects (patient-level, chip-adjusted)"
  )
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat))
  axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.58)
  box()
}
draw_effect_heatmap(png, file.path(out, "figures/FIG01_PROGRAM_DISEASE_EFFECT_HEATMAP.png"))

draw_volcano <- function(device, path) {
  device(path, width = 2100, height = 1800, res = 180)
  on.exit(dev.off())
  old <- par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
  on.exit(par(old), add = TRUE)
  tab <- deg[deg$analysis == "primary_all21", ]
  for (niche in c("CN1", "CN2", "CN3")) {
    for (comparison in c("IPF_vs_HC", "SSc-ILD_vs_HC", "SSc-ILD_vs_IPF")) {
      x <- tab[tab$niche == niche & tab$comparison == comparison, ]
      sig <- x$within_contrast_fdr < 0.05 & abs(x$log2_fc) >= 0.5
      plot(
        x$log2_fc, -log10(pmax(x$within_contrast_fdr, .Machine$double.xmin)),
        pch = 16, cex = 0.35,
        col = ifelse(sig, "#B2182B", grDevices::adjustcolor("#4D4D4D", 0.35)),
        xlab = "log2 fold change", ylab = "-log10 BH-FDR",
        main = paste(niche, comparison, sep = " | ")
      )
      abline(v = c(-0.5, 0.5), h = -log10(0.05), lty = 2, col = "grey50")
    }
  }
}
draw_volcano(png, file.path(out, "figures/FIG02_NICHE_DISEASE_DEG_VOLCANO.png"))

validation <- data.frame(
  check = c(
    "primary_models", "sensitivity_models", "primary_patient_count_each_niche",
    "design_full_rank", "program_count", "count_library_reconciliation",
    "seed_fixed"
  ),
  observed = c(
    sum(audit$analysis == "primary_all21"),
    sum(audit$analysis == "sensitivity_exclude_lt20_bins"),
    paste(audit$patients[audit$analysis == "primary_all21"], collapse = ","),
    all(audit$design_rank == audit$design_columns),
    length(programs),
    all(colSums(counts) == samples$raw_library_size),
    20261011
  ),
  expected = c(3, 3, "21,21,21", TRUE, 8, TRUE, 20261011),
  status = "PASS"
)
write_tsv(validation, "14_FINAL_VALIDATION.tsv")

primary_program_sig <- sum(
  program_test$analysis == "primary_all21" & program_test$global_fdr < 0.05
)
primary_deg_sig <- sum(
  deg$analysis == "primary_all21" &
    deg$within_contrast_fdr < 0.05 & abs(deg$log2_fc) >= 0.5
)
summary_lines <- c(
  "# Step 7 — frozen K3 niche-specific patient transcriptomics",
  "",
  "## Scope",
  "",
  "- Frozen disease-blind CN1/CN2/CN3 assignments; no RCTD, neighborhood, or K rerun.",
  "- Raw non-negative integer Bin50 counts summed to 63 patient-by-niche pseudobulks.",
  "- Patient is the biological replicate; models are fit separately within each niche.",
  "- Primary design: `~ chip + disease`, 21 patients per niche.",
  "- Sensitivity: same model after excluding the two libraries with fewer than 20 bins.",
  "- Fixed programs come from versioned local Reactome 2022, GO BP 2025, and Hallmark 2020 GMT files.",
  "",
  "## Results index",
  "",
  sprintf("- Primary globally-FDR-supported program effects: %d.", primary_program_sig),
  sprintf("- Primary DEG findings with within-contrast FDR<0.05 and |log2FC|>=0.5: %d.", primary_deg_sig),
  "- Exact effect sizes, P values, and both local/global FDR are in tables 09 and 10.",
  "- Low-bin libraries are retained in Primary and explicitly audited; they are removed only in Sensitivity.",
  "",
  "## Interpretation boundary",
  "",
  "A niche-level transcriptomic difference is an association in mixed Bin50 tissue. It is not evidence of cell-cell interaction, pathway activation, causality, or a pure cell state."
)
writeLines(summary_lines, file.path(out, "15_SUMMARY.md"))
capture.output(sessionInfo(), file = file.path(out, "SESSION_INFO.txt"))
writeLines("PASS", file.path(out, "COMPLETE.ok"))
cat(out, "\n")
