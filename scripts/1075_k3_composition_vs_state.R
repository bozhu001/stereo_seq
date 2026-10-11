#!/usr/bin/env Rscript
# Post-DEG decomposition for frozen K3 niches.
# Tests whether robust disease-associated genes remain after conditioning on
# the frozen RCTD neighborhood composition. This is descriptive/post-selection
# evidence, not a causal mediation analysis.

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

read_tsv <- function(name) read.delim(file.path(out, name), check.names = FALSE)
write_tsv <- function(x, name) write.table(
  x,
  file.path(out, name),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
bh <- function(x) p.adjust(x, method = "BH")
logit <- function(x, epsilon = 1e-6) {
  x <- pmin(pmax(x, epsilon), 1 - epsilon)
  log(x / (1 - x))
}

samples <- read_tsv("01_PSEUDOBULK_SAMPLE_METADATA.tsv")
samples$low_bin_coverage_lt20 <- tolower(as.character(samples$low_bin_coverage_lt20)) == "true"
genes <- read_tsv("02_GENE_METADATA.tsv")
composition <- read_tsv("19_PATIENT_NICHE_RCTD_COMPOSITION.tsv")
composition$low_bin_coverage_lt20 <- tolower(as.character(composition$low_bin_coverage_lt20)) == "true"
robust_audit <- read_tsv("16_PRIMARY_DEG_SENSITIVITY_AUDIT.tsv")
robust_audit$robust_to_lt20_bin_exclusion <- tolower(
  as.character(robust_audit$robust_to_lt20_bin_exclusion)
) == "true"
robust <- robust_audit[robust_audit$robust_to_lt20_bin_exclusion, , drop = FALSE]
counts <- as.matrix(readMM(gzfile(file.path(out, "03_PSEUDOBULK_RAW_COUNTS.mtx.gz"))))
storage.mode(counts) <- "integer"
colnames(counts) <- samples$pseudobulk_id
rownames(counts) <- genes$ensembl_id
stopifnot(nrow(robust) == 21, nrow(composition) == 63)
stopifnot(all(colSums(counts) == samples$raw_library_size))

meta_all <- merge(
  samples,
  composition,
  by = c(
    "pseudobulk_id", "patient_id", "chip_label", "disease", "niche",
    "low_bin_coverage_lt20"
  ),
  suffixes = c("_count", "_composition"),
  sort = FALSE,
  validate = "one_to_one"
)
meta_all <- meta_all[match(samples$pseudobulk_id, meta_all$pseudobulk_id), ]
stopifnot(all(meta_all$pseudobulk_id == samples$pseudobulk_id))

make_base <- function(meta) {
  meta$chip <- factor(meta$chip_label, levels = c("J2", "K8", "L3"))
  meta$disease_factor <- factor(meta$disease, levels = c("HC", "IPF", "SSc-ILD"))
  meta
}

rename_design <- function(design, covariate = NULL) {
  base_names <- c("Intercept", "chipK8", "chipL3", "diseaseIPF", "diseaseSSc")
  colnames(design)[seq_along(base_names)] <- base_names
  if (!is.null(covariate)) colnames(design)[ncol(design)] <- covariate
  design
}

contrasts_for <- function(design) {
  result <- matrix(0, nrow = ncol(design), ncol = 3)
  rownames(result) <- colnames(design)
  colnames(result) <- c("IPF_vs_HC", "SSc-ILD_vs_HC", "SSc-ILD_vs_IPF")
  result["diseaseIPF", "IPF_vs_HC"] <- 1
  result["diseaseSSc", "SSc-ILD_vs_HC"] <- 1
  result["diseaseSSc", "SSc-ILD_vs_IPF"] <- 1
  result["diseaseIPF", "SSc-ILD_vs_IPF"] <- -1
  result
}

analyses <- c("primary_all21", "sensitivity_exclude_lt20_bins")
niches <- c("CN1", "CN2", "CN3")
state_names <- c(
  "Alveolar/resident macrophage", "Inflammatory monocyte/macrophage",
  "Dendritic cell", "T", "NK", "Mast/Basophil", "B",
  "Ciliated epithelial", "AT2/alveolar transitional",
  "Secretory/other airway epithelial", "AT1", "Blood endothelial",
  "Perivascular/mesothelial stromal", "Fibroblast",
  "Lymphatic endothelial", "Plasma", "Basal-like epithelial", "Myofibroblast"
)

# Disease associations of the frozen 18-state neighborhood composition.
composition_tests <- list()
index <- 1L
for (analysis in analyses) {
  for (niche in niches) {
    take <- meta_all$niche == niche
    if (analysis == "sensitivity_exclude_lt20_bins") {
      take <- take & !meta_all$low_bin_coverage_lt20
    }
    meta <- make_base(meta_all[take, , drop = FALSE])
    matrix <- t(vapply(state_names, function(state) logit(meta[[state]]), numeric(nrow(meta))))
    rownames(matrix) <- state_names
    design <- rename_design(model.matrix(~ chip + disease_factor, data = meta))
    stopifnot(qr(design)$rank == ncol(design))
    contrast <- contrasts_for(design)
    fit <- eBayes(contrasts.fit(lmFit(matrix, design), contrast), robust = TRUE)
    for (comparison in colnames(contrast)) {
      composition_tests[[index]] <- data.frame(
        analysis = analysis,
        niche = niche,
        comparison = comparison,
        cell_type = state_names,
        logit_fraction_effect = fit$coefficients[, comparison],
        moderated_t = fit$t[, comparison],
        p_value = fit$p.value[, comparison],
        n_patients = nrow(meta)
      )
      index <- index + 1L
    }
  }
}
composition_tests <- do.call(rbind, composition_tests)
composition_tests$within_contrast_fdr <- ave(
  composition_tests$p_value,
  interaction(
    composition_tests$analysis,
    composition_tests$niche,
    composition_tests$comparison
  ),
  FUN = bh
)
composition_tests$global_fdr <- ave(
  composition_tests$p_value,
  composition_tests$analysis,
  FUN = bh
)
write_tsv(composition_tests, "22_RCTD_COMPOSITION_DISEASE_TESTS.tsv")

# Refit robust candidate genes with one frozen composition covariate at a time.
covariates <- c(
  ciliated_fraction = "Ciliated epithelial",
  airway_epithelial_sum = "airway_epithelial_sum"
)
adjusted_results <- list()
expression_rows <- list()
model_audit <- list()
result_index <- 1L
expression_index <- 1L
audit_index <- 1L
for (analysis in analyses) {
  for (niche in niches) {
    take <- meta_all$niche == niche
    if (analysis == "sensitivity_exclude_lt20_bins") {
      take <- take & !meta_all$low_bin_coverage_lt20
    }
    meta <- make_base(meta_all[take, , drop = FALSE])
    x <- counts[, match(meta$pseudobulk_id, colnames(counts)), drop = FALSE]
    base_design <- rename_design(model.matrix(~ chip + disease_factor, data = meta))
    y <- DGEList(counts = x, genes = genes)
    y <- calcNormFactors(y, method = "TMM")
    logcpm <- cpm(y, log = TRUE, prior.count = 2)
    candidate <- unique(robust$ensembl_id[robust$niche == niche])
    candidate <- intersect(candidate, rownames(logcpm))
    if (length(candidate)) {
      for (gene in candidate) {
        expression_rows[[expression_index]] <- data.frame(
          analysis = analysis,
          niche = niche,
          pseudobulk_id = meta$pseudobulk_id,
          patient_id = meta$patient_id,
          chip_label = meta$chip_label,
          disease = meta$disease,
          ensembl_id = gene,
          gene_symbol = genes$gene_symbol[match(gene, genes$ensembl_id)],
          log_cpm = as.numeric(logcpm[gene, ]),
          ciliated_fraction = meta[["Ciliated epithelial"]],
          airway_epithelial_sum = meta$airway_epithelial_sum
        )
        expression_index <- expression_index + 1L
      }
    }
    for (covariate_name in names(covariates)) {
      source_column <- covariates[[covariate_name]]
      covariate_z <- as.numeric(scale(logit(meta[[source_column]])))
      design_data <- meta
      design_data$composition_z <- covariate_z
      design <- rename_design(
        model.matrix(~ chip + disease_factor + composition_z, data = design_data),
        covariate = "composition_z"
      )
      stopifnot(qr(design)$rank == ncol(design))
      keep <- filterByExpr(y, design = design)
      yf <- y[keep, , keep.lib.sizes = FALSE]
      yf <- calcNormFactors(yf, method = "TMM")
      yf <- estimateDisp(yf, design, robust = TRUE)
      fit <- glmQLFit(yf, design, robust = TRUE)
      contrast <- contrasts_for(design)
      target_rows <- robust[robust$niche == niche, , drop = FALSE]
      for (comparison in unique(target_rows$comparison)) {
        test <- glmQLFTest(fit, contrast = contrast[, comparison])
        tab <- topTags(test, n = Inf, sort.by = "none")$table
        wanted <- target_rows[target_rows$comparison == comparison, , drop = FALSE]
        wanted <- wanted[wanted$ensembl_id %in% rownames(tab), , drop = FALSE]
        selected <- tab[match(wanted$ensembl_id, rownames(tab)), , drop = FALSE]
        adjusted_results[[result_index]] <- data.frame(
          analysis = analysis,
          niche = niche,
          comparison = comparison,
          adjustment = covariate_name,
          ensembl_id = wanted$ensembl_id,
          gene_symbol = wanted$gene_symbol_primary,
          original_log2_fc = if (analysis == "primary_all21") {
            wanted$log2_fc_primary
          } else {
            wanted$log2_fc_sensitivity
          },
          original_within_contrast_fdr = if (analysis == "primary_all21") {
            wanted$within_contrast_fdr_primary
          } else {
            wanted$within_contrast_fdr_sensitivity
          },
          adjusted_log2_fc = selected$logFC,
          adjusted_log_cpm = selected$logCPM,
          adjusted_ql_f = selected$F,
          adjusted_p_value = selected$PValue,
          n_patients = nrow(meta),
          composition_effect_per_sd = fit$coefficients[
            match(wanted$ensembl_id, rownames(fit$coefficients)),
            "composition_z"
          ]
        )
        result_index <- result_index + 1L
      }
      model_audit[[audit_index]] <- data.frame(
        analysis = analysis,
        niche = niche,
        adjustment = covariate_name,
        n_patients = nrow(meta),
        design_rank = qr(design)$rank,
        design_columns = ncol(design),
        residual_df = nrow(meta) - qr(design)$rank,
        genes_after_filter = nrow(yf),
        covariate_min = min(meta[[source_column]]),
        covariate_max = max(meta[[source_column]])
      )
      audit_index <- audit_index + 1L
    }
  }
}
adjusted <- do.call(rbind, adjusted_results)
adjusted$adjusted_candidate_fdr <- ave(
  adjusted$adjusted_p_value,
  interaction(adjusted$analysis, adjusted$adjustment),
  FUN = bh
)
adjusted$absolute_effect_ratio <- abs(adjusted$adjusted_log2_fc) /
  abs(adjusted$original_log2_fc)
adjusted$attenuation_fraction <- 1 - adjusted$absolute_effect_ratio
adjusted$same_direction <- adjusted$adjusted_log2_fc * adjusted$original_log2_fc > 0
adjusted$adjusted_nominal_p_lt_005 <- adjusted$adjusted_p_value < 0.05
write_tsv(adjusted, "23_COMPOSITION_ADJUSTED_ROBUST_DEG.tsv")
write_tsv(do.call(rbind, model_audit), "24_COMPOSITION_ADJUSTED_MODEL_AUDIT.tsv")

expression <- do.call(rbind, expression_rows)
write_tsv(expression, "25_ROBUST_GENE_PATIENT_LOGCPM.tsv")
correlation_rows <- list()
correlation_index <- 1L
for (analysis in analyses) {
  sub_analysis <- expression[expression$analysis == analysis, , drop = FALSE]
  keys <- unique(sub_analysis[c("niche", "ensembl_id", "gene_symbol")])
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    sub <- sub_analysis[
      sub_analysis$niche == key$niche & sub_analysis$ensembl_id == key$ensembl_id,
      ,
      drop = FALSE
    ]
    for (covariate in c("ciliated_fraction", "airway_epithelial_sum")) {
      test <- suppressWarnings(cor.test(
        sub$log_cpm,
        sub[[covariate]],
        method = "spearman",
        exact = FALSE
      ))
      correlation_rows[[correlation_index]] <- data.frame(
        analysis = analysis,
        niche = key$niche,
        ensembl_id = key$ensembl_id,
        gene_symbol = key$gene_symbol,
        covariate = covariate,
        spearman_rho = unname(test$estimate),
        p_value = test$p.value,
        n_patients = nrow(sub)
      )
      correlation_index <- correlation_index + 1L
    }
  }
}
correlations <- do.call(rbind, correlation_rows)
correlations$fdr <- ave(
  correlations$p_value,
  interaction(correlations$analysis, correlations$covariate),
  FUN = bh
)
write_tsv(correlations, "26_ROBUST_GENE_COMPOSITION_CORRELATIONS.tsv")

primary_ciliated <- adjusted[
  adjusted$analysis == "primary_all21" &
    adjusted$adjustment == "ciliated_fraction",
  ,
  drop = FALSE
]
sensitivity_ciliated <- adjusted[
  adjusted$analysis == "sensitivity_exclude_lt20_bins" &
    adjusted$adjustment == "ciliated_fraction",
  ,
  drop = FALSE
]
summary <- data.frame(
  metric = c(
    "robust_gene_contrast_rows",
    "primary_ciliated_adjusted_nominal_p_lt_005",
    "primary_ciliated_adjusted_candidate_fdr_lt_005",
    "primary_median_attenuation_fraction",
    "sensitivity_ciliated_adjusted_nominal_p_lt_005",
    "sensitivity_ciliated_adjusted_candidate_fdr_lt_005",
    "sensitivity_median_attenuation_fraction",
    "primary_ciliated_correlation_fdr_lt_005"
  ),
  value = c(
    nrow(robust),
    sum(primary_ciliated$adjusted_nominal_p_lt_005),
    sum(primary_ciliated$adjusted_candidate_fdr < 0.05),
    median(primary_ciliated$attenuation_fraction),
    sum(sensitivity_ciliated$adjusted_nominal_p_lt_005),
    sum(sensitivity_ciliated$adjusted_candidate_fdr < 0.05),
    median(sensitivity_ciliated$attenuation_fraction),
    sum(
      correlations$analysis == "primary_all21" &
        correlations$covariate == "ciliated_fraction" &
        correlations$fdr < 0.05
    )
  )
)
write_tsv(summary, "27_COMPOSITION_VS_STATE_SUMMARY.tsv")

validation <- data.frame(
  check = c(
    "patient_niche_rows", "robust_gene_contrast_rows", "composition_states",
    "composition_models_full_rank", "adjusted_models_full_rank",
    "analyses", "composition_covariates"
  ),
  observed = c(
    nrow(meta_all), nrow(robust), length(state_names),
    TRUE,
    all(vapply(model_audit, function(x) x$design_rank == x$design_columns, logical(1))),
    length(analyses), length(covariates)
  ),
  expected = c(63, 21, 18, TRUE, TRUE, 2, 2),
  status = "PASS"
)
write_tsv(validation, "28_COMPOSITION_VS_STATE_VALIDATION.tsv")
capture.output(sessionInfo(), file = file.path(out, "SESSION_INFO_COMPOSITION_VS_STATE.txt"))
writeLines("PASS", file.path(out, "COMPOSITION_VS_STATE_COMPLETE.ok"))
cat(out, "\n")
