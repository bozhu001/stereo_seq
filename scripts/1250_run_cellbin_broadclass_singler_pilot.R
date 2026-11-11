#!/usr/bin/env Rscript

# Small, gated SingleR pilot for broad CellBin annotation. This script reads
# only prepared sparse matrices and writes pilot diagnostics; it does not alter
# Joint V1 or expand to the full cohort.

suppressPackageStartupMessages({
  library(BiocParallel)
  library(Matrix)
  library(rhdf5)
  library(scuttle)
  library(SingleCellExperiment)
  library(SingleR)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "results", "cell bin", "CellBin_BroadClass_pilot_21samples")
work <- file.path(out, "work")
log_path <- file.path(out, "pilot_run.log")

log_line <- function(...) {
  value <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                 paste0(..., collapse = ""), sep = " | ")
  cat(value, "\n")
  cat(value, "\n", file = log_path, append = TRUE)
}

read_sparse_csc <- function(path) {
  shape <- as.integer(h5read(path, "matrix/shape"))
  matrix <- sparseMatrix(
    i = as.integer(h5read(path, "matrix/indices")) + 1L,
    p = as.integer(h5read(path, "matrix/indptr")),
    x = as.numeric(h5read(path, "matrix/data")),
    dims = shape,
    index1 = TRUE,
    giveCsparse = TRUE
  )
  rownames(matrix) <- as.character(h5read(path, "row_names"))
  colnames(matrix) <- as.character(h5read(path, "col_names"))
  matrix
}

prediction_frame <- function(prediction, ids) {
  scores <- as.data.frame(as.matrix(prediction$scores), check.names = FALSE)
  colnames(scores) <- paste0("score_", gsub("[^A-Za-z0-9]+", "_", colnames(scores)))
  data.frame(
    global_cellbin_id = ids,
    predicted_broad_class_unpruned = as.character(prediction$labels),
    delta_next = as.numeric(prediction$delta.next),
    predicted_broad_class_pruned = as.character(prediction$pruned.labels),
    scores,
    check.names = FALSE
  )
}

set.seed(20261111)
unlink(log_path)
log_line("START SingleR pilot; SingleR=", as.character(packageVersion("SingleR")))

reference_profiles <- read.delim(
  file.path(work, "reference_profiles.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
pilot_metadata <- read.delim(
  file.path(work, "pilot_query_metadata.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
reference_counts <- read_sparse_csc(file.path(work, "reference_profiles_counts.h5"))
pilot_counts <- read_sparse_csc(file.path(work, "pilot_query_counts.h5"))

stopifnot(
  identical(rownames(reference_counts), rownames(pilot_counts)),
  identical(colnames(reference_counts), reference_profiles$profile_id),
  identical(colnames(pilot_counts), pilot_metadata$global_cellbin_id),
  all(reference_counts@x >= 0),
  all(reference_counts@x == round(reference_counts@x)),
  all(pilot_counts@x >= 0),
  all(pilot_counts@x == round(pilot_counts@x))
)
log_line("Loaded sparse matrices; reference=", paste(dim(reference_counts), collapse = "x"),
         "; query=", paste(dim(pilot_counts), collapse = "x"))

reference_sce <- SingleCellExperiment(list(counts = reference_counts))
pilot_sce <- SingleCellExperiment(list(counts = pilot_counts))
reference_sce <- logNormCounts(reference_sce)
pilot_sce <- logNormCounts(pilot_sce)

bp <- SerialParam()
prediction <- SingleR(
  test = pilot_sce,
  ref = reference_sce,
  labels = reference_profiles$broad_class,
  genes = "de",
  de.method = "classic",
  fine.tune = TRUE,
  prune = TRUE,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  BPPARAM = bp
)
pilot_result <- prediction_frame(prediction, colnames(pilot_counts))
write.table(
  pilot_result,
  file.path(work, "pilot_singler_predictions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
log_line("Pilot prediction complete; pruned Unknown=",
         sum(is.na(pilot_result$predicted_broad_class_pruned)))

# Five-fold donor-grouped reference self-validation. Every donor appears in
# exactly one held-out fold; profiles from the same donor never enter both the
# training and test reference in a fold.
donors <- sort(unique(reference_profiles$donor_id))
donor_fold <- setNames((seq_along(donors) - 1L) %% 5L + 1L, donors)
reference_profiles$cv_fold <- donor_fold[reference_profiles$donor_id]
cv_rows <- list()
for (fold in seq_len(5L)) {
  test_index <- which(reference_profiles$cv_fold == fold)
  train_index <- which(reference_profiles$cv_fold != fold)
  fold_prediction <- SingleR(
    test = reference_sce[, test_index],
    ref = reference_sce[, train_index],
    labels = reference_profiles$broad_class[train_index],
    genes = "de",
    de.method = "classic",
    fine.tune = TRUE,
    prune = TRUE,
    assay.type.test = "logcounts",
    assay.type.ref = "logcounts",
    BPPARAM = bp
  )
  frame <- prediction_frame(
    fold_prediction,
    reference_profiles$profile_id[test_index]
  )
  frame$cv_fold <- fold
  frame$donor_id <- reference_profiles$donor_id[test_index]
  frame$reference18_state <- reference_profiles$reference18_state[test_index]
  frame$true_broad_class <- reference_profiles$broad_class[test_index]
  cv_rows[[fold]] <- frame
  log_line("Completed donor CV fold ", fold, "; n=", length(test_index))
}
cv <- do.call(rbind, cv_rows)
write.table(
  cv,
  file.path(work, "reference_donor_grouped_cv.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

session <- capture.output(sessionInfo())
writeLines(session, file.path(work, "pilot_R_sessionInfo.txt"))
log_line("COMPLETE; CV profiles=", nrow(cv))
