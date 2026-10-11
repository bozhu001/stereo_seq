suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(Matrix)
})
cat("edgeR=", as.character(packageVersion("edgeR")), "\n", sep = "")
cat("limma=", as.character(packageVersion("limma")), "\n", sep = "")
cat("Matrix=", as.character(packageVersion("Matrix")), "\n", sep = "")
