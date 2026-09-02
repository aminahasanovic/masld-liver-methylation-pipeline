# 02_feature_selection.R
# Feature selection helpers (run inside training folds only)

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(matrixStats)
  library(limma)
  library(Matrix)
})

# ---- Functions ----

select_top_variable_cpgs <- function(x_mat, n_top = 5000) {
  x <- Matrix::Matrix(x_mat, sparse = TRUE)
  col_means <- Matrix::colMeans(x)
  col_means_sq <- Matrix::colMeans(x^2)
  vars <- col_means_sq - (col_means^2)
  names(vars) <- colnames(x)
  names(sort(vars, decreasing = TRUE))[seq_len(min(n_top, length(vars)))]
}

select_diffmeth_cpgs <- function(x_mat, y, n_top = 5000, adj_p = 0.05) {
  # y is a factor (classes)
  design <- model.matrix(~ 0 + y)
  colnames(design) <- levels(y)

  fit <- limma::lmFit(t(x_mat), design)
  # One-vs-rest for each class, then rank by min adj.P.Val
  contrasts <- limma::makeContrasts(contrasts = colnames(design), levels = design)
  fit2 <- limma::contrasts.fit(fit, contrasts)
  fit2 <- limma::eBayes(fit2)

  tt <- limma::topTable(fit2, number = Inf, sort.by = "P")
  tt <- tt |>
    dplyr::filter(adj.P.Val <= adj_p)

  if (nrow(tt) == 0) {
    return(character())
  }

  tt <- tt |>
    dplyr::slice_head(n = min(n_top, nrow(tt)))

  rownames(tt)
}

# ---- Example usage (do not run globally) ----
# inp <- readRDS(run_file(out_paths$data_processed, "classifier_inputs"))
# split <- readRDS(run_file(out_paths$data_processed, "cv_split"))
# train_idx <- split$train_idx
# x_train <- x_cpg[train_idx, ]
# y_train <- meta$stage[train_idx]
# cpgs_var <- select_top_variable_cpgs(x_train, n_top = config$feature_selection$top_variable_n)
# cpgs_dm <- select_diffmeth_cpgs(x_train, y_train,
#   n_top = config$feature_selection$diffmeth_top_n,
#   adj_p = config$feature_selection$diffmeth_adj_p
# )
