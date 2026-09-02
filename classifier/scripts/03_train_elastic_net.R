# 03_train_elastic_net.R
# Elastic net training with caret/glmnet

source("scripts/00_setup_paths.R")
source("scripts/02_feature_selection.R")

suppressPackageStartupMessages({
  library(caret)
  library(dplyr)
  library(Matrix)
  library(glmnet)
})

# ---- Load prepared input ----
classifier_inputs_path <- run_file(out_paths$data_processed, "classifier_inputs")
cv_split_path <- run_file(out_paths$data_processed, "cv_split")

inp <- readRDS(check_file(classifier_inputs_path))
split <- readRDS(check_file(cv_split_path))

x_cpg <- inp$x_cpg
if (is.null(x_cpg)) {
  x_cpg <- readRDS(check_file(inp$x_cpg_path))
}
meta <- inp$meta
covar_cols <- inp$covar_cols %||% character()
categorical_covar_cols <- inp$categorical_covar_cols %||% covar_cols
numeric_covar_cols <- inp$numeric_covar_cols %||% character()

# ---- Defensive checks and alignment ----
if (!inherits(x_cpg, "Matrix")) {
  x_cpg <- Matrix::Matrix(x_cpg, sparse = TRUE)
}

if (is.null(rownames(x_cpg))) {
  stop("x_cpg has no rownames. Expected sample IDs as rownames.")
}

if (is.null(colnames(x_cpg))) {
  stop("x_cpg has no colnames. Expected CpG IDs as colnames.")
}

if (!"sample_id" %in% colnames(meta)) {
  stop("meta has no sample_id column.")
}

if (!"stage" %in% colnames(meta)) {
  stop("meta has no stage column.")
}

meta <- meta %>%
  filter(sample_id %in% rownames(x_cpg)) %>%
  arrange(match(sample_id, rownames(x_cpg)))

x_cpg <- x_cpg[meta$sample_id, , drop = FALSE]

stopifnot(
  nrow(x_cpg) == nrow(meta),
  identical(rownames(x_cpg), meta$sample_id)
)

meta$stage <- droplevels(meta$stage)

if (nlevels(meta$stage) < 2) {
  stop("Outcome stage has fewer than 2 classes.")
}

# Important:
# Do not use as.matrix() as fallback here.
# That densifies the full methylation matrix and causes the 1.5 GiB warning.
safe_subset <- function(mat, i, j) {
  mat[i, j, drop = FALSE]
}

cv_strategy <- config$modeling$cv_strategy %||% "kfold"

# ---- Define train/test rows ----
# For train_test:
#   feature selection is done only on the training subset.
# For kfold / LOSO:
#   feature selection is currently done once on all samples.
#   This is okay for a first implementation, but later you may want fold-wise FS.
if (cv_strategy == "train_test") {
  train_idx <- as.integer(split$train_idx)
  test_idx <- as.integer(split$test_idx)
  fs_idx <- train_idx
} else {
  train_idx <- seq_len(nrow(meta))
  test_idx <- integer(0)
  fs_idx <- train_idx
}

y_train <- droplevels(meta$stage[train_idx])
y_test <- if (length(test_idx) > 0) droplevels(meta$stage[test_idx]) else NULL

if (nlevels(y_train) < 2) {
  stop("Training outcome y_train has fewer than 2 classes.")
}

# ---- Cache covariate handling from training metadata only ----
if (length(covar_cols) > 0) {
  missing_covars <- setdiff(covar_cols, colnames(meta))

  if (length(missing_covars) > 0) {
    stop(
      paste0(
        "Missing covariate columns in metadata: ",
        paste(missing_covars, collapse = ", ")
      )
    )
  }
}

cov_levels <- NULL

if (length(categorical_covar_cols) > 0) {
  cov_levels <- lapply(meta[train_idx, categorical_covar_cols, drop = FALSE], function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- "Unknown"
    sort(unique(c(x, "Unknown")))
  })
}

numeric_impute_values <- numeric()

if (length(numeric_covar_cols) > 0) {
  numeric_impute_values <- vapply(numeric_covar_cols, function(nm) {
    x <- suppressWarnings(as.numeric(as.character(meta[[nm]][train_idx])))
    med <- stats::median(x, na.rm = TRUE)
    if (is.na(med)) 0 else med
  }, numeric(1))
}

# ---- Feature selection ----
message("Selecting top variable CpGs...")

selected_cpgs_raw <- select_top_variable_cpgs(
  x_mat = safe_subset(x_cpg, fs_idx, seq_len(ncol(x_cpg))),
  n_top = config$feature_selection$top_variable_n
)

# The feature selection function may return either numeric column indices
# or CpG names. Handle both safely.
if (is.numeric(selected_cpgs_raw)) {
  selected_cpgs <- colnames(x_cpg)[selected_cpgs_raw]
} else {
  selected_cpgs <- as.character(selected_cpgs_raw)
}

selected_cpgs <- intersect(selected_cpgs, colnames(x_cpg))

if (length(selected_cpgs) == 0) {
  stop("Feature selection returned no valid CpG column names.")
}

message("Selected CpGs: ", length(selected_cpgs))

# ---- Build design matrix ----
make_design <- function(x_mat, meta_df, include_covariates = TRUE) {
  x <- x_mat

  if (!inherits(x, "Matrix")) {
    x <- Matrix::Matrix(x, sparse = TRUE)
  }

  if (include_covariates && length(categorical_covar_cols) > 0) {
    cov_df <- meta_df |>
      dplyr::select(all_of(categorical_covar_cols)) |>
      dplyr::mutate(
        across(
          everything(),
          ~ {
            z <- as.character(.)
            z[is.na(z) | z == ""] <- "Unknown"
            z
          }
        )
      ) |>
      as.data.frame()

    for (nm in names(cov_levels)) {
      cov_df[[nm]][!cov_df[[nm]] %in% cov_levels[[nm]]] <- "Unknown"
      cov_df[[nm]] <- factor(cov_df[[nm]], levels = cov_levels[[nm]])
    }

    cov_mm <- Matrix::sparse.model.matrix(
      ~ . - 1,
      data = cov_df,
      na.action = stats::na.pass
    )

    if (nrow(cov_mm) != nrow(x)) {
      stop(
        paste0(
          "Covariate design matrix row count does not match CpG matrix.\n",
          "nrow(x) = ", nrow(x), "\n",
          "nrow(cov_mm) = ", nrow(cov_mm), "\n",
          "nrow(meta_df) = ", nrow(meta_df), "\n",
          "This usually means covariate factor levels or missing values were not handled consistently."
        )
      )
    }

    rownames(cov_mm) <- rownames(x)

    x <- Matrix::cbind2(x, cov_mm)
  }

  if (include_covariates && length(numeric_covar_cols) > 0) {
    numeric_df <- meta_df |>
      dplyr::select(all_of(numeric_covar_cols)) |>
      as.data.frame()

    for (nm in numeric_covar_cols) {
      numeric_df[[nm]] <- suppressWarnings(as.numeric(as.character(numeric_df[[nm]])))
      numeric_df[[nm]][is.na(numeric_df[[nm]])] <- numeric_impute_values[[nm]]
    }

    numeric_mm <- Matrix::Matrix(as.matrix(numeric_df), sparse = TRUE)
    rownames(numeric_mm) <- rownames(x)
    colnames(numeric_mm) <- numeric_covar_cols

    x <- Matrix::cbind2(x, numeric_mm)
  }

  x
}

message("Building training design matrix...")

x_train <- make_design(
  safe_subset(x_cpg, train_idx, selected_cpgs),
  meta[train_idx, , drop = FALSE]
)

stopifnot(nrow(x_train) == length(y_train))

message("Training samples: ", nrow(x_train))
message("Training features including covariates: ", ncol(x_train))
message("Training classes: ", paste(levels(y_train), collapse = ", "))

# ---- caret training control ----
if (cv_strategy == "kfold") {
  cv_index <- split$folds_train
  cv_index_out <- split$folds_out
} else if (cv_strategy == "loso") {
  cv_index <- split$loso_train
  cv_index_out <- split$loso_out
} else {
  cv_index <- NULL
  cv_index_out <- NULL
}

ctrl <- caret::trainControl(
  method = "cv",
  number = config$modeling$cv_folds,
  classProbs = TRUE,
  savePredictions = "final",
  index = cv_index,
  indexOut = cv_index_out
)

# ---- Tuning grid ----
alpha_grid <- as.numeric(unlist(config$elastic_net$alpha_grid))

if (length(alpha_grid) == 0 || any(is.na(alpha_grid))) {
  stop("config$elastic_net$alpha_grid is empty or contains NA.")
}

lambda_grid <- 10^seq(-4, 1, length.out = 20)

tune_grid <- expand.grid(
  alpha = alpha_grid,
  lambda = lambda_grid
)

message("Starting elastic net training...")

# ---- Train model ----
set.seed(config$modeling$seed)

fit <- caret::train(
  x = x_train,
  y = y_train,
  method = "glmnet",
  family = ifelse(nlevels(y_train) > 2, "multinomial", "binomial"),
  trControl = ctrl,
  tuneGrid = tune_grid
)

# ---- Optional held-out test prediction for train_test strategy ----
test_predictions <- NULL

if (length(test_idx) > 0) {
  message("Creating held-out test predictions...")

  x_test <- make_design(
    safe_subset(x_cpg, test_idx, selected_cpgs),
    meta[test_idx, , drop = FALSE]
  )

  test_predictions <- list(
    truth = y_test,
    class = predict(fit, newdata = x_test),
    prob = predict(fit, newdata = x_test, type = "prob")
  )
}

# ---- Save outputs ----
saveRDS(
  fit,
  file = run_file(out_paths$results_models, "elastic_net_model")
)

saveRDS(
  list(
    run_id = config$run_id,
    outcome_info = inp$outcome_info,
    covariate_set = inp$covariate_set,
    selected_cpgs = selected_cpgs,
    train_idx = train_idx,
    test_idx = test_idx,
    covar_cols = covar_cols,
    categorical_covar_cols = categorical_covar_cols,
    numeric_covar_cols = numeric_covar_cols,
    cov_levels = cov_levels,
    numeric_impute_values = numeric_impute_values,
    bestTune = fit$bestTune,
    test_predictions = test_predictions
  ),
  file = run_file(out_paths$results_features, "elastic_net_features")
)

message("Elastic net training complete.")
message("Run ID: ", config$run_id)
message("Saved model to: ", run_file(out_paths$results_models, "elastic_net_model"))
message("Saved feature info to: ", run_file(out_paths$results_features, "elastic_net_features"))
message("Best tuning parameters:")
print(fit$bestTune)
