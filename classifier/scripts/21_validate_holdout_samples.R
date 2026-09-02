# 21_validate_holdout_samples.R
# Predict optional validation-holdout samples with trained train/test models.

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(caret)
  library(dplyr)
  library(Matrix)
  library(readr)
  library(tibble)
  library(tidyr)
  library(xgboost)
})

parse_model_list <- function() {
  value <- Sys.getenv("CLASSIFIER_VALIDATION_MODELS", "elastic_net,random_forest,boosted_trees")
  models <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  models[nzchar(models)]
}

safe_subset <- function(mat, i, j) {
  mat[i, j, drop = FALSE]
}

read_existing <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

make_sparse_design <- function(x_mat, meta_df, feat) {
  x <- x_mat
  if (!inherits(x, "Matrix")) {
    x <- Matrix::Matrix(x, sparse = TRUE)
  }

  categorical_covar_cols <- feat$categorical_covar_cols %||% character()
  numeric_covar_cols <- feat$numeric_covar_cols %||% character()
  cov_levels <- feat$cov_levels %||% feat$categorical_levels %||% NULL
  numeric_impute_values <- feat$numeric_impute_values %||% numeric()

  if (length(categorical_covar_cols) > 0) {
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
    rownames(cov_mm) <- rownames(x)
    x <- Matrix::cbind2(x, cov_mm)
  }

  if (length(numeric_covar_cols) > 0) {
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

make_predictor_df <- function(x_mat, meta_df, feat) {
  x_df <- as.data.frame(as.matrix(x_mat), check.names = FALSE)

  categorical_covar_cols <- feat$categorical_covar_cols %||% character()
  numeric_covar_cols <- feat$numeric_covar_cols %||% character()
  categorical_levels <- feat$categorical_levels %||% feat$cov_levels %||% NULL
  numeric_impute_values <- feat$numeric_impute_values %||% numeric()

  if (length(categorical_covar_cols) > 0) {
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

    for (nm in names(categorical_levels)) {
      cov_df[[nm]][!cov_df[[nm]] %in% categorical_levels[[nm]]] <- "Unknown"
      cov_df[[nm]] <- factor(cov_df[[nm]], levels = categorical_levels[[nm]])
    }

    x_df <- dplyr::bind_cols(x_df, cov_df)
  }

  if (length(numeric_covar_cols) > 0) {
    numeric_df <- meta_df |>
      dplyr::select(all_of(numeric_covar_cols)) |>
      as.data.frame()

    for (nm in numeric_covar_cols) {
      numeric_df[[nm]] <- suppressWarnings(as.numeric(as.character(numeric_df[[nm]])))
      numeric_df[[nm]][is.na(numeric_df[[nm]])] <- numeric_impute_values[[nm]]
    }

    x_df <- dplyr::bind_cols(x_df, numeric_df)
  }

  x_df
}

align_xgb_matrix <- function(x_df, feature_names) {
  x_mm <- Matrix::sparse.model.matrix(~ . - 1, data = x_df)

  missing_cols <- setdiff(feature_names, colnames(x_mm))
  if (length(missing_cols) > 0) {
    zero_mm <- Matrix::Matrix(
      0,
      nrow = nrow(x_mm),
      ncol = length(missing_cols),
      sparse = TRUE
    )
    colnames(zero_mm) <- missing_cols
    rownames(zero_mm) <- rownames(x_mm)
    x_mm <- Matrix::cbind2(x_mm, zero_mm)
  }

  extra_cols <- setdiff(colnames(x_mm), feature_names)
  if (length(extra_cols) > 0) {
    x_mm <- x_mm[, setdiff(colnames(x_mm), extra_cols), drop = FALSE]
  }

  x_mm[, feature_names, drop = FALSE]
}

predict_xgb_prob <- function(booster, x_mat, class_levels) {
  raw_pred <- predict(booster, newdata = x_mat)

  if (length(class_levels) == 2) {
    prob <- data.frame(
      stats::setNames(list(1 - raw_pred, raw_pred), class_levels),
      check.names = FALSE
    )
  } else {
    prob_mat <- matrix(raw_pred, ncol = length(class_levels), byrow = FALSE)
    colnames(prob_mat) <- class_levels
    prob <- as.data.frame(prob_mat, check.names = FALSE)
  }

  prob
}

predict_one_model <- function(model_name, inp, x_cpg, holdout_meta, class_levels) {
  model_path <- run_file(out_paths$results_models, paste0(model_name, "_model"))
  feature_path <- run_file(out_paths$results_features, paste0(model_name, "_features"))

  fit <- read_existing(model_path)
  feat <- read_existing(feature_path)

  if (is.null(fit) || is.null(feat)) {
    warning("Skipping ", model_name, ": missing model or feature file for run_id ", config$run_id)
    return(NULL)
  }

  selected_cpgs <- intersect(feat$selected_cpgs, colnames(x_cpg))
  if (length(selected_cpgs) == 0) {
    stop("No valid selected CpGs for ", model_name, ".", call. = FALSE)
  }

  x_holdout_cpg <- safe_subset(x_cpg, holdout_meta$sample_id, selected_cpgs)

  if (model_name == "elastic_net") {
    x_holdout <- make_sparse_design(x_holdout_cpg, holdout_meta, feat)
    pred_class <- predict(fit, newdata = x_holdout)
    pred_prob <- predict(fit, newdata = x_holdout, type = "prob")
  } else if (model_name == "random_forest") {
    x_holdout <- make_predictor_df(x_holdout_cpg, holdout_meta, feat)
    pred_class <- predict(fit, newdata = x_holdout)
    pred_prob <- predict(fit, newdata = x_holdout, type = "prob")
  } else if (model_name == "boosted_trees") {
    x_holdout_df <- make_predictor_df(x_holdout_cpg, holdout_meta, feat)
    x_holdout <- align_xgb_matrix(x_holdout_df, fit$feature_names)
    pred_prob <- predict_xgb_prob(fit$booster, x_holdout, class_levels)
    pred_class <- class_levels[max.col(as.matrix(pred_prob), ties.method = "first")]
  } else {
    stop("Unsupported validation model: ", model_name, call. = FALSE)
  }

  prob_tbl <- tibble::as_tibble(pred_prob)
  missing_prob_cols <- setdiff(class_levels, colnames(prob_tbl))
  for (nm in missing_prob_cols) {
    prob_tbl[[nm]] <- NA_real_
  }
  prob_tbl <- prob_tbl[, class_levels, drop = FALSE]

  tibble::tibble(
    model = model_name,
    run_id = config$run_id,
    holdout = inp$holdout_info$name,
    sample_id = holdout_meta$sample_id,
    dataset = holdout_meta[[config$covariates$study_id_col]],
    source_disease_group = holdout_meta[[inp$outcome_info$source_column]],
    mapped_truth = factor(holdout_meta$stage, levels = class_levels),
    prediction = factor(pred_class, levels = class_levels)
  ) |>
    dplyr::bind_cols(prob_tbl)
}

inp <- readRDS(check_file(run_file(out_paths$data_processed, "classifier_inputs")))

if (is.null(inp$holdout_info) || is.null(inp$holdout_meta) || nrow(inp$holdout_meta) == 0) {
  stop(
    "No validation holdout stored in classifier inputs. Re-run 01_prepare_classifier_data.R with CLASSIFIER_HOLDOUT=...",
    call. = FALSE
  )
}

if ((config$modeling$cv_strategy %||% "") != "train_test") {
  stop("Holdout validation currently expects a train_test model run.", call. = FALSE)
}

x_cpg <- readRDS(check_file(inp$x_cpg_path))
if (!inherits(x_cpg, "Matrix")) {
  x_cpg <- Matrix::Matrix(x_cpg, sparse = TRUE)
}

holdout_meta <- inp$holdout_meta |>
  dplyr::filter(sample_id %in% rownames(x_cpg)) |>
  dplyr::arrange(match(sample_id, rownames(x_cpg)))

class_levels <- inp$outcome_info$levels
models <- parse_model_list()

predictions <- dplyr::bind_rows(lapply(models, predict_one_model,
  inp = inp,
  x_cpg = x_cpg,
  holdout_meta = holdout_meta,
  class_levels = class_levels
))

if (nrow(predictions) == 0) {
  stop("No holdout predictions were generated.", call. = FALSE)
}

summary_tbl <- predictions |>
  dplyr::count(model, prediction, name = "n") |>
  dplyr::group_by(model) |>
  dplyr::mutate(
    total = sum(n),
    fraction = n / total
  ) |>
  dplyr::ungroup() |>
  tidyr::complete(
    model,
    prediction = factor(class_levels, levels = class_levels),
    fill = list(n = 0, total = length(unique(predictions$sample_id)), fraction = 0)
  ) |>
  dplyr::arrange(model, prediction)

summary_by_source_tbl <- predictions |>
  dplyr::count(model, source_disease_group, mapped_truth, prediction, name = "n") |>
  dplyr::group_by(model, source_disease_group, mapped_truth) |>
  dplyr::mutate(
    total = sum(n),
    fraction = n / total
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(model, source_disease_group, mapped_truth, prediction)

prob_long <- predictions |>
  tidyr::pivot_longer(
    cols = all_of(class_levels),
    names_to = "class",
    values_to = "probability"
  )

prob_summary <- prob_long |>
  dplyr::group_by(model, class) |>
  dplyr::summarise(
    mean_probability = mean(probability, na.rm = TRUE),
    median_probability = stats::median(probability, na.rm = TRUE),
    q25_probability = stats::quantile(probability, 0.25, na.rm = TRUE),
    q75_probability = stats::quantile(probability, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(model, class)

prob_by_source_summary <- prob_long |>
  dplyr::group_by(model, source_disease_group, mapped_truth, class) |>
  dplyr::summarise(
    mean_probability = mean(probability, na.rm = TRUE),
    median_probability = stats::median(probability, na.rm = TRUE),
    q25_probability = stats::quantile(probability, 0.25, na.rm = TRUE),
    q75_probability = stats::quantile(probability, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(model, source_disease_group, mapped_truth, class)

readr::write_csv(
  predictions,
  run_file(out_paths$results_metrics, "holdout_predictions", "csv")
)
readr::write_csv(
  summary_tbl,
  run_file(out_paths$results_metrics, "holdout_prediction_summary", "csv")
)
readr::write_csv(
  summary_by_source_tbl,
  run_file(out_paths$results_metrics, "holdout_prediction_by_source_summary", "csv")
)
readr::write_csv(
  prob_summary,
  run_file(out_paths$results_metrics, "holdout_probability_summary", "csv")
)
readr::write_csv(
  prob_by_source_summary,
  run_file(out_paths$results_metrics, "holdout_probability_by_source_summary", "csv")
)

message("Holdout validation complete.")
message("Run ID: ", config$run_id)
message("Holdout: ", inp$holdout_info$name, " (", nrow(holdout_meta), " samples)")
message("Saved predictions to: ", run_file(out_paths$results_metrics, "holdout_predictions", "csv"))
print(summary_tbl, n = Inf)
print(summary_by_source_tbl, n = Inf)
