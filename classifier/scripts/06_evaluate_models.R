# 06_evaluate_models.R
# Compute evaluation metrics with balanced focus

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(caret)
  library(dplyr)
  library(yardstick)
  library(pROC)
  library(Matrix)
  library(tibble)
})

# ---- Load inputs ----
classifier_inputs_path <- run_file(out_paths$data_processed, "classifier_inputs")
cv_split_path <- run_file(out_paths$data_processed, "cv_split")
elastic_net_model_path <- run_file(out_paths$results_models, "elastic_net_model")
elastic_net_features_path <- run_file(out_paths$results_features, "elastic_net_features")

inp <- readRDS(check_file(classifier_inputs_path))
split <- readRDS(check_file(cv_split_path))
fit <- readRDS(check_file(elastic_net_model_path))
feat <- readRDS(check_file(elastic_net_features_path))

x_cpg <- inp$x_cpg
if (is.null(x_cpg)) {
  x_cpg <- readRDS(check_file(inp$x_cpg_path))
}
meta <- inp$meta
covar_cols <- feat$covar_cols %||% inp$covar_cols %||% character()
categorical_covar_cols <- feat$categorical_covar_cols %||% inp$categorical_covar_cols %||% covar_cols
numeric_covar_cols <- feat$numeric_covar_cols %||% inp$numeric_covar_cols %||% character()
cov_levels <- feat$cov_levels %||% NULL
numeric_impute_values <- feat$numeric_impute_values %||% numeric()

cv_strategy <- config$modeling$cv_strategy %||% "kfold"

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

# ---- Selected CpGs ----
selected_cpgs <- feat$selected_cpgs
selected_cpgs <- intersect(selected_cpgs, colnames(x_cpg))

if (length(selected_cpgs) == 0) {
  stop("No valid selected CpGs found in x_cpg.")
}

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

# ---- Helper: build design matrix same as in training ----
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
          "nrow(meta_df) = ", nrow(meta_df)
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

# ---- Helper: compute metrics ----
compute_metrics <- function(y_true, pred_class, pred_prob = NULL) {
  y_true <- droplevels(as.factor(y_true))
  pred_class <- factor(pred_class, levels = levels(y_true))

  cm <- caret::confusionMatrix(pred_class, y_true)

  results <- tibble::tibble(
    accuracy = yardstick::accuracy_vec(y_true, pred_class),
    balanced_accuracy = yardstick::bal_accuracy_vec(y_true, pred_class),
    mcc = yardstick::mcc_vec(y_true, pred_class),
    f1_macro = yardstick::f_meas_vec(y_true, pred_class, estimator = "macro")
  )

  by_class <- as.data.frame(cm$byClass)

  auc_tbl <- NULL

  if (!is.null(pred_prob)) {
    available_prob_cols <- intersect(levels(y_true), colnames(pred_prob))

    if (length(available_prob_cols) == nlevels(y_true)) {
      if (nlevels(y_true) == 2) {
        pos_level <- levels(y_true)[2]

        auc_val <- pROC::roc(
          response = y_true,
          predictor = pred_prob[[pos_level]],
          levels = levels(y_true),
          quiet = TRUE
        ) |>
          pROC::auc() |>
          as.numeric()

        auc_tbl <- tibble::tibble(
          auc = auc_val,
          type = "binary",
          positive_class = pos_level
        )

      } else {
        auc_vals <- vapply(levels(y_true), function(cls) {
          response_binary <- factor(
            y_true == cls,
            levels = c(FALSE, TRUE)
          )

          pROC::roc(
            response = response_binary,
            predictor = pred_prob[[cls]],
            levels = c(FALSE, TRUE),
            quiet = TRUE
          ) |>
            pROC::auc() |>
            as.numeric()
        }, numeric(1))

        auc_tbl <- tibble::tibble(
          class = names(auc_vals),
          auc = as.numeric(auc_vals),
          type = "one_vs_rest"
        ) |>
          bind_rows(
            tibble::tibble(
              class = "macro_mean",
              auc = mean(auc_vals, na.rm = TRUE),
              type = "ovr_macro"
            )
          )
      }
    }
  }

  list(
    metrics = results,
    confusion_matrix = cm,
    by_class = by_class,
    auc = auc_tbl
  )
}

# ---- Evaluation ----
if (cv_strategy == "train_test") {
  message("Evaluating held-out test set...")

  test_idx <- as.integer(split$test_idx)

  if (length(test_idx) == 0) {
    stop("cv_strategy is train_test, but split$test_idx is empty.")
  }

  x_test <- make_design(
    x_cpg[test_idx, selected_cpgs, drop = FALSE],
    meta[test_idx, , drop = FALSE]
  )

  y_test <- droplevels(meta$stage[test_idx])

  pred_class <- predict(fit, newdata = x_test)

  pred_prob <- tryCatch(
    predict(fit, newdata = x_test, type = "prob"),
    error = function(e) NULL
  )

  eval <- compute_metrics(
    y_true = y_test,
    pred_class = pred_class,
    pred_prob = pred_prob
  )

  predictions <- tibble::tibble(
    sample_id = meta$sample_id[test_idx],
    truth = y_test,
    prediction = pred_class
  )

  if (!is.null(pred_prob)) {
    predictions <- bind_cols(predictions, as_tibble(pred_prob))
  }

} else {
  message("Evaluating cross-validation predictions from caret fit$pred...")

  if (is.null(fit$pred) || nrow(fit$pred) == 0) {
    stop(
      "fit$pred is empty. Re-run 03_train_elastic_net.R with savePredictions = 'final' in trainControl()."
    )
  }

  pred_df <- fit$pred

  # If caret stored predictions for multiple tuning parameters, keep only bestTune.
  if (!is.null(fit$bestTune)) {
    for (nm in names(fit$bestTune)) {
      if (nm %in% colnames(pred_df)) {
        pred_df <- pred_df[pred_df[[nm]] == fit$bestTune[[nm]], , drop = FALSE]
      }
    }
  }

  if (!all(c("obs", "pred") %in% colnames(pred_df))) {
    stop("fit$pred must contain columns 'obs' and 'pred'.")
  }

  y_cv <- droplevels(as.factor(pred_df$obs))
  pred_class <- factor(pred_df$pred, levels = levels(y_cv))

  prob_cols <- intersect(levels(y_cv), colnames(pred_df))
  pred_prob <- NULL

  if (length(prob_cols) == nlevels(y_cv)) {
    pred_prob <- pred_df[, prob_cols, drop = FALSE]
  }

  eval <- compute_metrics(
    y_true = y_cv,
    pred_class = pred_class,
    pred_prob = pred_prob
  )

  predictions <- pred_df
}

# ---- Save outputs ----
saveRDS(
  list(
    run_id = config$run_id,
    outcome_info = inp$outcome_info,
    covariate_set = inp$covariate_set,
    cv_strategy = cv_strategy,
    bestTune = fit$bestTune,
    metrics = eval$metrics,
    confusion_matrix = eval$confusion_matrix,
    by_class = eval$by_class,
    auc = eval$auc,
    predictions = predictions
  ),
  file = run_file(out_paths$results_metrics, "elastic_net_metrics")
)

readr::write_csv(
  eval$metrics,
  run_file(out_paths$results_metrics, "elastic_net_metrics", "csv")
)

confusion_tbl <- as.data.frame(eval$confusion_matrix$table)
readr::write_csv(
  confusion_tbl,
  run_file(out_paths$results_metrics, "elastic_net_confusion_matrix", "csv")
)

by_class_tbl <- eval$by_class |>
  tibble::rownames_to_column("class")
readr::write_csv(
  by_class_tbl,
  run_file(out_paths$results_metrics, "elastic_net_by_class_metrics", "csv")
)

readr::write_csv(
  predictions,
  run_file(out_paths$results_metrics, "elastic_net_predictions", "csv")
)

if (!is.null(eval$auc)) {
  readr::write_csv(
    eval$auc,
    run_file(out_paths$results_metrics, "elastic_net_auc", "csv")
  )
}

yaml::write_yaml(
  inp$outcome_info,
  file = run_file(out_paths$results_metrics, "elastic_net_outcome_definition", "yaml")
)

message("Evaluation complete.")
message("Run ID: ", config$run_id)
message("Saved metrics to: ", run_file(out_paths$results_metrics, "elastic_net_metrics"))
print(eval$metrics)

if (!is.null(eval$auc)) {
  print(eval$auc)
}
