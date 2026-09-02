# 07_evaluate_loso_elastic_net.R
# Leave-one-study-out evaluation with fold-wise feature selection.

source("scripts/00_setup_paths.R")
source("scripts/02_feature_selection.R")

suppressPackageStartupMessages({
  library(caret)
  library(dplyr)
  library(Matrix)
  library(glmnet)
  library(pROC)
  library(tibble)
  library(yardstick)
})

if ((config$modeling$cv_strategy %||% "") != "loso") {
  stop("Set CLASSIFIER_CV_STRATEGY=loso or config$modeling$cv_strategy: loso before running this script.")
}

alpha_grid <- as.numeric(unlist(config$elastic_net$alpha_grid))
if (any(is.na(alpha_grid))) {
  stop("alpha_grid contains non-numeric values after parsing.", call. = FALSE)
}

analysis_mode <- as.character(config$analysis_mode %||% "")
if (identical(analysis_mode, "sparse_stability") && any(alpha_grid <= 0)) {
  stop(
    paste0(
      "Sparse-stability LOSO requires all elastic_net.alpha_grid values to be > 0; ",
      "ridge alpha = 0 is not valid for sparse feature-selection stability. ",
      "Loaded alpha grid: ",
      paste(alpha_grid, collapse = ", ")
    ),
    call. = FALSE
  )
}

# ---- Load prepared input ----
inp <- readRDS(check_file(run_file(out_paths$data_processed, "classifier_inputs")))
split <- readRDS(check_file(run_file(out_paths$data_processed, "cv_split")))

if (is.null(split$loso_train) || is.null(split$loso_out)) {
  stop("cv_split does not contain loso_train/loso_out. Re-run 01_prepare_classifier_data.R with cv_strategy = loso.")
}

x_cpg <- inp$x_cpg
if (is.null(x_cpg)) {
  x_cpg <- readRDS(check_file(inp$x_cpg_path))
}

meta <- inp$meta
covar_cols <- inp$covar_cols %||% character()
categorical_covar_cols <- inp$categorical_covar_cols %||% covar_cols
numeric_covar_cols <- inp$numeric_covar_cols %||% character()
tune_metric <- config$modeling$tune_metric %||% inp$tune_metric %||% "Accuracy"

if (!inherits(x_cpg, "Matrix")) {
  x_cpg <- Matrix::Matrix(x_cpg, sparse = TRUE)
}

meta <- meta |>
  dplyr::filter(sample_id %in% rownames(x_cpg)) |>
  dplyr::arrange(match(sample_id, rownames(x_cpg)))

x_cpg <- x_cpg[meta$sample_id, , drop = FALSE]
meta$stage <- droplevels(meta$stage)

stopifnot(
  nrow(x_cpg) == nrow(meta),
  identical(rownames(x_cpg), meta$sample_id)
)

safe_subset <- function(mat, i, j) {
  mat[i, j, drop = FALSE]
}

make_design <- function(x_mat, meta_df, categorical_levels, numeric_impute_values) {
  x <- x_mat

  if (!inherits(x, "Matrix")) {
    x <- Matrix::Matrix(x, sparse = TRUE)
  }

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

compute_metrics <- function(y_true, pred_class, pred_prob = NULL, class_levels = NULL) {
  if (is.null(class_levels)) {
    class_levels <- levels(droplevels(as.factor(y_true)))
  }

  y_true <- factor(y_true, levels = class_levels)
  pred_class <- factor(pred_class, levels = class_levels)

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
    available_prob_cols <- intersect(class_levels, colnames(pred_prob))

    if (length(available_prob_cols) == length(class_levels)) {
      if (length(class_levels) == 2) {
        pos_level <- class_levels[2]

        auc_val <- tryCatch(
          pROC::roc(
            response = y_true,
            predictor = pred_prob[[pos_level]],
            levels = class_levels,
            quiet = TRUE
          ) |>
            pROC::auc() |>
            as.numeric(),
          error = function(e) NA_real_
        )

        auc_tbl <- tibble::tibble(
          auc = auc_val,
          type = "binary",
          positive_class = pos_level
        )
      } else {
        auc_vals <- vapply(class_levels, function(cls) {
          response_binary <- factor(y_true == cls, levels = c(FALSE, TRUE))

          tryCatch(
            pROC::roc(
              response = response_binary,
              predictor = pred_prob[[cls]],
              levels = c(FALSE, TRUE),
              quiet = TRUE
            ) |>
              pROC::auc() |>
              as.numeric(),
            error = function(e) NA_real_
          )
        }, numeric(1))

        auc_tbl <- tibble::tibble(
          class = names(auc_vals),
          auc = as.numeric(auc_vals),
          type = "one_vs_rest"
        ) |>
          dplyr::bind_rows(
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

weighted_macro_f1 <- function(obs, pred, class_levels) {
  obs <- factor(obs, levels = class_levels)
  pred <- factor(pred, levels = class_levels)
  tab <- table(obs, pred)
  tp <- diag(tab)
  precision <- tp / colSums(tab)
  recall <- tp / rowSums(tab)
  f1 <- 2 * precision * recall / (precision + recall)
  f1[!is.finite(f1)] <- NA_real_
  mean(f1, na.rm = TRUE)
}

balanced_accuracy_manual <- function(obs, pred, class_levels) {
  obs <- factor(obs, levels = class_levels)
  pred <- factor(pred, levels = class_levels)
  tab <- table(obs, pred)
  recall <- diag(tab) / rowSums(tab)
  recall[!is.finite(recall)] <- NA_real_
  mean(recall, na.rm = TRUE)
}

weighted_summary <- function(data, lev = NULL, model = NULL) {
  if (is.null(lev)) {
    lev <- levels(data$obs)
  }

  obs <- factor(data$obs, levels = lev)
  pred <- factor(data$pred, levels = lev)

  c(
    Accuracy = mean(obs == pred, na.rm = TRUE),
    Balanced_Accuracy = balanced_accuracy_manual(obs, pred, lev),
    Macro_F1 = weighted_macro_f1(obs, pred, lev)
  )
}

lambda_grid <- 10^seq(-4, 1, length.out = 20)
tune_grid <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)
class_levels <- levels(meta$stage)

fold_predictions <- list()
fold_summaries <- list()
fold_features <- list()
fold_coefficients <- list()

extract_nonzero_coefficients <- function(fit, feature_names, selected_cpgs, class_levels) {
  coef_obj <- stats::coef(fit$finalModel, s = fit$bestTune$lambda)

  coef_one_class <- function(coef_mat, class_label) {
    coef_dense <- as.matrix(coef_mat)
    tibble::tibble(
      class = class_label,
      feature = rownames(coef_dense),
      coefficient = as.numeric(coef_dense[, 1])
    )
  }

  coef_tbl <- if (is.list(coef_obj)) {
    dplyr::bind_rows(lapply(names(coef_obj), function(class_label) {
      coef_one_class(coef_obj[[class_label]], class_label)
    }))
  } else {
    positive_class <- if (length(class_levels) == 2) class_levels[2] else NA_character_
    coef_one_class(coef_obj, positive_class)
  }

  coef_tbl |>
    dplyr::filter(coefficient != 0) |>
    dplyr::mutate(
      feature_type = dplyr::case_when(
        feature == "(Intercept)" ~ "intercept",
        feature %in% selected_cpgs ~ "CpG",
        feature %in% feature_names ~ "model_feature",
        TRUE ~ "unknown"
      ),
      abs_coefficient = abs(coefficient),
      alpha = fit$bestTune$alpha,
      lambda = fit$bestTune$lambda
    ) |>
    dplyr::select(class, feature, feature_type, coefficient, abs_coefficient, alpha, lambda)
}

message("Starting LOSO elastic net evaluation.")
message("Run ID: ", config$run_id)
message("Outcome: ", inp$outcome_info$name)
message("Covariate set: ", inp$covariate_set)
message("Tune metric: ", tune_metric)
message("Alpha grid: ", paste(alpha_grid, collapse = ", "))

for (fold_name in names(split$loso_out)) {
  message("LOSO fold held out: ", fold_name)

  train_idx <- as.integer(split$loso_train[[fold_name]])
  test_idx <- as.integer(split$loso_out[[fold_name]])

  y_train <- droplevels(meta$stage[train_idx])
  y_test <- factor(meta$stage[test_idx], levels = class_levels)

  if (nlevels(y_train) < 2) {
    stop("Training fold for ", fold_name, " has fewer than 2 classes.")
  }

  missing_train_classes <- setdiff(class_levels, levels(y_train))
  if (length(missing_train_classes) > 0) {
    stop(
      "Training fold for ", fold_name,
      " is missing class(es): ", paste(missing_train_classes, collapse = ", ")
    )
  }

  selected_cpgs <- select_top_variable_cpgs(
    x_mat = safe_subset(x_cpg, train_idx, seq_len(ncol(x_cpg))),
    n_top = config$feature_selection$top_variable_n
  )

  selected_cpgs <- intersect(as.character(selected_cpgs), colnames(x_cpg))

  categorical_levels <- NULL
  if (length(categorical_covar_cols) > 0) {
    categorical_levels <- lapply(meta[train_idx, categorical_covar_cols, drop = FALSE], function(x) {
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

  x_train <- make_design(
    safe_subset(x_cpg, train_idx, selected_cpgs),
    meta[train_idx, , drop = FALSE],
    categorical_levels,
    numeric_impute_values
  )

  x_test <- make_design(
    safe_subset(x_cpg, test_idx, selected_cpgs),
    meta[test_idx, , drop = FALSE],
    categorical_levels,
    numeric_impute_values
  )

  min_class_n <- min(table(y_train))
  inner_folds <- min(5, as.integer(min_class_n))

  if (inner_folds < 2) {
    stop("Training fold for ", fold_name, " has a class with fewer than 2 samples.")
  }

  ctrl <- caret::trainControl(
    method = "cv",
    number = inner_folds,
    classProbs = TRUE,
    savePredictions = "final",
    summaryFunction = weighted_summary
  )

  set.seed(config$modeling$seed)
  fit <- caret::train(
    x = x_train,
    y = y_train,
    method = "glmnet",
    family = ifelse(length(class_levels) > 2, "multinomial", "binomial"),
    trControl = ctrl,
    tuneGrid = tune_grid,
    metric = tune_metric
  )

  pred_class <- predict(fit, newdata = x_test)
  pred_prob <- predict(fit, newdata = x_test, type = "prob")

  pred_tbl <- tibble::tibble(
    outer_fold = fold_name,
    sample_id = meta$sample_id[test_idx],
    dataset = meta[[config$covariates$study_id_col]][test_idx],
    truth = y_test,
    prediction = factor(pred_class, levels = class_levels)
  ) |>
    dplyr::bind_cols(tibble::as_tibble(pred_prob))

  fold_eval <- compute_metrics(
    y_true = y_test,
    pred_class = pred_class,
    pred_prob = pred_prob,
    class_levels = class_levels
  )

  fold_predictions[[fold_name]] <- pred_tbl
  fold_features[[fold_name]] <- selected_cpgs
  fold_coefficients[[fold_name]] <- extract_nonzero_coefficients(
    fit = fit,
    feature_names = colnames(x_train),
    selected_cpgs = selected_cpgs,
    class_levels = class_levels
  ) |>
    dplyr::mutate(
      outer_fold = fold_name,
      train_n = length(train_idx),
      test_n = length(test_idx),
      .before = 1
    )
  fold_summaries[[fold_name]] <- dplyr::bind_cols(
    tibble::tibble(
      outer_fold = fold_name,
      train_n = length(train_idx),
      test_n = length(test_idx),
      selected_cpgs_n = length(selected_cpgs),
      inner_folds = inner_folds,
      tune_metric = tune_metric,
      min_class_n = min_class_n,
      alpha = fit$bestTune$alpha,
      lambda = fit$bestTune$lambda
    ),
    fold_eval$metrics
  )
}

predictions <- dplyr::bind_rows(fold_predictions)
fold_metrics <- dplyr::bind_rows(fold_summaries)

pooled_eval <- compute_metrics(
  y_true = predictions$truth,
  pred_class = predictions$prediction,
  pred_prob = predictions[, class_levels, drop = FALSE],
  class_levels = class_levels
)

feature_frequency <- table(unlist(fold_features))
feature_frequency <- tibble::tibble(
  cpg = names(feature_frequency),
  selected_folds = as.integer(feature_frequency)
) |>
  dplyr::arrange(dplyr::desc(selected_folds), cpg)

nonzero_coefficients <- dplyr::bind_rows(fold_coefficients)

nonzero_feature_frequency <- nonzero_coefficients |>
  dplyr::filter(feature_type == "CpG") |>
  dplyr::group_by(feature) |>
  dplyr::summarise(
    # Multinomial CpGs are selected for a fold if any class-specific coefficient is non-zero.
    nonzero_folds = dplyr::n_distinct(outer_fold),
    nonzero_fold_class_terms = dplyr::n(),
    nonzero_classes = dplyr::n_distinct(class),
    mean_abs_coefficient = mean(abs_coefficient),
    max_abs_coefficient = max(abs_coefficient),
    mean_coefficient = mean(coefficient),
    folds = paste(sort(unique(outer_fold)), collapse = ";"),
    classes = paste(sort(unique(class)), collapse = ";"),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(nonzero_folds),
    dplyr::desc(mean_abs_coefficient),
    feature
  )

saveRDS(
  list(
    run_id = config$run_id,
    outcome_info = inp$outcome_info,
    covariate_set = inp$covariate_set,
    cv_strategy = "loso",
    metrics = pooled_eval$metrics,
    confusion_matrix = pooled_eval$confusion_matrix,
    by_class = pooled_eval$by_class,
    auc = pooled_eval$auc,
    fold_metrics = fold_metrics,
    predictions = predictions,
    tune_metric = tune_metric,
    alpha_grid = alpha_grid,
    selected_features_by_fold = fold_features,
    feature_frequency = feature_frequency,
    nonzero_coefficients = nonzero_coefficients,
    nonzero_feature_frequency = nonzero_feature_frequency
  ),
  file = run_file(out_paths$results_metrics, "elastic_net_loso_metrics")
)

readr::write_csv(
  pooled_eval$metrics,
  run_file(out_paths$results_metrics, "elastic_net_loso_metrics", "csv")
)

readr::write_csv(
  as.data.frame(pooled_eval$confusion_matrix$table),
  run_file(out_paths$results_metrics, "elastic_net_loso_confusion_matrix", "csv")
)

readr::write_csv(
  pooled_eval$by_class |>
    tibble::rownames_to_column("class"),
  run_file(out_paths$results_metrics, "elastic_net_loso_by_class_metrics", "csv")
)

if (!is.null(pooled_eval$auc)) {
  readr::write_csv(
    pooled_eval$auc,
    run_file(out_paths$results_metrics, "elastic_net_loso_auc", "csv")
  )
}

readr::write_csv(
  predictions,
  run_file(out_paths$results_metrics, "elastic_net_loso_predictions", "csv")
)

readr::write_csv(
  fold_metrics,
  run_file(out_paths$results_metrics, "elastic_net_loso_fold_metrics", "csv")
)

readr::write_csv(
  feature_frequency,
  run_file(out_paths$results_features, "elastic_net_loso_feature_frequency", "csv")
)

readr::write_csv(
  nonzero_coefficients,
  run_file(out_paths$results_features, "elastic_net_loso_nonzero_coefficients", "csv")
)

readr::write_csv(
  nonzero_feature_frequency,
  run_file(out_paths$results_features, "elastic_net_loso_nonzero_feature_frequency", "csv")
)

message("LOSO evaluation complete.")
message("Saved pooled metrics to: ", run_file(out_paths$results_metrics, "elastic_net_loso_metrics"))
message("Saved non-zero coefficients to: ", run_file(out_paths$results_features, "elastic_net_loso_nonzero_coefficients", "csv"))
print(pooled_eval$metrics)

if (!is.null(pooled_eval$auc)) {
  print(pooled_eval$auc)
}
