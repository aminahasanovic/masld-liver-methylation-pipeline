# 05_train_boosted_trees.R
# Run-specific XGBoost baseline for train/test or LOSO evaluation.

source("scripts/00_setup_paths.R")
source("scripts/02_feature_selection.R")

suppressPackageStartupMessages({
  library(caret)
  library(dplyr)
  library(Matrix)
  library(pROC)
  library(readr)
  library(tibble)
  library(xgboost)
  library(yardstick)
})

# ---- Runtime settings ----
parse_int_env <- function(name, default) {
  value <- Sys.getenv(name, "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }

  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }

  parsed
}

parse_int_vector_env <- function(name, default) {
  value <- Sys.getenv(name, "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }

  parsed <- suppressWarnings(as.integer(strsplit(value, ",", fixed = TRUE)[[1]]))
  parsed <- parsed[!is.na(parsed) & parsed > 0]

  if (length(parsed) == 0) {
    stop(name, " must contain at least one positive integer.", call. = FALSE)
  }

  unique(parsed)
}

parse_num_vector_env <- function(name, default) {
  value <- Sys.getenv(name, "")
  if (!nzchar(value)) {
    return(as.numeric(default))
  }

  parsed <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1]]))
  parsed <- parsed[!is.na(parsed)]

  if (length(parsed) == 0) {
    stop(name, " must contain at least one numeric value.", call. = FALSE)
  }

  unique(parsed)
}

xgb_cv_folds <- parse_int_env("CLASSIFIER_XGB_CV_FOLDS", 3)
xgb_threads <- parse_int_env("CLASSIFIER_XGB_THREADS", 2)
xgb_top_variable_n <- parse_int_env(
  "CLASSIFIER_XGB_TOP_VARIABLE_N",
  config$feature_selection$top_variable_n
)
xgb_tune_metric <- Sys.getenv("CLASSIFIER_XGB_TUNE_METRIC", "Balanced_Accuracy")
xgb_nrounds <- parse_int_vector_env("CLASSIFIER_XGB_NROUNDS", c(50L, 100L))
xgb_max_depth <- parse_int_vector_env("CLASSIFIER_XGB_MAX_DEPTH", c(2L, 3L))
xgb_eta <- parse_num_vector_env("CLASSIFIER_XGB_ETA", c(0.05, 0.1))
xgb_gamma <- parse_num_vector_env("CLASSIFIER_XGB_GAMMA", 0)
xgb_colsample <- parse_num_vector_env("CLASSIFIER_XGB_COLSAMPLE_BYTREE", 0.8)
xgb_min_child <- parse_num_vector_env("CLASSIFIER_XGB_MIN_CHILD_WEIGHT", c(1, 5))
xgb_subsample <- parse_num_vector_env("CLASSIFIER_XGB_SUBSAMPLE", 0.8)
save_fold_models <- tolower(Sys.getenv("CLASSIFIER_XGB_SAVE_FOLD_MODELS", "false")) %in%
  c("1", "true", "yes")

if (!xgb_tune_metric %in% c("Accuracy", "Balanced_Accuracy", "Macro_F1")) {
  stop(
    "CLASSIFIER_XGB_TUNE_METRIC must be Accuracy, Balanced_Accuracy, or Macro_F1.",
    call. = FALSE
  )
}

xgb_grid <- expand.grid(
  nrounds = xgb_nrounds,
  max_depth = xgb_max_depth,
  eta = xgb_eta,
  gamma = xgb_gamma,
  colsample_bytree = xgb_colsample,
  min_child_weight = xgb_min_child,
  subsample = xgb_subsample
)

safe_subset <- function(mat, i, j) {
  mat[i, j, drop = FALSE]
}

balanced_accuracy_manual <- function(obs, pred, class_levels) {
  obs <- factor(obs, levels = class_levels)
  pred <- factor(pred, levels = class_levels)
  tab <- table(obs, pred)
  recall <- diag(tab) / rowSums(tab)
  recall[!is.finite(recall)] <- NA_real_
  mean(recall, na.rm = TRUE)
}

macro_f1_manual <- function(obs, pred, class_levels) {
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

xgb_summary <- function(data, lev = NULL, model = NULL) {
  if (is.null(lev)) {
    lev <- levels(data$obs)
  }

  obs <- factor(data$obs, levels = lev)
  pred <- factor(data$pred, levels = lev)

  c(
    Accuracy = mean(obs == pred, na.rm = TRUE),
    Balanced_Accuracy = balanced_accuracy_manual(obs, pred, lev),
    Macro_F1 = macro_f1_manual(obs, pred, lev)
  )
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
    balanced_accuracy = suppressWarnings(yardstick::bal_accuracy_vec(y_true, pred_class)),
    macro_recall_present = balanced_accuracy_manual(y_true, pred_class, class_levels),
    mcc = yardstick::mcc_vec(y_true, pred_class),
    f1_macro = suppressWarnings(yardstick::f_meas_vec(y_true, pred_class, estimator = "macro")),
    f1_macro_present = macro_f1_manual(y_true, pred_class, class_levels)
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

# ---- Load prepared input ----
inp <- readRDS(check_file(run_file(out_paths$data_processed, "classifier_inputs")))
split <- readRDS(check_file(run_file(out_paths$data_processed, "cv_split")))

x_cpg <- inp$x_cpg
if (is.null(x_cpg)) {
  x_cpg <- readRDS(check_file(inp$x_cpg_path))
}

meta <- inp$meta
covar_cols <- inp$covar_cols %||% character()
categorical_covar_cols <- inp$categorical_covar_cols %||% covar_cols
numeric_covar_cols <- inp$numeric_covar_cols %||% character()
cv_strategy <- config$modeling$cv_strategy %||% "train_test"

if (!inherits(x_cpg, "Matrix")) {
  x_cpg <- Matrix::Matrix(x_cpg, sparse = TRUE)
}

meta <- meta |>
  dplyr::filter(sample_id %in% rownames(x_cpg)) |>
  dplyr::arrange(match(sample_id, rownames(x_cpg)))

x_cpg <- x_cpg[meta$sample_id, , drop = FALSE]
meta$stage <- droplevels(meta$stage)
class_levels <- levels(meta$stage)

stopifnot(
  nrow(x_cpg) == nrow(meta),
  identical(rownames(x_cpg), meta$sample_id)
)

if (length(covar_cols) > 0) {
  missing_covars <- setdiff(covar_cols, colnames(meta))

  if (length(missing_covars) > 0) {
    stop(
      "Missing covariate columns in metadata: ",
      paste(missing_covars, collapse = ", "),
      call. = FALSE
    )
  }
}

make_predictor_df <- function(x_mat, meta_df, categorical_levels, numeric_impute_values) {
  x_df <- as.data.frame(as.matrix(x_mat), check.names = FALSE)

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

  if (anyDuplicated(colnames(x_df)) > 0) {
    stop("Duplicate predictor names after adding covariates.", call. = FALSE)
  }

  x_df
}

encode_predictors <- function(train_df, test_df) {
  train_mm <- Matrix::sparse.model.matrix(~ . - 1, data = train_df)
  test_mm <- Matrix::sparse.model.matrix(~ . - 1, data = test_df)

  missing_in_test <- setdiff(colnames(train_mm), colnames(test_mm))
  if (length(missing_in_test) > 0) {
    zero_mm <- Matrix::Matrix(
      0,
      nrow = nrow(test_mm),
      ncol = length(missing_in_test),
      sparse = TRUE
    )
    colnames(zero_mm) <- missing_in_test
    rownames(zero_mm) <- rownames(test_mm)
    test_mm <- Matrix::cbind2(test_mm, zero_mm)
  }

  extra_in_test <- setdiff(colnames(test_mm), colnames(train_mm))
  if (length(extra_in_test) > 0) {
    test_mm <- test_mm[, setdiff(colnames(test_mm), extra_in_test), drop = FALSE]
  }

  test_mm <- test_mm[, colnames(train_mm), drop = FALSE]

  list(train = train_mm, test = test_mm)
}

prepare_xgb_labels <- function(y) {
  as.integer(factor(y, levels = class_levels)) - 1L
}

make_xgb_params <- function(grid_row) {
  params <- list(
    max_depth = as.integer(grid_row$max_depth),
    eta = as.numeric(grid_row$eta),
    gamma = as.numeric(grid_row$gamma),
    colsample_bytree = as.numeric(grid_row$colsample_bytree),
    min_child_weight = as.numeric(grid_row$min_child_weight),
    subsample = as.numeric(grid_row$subsample),
    nthread = xgb_threads
  )

  if (length(class_levels) == 2) {
    params$objective <- "binary:logistic"
    params$eval_metric <- "logloss"
  } else {
    params$objective <- "multi:softprob"
    params$num_class <- length(class_levels)
    params$eval_metric <- "mlogloss"
  }

  params
}

predict_xgb_prob <- function(booster, x_mat) {
  raw_pred <- predict(booster, newdata = x_mat)

  if (length(class_levels) == 2) {
    prob <- data.frame(
      stats::setNames(
        list(1 - raw_pred, raw_pred),
        class_levels
      ),
      check.names = FALSE
    )
  } else {
    # xgboost returns multi:softprob values grouped by class in the R API.
    prob_mat <- matrix(raw_pred, ncol = length(class_levels), byrow = FALSE)
    colnames(prob_mat) <- class_levels
    prob <- as.data.frame(prob_mat, check.names = FALSE)
  }

  prob
}

predict_xgb_class <- function(prob_df) {
  class_levels[max.col(as.matrix(prob_df), ties.method = "first")]
}

score_metric <- function(y_true, pred_class) {
  if (xgb_tune_metric == "Accuracy") {
    return(mean(factor(pred_class, levels = class_levels) == factor(y_true, levels = class_levels)))
  }

  if (xgb_tune_metric == "Macro_F1") {
    return(macro_f1_manual(y_true, pred_class, class_levels))
  }

  balanced_accuracy_manual(y_true, pred_class, class_levels)
}

is_covariate_feature <- function(feature) {
  if (length(covar_cols) == 0) {
    return(rep(FALSE, length(feature)))
  }

  vapply(feature, function(x) {
    any(vapply(covar_cols, function(prefix) startsWith(x, prefix), logical(1)))
  }, logical(1))
}

fit_xgb_fold <- function(train_idx, test_idx, fold_name = NA_character_, seed_offset = 0L) {
  y_train <- droplevels(meta$stage[train_idx])
  y_test <- factor(meta$stage[test_idx], levels = class_levels)

  if (nlevels(y_train) < 2) {
    stop("Training fold has fewer than 2 classes.", call. = FALSE)
  }

  missing_train_classes <- setdiff(class_levels, levels(y_train))
  if (length(missing_train_classes) > 0) {
    stop(
      "Training fold is missing class(es): ",
      paste(missing_train_classes, collapse = ", "),
      call. = FALSE
    )
  }

  selected_cpgs <- select_top_variable_cpgs(
    x_mat = safe_subset(x_cpg, train_idx, seq_len(ncol(x_cpg))),
    n_top = xgb_top_variable_n
  )
  selected_cpgs <- intersect(as.character(selected_cpgs), colnames(x_cpg))

  if (length(selected_cpgs) == 0) {
    stop("Feature selection returned no valid CpGs.", call. = FALSE)
  }

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

  x_train_df <- make_predictor_df(
    safe_subset(x_cpg, train_idx, selected_cpgs),
    meta[train_idx, , drop = FALSE],
    categorical_levels,
    numeric_impute_values
  )
  x_test_df <- make_predictor_df(
    safe_subset(x_cpg, test_idx, selected_cpgs),
    meta[test_idx, , drop = FALSE],
    categorical_levels,
    numeric_impute_values
  )

  encoded <- encode_predictors(x_train_df, x_test_df)
  x_train <- encoded$train
  x_test <- encoded$test

  min_class_n <- min(table(y_train))
  inner_folds <- min(xgb_cv_folds, as.integer(min_class_n))

  if (inner_folds < 2) {
    stop("Training fold has a class with fewer than 2 samples.", call. = FALSE)
  }

  set.seed(config$modeling$seed + seed_offset)
  inner_out <- caret::createFolds(y_train, k = inner_folds, list = TRUE)

  tune_results <- lapply(seq_len(nrow(xgb_grid)), function(grid_i) {
    grid_row <- xgb_grid[grid_i, , drop = FALSE]
    fold_scores <- vapply(seq_along(inner_out), function(inner_i) {
      valid_local <- inner_out[[inner_i]]
      train_local <- setdiff(seq_along(y_train), valid_local)

      dtrain <- xgboost::xgb.DMatrix(
        data = x_train[train_local, , drop = FALSE],
        label = prepare_xgb_labels(y_train[train_local])
      )

      booster <- xgboost::xgb.train(
        params = make_xgb_params(grid_row),
        data = dtrain,
        nrounds = as.integer(grid_row$nrounds),
        verbose = 0
      )

      prob_valid <- predict_xgb_prob(
        booster,
        x_train[valid_local, , drop = FALSE]
      )
      pred_valid <- predict_xgb_class(prob_valid)
      score_metric(y_train[valid_local], pred_valid)
    }, numeric(1))

    dplyr::bind_cols(
      grid_row,
      tibble::tibble(
        mean_metric = mean(fold_scores, na.rm = TRUE),
        sd_metric = stats::sd(fold_scores, na.rm = TRUE)
      )
    )
  }) |>
    dplyr::bind_rows() |>
    dplyr::arrange(dplyr::desc(mean_metric), nrounds, max_depth)

  if (nrow(tune_results) == 0 || all(is.na(tune_results$mean_metric))) {
    stop("All XGBoost tuning results are NA.", call. = FALSE)
  }

  best_tune <- tune_results |>
    dplyr::slice(1) |>
    dplyr::select(-mean_metric, -sd_metric)

  dtrain_final <- xgboost::xgb.DMatrix(
    data = x_train,
    label = prepare_xgb_labels(y_train)
  )

  booster <- xgboost::xgb.train(
    params = make_xgb_params(best_tune),
    data = dtrain_final,
    nrounds = as.integer(best_tune$nrounds),
    verbose = 0
  )

  fit <- list(
    booster = booster,
    bestTune = best_tune,
    tune_results = tune_results,
    feature_names = colnames(x_train)
  )

  pred_prob <- predict_xgb_prob(booster, x_test)
  pred_class <- predict_xgb_class(pred_prob)

  eval <- compute_metrics(
    y_true = y_test,
    pred_class = pred_class,
    pred_prob = pred_prob,
    class_levels = class_levels
  )

  predictions <- tibble::tibble(
    sample_id = meta$sample_id[test_idx],
    dataset = meta[[config$covariates$study_id_col]][test_idx],
    truth = y_test,
    prediction = factor(pred_class, levels = class_levels)
  ) |>
    dplyr::bind_cols(tibble::as_tibble(pred_prob))

  if (!is.na(fold_name)) {
    predictions <- predictions |>
      dplyr::mutate(outer_fold = fold_name, .before = 1)
  }

  importance_tbl <- tryCatch(
    {
      xgboost::xgb.importance(
        feature_names = colnames(x_train),
        model = booster
      ) |>
        tibble::as_tibble() |>
        dplyr::rename(feature = Feature) |>
        dplyr::mutate(
          feature_type = dplyr::case_when(
            feature %in% selected_cpgs ~ "CpG",
            is_covariate_feature(feature) ~ "covariate",
            TRUE ~ "other"
          ),
          .after = feature
        )
    },
    error = function(e) {
      warning("Could not extract XGBoost variable importance: ", conditionMessage(e))
      tibble::tibble()
    }
  )

  if (!is.na(fold_name) && nrow(importance_tbl) > 0) {
    importance_tbl <- importance_tbl |>
      dplyr::mutate(outer_fold = fold_name, .before = 1)
  }

  summary <- dplyr::bind_cols(
    tibble::tibble(
      outer_fold = fold_name,
      train_n = length(train_idx),
      test_n = length(test_idx),
      selected_cpgs_n = length(selected_cpgs),
      predictor_n = ncol(x_train),
      inner_folds = inner_folds,
      tune_metric = xgb_tune_metric,
      tune_metric_value = tune_results$mean_metric[1],
      min_class_n = min_class_n,
      train_class_counts = paste(names(table(y_train)), as.integer(table(y_train)), sep = "=", collapse = ";"),
      test_class_counts = paste(names(table(y_test)), as.integer(table(y_test)), sep = "=", collapse = ";")
    ),
    best_tune,
    eval$metrics
  )

  list(
    fit = fit,
    selected_cpgs = selected_cpgs,
    categorical_levels = categorical_levels,
    numeric_impute_values = numeric_impute_values,
    predictions = predictions,
    eval = eval,
    importance = importance_tbl,
    summary = summary
  )
}

message("Starting XGBoost boosted-tree evaluation.")
message("Run ID: ", config$run_id)
message("Outcome: ", inp$outcome_info$name)
message("Covariate set: ", inp$covariate_set)
message("CV strategy: ", cv_strategy)
message("Tune metric: ", xgb_tune_metric)
message("Top variable CpGs: ", xgb_top_variable_n)
message("Grid rows: ", nrow(xgb_grid))

if (cv_strategy == "train_test") {
  train_idx <- as.integer(split$train_idx)
  test_idx <- as.integer(split$test_idx)

  if (length(train_idx) == 0 || length(test_idx) == 0) {
    stop("train_test split must contain non-empty train_idx and test_idx.", call. = FALSE)
  }

  result <- fit_xgb_fold(train_idx, test_idx)

  saveRDS(
    result$fit,
    file = run_file(out_paths$results_models, "boosted_trees_model")
  )

  saveRDS(
    list(
      run_id = config$run_id,
      outcome_info = inp$outcome_info,
      covariate_set = inp$covariate_set,
      cv_strategy = "train_test",
      selected_cpgs = result$selected_cpgs,
      train_idx = train_idx,
      test_idx = test_idx,
      covar_cols = covar_cols,
      categorical_covar_cols = categorical_covar_cols,
      numeric_covar_cols = numeric_covar_cols,
      categorical_levels = result$categorical_levels,
      numeric_impute_values = result$numeric_impute_values,
      xgb_tune_metric = xgb_tune_metric,
      bestTune = result$fit$bestTune,
      variable_importance = result$importance
    ),
    file = run_file(out_paths$results_features, "boosted_trees_features")
  )

  saveRDS(
    list(
      run_id = config$run_id,
      outcome_info = inp$outcome_info,
      covariate_set = inp$covariate_set,
      cv_strategy = "train_test",
      bestTune = result$fit$bestTune,
      metrics = result$eval$metrics,
      confusion_matrix = result$eval$confusion_matrix,
      by_class = result$eval$by_class,
      auc = result$eval$auc,
      predictions = result$predictions
    ),
    file = run_file(out_paths$results_metrics, "boosted_trees_metrics")
  )

  readr::write_csv(
    result$eval$metrics,
    run_file(out_paths$results_metrics, "boosted_trees_metrics", "csv")
  )
  readr::write_csv(
    as.data.frame(result$eval$confusion_matrix$table),
    run_file(out_paths$results_metrics, "boosted_trees_confusion_matrix", "csv")
  )
  readr::write_csv(
    result$eval$by_class |>
      tibble::rownames_to_column("class"),
    run_file(out_paths$results_metrics, "boosted_trees_by_class_metrics", "csv")
  )
  if (!is.null(result$eval$auc)) {
    readr::write_csv(
      result$eval$auc,
      run_file(out_paths$results_metrics, "boosted_trees_auc", "csv")
    )
  }
  readr::write_csv(
    result$predictions,
    run_file(out_paths$results_metrics, "boosted_trees_predictions", "csv")
  )
  if (nrow(result$importance) > 0) {
    readr::write_csv(
      result$importance,
      run_file(out_paths$results_features, "boosted_trees_variable_importance", "csv")
    )
  }
  yaml::write_yaml(
    inp$outcome_info,
    file = run_file(out_paths$results_metrics, "boosted_trees_outcome_definition", "yaml")
  )

  message("XGBoost train/test evaluation complete.")
  print(result$eval$metrics)
  if (!is.null(result$eval$auc)) print(result$eval$auc)

} else if (cv_strategy == "loso") {
  if (is.null(split$loso_train) || is.null(split$loso_out)) {
    stop("cv_split does not contain loso_train/loso_out. Re-run 01_prepare_classifier_data.R with LOSO.", call. = FALSE)
  }

  fold_results <- vector("list", length(split$loso_out))
  names(fold_results) <- names(split$loso_out)

  for (fold_counter in seq_along(split$loso_out)) {
    fold_name <- names(split$loso_out)[[fold_counter]]
    message("LOSO fold held out: ", fold_name)
    fold_results[[fold_name]] <- fit_xgb_fold(
      train_idx = as.integer(split$loso_train[[fold_name]]),
      test_idx = as.integer(split$loso_out[[fold_name]]),
      fold_name = fold_name,
      seed_offset = fold_counter
    )
  }

  predictions <- dplyr::bind_rows(lapply(fold_results, `[[`, "predictions"))
  fold_metrics <- dplyr::bind_rows(lapply(fold_results, `[[`, "summary"))
  importance_all <- dplyr::bind_rows(lapply(fold_results, `[[`, "importance"))

  pooled_eval <- compute_metrics(
    y_true = predictions$truth,
    pred_class = predictions$prediction,
    pred_prob = predictions[, class_levels, drop = FALSE],
    class_levels = class_levels
  )

  fold_features <- lapply(fold_results, `[[`, "selected_cpgs")
  feature_frequency <- table(unlist(fold_features))
  feature_frequency <- tibble::tibble(
    cpg = names(feature_frequency),
    selected_folds = as.integer(feature_frequency)
  ) |>
    dplyr::arrange(dplyr::desc(selected_folds), cpg)

  importance_summary <- tibble::tibble()
  if (nrow(importance_all) > 0) {
    importance_value_cols <- setdiff(
      colnames(importance_all),
      c("outer_fold", "feature", "feature_type")
    )

    importance_summary <- importance_all |>
      dplyr::mutate(
        importance_mean_row = rowMeans(
          dplyr::across(all_of(importance_value_cols)),
          na.rm = TRUE
        )
      ) |>
      dplyr::group_by(feature, feature_type) |>
      dplyr::summarise(
        importance_folds = dplyr::n_distinct(outer_fold),
        mean_importance = mean(importance_mean_row, na.rm = TRUE),
        max_importance = max(importance_mean_row, na.rm = TRUE),
        folds = paste(sort(unique(outer_fold)), collapse = ";"),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(importance_folds), dplyr::desc(mean_importance), feature)
  }

  metrics_object <- list(
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
    xgb_tune_metric = xgb_tune_metric,
    xgb_grid = xgb_grid,
    selected_features_by_fold = fold_features,
    feature_frequency = feature_frequency,
    variable_importance = importance_all,
    variable_importance_summary = importance_summary
  )

  if (save_fold_models) {
    metrics_object$fold_models <- lapply(fold_results, `[[`, "fit")
  }

  saveRDS(
    metrics_object,
    file = run_file(out_paths$results_metrics, "boosted_trees_loso_metrics")
  )

  readr::write_csv(
    pooled_eval$metrics,
    run_file(out_paths$results_metrics, "boosted_trees_loso_metrics", "csv")
  )
  readr::write_csv(
    as.data.frame(pooled_eval$confusion_matrix$table),
    run_file(out_paths$results_metrics, "boosted_trees_loso_confusion_matrix", "csv")
  )
  readr::write_csv(
    pooled_eval$by_class |>
      tibble::rownames_to_column("class"),
    run_file(out_paths$results_metrics, "boosted_trees_loso_by_class_metrics", "csv")
  )
  if (!is.null(pooled_eval$auc)) {
    readr::write_csv(
      pooled_eval$auc,
      run_file(out_paths$results_metrics, "boosted_trees_loso_auc", "csv")
    )
  }
  readr::write_csv(
    predictions,
    run_file(out_paths$results_metrics, "boosted_trees_loso_predictions", "csv")
  )
  readr::write_csv(
    fold_metrics,
    run_file(out_paths$results_metrics, "boosted_trees_loso_fold_metrics", "csv")
  )
  readr::write_csv(
    feature_frequency,
    run_file(out_paths$results_features, "boosted_trees_loso_feature_frequency", "csv")
  )
  if (nrow(importance_all) > 0) {
    readr::write_csv(
      importance_all,
      run_file(out_paths$results_features, "boosted_trees_loso_variable_importance", "csv")
    )
  }
  if (nrow(importance_summary) > 0) {
    readr::write_csv(
      importance_summary,
      run_file(out_paths$results_features, "boosted_trees_loso_variable_importance_summary", "csv")
    )
  }
  yaml::write_yaml(
    inp$outcome_info,
    file = run_file(out_paths$results_metrics, "boosted_trees_loso_outcome_definition", "yaml")
  )

  message("XGBoost LOSO evaluation complete.")
  print(pooled_eval$metrics)
  if (!is.null(pooled_eval$auc)) print(pooled_eval$auc)

} else {
  stop("05_train_boosted_trees.R supports cv_strategy train_test or loso.", call. = FALSE)
}
