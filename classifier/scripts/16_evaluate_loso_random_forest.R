# 16_evaluate_loso_random_forest.R
# Leave-one-study-out Random Forest evaluation with fold-wise feature selection.

source("scripts/00_setup_paths.R")
source("scripts/02_feature_selection.R")

suppressPackageStartupMessages({
  library(caret)
  library(dplyr)
  library(Matrix)
  library(pROC)
  library(ranger)
  library(readr)
  library(tibble)
  library(yardstick)
})

if ((config$modeling$cv_strategy %||% "") != "loso") {
  stop(
    "Set CLASSIFIER_CV_STRATEGY=loso before running this script.",
    call. = FALSE
  )
}

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

parse_int_vector_env <- function(name, default = NULL) {
  value <- Sys.getenv(name, "")
  if (!nzchar(value)) {
    return(default)
  }

  parsed <- suppressWarnings(as.integer(strsplit(value, ",", fixed = TRUE)[[1]]))
  parsed <- parsed[!is.na(parsed) & parsed > 0]

  if (length(parsed) == 0) {
    stop(name, " must contain at least one positive integer.", call. = FALSE)
  }

  unique(parsed)
}

rf_num_trees <- parse_int_env("CLASSIFIER_RF_NUM_TREES", 500)
rf_threads <- parse_int_env("CLASSIFIER_RF_THREADS", 2)
rf_cv_folds <- parse_int_env("CLASSIFIER_RF_CV_FOLDS", 5)
rf_importance <- Sys.getenv("CLASSIFIER_RF_IMPORTANCE", "impurity")
rf_tune_metric <- Sys.getenv("CLASSIFIER_RF_TUNE_METRIC", "Balanced_Accuracy")
rf_mtry_grid_env <- parse_int_vector_env("CLASSIFIER_RF_MTRY_GRID", NULL)
rf_min_node_size_grid <- parse_int_vector_env(
  "CLASSIFIER_RF_MIN_NODE_SIZE_GRID",
  c(1L, 5L, 10L)
)
save_fold_models <- tolower(Sys.getenv("CLASSIFIER_RF_SAVE_FOLD_MODELS", "false")) %in%
  c("1", "true", "yes")

if (!rf_tune_metric %in% c("Accuracy", "Balanced_Accuracy", "Macro_F1")) {
  stop(
    "CLASSIFIER_RF_TUNE_METRIC must be Accuracy, Balanced_Accuracy, or Macro_F1.",
    call. = FALSE
  )
}

build_mtry_grid <- function(n_features) {
  if (!is.null(rf_mtry_grid_env)) {
    return(unique(pmax(1L, pmin(n_features, rf_mtry_grid_env))))
  }

  unique(pmax(
    1L,
    pmin(
      n_features,
      as.integer(round(c(sqrt(n_features), n_features / 20, n_features / 10)))
    )
  ))
}

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

rf_summary <- function(data, lev = NULL, model = NULL) {
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

if (is.null(split$loso_train) || is.null(split$loso_out)) {
  stop("cv_split does not contain loso_train/loso_out. Re-run 01_prepare_classifier_data.R with LOSO.", call. = FALSE)
}

x_cpg <- inp$x_cpg
if (is.null(x_cpg)) {
  x_cpg <- readRDS(check_file(inp$x_cpg_path))
}

meta <- inp$meta
covar_cols <- inp$covar_cols %||% character()
categorical_covar_cols <- inp$categorical_covar_cols %||% covar_cols
numeric_covar_cols <- inp$numeric_covar_cols %||% character()

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

message("Starting LOSO Random Forest evaluation.")
message("Run ID: ", config$run_id)
message("Outcome: ", inp$outcome_info$name)
message("Covariate set: ", inp$covariate_set)
message("Tune metric: ", rf_tune_metric)
message("Trees: ", rf_num_trees)

fold_predictions <- list()
fold_summaries <- list()
fold_features <- list()
fold_importance <- list()
fold_models <- list()

for (fold_counter in seq_along(split$loso_out)) {
  fold_name <- names(split$loso_out)[[fold_counter]]
  message("LOSO fold held out: ", fold_name)

  train_idx <- as.integer(split$loso_train[[fold_name]])
  test_idx <- as.integer(split$loso_out[[fold_name]])

  y_train <- droplevels(meta$stage[train_idx])
  y_test <- factor(meta$stage[test_idx], levels = class_levels)

  if (nlevels(y_train) < 2) {
    stop("Training fold for ", fold_name, " has fewer than 2 classes.", call. = FALSE)
  }

  missing_train_classes <- setdiff(class_levels, levels(y_train))
  if (length(missing_train_classes) > 0) {
    stop(
      "Training fold for ", fold_name,
      " is missing class(es): ", paste(missing_train_classes, collapse = ", "),
      call. = FALSE
    )
  }

  selected_cpgs <- select_top_variable_cpgs(
    x_mat = safe_subset(x_cpg, train_idx, seq_len(ncol(x_cpg))),
    n_top = config$feature_selection$top_variable_n
  )
  selected_cpgs <- intersect(as.character(selected_cpgs), colnames(x_cpg))

  if (length(selected_cpgs) == 0) {
    stop("Feature selection returned no valid CpGs for fold ", fold_name, ".", call. = FALSE)
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

  x_train <- make_predictor_df(
    safe_subset(x_cpg, train_idx, selected_cpgs),
    meta[train_idx, , drop = FALSE],
    categorical_levels,
    numeric_impute_values
  )
  x_test <- make_predictor_df(
    safe_subset(x_cpg, test_idx, selected_cpgs),
    meta[test_idx, , drop = FALSE],
    categorical_levels,
    numeric_impute_values
  )

  min_class_n <- min(table(y_train))
  inner_folds <- min(rf_cv_folds, as.integer(min_class_n))

  if (inner_folds < 2) {
    stop("Training fold for ", fold_name, " has a class with fewer than 2 samples.", call. = FALSE)
  }

  rf_grid <- expand.grid(
    mtry = build_mtry_grid(ncol(x_train)),
    splitrule = "gini",
    min.node.size = unique(pmax(1L, pmin(nrow(x_train), rf_min_node_size_grid)))
  )

  ctrl <- caret::trainControl(
    method = "cv",
    number = inner_folds,
    classProbs = TRUE,
    savePredictions = "final",
    summaryFunction = rf_summary
  )

  set.seed(config$modeling$seed + fold_counter)
  fit <- caret::train(
    x = x_train,
    y = y_train,
    method = "ranger",
    trControl = ctrl,
    tuneGrid = rf_grid,
    metric = rf_tune_metric,
    num.trees = rf_num_trees,
    importance = rf_importance,
    num.threads = rf_threads
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

  importance_tbl <- tryCatch(
    {
      caret::varImp(fit, scale = FALSE)$importance |>
        as.data.frame() |>
        tibble::rownames_to_column("feature") |>
        dplyr::mutate(
          outer_fold = fold_name,
          feature_type = dplyr::case_when(
            feature %in% selected_cpgs ~ "CpG",
            feature %in% covar_cols ~ "covariate",
            TRUE ~ "other"
          ),
          .before = 1
        )
    },
    error = function(e) {
      warning("Could not extract RF variable importance for ", fold_name, ": ", conditionMessage(e))
      tibble::tibble()
    }
  )

  fold_predictions[[fold_name]] <- pred_tbl
  fold_features[[fold_name]] <- selected_cpgs
  fold_importance[[fold_name]] <- importance_tbl
  if (save_fold_models) {
    fold_models[[fold_name]] <- fit
  }

  fold_summaries[[fold_name]] <- dplyr::bind_cols(
    tibble::tibble(
      outer_fold = fold_name,
      train_n = length(train_idx),
      test_n = length(test_idx),
      selected_cpgs_n = length(selected_cpgs),
      predictor_n = ncol(x_train),
      inner_folds = inner_folds,
      tune_metric = rf_tune_metric,
      rf_num_trees = rf_num_trees,
      min_class_n = min_class_n,
      mtry = fit$bestTune$mtry,
      splitrule = fit$bestTune$splitrule,
      min.node.size = fit$bestTune$min.node.size,
      train_class_counts = paste(names(table(y_train)), as.integer(table(y_train)), sep = "=", collapse = ";"),
      test_class_counts = paste(names(table(y_test)), as.integer(table(y_test)), sep = "=", collapse = ";")
    ),
    fold_eval$metrics
  )
}

predictions <- dplyr::bind_rows(fold_predictions)
fold_metrics <- dplyr::bind_rows(fold_summaries)
importance_all <- dplyr::bind_rows(fold_importance)

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
  rf_tune_metric = rf_tune_metric,
  rf_num_trees = rf_num_trees,
  rf_importance = rf_importance,
  selected_features_by_fold = fold_features,
  feature_frequency = feature_frequency,
  variable_importance = importance_all,
  variable_importance_summary = importance_summary
)

if (save_fold_models) {
  metrics_object$fold_models <- fold_models
}

saveRDS(
  metrics_object,
  file = run_file(out_paths$results_metrics, "random_forest_loso_metrics")
)

readr::write_csv(
  pooled_eval$metrics,
  run_file(out_paths$results_metrics, "random_forest_loso_metrics", "csv")
)

readr::write_csv(
  as.data.frame(pooled_eval$confusion_matrix$table),
  run_file(out_paths$results_metrics, "random_forest_loso_confusion_matrix", "csv")
)

readr::write_csv(
  pooled_eval$by_class |>
    tibble::rownames_to_column("class"),
  run_file(out_paths$results_metrics, "random_forest_loso_by_class_metrics", "csv")
)

if (!is.null(pooled_eval$auc)) {
  readr::write_csv(
    pooled_eval$auc,
    run_file(out_paths$results_metrics, "random_forest_loso_auc", "csv")
  )
}

readr::write_csv(
  predictions,
  run_file(out_paths$results_metrics, "random_forest_loso_predictions", "csv")
)

readr::write_csv(
  fold_metrics,
  run_file(out_paths$results_metrics, "random_forest_loso_fold_metrics", "csv")
)

readr::write_csv(
  feature_frequency,
  run_file(out_paths$results_features, "random_forest_loso_feature_frequency", "csv")
)

if (nrow(importance_all) > 0) {
  readr::write_csv(
    importance_all,
    run_file(out_paths$results_features, "random_forest_loso_variable_importance", "csv")
  )
}

if (nrow(importance_summary) > 0) {
  readr::write_csv(
    importance_summary,
    run_file(out_paths$results_features, "random_forest_loso_variable_importance_summary", "csv")
  )
}

yaml::write_yaml(
  inp$outcome_info,
  file = run_file(out_paths$results_metrics, "random_forest_loso_outcome_definition", "yaml")
)

message("LOSO Random Forest evaluation complete.")
message("Saved pooled metrics to: ", run_file(out_paths$results_metrics, "random_forest_loso_metrics"))
message("Saved fold metrics to: ", run_file(out_paths$results_metrics, "random_forest_loso_fold_metrics", "csv"))
print(pooled_eval$metrics)

if (!is.null(pooled_eval$auc)) {
  print(pooled_eval$auc)
}
