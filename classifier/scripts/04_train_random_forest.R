# 04_train_random_forest.R
# Train a run-specific Random Forest baseline with ranger.
#
# This script is intended for train/test runs. LOSO evaluation uses
# scripts/16_evaluate_loso_random_forest.R because feature selection and
# tuning must be repeated inside each held-out-study fold.

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

if ((config$modeling$cv_strategy %||% "") != "train_test") {
  stop(
    "04_train_random_forest.R is for train_test runs. ",
    "Use scripts/16_evaluate_loso_random_forest.R for LOSO.",
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

train_idx <- as.integer(split$train_idx)
test_idx <- as.integer(split$test_idx)

if (length(train_idx) == 0 || length(test_idx) == 0) {
  stop("train_test split must contain non-empty train_idx and test_idx.", call. = FALSE)
}

y_train <- droplevels(meta$stage[train_idx])
y_test <- factor(meta$stage[test_idx], levels = class_levels)

missing_train_classes <- setdiff(class_levels, levels(y_train))
if (length(missing_train_classes) > 0) {
  stop(
    "Training split is missing class(es): ",
    paste(missing_train_classes, collapse = ", "),
    call. = FALSE
  )
}

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

message("Selecting top variable CpGs on training samples only...")
selected_cpgs <- select_top_variable_cpgs(
  x_mat = safe_subset(x_cpg, train_idx, seq_len(ncol(x_cpg))),
  n_top = config$feature_selection$top_variable_n
)
selected_cpgs <- intersect(as.character(selected_cpgs), colnames(x_cpg))

if (length(selected_cpgs) == 0) {
  stop("Feature selection returned no valid CpG column names.", call. = FALSE)
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

make_predictor_df <- function(x_mat, meta_df) {
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

message("Building Random Forest training and test matrices...")
x_train <- make_predictor_df(
  safe_subset(x_cpg, train_idx, selected_cpgs),
  meta[train_idx, , drop = FALSE]
)
x_test <- make_predictor_df(
  safe_subset(x_cpg, test_idx, selected_cpgs),
  meta[test_idx, , drop = FALSE]
)

min_class_n <- min(table(y_train))
inner_folds <- min(rf_cv_folds, as.integer(min_class_n))

if (inner_folds < 2) {
  stop("Training split has a class with fewer than 2 samples.", call. = FALSE)
}

mtry_grid <- build_mtry_grid(ncol(x_train))
rf_grid <- expand.grid(
  mtry = mtry_grid,
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

message("Starting Random Forest training...")
message("Run ID: ", config$run_id)
message("Outcome: ", inp$outcome_info$name)
message("Covariate set: ", inp$covariate_set)
message("Selected CpGs: ", length(selected_cpgs))
message("Predictor columns including covariates: ", ncol(x_train))
message("Trees: ", rf_num_trees)
message("Inner CV folds: ", inner_folds)
message("Tune metric: ", rf_tune_metric)
message("Mtry grid: ", paste(unique(rf_grid$mtry), collapse = ", "))
message("Min node size grid: ", paste(unique(rf_grid$min.node.size), collapse = ", "))

set.seed(config$modeling$seed)
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

message("Creating held-out test predictions...")
pred_class <- predict(fit, newdata = x_test)
pred_prob <- predict(fit, newdata = x_test, type = "prob")

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

importance_tbl <- tryCatch(
  {
    caret::varImp(fit, scale = FALSE)$importance |>
      as.data.frame() |>
      tibble::rownames_to_column("feature") |>
      dplyr::mutate(
        feature_type = dplyr::case_when(
          feature %in% selected_cpgs ~ "CpG",
          feature %in% covar_cols ~ "covariate",
          TRUE ~ "other"
        ),
        .after = feature
      ) |>
      dplyr::arrange(dplyr::desc(dplyr::coalesce(.data$Overall, 0)))
  },
  error = function(e) {
    warning("Could not extract Random Forest variable importance: ", conditionMessage(e))
    tibble::tibble()
  }
)

# ---- Save outputs ----
saveRDS(
  fit,
  file = run_file(out_paths$results_models, "random_forest_model")
)

saveRDS(
  list(
    run_id = config$run_id,
    outcome_info = inp$outcome_info,
    covariate_set = inp$covariate_set,
    cv_strategy = "train_test",
    selected_cpgs = selected_cpgs,
    train_idx = train_idx,
    test_idx = test_idx,
    covar_cols = covar_cols,
    categorical_covar_cols = categorical_covar_cols,
    numeric_covar_cols = numeric_covar_cols,
    categorical_levels = categorical_levels,
    numeric_impute_values = numeric_impute_values,
    rf_num_trees = rf_num_trees,
    rf_threads = rf_threads,
    rf_importance = rf_importance,
    rf_tune_metric = rf_tune_metric,
    bestTune = fit$bestTune,
    variable_importance = importance_tbl
  ),
  file = run_file(out_paths$results_features, "random_forest_features")
)

saveRDS(
  list(
    run_id = config$run_id,
    outcome_info = inp$outcome_info,
    covariate_set = inp$covariate_set,
    cv_strategy = "train_test",
    bestTune = fit$bestTune,
    metrics = eval$metrics,
    confusion_matrix = eval$confusion_matrix,
    by_class = eval$by_class,
    auc = eval$auc,
    predictions = predictions
  ),
  file = run_file(out_paths$results_metrics, "random_forest_metrics")
)

readr::write_csv(
  eval$metrics,
  run_file(out_paths$results_metrics, "random_forest_metrics", "csv")
)

readr::write_csv(
  as.data.frame(eval$confusion_matrix$table),
  run_file(out_paths$results_metrics, "random_forest_confusion_matrix", "csv")
)

readr::write_csv(
  eval$by_class |>
    tibble::rownames_to_column("class"),
  run_file(out_paths$results_metrics, "random_forest_by_class_metrics", "csv")
)

if (!is.null(eval$auc)) {
  readr::write_csv(
    eval$auc,
    run_file(out_paths$results_metrics, "random_forest_auc", "csv")
  )
}

readr::write_csv(
  predictions,
  run_file(out_paths$results_metrics, "random_forest_predictions", "csv")
)

if (nrow(importance_tbl) > 0) {
  readr::write_csv(
    importance_tbl,
    run_file(out_paths$results_features, "random_forest_variable_importance", "csv")
  )
}

yaml::write_yaml(
  inp$outcome_info,
  file = run_file(out_paths$results_metrics, "random_forest_outcome_definition", "yaml")
)

message("Random Forest train/test evaluation complete.")
message("Saved model to: ", run_file(out_paths$results_models, "random_forest_model"))
message("Saved metrics to: ", run_file(out_paths$results_metrics, "random_forest_metrics"))
message("Saved feature info to: ", run_file(out_paths$results_features, "random_forest_features"))
print(eval$metrics)

if (!is.null(eval$auc)) {
  print(eval$auc)
}
