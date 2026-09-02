# 00_setup_paths.R
# Setup paths, config, and common helpers

suppressPackageStartupMessages({
  library(fs)
  library(yaml)
  library(glue)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    y
  } else {
    x
  }
}

sanitize_id <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

is_absolute_path <- function(path_in) {
  grepl("^(/|[A-Za-z]:[/\\\\])", path_in)
}

resolve_existing_config <- function(path_in) {
  if (fs::file_exists(path_in)) {
    return(path_in)
  }

  alt_path <- fs::path("classifier", path_in)
  if (fs::file_exists(alt_path)) {
    return(alt_path)
  }

  stop(glue("Missing classifier config: {path_in}"), call. = FALSE)
}

resolve_path <- function(path_in, base_dir, must_work = FALSE) {
  if (is.null(path_in) || identical(path_in, "")) {
    return(path_in)
  }
  if (is_absolute_path(path_in)) {
    return(normalizePath(path_in, mustWork = must_work))
  }
  normalizePath(fs::path(base_dir, path_in), mustWork = must_work)
}

config_path <- resolve_existing_config(Sys.getenv("CLASSIFIER_CONFIG", "config/classifier_config.yaml"))
config <- yaml::read_yaml(config_path)
config_dir <- dirname(normalizePath(config_path, mustWork = TRUE))
classifier_root <- normalizePath(file.path(config_dir, ".."), mustWork = TRUE)

env_outcome <- Sys.getenv("CLASSIFIER_OUTCOME", "")
if (nzchar(env_outcome)) {
  config$outcome$active <- env_outcome
}

env_covariate_set <- Sys.getenv("CLASSIFIER_COVARIATE_SET", "")
if (nzchar(env_covariate_set)) {
  config$modeling$covariate_set <- env_covariate_set
}

env_cv_strategy <- Sys.getenv("CLASSIFIER_CV_STRATEGY", "")
if (nzchar(env_cv_strategy)) {
  config$modeling$cv_strategy <- env_cv_strategy
}

env_tune_metric <- Sys.getenv("CLASSIFIER_TUNE_METRIC", "")
if (nzchar(env_tune_metric)) {
  config$modeling$tune_metric <- env_tune_metric
}

env_alpha_grid <- Sys.getenv("CLASSIFIER_ALPHA_GRID", "")
if (nzchar(env_alpha_grid)) {
  config$elastic_net$alpha_grid <- strsplit(env_alpha_grid, ",", fixed = TRUE)[[1]]
}

env_run_id <- Sys.getenv("CLASSIFIER_RUN_ID", "")
if (nzchar(env_run_id)) {
  config$run_id <- env_run_id
}

env_holdout <- Sys.getenv("CLASSIFIER_HOLDOUT", "")
if (nzchar(env_holdout)) {
  config$validation_holdout$active <- env_holdout
}

project_root <- config$project_root
if (is.null(project_root) || identical(project_root, "")) {
  project_root <- classifier_root
} else {
  project_root <- resolve_path(project_root, classifier_root, must_work = FALSE)
}
config$project_root <- project_root

input_keys <- c("beta_matrix_path", "metadata_path", "shared_cpgs_path")
for (input_key in intersect(input_keys, names(config$inputs))) {
  config$inputs[[input_key]] <- resolve_path(config$inputs[[input_key]], project_root, must_work = FALSE)
}

active_outcome <- config$outcome$active %||% config$outcome$outcome_variant %||% "disease_group_4"
active_covariate_set <- config$modeling$covariate_set %||% "all_covariates"
active_cv_strategy <- config$modeling$cv_strategy %||% "kfold"
active_holdout <- config$validation_holdout$active %||% ""

if (is.null(config$run_id) || identical(config$run_id, "")) {
  run_parts <- c(active_outcome, active_covariate_set, active_cv_strategy)
  if (nzchar(active_holdout)) {
    run_parts <- c(run_parts, paste0(active_holdout, "holdout"))
  }
  config$run_id <- sanitize_id(paste(run_parts, collapse = "_"))
} else {
  config$run_id <- sanitize_id(config$run_id)
}

config$active_outcome <- active_outcome
config$active_covariate_set <- active_covariate_set
config$active_cv_strategy <- active_cv_strategy
config$active_holdout <- active_holdout

# Output paths
out_paths <- list(
  data_processed = path(project_root, config$outputs$data_processed),
  data_intermediate = path(project_root, config$outputs$data_intermediate),
  results_models = path(project_root, config$outputs$results_models),
  results_metrics = path(project_root, config$outputs$results_metrics),
  results_features = path(project_root, config$outputs$results_features),
  results_plots = path(project_root, config$outputs$results_plots),
  logs = path(project_root, config$outputs$logs)
)

# Ensure directories exist
purrr::walk(out_paths, fs::dir_create)

# Helpers
check_file <- function(path_in) {
  if (is.null(path_in) || path_in == "" || !fs::file_exists(path_in)) {
    stop(glue("Missing file: {path_in}"))
  }
  path_in
}

run_file <- function(base_dir, stem, ext = "rds") {
  file.path(base_dir, paste0(stem, "_", config$run_id, ".", ext))
}

set.seed(config$modeling$seed)
