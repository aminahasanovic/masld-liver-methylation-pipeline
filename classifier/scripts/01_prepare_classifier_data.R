# 01_prepare_classifier_data.R
# Load beta matrix + metadata, harmonize CpGs, build modeling dataset

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(janitor)
  library(caret)
  library(Matrix)
  library(glue)
})

# ---- Load inputs ----
beta_path <- check_file(config$inputs$beta_matrix_path)
meta_path <- check_file(config$inputs$metadata_path)
shared_cpgs_path <- config$inputs$shared_cpgs_path
sample_id_col <- config$inputs$sample_id_col %||% "Sample_Name"

# Expectation:
# beta matrix = CpGs x Samples
# rownames(beta) = CpG IDs
# colnames(beta) = Sample IDs
beta <- readRDS(beta_path)
meta <- readRDS(meta_path)

# Phenotype tables produced before the cohort was renamed carry the previous
# label for the ITEN dataset. Relabel it here, where the study column enters the
# classifier chain, so study names in the candidate outputs match the thesis.
study_id_col <- config$covariates$study_id_col %||% "Dataset"
if (study_id_col %in% colnames(meta)) {
  old_label <- as.character(meta[[study_id_col]]) == "Kim"
  if (any(old_label, na.rm = TRUE)) {
    message(
      "Phenotype table uses the previous label 'Kim' for the ITEN cohort; relabelling ",
      sum(old_label, na.rm = TRUE), " samples."
    )
    meta[[study_id_col]] <- as.character(meta[[study_id_col]])
    meta[[study_id_col]][old_label] <- "ITEN"
  }
}

# ---- Basic checks ----
if (is.null(rownames(beta))) {
  stop("Beta matrix has no rownames. Expected CpG IDs as rownames.")
}

if (is.null(colnames(beta))) {
  stop("Beta matrix has no colnames. Expected sample IDs as colnames.")
}

if (!sample_id_col %in% colnames(meta)) {
  stop(glue::glue(
    "Configured sample_id_col '{sample_id_col}' not found in metadata. Available columns: {paste(colnames(meta), collapse = ', ')}"
  ))
}

if (anyDuplicated(colnames(beta)) > 0) {
  stop("Duplicate sample IDs found in beta matrix colnames.")
}

# ---- Optional shared CpG list ----
if (!is.null(shared_cpgs_path) && shared_cpgs_path != "") {
  shared_cpgs <- readr::read_lines(check_file(shared_cpgs_path))
} else {
  shared_cpgs <- rownames(beta)
}

# ---- Harmonize CpGs ----
common_cpgs <- intersect(rownames(beta), shared_cpgs)

if (length(common_cpgs) == 0) {
  stop("No shared CpGs found between beta matrix and provided shared_cpgs list.")
}

beta <- beta[common_cpgs, , drop = FALSE]

# ---- Match samples ----
# Normalize sample IDs in metadata to allow matching if needed
if (!"normalized_sample_id" %in% colnames(meta)) {
  meta <- meta %>%
    mutate(normalized_sample_id = as.character(.data[[sample_id_col]]))
}

if (anyDuplicated(meta$normalized_sample_id) > 0) {
  dup_ids <- unique(meta$normalized_sample_id[duplicated(meta$normalized_sample_id)])
  stop(glue::glue(
    "Duplicate sample IDs found in metadata normalized_sample_id. Examples: {paste(head(dup_ids, 10), collapse = ', ')}"
  ))
}

common_samples <- intersect(colnames(beta), meta$normalized_sample_id)

if (length(common_samples) == 0) {
  stop("No overlapping samples between beta matrix and metadata. Check sample_id_col and sample names.")
}

beta <- beta[, common_samples, drop = FALSE]

meta <- meta %>%
  filter(normalized_sample_id %in% common_samples) %>%
  mutate(sample_id = normalized_sample_id) %>%
  arrange(match(sample_id, common_samples))

stopifnot(identical(colnames(beta), meta$sample_id))

# ---- Outcome ----
outcome_name <- config$active_outcome
outcome_variants <- config$outcome$variants

if (is.null(outcome_variants) || !outcome_name %in% names(outcome_variants)) {
  stop(glue::glue(
    "Outcome variant '{outcome_name}' is not defined in config$outcome$variants."
  ))
}

outcome_variant <- outcome_variants[[outcome_name]]
stage_col <- outcome_variant$source_column %||%
  config$outcome$source_column %||%
  config$outcome$stage_column

if (!stage_col %in% colnames(meta)) {
  stop(glue::glue(
    "Configured outcome source column '{stage_col}' not found in metadata. Available columns: {paste(colnames(meta), collapse = ', ')}"
  ))
}

map_outcome <- function(raw_stage, mapping) {
  raw_stage <- as.character(raw_stage)
  mapped_stage <- rep(NA_character_, length(raw_stage))

  for (target_stage in names(mapping)) {
    source_values <- unname(unlist(mapping[[target_stage]]))
    mapped_stage[raw_stage %in% source_values] <- target_stage
  }

  mapped_stage
}

stage_raw <- as.character(meta[[stage_col]])
stage_mapped <- map_outcome(stage_raw, outcome_variant$mapping)

unmapped_stages <- sort(unique(stage_raw[!is.na(stage_raw) & is.na(stage_mapped)]))

if (length(unmapped_stages) > 0) {
  warning(glue::glue(
    "Outcome variant '{outcome_name}' does not map these source values and they will be removed: {paste(unmapped_stages, collapse = ', ')}"
  ))
}

if (!is.null(outcome_variant$levels) && length(outcome_variant$levels) > 0) {
  meta$stage <- factor(stage_mapped, levels = outcome_variant$levels)
} else {
  meta$stage <- factor(stage_mapped)
}

outcome_info <- list(
  name = outcome_name,
  source_column = stage_col,
  description = outcome_variant$description %||% "",
  levels = levels(meta$stage),
  mapping = outcome_variant$mapping
)

# Remove samples with undefined outcome
if (any(is.na(meta$stage))) {
  n_missing_stage <- sum(is.na(meta$stage))
  warning(glue::glue("Removing {n_missing_stage} samples with NA stage after factor level assignment."))
  meta <- meta %>% filter(!is.na(stage))
  beta <- beta[, meta$sample_id, drop = FALSE]
}

meta$stage <- droplevels(meta$stage)

if (nlevels(meta$stage) < 2) {
  stop("Outcome stage has fewer than 2 classes after filtering.")
}

meta_shared <- meta
beta_shared <- beta

# ---- Optional validation holdout ----
holdout_info <- NULL
holdout_meta <- meta[0, , drop = FALSE]

holdout_name <- config$active_holdout %||% ""
if (nzchar(holdout_name)) {
  holdout_defs <- config$validation_holdout$definitions %||% list()

  if (!holdout_name %in% names(holdout_defs)) {
    stop(glue::glue(
      "Validation holdout '{holdout_name}' is not defined in config$validation_holdout$definitions."
    ))
  }

  holdout_def <- holdout_defs[[holdout_name]]
  holdout_dataset_col <- holdout_def$dataset_col %||% config$covariates$study_id_col
  holdout_source_col <- holdout_def$source_column %||% stage_col

  missing_holdout_cols <- setdiff(
    c(holdout_dataset_col, holdout_source_col),
    colnames(meta)
  )

  if (length(missing_holdout_cols) > 0) {
    stop(glue::glue(
      "Validation holdout '{holdout_name}' needs missing metadata columns: {paste(missing_holdout_cols, collapse = ', ')}"
    ))
  }

  holdout_dataset_values <- as.character(unlist(holdout_def$dataset_values))
  holdout_source_values <- as.character(unlist(holdout_def$source_values %||% character()))

  holdout_mask <- as.character(meta[[holdout_dataset_col]]) %in% holdout_dataset_values
  if (length(holdout_source_values) > 0) {
    holdout_mask <- holdout_mask &
      as.character(meta[[holdout_source_col]]) %in% holdout_source_values
  }

  holdout_meta <- meta[holdout_mask, , drop = FALSE]

  if (nrow(holdout_meta) == 0) {
    stop(glue::glue(
      "Validation holdout '{holdout_name}' matched 0 samples. Check dataset/source values."
    ))
  }

  meta <- meta[!holdout_mask, , drop = FALSE]
  beta <- beta[, meta$sample_id, drop = FALSE]
  meta$stage <- droplevels(meta$stage)

  if (nlevels(meta$stage) < 2) {
    stop("Outcome stage has fewer than 2 classes after validation holdout removal.")
  }

  missing_classes_after_holdout <- setdiff(outcome_info$levels, levels(meta$stage))
  if (length(missing_classes_after_holdout) > 0) {
    stop(glue::glue(
      "Validation holdout '{holdout_name}' removed all training samples for class(es): {paste(missing_classes_after_holdout, collapse = ', ')}"
    ))
  }

  holdout_info <- list(
    name = holdout_name,
    description = holdout_def$description %||% "",
    dataset_col = holdout_dataset_col,
    dataset_values = holdout_dataset_values,
    source_column = holdout_source_col,
    source_values = holdout_source_values,
    n = nrow(holdout_meta),
    sample_ids = holdout_meta$sample_id,
    mapped_stage_counts = as.list(table(holdout_meta$stage)),
    training_stage_counts_after_removal = as.list(table(meta$stage))
  )
}

# ---- Covariates ----
covariate_set <- config$active_covariate_set
covariate_keys <- switch(
  covariate_set,
  cpg_only = character(),
  clinical = c("age", "sex"),
  study_array = c("study", "array"),
  all_covariates = c("study", "array", "age", "sex"),
  stop(glue::glue("Unknown modeling covariate_set: {covariate_set}"))
)

covariate_lookup <- list(
  study = config$covariates$study_id_col,
  array = config$covariates$array_type_col,
  age = config$covariates$age_col,
  sex = config$covariates$sex_col
)

covar_cols <- unname(unlist(covariate_lookup[covariate_keys]))
covar_cols <- covar_cols[!is.na(covar_cols) & covar_cols != ""]

numeric_covar_cols <- unname(unlist(covariate_lookup[intersect(covariate_keys, "age")]))
numeric_covar_cols <- numeric_covar_cols[!is.na(numeric_covar_cols) & numeric_covar_cols != ""]

categorical_keys <- setdiff(covariate_keys, "age")
categorical_covar_cols <- unname(unlist(covariate_lookup[categorical_keys]))
categorical_covar_cols <- categorical_covar_cols[!is.na(categorical_covar_cols) & categorical_covar_cols != ""]

missing_covars <- setdiff(covar_cols, colnames(meta))

if (length(missing_covars) > 0) {
  stop(glue::glue(
    "Missing covariate columns in metadata: {paste(missing_covars, collapse = ', ')}"
  ))
}

if (length(categorical_covar_cols) > 0) {
  meta <- meta %>%
    mutate(across(all_of(categorical_covar_cols), ~ as.factor(.)))
}

if (length(numeric_covar_cols) > 0) {
  meta <- meta %>%
    mutate(across(all_of(numeric_covar_cols), ~ suppressWarnings(as.numeric(as.character(.)))))
}

# ---- Create modeling frame ----
# Transpose to Samples x CpGs
beta_t <- t(beta)

# Use sparse matrix for high-dimensional modeling
x_cpg <- Matrix::Matrix(beta_t, sparse = TRUE)

# Make sure dimnames are preserved
rownames(x_cpg) <- rownames(beta_t)
colnames(x_cpg) <- colnames(beta_t)

# ---- Final sample synchronization BEFORE saving ----
meta <- meta %>%
  filter(sample_id %in% rownames(x_cpg)) %>%
  arrange(match(sample_id, rownames(x_cpg)))

x_cpg <- x_cpg[meta$sample_id, , drop = FALSE]

stopifnot(
  nrow(x_cpg) == nrow(meta),
  identical(rownames(x_cpg), meta$sample_id)
)

# ---- Save/reuse shared CpG matrix ----
x_cpg_shared_path <- file.path(out_paths$data_processed, "classifier_x_cpg_shared.rds")

source_input_times <- fs::file_info(c(beta_path, meta_path))$modification_time
newest_source_input <- max(source_input_times, na.rm = TRUE)
shared_cache_is_current <- fs::file_exists(x_cpg_shared_path) &&
  fs::file_info(x_cpg_shared_path)$modification_time >= newest_source_input

if (shared_cache_is_current) {
  message("Reusing current shared CpG matrix cache: ", x_cpg_shared_path)
} else {
  message("Building shared CpG matrix cache: ", x_cpg_shared_path)
  beta_shared_t <- t(beta_shared)
  x_cpg_shared <- Matrix::Matrix(beta_shared_t, sparse = TRUE)
  rownames(x_cpg_shared) <- rownames(beta_shared_t)
  colnames(x_cpg_shared) <- colnames(beta_shared_t)
  x_cpg_shared <- x_cpg_shared[meta_shared$sample_id, , drop = FALSE]
  saveRDS(x_cpg_shared, file = x_cpg_shared_path)
}

# ---- Save modeling input ----
saveRDS(
  list(
    x_cpg_path = x_cpg_shared_path,
    x_cpg_dim = dim(x_cpg),
    meta = meta,
    sample_id_col = sample_id_col,
    covar_cols = covar_cols,
    categorical_covar_cols = categorical_covar_cols,
    numeric_covar_cols = numeric_covar_cols,
    covariate_set = covariate_set,
    outcome_info = outcome_info,
    tune_metric = config$modeling$tune_metric %||% "Accuracy",
    holdout_info = holdout_info,
    holdout_meta = holdout_meta,
    run_id = config$run_id
  ),
  file = run_file(out_paths$data_processed, "classifier_inputs")
)

yaml::write_yaml(
  outcome_info,
  file = run_file(out_paths$data_processed, "outcome_definition", "yaml")
)

# ---- Split / CV scaffold ----
set.seed(config$modeling$seed)

cv_strategy <- config$modeling$cv_strategy %||% "kfold"

split <- list()

if (cv_strategy == "train_test") {
  train_idx <- caret::createDataPartition(
    meta$stage,
    p = 1 - config$modeling$test_split,
    list = FALSE
  )

  train_idx <- as.integer(train_idx)

  split$train_idx <- train_idx
  split$test_idx <- setdiff(seq_len(nrow(meta)), train_idx)

} else if (cv_strategy == "kfold") {
  folds_out <- caret::createFolds(
    meta$stage,
    k = config$modeling$cv_folds,
    list = TRUE
  )

  # caret::trainControl(index = ...) expects TRAINING indices.
  # caret::createFolds() returns held-out/test fold indices.
  folds_train <- lapply(folds_out, function(test_idx) {
    setdiff(seq_len(nrow(meta)), test_idx)
  })

  split$folds_out <- folds_out
  split$folds_train <- folds_train

} else if (cv_strategy == "loso") {
  loso_col <- config$modeling$loso_group_col %||% config$covariates$study_id_col

  if (!loso_col %in% colnames(meta)) {
    stop(glue::glue("LOSO group column '{loso_col}' not found in metadata."))
  }

  groups <- unique(meta[[loso_col]])

  loso_out <- lapply(groups, function(g) {
    which(meta[[loso_col]] == g)
  })

  names(loso_out) <- as.character(groups)

  loso_train <- lapply(loso_out, function(test_idx) {
    setdiff(seq_len(nrow(meta)), test_idx)
  })

  split$loso_out <- loso_out
  split$loso_train <- loso_train

} else {
  stop(glue::glue("Unknown cv_strategy: {cv_strategy}"))
}

saveRDS(
  split,
  file = run_file(out_paths$data_processed, "cv_split")
)

message("Prepared classifier inputs and CV scaffold.")
message("Run ID: ", config$run_id)
message("Outcome: ", outcome_info$name, " - ", outcome_info$description)
message("Tune metric: ", config$modeling$tune_metric %||% "Accuracy")
if (!is.null(holdout_info)) {
  message("Validation holdout: ", holdout_info$name, " (", holdout_info$n, " samples removed before training)")
}
message("Saved classifier inputs to: ", run_file(out_paths$data_processed, "classifier_inputs"))
message("Shared CpG matrix: ", x_cpg_shared_path)
message("Saved CV split to: ", run_file(out_paths$data_processed, "cv_split"))
message("Samples: ", nrow(meta))
message("CpGs: ", ncol(x_cpg))
message("Classes: ", paste(levels(meta$stage), collapse = ", "))
