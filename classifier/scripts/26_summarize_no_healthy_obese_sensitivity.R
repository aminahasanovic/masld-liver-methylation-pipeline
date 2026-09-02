# 26_summarize_no_healthy_obese_sensitivity.R
# Summarize Elastic Net performance with Healthy_Obese source samples excluded.

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(glue)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

message("Summarizing Healthy_Obese exclusion sensitivity.")

outcome_pairs <- tibble::tribble(
  ~analysis_label, ~baseline_outcome, ~sensitivity_outcome,
  "4-class", "disease_group_4", "disease_group_4_no_healthy_obese",
  "3-class", "disease_group_3", "disease_group_3_no_healthy_obese",
  "Binary", "binary_healthy_vs_disease", "binary_healthy_vs_disease_no_healthy_obese"
)

run_specs <- tibble::tibble(strategy = c("train_test", "loso")) |>
  tidyr::crossing(outcome_pairs)

make_run_id <- function(outcome, strategy) {
  paste(c(outcome, "cpg_only", strategy), collapse = "_")
}

read_auc_value <- function(run_id, strategy) {
  auc_stem <- if (identical(strategy, "loso")) {
    "elastic_net_loso_auc"
  } else {
    "elastic_net_auc"
  }

  auc_path <- file.path(out_paths$results_metrics, paste0(auc_stem, "_", run_id, ".csv"))
  if (!file.exists(auc_path)) {
    return(NA_real_)
  }

  auc_tbl <- readr::read_csv(auc_path, show_col_types = FALSE)

  if ("class" %in% names(auc_tbl) && "macro_mean" %in% auc_tbl$class) {
    return(auc_tbl$auc[auc_tbl$class == "macro_mean"][1])
  }

  if ("type" %in% names(auc_tbl) && "ovr_macro" %in% auc_tbl$type) {
    return(auc_tbl$auc[auc_tbl$type == "ovr_macro"][1])
  }

  if ("auc" %in% names(auc_tbl)) {
    return(auc_tbl$auc[1])
  }

  NA_real_
}

read_run_metrics <- function(outcome, analysis_label, strategy, analysis_set) {
  run_id <- make_run_id(outcome, strategy)
  metric_stem <- if (identical(strategy, "loso")) {
    "elastic_net_loso_metrics"
  } else {
    "elastic_net_metrics"
  }

  metric_path <- file.path(out_paths$results_metrics, paste0(metric_stem, "_", run_id, ".csv"))
  input_path <- file.path(out_paths$data_processed, paste0("classifier_inputs_", run_id, ".rds"))

  if (!file.exists(metric_path)) {
    return(tibble::tibble(
      analysis_label = analysis_label,
      analysis_set = analysis_set,
      outcome = outcome,
      strategy = strategy,
      run_id = run_id,
      available = FALSE,
      missing_file = fs::path_rel(metric_path, start = project_root)
    ))
  }

  metrics <- readr::read_csv(metric_path, show_col_types = FALSE)

  sample_n <- NA_integer_
  class_counts <- NA_character_
  if (file.exists(input_path)) {
    inp <- readRDS(input_path)
    sample_n <- nrow(inp$meta)
    class_counts <- paste(
      names(table(inp$meta$stage)),
      as.integer(table(inp$meta$stage)),
      sep = "=",
      collapse = "; "
    )
  }

  metrics |>
    dplyr::mutate(
      auc = read_auc_value(run_id, strategy),
      analysis_label = analysis_label,
      analysis_set = analysis_set,
      outcome = outcome,
      strategy = strategy,
      run_id = run_id,
      available = TRUE,
      missing_file = NA_character_,
      sample_n = sample_n,
      class_counts = class_counts,
      .before = 1
    )
}

all_runs <- purrr::pmap_dfr(run_specs, function(strategy, analysis_label,
                                                baseline_outcome, sensitivity_outcome) {
  dplyr::bind_rows(
    read_run_metrics(
      outcome = baseline_outcome,
      analysis_label = analysis_label,
      strategy = strategy,
      analysis_set = "baseline"
    ),
    read_run_metrics(
      outcome = sensitivity_outcome,
      analysis_label = analysis_label,
      strategy = strategy,
      analysis_set = "no_healthy_obese"
    )
  )
})

metric_cols <- intersect(
  c("accuracy", "balanced_accuracy", "mcc", "f1_macro", "auc"),
  colnames(all_runs)
)

comparison <- all_runs |>
  dplyr::filter(available) |>
  dplyr::select(
    analysis_label,
    strategy,
    analysis_set,
    outcome,
    run_id,
    sample_n,
    class_counts,
    dplyr::all_of(metric_cols)
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(metric_cols),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::select(
    analysis_label,
    strategy,
    metric,
    analysis_set,
    value,
    outcome,
    run_id,
    sample_n,
    class_counts
  ) |>
  tidyr::pivot_wider(
    names_from = analysis_set,
    values_from = c(value, outcome, run_id, sample_n, class_counts),
    names_sep = "__"
  )

required_numeric_cols <- c(
  "value__baseline",
  "value__no_healthy_obese",
  "sample_n__baseline",
  "sample_n__no_healthy_obese"
)
required_character_cols <- c(
  "outcome__baseline",
  "outcome__no_healthy_obese",
  "run_id__baseline",
  "run_id__no_healthy_obese",
  "class_counts__baseline",
  "class_counts__no_healthy_obese"
)

for (nm in required_numeric_cols) {
  if (!nm %in% names(comparison)) {
    comparison[[nm]] <- NA_real_
  }
}

for (nm in required_character_cols) {
  if (!nm %in% names(comparison)) {
    comparison[[nm]] <- NA_character_
  }
}

comparison <- comparison |>
  dplyr::mutate(
    delta_no_healthy_obese_minus_baseline =
      .data$value__no_healthy_obese - .data$value__baseline
  ) |>
  dplyr::arrange(analysis_label, strategy, metric)

source_meta <- readRDS(check_file(config$inputs$metadata_path))
source_col <- config$outcome$source_column %||% "DiseaseGroup"
dataset_col <- config$covariates$study_id_col %||% "Dataset"

removed_summary <- source_meta |>
  dplyr::filter(.data[[source_col]] == "Healthy_Obese") |>
  dplyr::count(.data[[dataset_col]], name = "n_removed") |>
  dplyr::arrange(dplyr::desc(n_removed))

remaining_summary <- source_meta |>
  dplyr::filter(.data[[source_col]] != "Healthy_Obese") |>
  dplyr::count(.data[[source_col]], name = "n_remaining") |>
  dplyr::arrange(.data[[source_col]])

run_metrics_path <- file.path(
  out_paths$results_metrics,
  "no_healthy_obese_elastic_net_run_metrics.csv"
)
comparison_path <- file.path(
  out_paths$results_metrics,
  "no_healthy_obese_elastic_net_performance_comparison.csv"
)
removed_path <- file.path(
  out_paths$results_metrics,
  "no_healthy_obese_removed_sample_summary.csv"
)
remaining_path <- file.path(
  out_paths$results_metrics,
  "no_healthy_obese_remaining_source_group_summary.csv"
)
summary_path <- file.path(
  out_paths$results_metrics,
  "no_healthy_obese_elastic_net_summary.md"
)

readr::write_csv(all_runs, run_metrics_path)
readr::write_csv(comparison, comparison_path)
readr::write_csv(removed_summary, removed_path)
readr::write_csv(remaining_summary, remaining_path)

format_num <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.3f", x))
}

capture_tbl <- function(x) {
  paste(capture.output(print(tibble::as_tibble(x), n = Inf)), collapse = "\n")
}

comparison_md <- comparison |>
  dplyr::filter(metric %in% c("balanced_accuracy", "f1_macro", "mcc", "auc")) |>
  dplyr::mutate(
    baseline = format_num(.data$value__baseline),
    no_healthy_obese = format_num(.data$value__no_healthy_obese),
    delta = format_num(.data$delta_no_healthy_obese_minus_baseline)
  ) |>
  dplyr::select(
    analysis_label,
    strategy,
    metric,
    baseline,
    no_healthy_obese,
    delta
  )

missing_runs <- all_runs |>
  dplyr::filter(!available) |>
  dplyr::select(analysis_label, analysis_set, strategy, run_id, missing_file)

md <- c(
  "# Healthy_Obese Exclusion Sensitivity",
  "",
  glue::glue(
    "Removed source group: `Healthy_Obese` ({sum(removed_summary$n_removed)} samples)."
  ),
  "",
  "## Removed Samples By Dataset",
  "",
  capture_tbl(removed_summary),
  "",
  "## Remaining Source Groups",
  "",
  capture_tbl(remaining_summary),
  "",
  "## Performance Deltas",
  "",
  "Delta = no Healthy_Obese minus baseline.",
  "",
  capture_tbl(comparison_md),
  "",
  "## Output Files",
  "",
  glue::glue("- `{fs::path_rel(run_metrics_path, start = project_root)}`"),
  glue::glue("- `{fs::path_rel(comparison_path, start = project_root)}`"),
  glue::glue("- `{fs::path_rel(removed_path, start = project_root)}`"),
  glue::glue("- `{fs::path_rel(remaining_path, start = project_root)}`")
)

if (nrow(missing_runs) > 0) {
  md <- c(
    md,
    "",
    "## Missing Runs",
    "",
    capture_tbl(missing_runs)
  )
}

writeLines(md, summary_path)

message("Saved run metrics to: ", run_metrics_path)
message("Saved comparison to: ", comparison_path)
message("Saved summary to: ", summary_path)
