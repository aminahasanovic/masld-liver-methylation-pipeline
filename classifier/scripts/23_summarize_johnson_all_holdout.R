# 23_summarize_johnson_all_holdout.R
# Compact summary for the strict all-Johnson holdout.

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

holdout_id <- "johnson_allholdout"
outcomes <- c("disease_group_4", "disease_group_5_obesity_split")
models <- c("elastic_net", "random_forest", "boosted_trees")
model_labels <- c(
  elastic_net = "Elastic Net",
  random_forest = "Random Forest",
  boosted_trees = "Boosted Trees"
)
outcome_labels <- c(
  disease_group_4 = "4-class",
  disease_group_5_obesity_split = "5-class"
)

run_id_for <- function(outcome) {
  paste(outcome, "cpg_only", "train_test", holdout_id, sep = "_")
}

read_metric <- function(outcome, model) {
  outcome_value <- outcome
  run_id <- run_id_for(outcome)
  path <- file.path(out_paths$results_metrics, paste0(model, "_metrics_", run_id, ".csv"))
  if (!file.exists(path)) {
    warning("Missing metric file: ", path)
    return(tibble::tibble())
  }

  readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::mutate(
      outcome = outcome_value,
      outcome_label = unname(outcome_labels[[outcome_value]]),
      model = model,
      model_label = unname(model_labels[[model]]),
      run_id = run_id,
      .before = 1
    )
}

read_by_source <- function(outcome) {
  outcome_value <- outcome
  run_id <- run_id_for(outcome)
  path <- file.path(out_paths$results_metrics, paste0("holdout_prediction_by_source_summary_", run_id, ".csv"))
  if (!file.exists(path)) {
    warning("Missing holdout by-source file: ", path)
    return(tibble::tibble())
  }

  readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::mutate(
      outcome = outcome_value,
      outcome_label = unname(outcome_labels[[outcome_value]]),
      run_id = run_id,
      .before = 1
    )
}

read_prob_by_source <- function(outcome) {
  outcome_value <- outcome
  run_id <- run_id_for(outcome)
  path <- file.path(out_paths$results_metrics, paste0("holdout_probability_by_source_summary_", run_id, ".csv"))
  if (!file.exists(path)) {
    warning("Missing holdout probability by-source file: ", path)
    return(tibble::tibble())
  }

  readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::mutate(
      outcome = outcome_value,
      outcome_label = unname(outcome_labels[[outcome_value]]),
      run_id = run_id,
      .before = 1
    )
}

metrics <- dplyr::bind_rows(lapply(outcomes, function(outcome) {
  dplyr::bind_rows(lapply(models, function(model) read_metric(outcome, model)))
})) |>
  dplyr::select(
    outcome,
    outcome_label,
    model,
    model_label,
    run_id,
    accuracy,
    balanced_accuracy,
    dplyr::any_of(c("macro_recall_present")),
    mcc,
    f1_macro,
    dplyr::any_of(c("f1_macro_present"))
  ) |>
  dplyr::arrange(outcome, model)

by_source <- dplyr::bind_rows(lapply(outcomes, read_by_source)) |>
  dplyr::arrange(outcome, model, source_disease_group, prediction)

prob_by_source <- dplyr::bind_rows(lapply(outcomes, read_prob_by_source)) |>
  dplyr::arrange(outcome, model, source_disease_group, class)

healthy_obese_predictions <- by_source |>
  dplyr::filter(source_disease_group == "Healthy_Obese") |>
  dplyr::select(outcome_label, model, prediction, n, total, fraction)

advanced_predictions <- by_source |>
  dplyr::filter(source_disease_group == "Advanced_Fibrosis") |>
  dplyr::select(outcome_label, model, prediction, n, total, fraction)

readr::write_csv(
  metrics,
  file.path(out_paths$results_metrics, "johnson_all_holdout_internal_metrics.csv")
)
readr::write_csv(
  by_source,
  file.path(out_paths$results_metrics, "johnson_all_holdout_prediction_by_source_summary.csv")
)
readr::write_csv(
  prob_by_source,
  file.path(out_paths$results_metrics, "johnson_all_holdout_probability_by_source_summary.csv")
)

report_lines <- c(
  "# Johnson-all holdout summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Holdout definition: all `Dataset == Johnson` samples removed before model training and predicted afterward.",
  "",
  "## Internal non-Johnson train/test metrics",
  "",
  paste(capture.output(print(metrics, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "## Johnson Healthy_Obese predictions",
  "",
  paste(capture.output(print(healthy_obese_predictions, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "## Johnson Advanced_Fibrosis predictions",
  "",
  paste(capture.output(print(advanced_predictions, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "## Output files",
  "",
  "- `results/metrics/johnson_all_holdout_internal_metrics.csv`",
  "- `results/metrics/johnson_all_holdout_prediction_by_source_summary.csv`",
  "- `results/metrics/johnson_all_holdout_probability_by_source_summary.csv`"
)

writeLines(
  report_lines,
  file.path(out_paths$results_metrics, "johnson_all_holdout_summary.md")
)

message("Johnson-all holdout summary complete.")
message("Saved summary to: ", file.path(out_paths$results_metrics, "johnson_all_holdout_summary.md"))
print(metrics, n = Inf)
print(healthy_obese_predictions, n = Inf)
print(advanced_predictions, n = Inf)
