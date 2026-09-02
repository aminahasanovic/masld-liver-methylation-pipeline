# 22_summarize_obesity_split_sensitivity.R
# Compact summary for the Healthy versus Healthy_Obese sensitivity analyses.

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

metric_path <- file.path(out_paths$results_metrics, "model_comparison_three_classifiers.csv")
if (!file.exists(metric_path)) {
  stop("Missing model comparison table. Run scripts/17_compare_classifier_models.R first.", call. = FALSE)
}

metrics <- readr::read_csv(metric_path, show_col_types = FALSE)

focused_metrics <- metrics |>
  dplyr::filter(outcome %in% c("disease_group_5_obesity_split", "disease_group_4")) |>
  dplyr::select(
    outcome,
    outcome_label,
    strategy,
    model_label,
    accuracy,
    balanced_accuracy,
    f1_macro,
    mcc,
    auc
  ) |>
  dplyr::arrange(outcome, strategy, model_label)

read_optional_csv <- function(filename) {
  path <- file.path(out_paths$results_metrics, filename)
  if (!file.exists(path)) {
    warning("Missing optional file: ", path)
    return(tibble::tibble())
  }
  readr::read_csv(path, show_col_types = FALSE)
}

holdout_4 <- read_optional_csv(
  "holdout_prediction_summary_disease_group_4_cpg_only_train_test_johnson_nonfibroticholdout.csv"
) |>
  dplyr::mutate(outcome = "disease_group_4", .before = 1)

holdout_5 <- read_optional_csv(
  "holdout_prediction_summary_disease_group_5_obesity_split_cpg_only_train_test_johnson_nonfibroticholdout.csv"
) |>
  dplyr::mutate(outcome = "disease_group_5_obesity_split", .before = 1)

holdout_predictions <- dplyr::bind_rows(holdout_4, holdout_5) |>
  dplyr::arrange(outcome, model, prediction)

holdout_prob_5 <- read_optional_csv(
  "holdout_probability_summary_disease_group_5_obesity_split_cpg_only_train_test_johnson_nonfibroticholdout.csv"
) |>
  dplyr::mutate(outcome = "disease_group_5_obesity_split", .before = 1)

readr::write_csv(
  focused_metrics,
  file.path(out_paths$results_metrics, "obesity_split_sensitivity_model_metrics.csv")
)
readr::write_csv(
  holdout_predictions,
  file.path(out_paths$results_metrics, "obesity_split_sensitivity_johnson_holdout_predictions.csv")
)
readr::write_csv(
  holdout_prob_5,
  file.path(out_paths$results_metrics, "obesity_split_sensitivity_johnson_holdout_probabilities_5class.csv")
)

fmt <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.3f", x))
}

best_train_test <- focused_metrics |>
  dplyr::filter(strategy == "Train/test") |>
  dplyr::group_by(outcome_label) |>
  dplyr::slice_max(order_by = balanced_accuracy, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

best_loso <- focused_metrics |>
  dplyr::filter(strategy == "LOSO") |>
  dplyr::group_by(outcome_label) |>
  dplyr::slice_max(order_by = balanced_accuracy, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

holdout_advanced <- holdout_predictions |>
  dplyr::filter(prediction == "Advanced_Fibrosis") |>
  dplyr::select(outcome, model, n, total, fraction)

report_lines <- c(
  "# Obesity-split classifier sensitivity",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Best train/test balanced accuracy",
  "",
  paste(capture.output(print(best_train_test, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "## Best LOSO balanced accuracy",
  "",
  paste(capture.output(print(best_loso, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "## Johnson non-fibrotic holdout",
  "",
  "Holdout definition: `Dataset == Johnson` and source `DiseaseGroup == Healthy_Obese`; samples were removed before training and predicted afterward.",
  "",
  paste(capture.output(print(holdout_advanced, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "Interpretation note: if only Johnson non-fibrotic samples are held out, Johnson advanced-fibrosis samples remain in training. This is useful for testing the disputed label, but it can also expose Johnson-specific dataset signal. A stricter follow-up is to hold out all Johnson samples and inspect whether non-fibrotic and advanced Johnson samples separate without any Johnson data in training.",
  "",
  "## Output files",
  "",
  "- `results/metrics/obesity_split_sensitivity_model_metrics.csv`",
  "- `results/metrics/obesity_split_sensitivity_johnson_holdout_predictions.csv`",
  "- `results/metrics/obesity_split_sensitivity_johnson_holdout_probabilities_5class.csv`"
)

writeLines(
  report_lines,
  file.path(out_paths$results_metrics, "obesity_split_sensitivity_summary.md")
)

message("Obesity-split sensitivity summary complete.")
message("Saved summary to: ", file.path(out_paths$results_metrics, "obesity_split_sensitivity_summary.md"))
print(best_train_test, n = Inf)
print(best_loso, n = Inf)
print(holdout_advanced, n = Inf)
