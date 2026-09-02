# 17_compare_classifier_models.R
# Compact Elastic Net, Random Forest, and Boosted Trees comparison for the main CpG-only runs.

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

outcomes <- c(
  "disease_group_5_obesity_split",
  "disease_group_4",
  "disease_group_3",
  "binary_healthy_vs_disease"
)

outcome_labels <- c(
  disease_group_5_obesity_split = "5-class obesity split",
  disease_group_4 = "4-class",
  disease_group_3 = "3-class",
  binary_healthy_vs_disease = "Binary"
)

model_specs <- tibble::tribble(
  ~model, ~model_label, ~strategy, ~metric_stem, ~confusion_stem, ~auc_stem, ~fold_stem,
  "elastic_net", "Elastic Net", "train_test", "elastic_net_metrics", "elastic_net_confusion_matrix", "elastic_net_auc", NA_character_,
  "elastic_net", "Elastic Net", "loso", "elastic_net_loso_metrics", "elastic_net_loso_confusion_matrix", "elastic_net_loso_auc", "elastic_net_loso_fold_metrics",
  "random_forest", "Random Forest", "train_test", "random_forest_metrics", "random_forest_confusion_matrix", "random_forest_auc", NA_character_,
  "random_forest", "Random Forest", "loso", "random_forest_loso_metrics", "random_forest_loso_confusion_matrix", "random_forest_loso_auc", "random_forest_loso_fold_metrics",
  "boosted_trees", "Boosted Trees", "train_test", "boosted_trees_metrics", "boosted_trees_confusion_matrix", "boosted_trees_auc", NA_character_,
  "boosted_trees", "Boosted Trees", "loso", "boosted_trees_loso_metrics", "boosted_trees_loso_confusion_matrix", "boosted_trees_loso_auc", "boosted_trees_loso_fold_metrics"
)

model_label_levels <- c("Elastic Net", "Random Forest", "Boosted Trees")

make_run_id <- function(outcome, strategy) {
  paste(outcome, "cpg_only", strategy, sep = "_")
}

read_metric <- function(outcome, spec) {
  outcome_value <- outcome
  outcome_label_value <- unname(outcome_labels[[outcome_value]])
  run_id <- make_run_id(outcome_value, spec$strategy)
  path_in <- file.path(out_paths$results_metrics, paste0(spec$metric_stem, "_", run_id, ".csv"))

  if (!file.exists(path_in)) {
    warning("Missing metric file: ", path_in)
    return(tibble::tibble())
  }

  metrics <- readr::read_csv(path_in, show_col_types = FALSE)

  auc_path <- file.path(out_paths$results_metrics, paste0(spec$auc_stem, "_", run_id, ".csv"))
  auc_value <- NA_real_

  if (file.exists(auc_path)) {
    auc_tbl <- readr::read_csv(auc_path, show_col_types = FALSE)

    if ("class" %in% names(auc_tbl) && "macro_mean" %in% auc_tbl$class) {
      auc_value <- auc_tbl$auc[auc_tbl$class == "macro_mean"][1]
    } else if ("auc" %in% names(auc_tbl)) {
      auc_value <- auc_tbl$auc[1]
    }
  }

  metrics |>
    dplyr::mutate(
      model = spec$model,
      model_label = spec$model_label,
      outcome = outcome_value,
      outcome_label = outcome_label_value,
      strategy = spec$strategy,
      run_id = run_id,
      auc = auc_value,
      .before = 1
    )
}

read_confusion <- function(outcome, spec) {
  outcome_value <- outcome
  outcome_label_value <- unname(outcome_labels[[outcome_value]])
  run_id <- make_run_id(outcome_value, spec$strategy)
  path_in <- file.path(out_paths$results_metrics, paste0(spec$confusion_stem, "_", run_id, ".csv"))

  if (!file.exists(path_in)) {
    warning("Missing confusion matrix file: ", path_in)
    return(tibble::tibble())
  }

  readr::read_csv(path_in, show_col_types = FALSE) |>
    dplyr::mutate(
      model = spec$model,
      model_label = spec$model_label,
      outcome = outcome_value,
      outcome_label = outcome_label_value,
      strategy = spec$strategy,
      run_id = run_id,
      .before = 1
    )
}

read_fold_metrics <- function(outcome, spec) {
  if (is.na(spec$fold_stem)) {
    return(tibble::tibble())
  }

  outcome_value <- outcome
  outcome_label_value <- unname(outcome_labels[[outcome_value]])
  run_id <- make_run_id(outcome_value, spec$strategy)
  path_in <- file.path(out_paths$results_metrics, paste0(spec$fold_stem, "_", run_id, ".csv"))

  if (!file.exists(path_in)) {
    warning("Missing fold metric file: ", path_in)
    return(tibble::tibble())
  }

  readr::read_csv(path_in, show_col_types = FALSE) |>
    dplyr::mutate(
      model = spec$model,
      model_label = spec$model_label,
      outcome = outcome_value,
      outcome_label = outcome_label_value,
      strategy = spec$strategy,
      run_id = run_id,
      .before = 1
    )
}

metrics <- purrr::map_dfr(outcomes, function(outcome) {
  purrr::pmap_dfr(model_specs, function(...) {
    read_metric(outcome, tibble::tibble(...))
  })
}) |>
  dplyr::mutate(
    model_label = factor(model_label, levels = model_label_levels),
    outcome = factor(outcome, levels = outcomes),
    outcome_label = factor(outcome_label, levels = unname(outcome_labels[outcomes])),
    strategy = factor(strategy, levels = c("train_test", "loso"), labels = c("Train/test", "LOSO"))
  ) |>
  dplyr::arrange(outcome, strategy, model_label)

readr::write_csv(
  metrics,
  file.path(out_paths$results_metrics, "model_comparison_elastic_net_random_forest.csv")
)

readr::write_csv(
  metrics,
  file.path(out_paths$results_metrics, "model_comparison_three_classifiers.csv")
)

metric_cols <- intersect(
  c("accuracy", "balanced_accuracy", "mcc", "f1_macro", "auc"),
  colnames(metrics)
)

metrics_long <- metrics |>
  dplyr::select(model_label, outcome_label, strategy, all_of(metric_cols)) |>
  tidyr::pivot_longer(
    cols = all_of(metric_cols),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("accuracy", "balanced_accuracy", "f1_macro", "mcc", "auc"),
      labels = c("Accuracy", "Balanced accuracy", "Macro F1", "MCC", "AUC")
    )
  )

p_metrics <- ggplot(metrics_long, aes(x = model_label, y = value, fill = model_label)) +
  geom_col(width = 0.72, na.rm = TRUE) +
  facet_grid(metric ~ outcome_label + strategy, scales = "free_y") +
  labs(
    x = NULL,
    y = NULL,
    title = "Classifier Model Comparison",
    subtitle = "CpG-only models; feature selection performed inside the training data"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(out_paths$results_plots, "model_comparison_elastic_net_vs_random_forest_metrics.png"),
  plot = p_metrics,
  width = 13,
  height = 9,
  dpi = 300
)

ggsave(
  filename = file.path(out_paths$results_plots, "model_comparison_three_classifiers_metrics.png"),
  plot = p_metrics,
  width = 13,
  height = 9,
  dpi = 300
)

fold_metrics <- purrr::map_dfr(outcomes, function(outcome) {
  purrr::pmap_dfr(model_specs |> dplyr::filter(strategy == "loso"), function(...) {
    read_fold_metrics(outcome, tibble::tibble(...))
  })
}) |>
  dplyr::mutate(
    model_label = factor(model_label, levels = model_label_levels),
    outcome = factor(outcome, levels = outcomes),
    outcome_label = factor(outcome_label, levels = unname(outcome_labels[outcomes]))
  )

if (nrow(fold_metrics) > 0) {
  readr::write_csv(
    fold_metrics,
    file.path(out_paths$results_metrics, "model_comparison_loso_fold_metrics_elastic_net_random_forest.csv")
  )

  readr::write_csv(
    fold_metrics,
    file.path(out_paths$results_metrics, "model_comparison_loso_fold_metrics_three_classifiers.csv")
  )

  fold_metric_cols <- intersect(
    c("accuracy", "balanced_accuracy", "macro_recall_present", "mcc", "f1_macro", "f1_macro_present"),
    colnames(fold_metrics)
  )

  fold_long <- fold_metrics |>
    dplyr::select(model_label, outcome_label, outer_fold, all_of(fold_metric_cols)) |>
    tidyr::pivot_longer(
      cols = all_of(fold_metric_cols),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::mutate(
      metric = factor(
        metric,
        levels = c("accuracy", "balanced_accuracy", "macro_recall_present", "f1_macro", "f1_macro_present", "mcc"),
        labels = c("Accuracy", "Balanced accuracy", "Macro recall present", "Macro F1", "Macro F1 present", "MCC")
      )
    )

  fold_line_data <- fold_long |>
    dplyr::group_by(model_label, outcome_label, metric) |>
    dplyr::filter(sum(!is.na(value)) >= 2) |>
    dplyr::ungroup()

  p_folds <- ggplot(fold_long, aes(x = outer_fold, y = value, color = model_label, group = model_label)) +
    geom_point(size = 2, na.rm = TRUE) +
    geom_line(data = fold_line_data, na.rm = TRUE) +
    facet_grid(metric ~ outcome_label, scales = "free_y") +
    labs(
      x = "Held-out study",
      y = NULL,
      color = NULL,
      title = "LOSO Fold Performance",
      subtitle = "Fold-level values are sensitive to which classes are present in the held-out study"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(out_paths$results_plots, "model_comparison_loso_fold_metrics_elastic_net_vs_random_forest.png"),
    plot = p_folds,
    width = 12,
    height = 10,
    dpi = 300
  )

  ggsave(
    filename = file.path(out_paths$results_plots, "model_comparison_loso_fold_metrics_three_classifiers.png"),
    plot = p_folds,
    width = 12,
    height = 10,
    dpi = 300
  )
}

confusion <- purrr::map_dfr(outcomes, function(outcome) {
  purrr::pmap_dfr(model_specs, function(...) {
    read_confusion(outcome, tibble::tibble(...))
  })
}) |>
  dplyr::mutate(
    model_label = factor(model_label, levels = model_label_levels),
    outcome = factor(outcome, levels = outcomes),
    outcome_label = factor(outcome_label, levels = unname(outcome_labels[outcomes])),
    strategy = factor(strategy, levels = c("train_test", "loso"), labels = c("Train/test", "LOSO"))
  ) |>
  dplyr::group_by(run_id, model_label, outcome_label, strategy, Reference) |>
  dplyr::mutate(reference_fraction = Freq / sum(Freq)) |>
  dplyr::ungroup()

readr::write_csv(
  confusion,
  file.path(out_paths$results_metrics, "model_comparison_confusion_matrices_elastic_net_random_forest.csv")
)

readr::write_csv(
  confusion,
  file.path(out_paths$results_metrics, "model_comparison_confusion_matrices_three_classifiers.csv")
)

plot_confusion <- function(strategy_label, filename) {
  plot_data <- confusion |>
    dplyr::filter(strategy == strategy_label)

  if (nrow(plot_data) == 0) {
    return(invisible(NULL))
  }

  p <- ggplot(plot_data, aes(x = Reference, y = Prediction, fill = reference_fraction)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = Freq), size = 2.4) +
    facet_grid(outcome_label ~ model_label, scales = "free", space = "free") +
    scale_fill_viridis_c(option = "C", limits = c(0, 1)) +
    labs(
      x = "True class",
      y = "Predicted class",
      fill = "Fraction\nwithin true class",
      title = paste(strategy_label, "Confusion Matrices")
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(
    filename = file.path(out_paths$results_plots, filename),
    plot = p,
    width = 11,
    height = 9,
    dpi = 300
  )
}

plot_confusion("Train/test", "model_comparison_confusion_heatmaps_train_test.png")
plot_confusion("LOSO", "model_comparison_confusion_heatmaps_loso.png")

best_rows <- metrics |>
  dplyr::group_by(outcome_label, strategy) |>
  dplyr::slice_max(order_by = balanced_accuracy, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(outcome_label, strategy, model_label, balanced_accuracy, accuracy, mcc, f1_macro, auc)

report_lines <- c(
"# Classifier model comparison",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "All compared runs use `cpg_only` predictors and run-specific inputs.",
  "For RF and Boosted Trees, CpG variance filtering is performed on the training data only; for LOSO this is repeated within each held-out-study fold.",
  "",
  "## Best model per outcome/CV by balanced accuracy",
  "",
  paste(capture.output(print(best_rows, n = Inf, width = Inf)), collapse = "\n"),
  "",
  "## Output files",
  "",
  "- `results/metrics/model_comparison_elastic_net_random_forest.csv`",
  "- `results/metrics/model_comparison_three_classifiers.csv`",
  "- `results/metrics/model_comparison_loso_fold_metrics_elastic_net_random_forest.csv`",
  "- `results/metrics/model_comparison_loso_fold_metrics_three_classifiers.csv`",
  "- `results/metrics/model_comparison_confusion_matrices_elastic_net_random_forest.csv`",
  "- `results/metrics/model_comparison_confusion_matrices_three_classifiers.csv`",
  "- `results/plots/model_comparison_elastic_net_vs_random_forest_metrics.png`",
  "- `results/plots/model_comparison_three_classifiers_metrics.png`",
  "- `results/plots/model_comparison_loso_fold_metrics_elastic_net_vs_random_forest.png`",
  "- `results/plots/model_comparison_loso_fold_metrics_three_classifiers.png`",
  "- `results/plots/model_comparison_confusion_heatmaps_train_test.png`",
  "- `results/plots/model_comparison_confusion_heatmaps_loso.png`"
)

writeLines(
  report_lines,
  con = file.path(out_paths$results_metrics, "model_comparison_elastic_net_random_forest.md")
)

writeLines(
  report_lines,
  con = file.path(out_paths$results_metrics, "model_comparison_three_classifiers.md")
)

message("Model comparison complete.")
message("Saved comparison table to: ", file.path(out_paths$results_metrics, "model_comparison_three_classifiers.csv"))
message("Saved comparison plots to: ", out_paths$results_plots)
print(metrics |> dplyr::select(outcome_label, strategy, model_label, accuracy, balanced_accuracy, mcc, f1_macro, auc), n = Inf)
