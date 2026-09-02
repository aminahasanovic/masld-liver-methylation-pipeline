# 08_plot_elastic_net_results.R
# Compact plots for Elastic Net LOSO results.

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

metric_files <- list.files(
  out_paths$results_metrics,
  pattern = "^elastic_net_loso_metrics_.*\\.rds$",
  full.names = TRUE
)

if (length(metric_files) == 0) {
  stop("No LOSO elastic net metric RDS files found.", call. = FALSE)
}

read_metric_file <- function(path_in) {
  obj <- readRDS(path_in)
  run_id <- stringr::str_remove(basename(path_in), "^elastic_net_loso_metrics_")
  run_id <- stringr::str_remove(run_id, "\\.rds$")

  auc_val <- NA_real_
  if (!is.null(obj$auc)) {
    if ("class" %in% names(obj$auc) && "macro_mean" %in% obj$auc$class) {
      auc_val <- obj$auc$auc[obj$auc$class == "macro_mean"][1]
    } else if ("auc" %in% names(obj$auc)) {
      auc_val <- obj$auc$auc[1]
    }
  }

  tibble::tibble(
    run_id = run_id,
    outcome = obj$outcome_info$name %||% NA_character_,
    covariate_set = obj$covariate_set %||% NA_character_,
    tune_metric = obj$tune_metric %||% NA_character_,
    accuracy = obj$metrics$accuracy[1],
    balanced_accuracy = obj$metrics$balanced_accuracy[1],
    mcc = obj$metrics$mcc[1],
    f1_macro = obj$metrics$f1_macro[1],
    auc = auc_val
  )
}

metrics <- purrr::map_dfr(metric_files, read_metric_file) |>
  dplyr::filter(covariate_set == "cpg_only") |>
  dplyr::mutate(
    outcome = factor(
      outcome,
      levels = c(
        "disease_group_5_obesity_split",
        "disease_group_4",
        "disease_group_4_no_healthy_obese",
        "disease_group_3",
        "disease_group_3_no_healthy_obese",
        "binary_healthy_vs_disease",
        "binary_healthy_vs_disease_no_healthy_obese"
      )
    )
  ) |>
  dplyr::arrange(outcome)

readr::write_csv(
  metrics,
  file.path(out_paths$results_metrics, "elastic_net_loso_metric_summary.csv")
)

metric_long <- metrics |>
  dplyr::select(run_id, outcome, balanced_accuracy, f1_macro, mcc, auc) |>
  tidyr::pivot_longer(
    cols = c(balanced_accuracy, f1_macro, mcc, auc),
    names_to = "metric",
    values_to = "value"
  )

p_metrics <- ggplot(metric_long, aes(x = outcome, y = value, fill = outcome)) +
  geom_col(width = 0.72, na.rm = TRUE) +
  facet_wrap(~metric, scales = "free_y") +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = NULL, y = NULL, title = "Elastic Net LOSO Performance by Outcome") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(out_paths$results_plots, "elastic_net_loso_metric_summary.png"),
  plot = p_metrics,
  width = 11,
  height = 8,
  dpi = 300
)

fold_metrics <- purrr::map_dfr(metrics$run_id, function(run_id) {
  path_in <- file.path(out_paths$results_metrics, paste0("elastic_net_loso_fold_metrics_", run_id, ".csv"))
  if (!file.exists(path_in)) return(tibble::tibble())
  readr::read_csv(path_in, show_col_types = FALSE) |>
    dplyr::mutate(run_id = run_id)
}) |>
  dplyr::left_join(metrics |> dplyr::select(run_id, outcome), by = "run_id")

if (nrow(fold_metrics) > 0) {
  fold_long <- fold_metrics |>
    dplyr::select(run_id, outcome, outer_fold, accuracy, balanced_accuracy, f1_macro, mcc) |>
    tidyr::pivot_longer(
      cols = c(accuracy, balanced_accuracy, f1_macro, mcc),
      names_to = "metric",
      values_to = "value"
    )

  p_folds <- ggplot(fold_long, aes(x = outer_fold, y = value, group = 1)) +
    geom_point(size = 2, na.rm = TRUE) +
    geom_line(na.rm = TRUE) +
    facet_grid(metric ~ outcome, scales = "free_y") +
    coord_cartesian(ylim = c(0, NA)) +
    labs(x = "Held-out study", y = NULL, title = "LOSO Fold Metrics") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(out_paths$results_plots, "elastic_net_loso_fold_metrics.png"),
    plot = p_folds,
    width = 11,
    height = 8,
    dpi = 300
  )

  tune_long <- fold_metrics |>
    dplyr::mutate(log10_lambda = log10(lambda)) |>
    dplyr::select(run_id, outcome, outer_fold, alpha, log10_lambda) |>
    tidyr::pivot_longer(
      cols = c(alpha, log10_lambda),
      names_to = "parameter",
      values_to = "value"
    )

  p_tuning <- ggplot(tune_long, aes(x = outer_fold, y = value, group = 1)) +
    geom_point(size = 2, na.rm = TRUE) +
    geom_line(na.rm = TRUE) +
    facet_grid(parameter ~ outcome, scales = "free_y") +
    labs(x = "Held-out study", y = NULL, title = "Selected Elastic Net Parameters by LOSO Fold") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(out_paths$results_plots, "elastic_net_loso_tuning_by_fold.png"),
    plot = p_tuning,
    width = 11,
    height = 7,
    dpi = 300
  )
}

nonzero_summary <- purrr::map_dfr(metrics$run_id, function(run_id) {
  freq_path <- file.path(out_paths$results_features, paste0("elastic_net_loso_nonzero_feature_frequency_", run_id, ".csv"))
  coef_path <- file.path(out_paths$results_features, paste0("elastic_net_loso_nonzero_coefficients_", run_id, ".csv"))
  if (!file.exists(freq_path) || !file.exists(coef_path)) return(tibble::tibble())

  freq <- readr::read_csv(freq_path, show_col_types = FALSE)
  coef <- readr::read_csv(coef_path, show_col_types = FALSE)

  tibble::tibble(
    run_id = run_id,
    unique_nonzero_cpgs = nrow(freq),
    nonzero_cpgs_all_folds = sum(freq$nonzero_folds == 5),
    nonzero_cpgs_at_least_4_folds = sum(freq$nonzero_folds >= 4),
    nonzero_cpg_terms = sum(coef$feature_type == "CpG")
  )
}) |>
  dplyr::left_join(metrics |> dplyr::select(run_id, outcome), by = "run_id")

if (nrow(nonzero_summary) > 0) {
  readr::write_csv(
    nonzero_summary,
    file.path(out_paths$results_metrics, "elastic_net_loso_nonzero_cpg_summary.csv")
  )

  nonzero_long <- nonzero_summary |>
    tidyr::pivot_longer(
      cols = c(unique_nonzero_cpgs, nonzero_cpgs_all_folds, nonzero_cpgs_at_least_4_folds),
      names_to = "metric",
      values_to = "n"
    )

  p_nonzero <- ggplot(nonzero_long, aes(x = outcome, y = n, fill = outcome)) +
    geom_col(width = 0.72) +
    facet_wrap(~metric, scales = "free_y") +
    labs(x = NULL, y = "CpGs", title = "Non-zero Elastic Net CpG Counts") +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(out_paths$results_plots, "elastic_net_loso_nonzero_cpg_summary.png"),
    plot = p_nonzero,
    width = 11,
    height = 8,
    dpi = 300
  )
}

confusion_tbl <- purrr::map_dfr(metrics$run_id, function(run_id) {
  path_in <- file.path(out_paths$results_metrics, paste0("elastic_net_loso_confusion_matrix_", run_id, ".csv"))
  if (!file.exists(path_in)) return(tibble::tibble())
  readr::read_csv(path_in, show_col_types = FALSE) |>
    dplyr::mutate(run_id = run_id)
}) |>
  dplyr::left_join(metrics |> dplyr::select(run_id, outcome), by = "run_id") |>
  dplyr::group_by(run_id, Reference) |>
  dplyr::mutate(reference_fraction = Freq / sum(Freq)) |>
  dplyr::ungroup()

if (nrow(confusion_tbl) > 0) {
  p_confusion <- ggplot(confusion_tbl, aes(x = Reference, y = Prediction, fill = reference_fraction)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = Freq), size = 2.6) +
    facet_wrap(~run_id, scales = "free", ncol = 2) +
    scale_fill_viridis_c(option = "C", limits = c(0, 1)) +
    labs(x = "True class", y = "Predicted class", fill = "Fraction\nwithin true class", title = "LOSO Confusion Matrices") +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(
    filename = file.path(out_paths$results_plots, "elastic_net_loso_confusion_heatmaps.png"),
    plot = p_confusion,
    width = 13,
    height = 12,
    dpi = 300
  )
}

message("Elastic Net LOSO plots written to: ", out_paths$results_plots)
