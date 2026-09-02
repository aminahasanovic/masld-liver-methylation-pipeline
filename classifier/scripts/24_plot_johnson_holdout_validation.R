# 24_plot_johnson_holdout_validation.R
# Summarize and plot Johnson holdout validation after cohort updates.

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(tidyr)
})

outcomes <- c("disease_group_4", "disease_group_5_obesity_split")
outcome_labels <- c(
  disease_group_4 = "4-class",
  disease_group_5_obesity_split = "5-class"
)
outcome_levels <- list(
  disease_group_4 = c("Healthy", "MASL_MASH", "Mild_Fibrosis", "Advanced_Fibrosis"),
  disease_group_5_obesity_split = c("Healthy", "Healthy_Obese", "MASL_MASH", "Mild_Fibrosis", "Advanced_Fibrosis")
)

holdouts <- c("johnson_nonfibrotic", "johnson_all")
holdout_labels <- c(
  johnson_nonfibrotic = "Johnson Healthy_Obese only",
  johnson_all = "All Johnson"
)
holdout_suffix <- c(
  johnson_nonfibrotic = "johnson_nonfibroticholdout",
  johnson_all = "johnson_allholdout"
)

models <- c("elastic_net", "random_forest", "boosted_trees")
model_labels <- c(
  elastic_net = "Elastic Net",
  random_forest = "Random Forest",
  boosted_trees = "Boosted Trees"
)
model_label_levels <- unname(model_labels[models])

run_id_for <- function(outcome, holdout) {
  paste(outcome, "cpg_only", "train_test", unname(holdout_suffix[[holdout]]), sep = "_")
}

read_optional_csv <- function(path) {
  if (!file.exists(path)) {
    warning("Missing file: ", path)
    return(tibble::tibble())
  }

  readr::read_csv(path, show_col_types = FALSE)
}

read_metric <- function(outcome, holdout, model) {
  outcome_value <- outcome
  holdout_value <- holdout
  model_value <- model
  run_id <- run_id_for(outcome_value, holdout_value)
  path <- file.path(out_paths$results_metrics, paste0(model, "_metrics_", run_id, ".csv"))

  read_optional_csv(path) |>
    dplyr::mutate(
      holdout = holdout_value,
      holdout_label = unname(holdout_labels[[holdout_value]]),
      outcome = outcome_value,
      outcome_label = unname(outcome_labels[[outcome_value]]),
      model = model_value,
      model_label = unname(model_labels[[model_value]]),
      run_id = run_id,
      .before = 1
    )
}

read_holdout_predictions <- function(outcome, holdout) {
  outcome_value <- outcome
  holdout_value <- holdout
  run_id <- run_id_for(outcome_value, holdout_value)
  path <- file.path(out_paths$results_metrics, paste0("holdout_prediction_by_source_summary_", run_id, ".csv"))

  read_optional_csv(path) |>
    dplyr::mutate(
      holdout = holdout_value,
      holdout_label = unname(holdout_labels[[holdout_value]]),
      outcome = outcome_value,
      outcome_label = unname(outcome_labels[[outcome_value]]),
      run_id = run_id,
      model_label = unname(model_labels[as.character(model)]),
      source_label = dplyr::case_when(
        source_disease_group == "Healthy_Obese" ~ "Johnson Healthy_Obese",
        source_disease_group == "Advanced_Fibrosis" ~ "Johnson Advanced_Fibrosis",
        TRUE ~ as.character(source_disease_group)
      ),
      .before = 1
    )
}

read_holdout_probabilities <- function(outcome, holdout) {
  outcome_value <- outcome
  holdout_value <- holdout
  run_id <- run_id_for(outcome_value, holdout_value)
  path <- file.path(out_paths$results_metrics, paste0("holdout_probability_by_source_summary_", run_id, ".csv"))

  read_optional_csv(path) |>
    dplyr::mutate(
      holdout = holdout_value,
      holdout_label = unname(holdout_labels[[holdout_value]]),
      outcome = outcome_value,
      outcome_label = unname(outcome_labels[[outcome_value]]),
      run_id = run_id,
      model_label = unname(model_labels[as.character(model)]),
      source_label = dplyr::case_when(
        source_disease_group == "Healthy_Obese" ~ "Johnson Healthy_Obese",
        source_disease_group == "Advanced_Fibrosis" ~ "Johnson Advanced_Fibrosis",
        TRUE ~ as.character(source_disease_group)
      ),
      .before = 1
    )
}

internal_metrics <- dplyr::bind_rows(lapply(outcomes, function(outcome) {
  dplyr::bind_rows(lapply(holdouts, function(holdout) {
    dplyr::bind_rows(lapply(models, function(model) {
      read_metric(outcome, holdout, model)
    }))
  }))
})) |>
  dplyr::select(
    holdout,
    holdout_label,
    outcome,
    outcome_label,
    model,
    model_label,
    run_id,
    accuracy,
    balanced_accuracy,
    dplyr::any_of("macro_recall_present"),
    mcc,
    f1_macro,
    dplyr::any_of("f1_macro_present")
  ) |>
  dplyr::arrange(holdout, outcome, model)

prediction_summary <- dplyr::bind_rows(lapply(outcomes, function(outcome) {
  dplyr::bind_rows(lapply(holdouts, function(holdout) {
    read_holdout_predictions(outcome, holdout)
  }))
})) |>
  dplyr::mutate(
    holdout_label = factor(holdout_label, levels = unname(holdout_labels[holdouts])),
    outcome_label = factor(outcome_label, levels = unname(outcome_labels[outcomes])),
    model_label = factor(model_label, levels = model_label_levels),
    source_label = factor(source_label, levels = c("Johnson Healthy_Obese", "Johnson Advanced_Fibrosis"))
  ) |>
  dplyr::arrange(holdout_label, outcome_label, source_label, model_label, prediction)

probability_summary <- dplyr::bind_rows(lapply(outcomes, function(outcome) {
  dplyr::bind_rows(lapply(holdouts, function(holdout) {
    read_holdout_probabilities(outcome, holdout)
  }))
})) |>
  dplyr::mutate(
    holdout_label = factor(holdout_label, levels = unname(holdout_labels[holdouts])),
    outcome_label = factor(outcome_label, levels = unname(outcome_labels[outcomes])),
    model_label = factor(model_label, levels = model_label_levels),
    source_label = factor(source_label, levels = c("Johnson Healthy_Obese", "Johnson Advanced_Fibrosis"))
  ) |>
  dplyr::arrange(holdout_label, outcome_label, source_label, model_label, class)

top_predictions <- prediction_summary |>
  dplyr::group_by(holdout_label, outcome_label, source_label, model_label) |>
  dplyr::slice_max(order_by = fraction, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(
    holdout_label,
    outcome_label,
    source_label,
    model_label,
    top_prediction = prediction,
    n,
    total,
    fraction
  )

dataset_counts <- readRDS(config$inputs$metadata_path) |>
  dplyr::count(Dataset, DiseaseGroup, name = "n") |>
  dplyr::arrange(Dataset, DiseaseGroup)

capture_tbl <- function(x) {
  paste(capture.output(print(tibble::as_tibble(x), n = Inf, width = Inf)), collapse = "\n")
}

readr::write_csv(
  internal_metrics,
  file.path(out_paths$results_metrics, "johnson_validation_after_kim_internal_metrics.csv")
)
readr::write_csv(
  prediction_summary,
  file.path(out_paths$results_metrics, "johnson_validation_after_kim_prediction_summary.csv")
)
readr::write_csv(
  probability_summary,
  file.path(out_paths$results_metrics, "johnson_validation_after_kim_probability_summary.csv")
)
readr::write_csv(
  top_predictions,
  file.path(out_paths$results_metrics, "johnson_validation_after_kim_top_predictions.csv")
)
readr::write_csv(
  dataset_counts,
  file.path(out_paths$results_metrics, "johnson_validation_after_kim_dataset_counts.csv")
)

percent_label <- function(x) {
  paste0(round(100 * x), "%")
}

clean_class_label <- function(x) {
  dplyr::recode(
    as.character(x),
    Healthy = "Healthy",
    Healthy_Obese = "Healthy obese",
    MASL_MASH = "MASL/MASH",
    Mild_Fibrosis = "Mild fibrosis",
    Advanced_Fibrosis = "Advanced fibrosis",
    .default = gsub("_", " ", as.character(x))
  )
}

plot_heatmap <- function(data, filename, height = 8) {
  p <- ggplot(
    data |> dplyr::mutate(prediction_label = clean_class_label(prediction)),
    aes(x = prediction_label, y = model_label, fill = fraction)
  ) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = percent_label(fraction)), size = 3) +
    facet_grid(source_label ~ outcome_label, scales = "free_x", space = "free_x") +
    scale_fill_gradient(
      low = "#F7FBFF",
      high = "#08519C",
      limits = c(0, 1),
      labels = percent_label
    ) +
    labs(
      x = "Predicted class",
      y = "Model",
      fill = "Fraction"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggplot2::ggsave(
    filename = file.path(out_paths$results_plots, filename),
    plot = p,
    width = 12,
    height = height,
    dpi = 300
  )

  ggplot2::ggsave(
    filename = file.path(out_paths$results_plots, sub("\\.png$", ".pdf", filename)),
    plot = p,
    width = 12,
    height = height
  )
}

plot_heatmap(
  data = prediction_summary |> dplyr::filter(holdout == "johnson_nonfibrotic"),
  filename = "johnson_nonfibrotic_holdout_predictions_after_kim.png",
  height = 5.5
)

plot_heatmap(
  data = prediction_summary |> dplyr::filter(holdout == "johnson_all"),
  filename = "johnson_all_holdout_predictions_after_kim.png",
  height = 7
)

f27a_class_levels <- c("Healthy", "MASL/MASH", "Mild fibrosis", "Advanced fibrosis")

p_all_johnson_4class <- prediction_summary |>
  dplyr::filter(holdout == "johnson_all", outcome == "disease_group_4") |>
  dplyr::mutate(
    prediction_label = factor(
      clean_class_label(prediction),
      levels = f27a_class_levels
    )
  ) |>
  tidyr::complete(
    source_label,
    model_label,
    prediction_label = factor(f27a_class_levels, levels = f27a_class_levels),
    fill = list(fraction = 0)
  ) |>
  ggplot(aes(x = prediction_label, y = model_label, fill = fraction)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = percent_label(fraction)), size = 3.2) +
  facet_wrap(~source_label, nrow = 1) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#08519C",
    limits = c(0, 1),
    labels = percent_label
  ) +
  labs(
    x = "Predicted class",
    y = "Model",
    fill = "Fraction"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = file.path(out_paths$results_plots, "johnson_all_holdout_4class_predictions_after_kim.png"),
  plot = p_all_johnson_4class,
  width = 8.8,
  height = 4.6,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(out_paths$results_plots, "johnson_all_holdout_4class_predictions_after_kim.pdf"),
  plot = p_all_johnson_4class,
  width = 8.8,
  height = 4.6
)

metrics_long <- internal_metrics |>
  dplyr::select(holdout_label, outcome_label, model_label, accuracy, balanced_accuracy, f1_macro, mcc) |>
  tidyr::pivot_longer(
    cols = c(accuracy, balanced_accuracy, f1_macro, mcc),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("accuracy", "balanced_accuracy", "f1_macro", "mcc"),
      labels = c("Accuracy", "Balanced accuracy", "Macro F1", "MCC")
    ),
    holdout_label = factor(holdout_label, levels = unname(holdout_labels[holdouts])),
    outcome_label = factor(outcome_label, levels = unname(outcome_labels[outcomes])),
    model_label = factor(model_label, levels = model_label_levels)
  )

p_metrics <- ggplot(metrics_long, aes(x = model_label, y = value, fill = model_label)) +
  geom_col(width = 0.72, na.rm = TRUE) +
  facet_grid(metric ~ holdout_label + outcome_label, scales = "free_y") +
  labs(
    x = "Model",
    y = "Metric value",
    fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = file.path(out_paths$results_plots, "johnson_holdout_internal_metrics_after_kim.png"),
  plot = p_metrics,
  width = 13,
  height = 8,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(out_paths$results_plots, "johnson_holdout_internal_metrics_after_kim.pdf"),
  plot = p_metrics,
  width = 13,
  height = 8
)

report_lines <- c(
  "# Johnson validation after ITEN cohort update",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "The input metadata now contains the ITEN cohort, adding 106 samples, including 17 MASL and 58 MASH samples.",
  "",
  "## Dataset counts",
  "",
  capture_tbl(dataset_counts),
  "",
  "## Internal non-holdout metrics",
  "",
  capture_tbl(internal_metrics),
  "",
  "## Top Johnson holdout predictions",
  "",
  capture_tbl(top_predictions),
  "",
  "## Interpretation notes",
  "",
  "- In the Johnson Healthy_Obese-only holdout, Johnson Advanced_Fibrosis remains in training, so strong Advanced_Fibrosis predictions can still reflect Johnson-specific cohort signal.",
  "- In the all-Johnson holdout, no Johnson samples are used for training; this is stricter for external-like validation but tests transfer to a cohort with strong domain shift.",
  "- The ITEN cohort improves MASL/MASH representation in the non-Johnson training data, so shifts toward MASL_MASH in the all-Johnson holdout are especially informative.",
  "",
  "## Output files",
  "",
  "- `results/metrics/johnson_validation_after_kim_internal_metrics.csv`",
  "- `results/metrics/johnson_validation_after_kim_prediction_summary.csv`",
  "- `results/metrics/johnson_validation_after_kim_probability_summary.csv`",
  "- `results/metrics/johnson_validation_after_kim_top_predictions.csv`",
  "- `results/plots/johnson_nonfibrotic_holdout_predictions_after_kim.png`",
  "- `results/plots/johnson_all_holdout_predictions_after_kim.png`",
  "- `results/plots/johnson_holdout_internal_metrics_after_kim.png`"
)

writeLines(
  report_lines,
  file.path(out_paths$results_metrics, "johnson_validation_after_kim_summary.md")
)

message("Johnson validation plots and tables complete.")
message("Saved summary to: ", file.path(out_paths$results_metrics, "johnson_validation_after_kim_summary.md"))
message("Saved plots to: ", out_paths$results_plots)
print(top_predictions, n = Inf)
