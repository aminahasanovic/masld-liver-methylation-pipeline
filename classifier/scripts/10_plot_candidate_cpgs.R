# 10_plot_candidate_cpgs.R
# Prioritize annotated Elastic Net CpGs and plot beta-value distributions.

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

top_n <- as.integer(Sys.getenv("CLASSIFIER_CPG_PLOT_TOP_N", "24"))

if (is.na(top_n) || top_n < 1) {
  stop("CLASSIFIER_CPG_PLOT_TOP_N must be a positive integer.")
}

annotation_unique_path <- file.path(
  out_paths$results_features,
  "elastic_net_cpg_annotation_priority_unique.csv"
)

annotation_long_path <- file.path(
  out_paths$results_features,
  "elastic_net_cpg_annotation_priority_long.csv"
)

if (!file.exists(annotation_unique_path) || !file.exists(annotation_long_path)) {
  stop("Missing CpG annotation outputs. Run scripts/09_annotate_elastic_net_cpgs.R first.")
}

beta_path <- check_file(config$inputs$beta_matrix_path)
meta_path <- check_file(config$inputs$metadata_path)
sample_id_col <- config$inputs$sample_id_col %||% "Sample_Name"

message("Prioritizing and plotting annotated CpGs.")
message("Top CpGs to plot: ", top_n)

anno_unique <- readr::read_csv(annotation_unique_path, show_col_types = FALSE)
anno_long <- readr::read_csv(annotation_long_path, show_col_types = FALSE)

prioritized <- anno_unique |>
  dplyr::mutate(
    has_disease_group_4 = stringr::str_detect(evidence_outcomes, "disease_group_4"),
    has_disease_group_3 = stringr::str_detect(evidence_outcomes, "disease_group_3"),
    has_binary = stringr::str_detect(evidence_outcomes, "binary_healthy_vs_disease"),
    all5_any_run = !is.na(all5_runs) & all5_runs != "",
    priority_label = dplyr::case_when(
      has_disease_group_4 & has_binary & all5_any_run ~ "stage_and_binary_all5",
      has_disease_group_4 & has_binary ~ "stage_and_binary",
      has_disease_group_4 & has_disease_group_3 & all5_any_run ~ "stage_replicated_all5",
      has_disease_group_4 & all5_any_run ~ "stage_all5",
      has_binary & all5_any_run ~ "binary_all5",
      TRUE ~ "stable_secondary"
    ),
    priority_label_rank = dplyr::case_when(
      priority_label == "stage_and_binary_all5" ~ 1L,
      priority_label == "stage_and_binary" ~ 2L,
      priority_label == "stage_replicated_all5" ~ 3L,
      priority_label == "stage_all5" ~ 4L,
      priority_label == "binary_all5" ~ 5L,
      TRUE ~ 6L
    )
  ) |>
  dplyr::arrange(
    priority_label_rank,
    best_priority_rank,
    dplyr::desc(best_nonzero_folds),
    dplyr::desc(max_mean_abs_coefficient),
    cpg
  )

top_cpgs <- prioritized |>
  dplyr::slice_head(n = top_n) |>
  dplyr::pull(cpg)

beta <- readRDS(beta_path)
meta <- readRDS(meta_path)

all_candidate_cpgs <- intersect(prioritized$cpg, rownames(beta))
missing_cpgs <- setdiff(top_cpgs, rownames(beta))

if (length(missing_cpgs) > 0) {
  warning(
    "Dropping top CpGs missing from beta matrix: ",
    paste(missing_cpgs, collapse = ", ")
  )
}

top_cpgs <- intersect(top_cpgs, rownames(beta))

if (length(top_cpgs) == 0) {
  stop("None of the prioritized CpGs are present in the beta matrix.")
}

if (!sample_id_col %in% colnames(meta)) {
  stop("Sample ID column missing from metadata: ", sample_id_col)
}

common_samples <- intersect(colnames(beta), as.character(meta[[sample_id_col]]))
meta <- meta |>
  dplyr::filter(.data[[sample_id_col]] %in% common_samples) |>
  dplyr::mutate(sample_id = as.character(.data[[sample_id_col]])) |>
  dplyr::arrange(match(sample_id, common_samples))

beta_sub <- beta[top_cpgs, meta$sample_id, drop = FALSE]
beta_all_candidates <- beta[all_candidate_cpgs, meta$sample_id, drop = FALSE]

candidate_beta_df <- as.data.frame(t(beta_all_candidates)) |>
  tibble::rownames_to_column("sample_id") |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(all_candidate_cpgs),
    names_to = "cpg",
    values_to = "beta_value"
  ) |>
  dplyr::left_join(
    meta |>
      dplyr::select(
        sample_id,
        Dataset,
        Array,
        DiseaseGroup,
        Progression3,
        Age,
        Sex
      ),
    by = "sample_id"
  )

safe_r2 <- function(df, formula_in, adjusted = FALSE) {
  fit <- tryCatch(
    stats::lm(formula_in, data = df),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(NA_real_)
  }

  fit_summary <- summary(fit)

  if (adjusted) {
    return(unname(fit_summary$adj.r.squared))
  }

  unname(fit_summary$r.squared)
}

signal_audit <- candidate_beta_df |>
  dplyr::mutate(
    DiseaseGroup = factor(DiseaseGroup),
    Dataset = factor(Dataset),
    Array = factor(Array)
  ) |>
  dplyr::group_by(cpg) |>
  dplyr::group_modify(function(.x, .y) {
    r2_disease_group <- safe_r2(.x, beta_value ~ DiseaseGroup)
    r2_dataset <- safe_r2(.x, beta_value ~ Dataset)
    r2_array <- safe_r2(.x, beta_value ~ Array)
    r2_dataset_disease <- safe_r2(.x, beta_value ~ Dataset + DiseaseGroup)

    tibble::tibble(
      n_samples = sum(!is.na(.x$beta_value)),
      r2_disease_group = r2_disease_group,
      r2_dataset = r2_dataset,
      r2_array = r2_array,
      r2_dataset_plus_disease_group = r2_dataset_disease,
      disease_r2_after_dataset = r2_dataset_disease - r2_dataset,
      dataset_r2_after_disease = r2_dataset_disease - r2_disease_group,
      likely_dataset_sensitive = !is.na(r2_dataset) &
        !is.na(r2_disease_group) &
        r2_dataset > r2_disease_group
    )
  }) |>
  dplyr::ungroup()

prioritized <- prioritized |>
  dplyr::left_join(signal_audit, by = "cpg") |>
  dplyr::arrange(
    priority_label_rank,
    best_priority_rank,
    dplyr::desc(best_nonzero_folds),
    dplyr::desc(max_mean_abs_coefficient),
    cpg
  )

plot_df <- as.data.frame(t(beta_sub)) |>
  tibble::rownames_to_column("sample_id") |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(top_cpgs),
    names_to = "cpg",
    values_to = "beta_value"
  ) |>
  dplyr::left_join(
    meta |>
      dplyr::select(
        sample_id,
        Dataset,
        Array,
        DiseaseGroup,
        Progression3,
        Age,
        Sex
      ),
    by = "sample_id"
  ) |>
  dplyr::left_join(
    prioritized |>
      dplyr::select(
        cpg,
        priority_label,
        evidence_outcomes,
        best_nonzero_folds,
        max_mean_abs_coefficient,
        chr,
        pos,
        UCSC_RefGene_Name,
        Relation_to_Island
      ),
    by = "cpg"
  ) |>
  dplyr::mutate(
    DiseaseGroup = factor(
      DiseaseGroup,
      levels = c(
        "Healthy",
        "Healthy_Obese",
        "MASL",
        "MASH",
        "Mild_Fibrosis",
        "Advanced_Fibrosis"
      )
    ),
    Dataset = factor(Dataset),
    Array = factor(Array),
    cpg_label = dplyr::if_else(
      !is.na(UCSC_RefGene_Name) & UCSC_RefGene_Name != "",
      paste0(cpg, "\n", stringr::str_trunc(UCSC_RefGene_Name, 28)),
      cpg
    )
  )

candidate_path <- file.path(out_paths$results_features, "elastic_net_cpg_candidate_priority.csv")
signal_audit_path <- file.path(out_paths$results_features, "elastic_net_cpg_candidate_signal_audit.csv")
plot_data_path <- file.path(out_paths$results_features, "elastic_net_top_cpg_beta_values.csv")

readr::write_csv(prioritized, candidate_path)
readr::write_csv(signal_audit, signal_audit_path)
readr::write_csv(plot_df, plot_data_path)

disease_group_palette <- c(
  Healthy = "#2C7BB6",
  Healthy_Obese = "#74ADD1",
  MASL = "#FDAE61",
  MASH = "#F46D43",
  Mild_Fibrosis = "#B2182B",
  Advanced_Fibrosis = "#5E3C99"
)

dataset_palette <- c(
  Ahrens = "#1B9E77",
  Horvath = "#D95F02",
  Johnson = "#7570B3",
  Murphy = "#E7298A",
  VanDijck = "#66A61E",
  ITEN = "#A6761D"
)

theme_cpg <- function() {
  ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 8, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1),
      legend.position = "bottom"
    )
}

p_group <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = DiseaseGroup, y = beta_value, color = Dataset)
) +
  ggplot2::geom_boxplot(
    ggplot2::aes(group = DiseaseGroup),
    outlier.shape = NA,
    color = "grey35",
    fill = "grey95",
    linewidth = 0.25
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.18, height = 0, seed = 123),
    size = 0.8,
    alpha = 0.65
  ) +
  ggplot2::facet_wrap(~ cpg_label, scales = "free_y", ncol = 4) +
  ggplot2::scale_color_manual(values = dataset_palette, drop = FALSE) +
  ggplot2::labs(
    x = NULL,
    y = "Beta value",
    color = "Dataset",
    title = "Prioritized Elastic Net CpGs by DiseaseGroup",
    subtitle = "Points are samples; color marks cohort/dataset to reveal possible cohort-driven signals."
  ) +
  theme_cpg()

p_dataset <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = Dataset, y = beta_value, color = DiseaseGroup)
) +
  ggplot2::geom_boxplot(
    ggplot2::aes(group = Dataset),
    outlier.shape = NA,
    color = "grey35",
    fill = "grey95",
    linewidth = 0.25
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.18, height = 0, seed = 123),
    size = 0.8,
    alpha = 0.65
  ) +
  ggplot2::facet_wrap(~ cpg_label, scales = "free_y", ncol = 4) +
  ggplot2::scale_color_manual(values = disease_group_palette, drop = FALSE) +
  ggplot2::labs(
    x = NULL,
    y = "Beta value",
    color = "DiseaseGroup",
    title = "Prioritized Elastic Net CpGs by Dataset",
    subtitle = "If separation follows dataset more than biology, this supports a cross-study generalization problem."
  ) +
  theme_cpg()

mean_df <- plot_df |>
  dplyr::group_by(cpg, cpg_label, Dataset, DiseaseGroup) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_beta = mean(beta_value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(n >= 2)

p_mean <- ggplot2::ggplot(
  mean_df,
  ggplot2::aes(x = DiseaseGroup, y = Dataset, fill = mean_beta)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.2) +
  ggplot2::facet_wrap(~ cpg_label, ncol = 4) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0.5,
    limits = c(0, 1),
    oob = scales::squish
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "Mean beta",
    title = "Mean beta values by Dataset and DiseaseGroup",
    subtitle = "Cells with fewer than two samples are omitted."
  ) +
  theme_cpg()

group_plot_path <- file.path(out_paths$results_plots, "elastic_net_top_cpg_beta_by_disease_group.png")
dataset_plot_path <- file.path(out_paths$results_plots, "elastic_net_top_cpg_beta_by_dataset.png")
mean_plot_path <- file.path(out_paths$results_plots, "elastic_net_top_cpg_mean_beta_heatmap.png")

ggplot2::ggsave(group_plot_path, p_group, width = 14, height = 12, dpi = 300)
ggplot2::ggsave(dataset_plot_path, p_dataset, width = 14, height = 12, dpi = 300)
ggplot2::ggsave(mean_plot_path, p_mean, width = 14, height = 12, dpi = 300)

message("Candidate priority table: ", candidate_path)
message("Candidate signal audit: ", signal_audit_path)
message("Top-CpG plot data: ", plot_data_path)
message("DiseaseGroup beta plot: ", group_plot_path)
message("Dataset beta plot: ", dataset_plot_path)
message("Mean beta heatmap: ", mean_plot_path)

message("Top prioritized CpGs:")
print(
  prioritized |>
    dplyr::select(
      cpg,
      priority_label,
      evidence_outcomes,
      best_nonzero_folds,
      max_mean_abs_coefficient,
      chr,
      pos,
      UCSC_RefGene_Name,
      Relation_to_Island
    ) |>
    utils::head(top_n)
)
