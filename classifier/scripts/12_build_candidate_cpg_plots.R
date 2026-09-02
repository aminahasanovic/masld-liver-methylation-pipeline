source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

top_n <- as.integer(Sys.getenv("CLASSIFIER_CANDIDATE_TOP_N", "15"))
top_cpg_override <- Sys.getenv("CLASSIFIER_CANDIDATE_TOP_CPGS", "")
plot_suffix <- sanitize_id(Sys.getenv("CLASSIFIER_CANDIDATE_PLOT_SUFFIX", "top15"))

if (is.na(top_n) || top_n < 1) {
  stop("CLASSIFIER_CANDIDATE_TOP_N must be a positive integer.")
}

meta_results_path <- file.path(out_paths$results_features, "elastic_net_candidate_cpg_meta_results.csv")
dataset_effects_path <- file.path(out_paths$results_features, "elastic_net_candidate_cpg_meta_dataset_effects.csv")

if (!file.exists(meta_results_path) || !file.exists(dataset_effects_path)) {
  stop("Missing meta-analysis outputs. Run scripts/11_meta_analyze_candidate_cpgs.R first.")
}

beta_path <- check_file(config$inputs$beta_matrix_path)
meta_path <- check_file(config$inputs$metadata_path)
sample_id_col <- config$inputs$sample_id_col %||% "Sample_Name"

message("Building CpG candidate summary.")
message("Top CpGs to plot: ", top_n)

meta_results <- readr::read_csv(meta_results_path, show_col_types = FALSE)
dataset_effects <- readr::read_csv(dataset_effects_path, show_col_types = FALSE)

robust_results <- meta_results |>
  dplyr::filter(evidence_tier != "not_robust")

if (nrow(robust_results) == 0) {
  stop("No robust CpG/contrast rows found in meta-analysis results.")
}

collapse_unique <- function(x, sep = ";") {
  x <- as.character(x)
  x <- unlist(strsplit(x, sep, fixed = TRUE), use.names = FALSE)
  x <- x[!is.na(x) & nzchar(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(sort(unique(x)), collapse = sep)
}

first_non_empty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[[1]]
}

model_priority_score <- function(priority_label) {
  dplyr::case_when(
    priority_label == "stage_and_binary_all5" ~ 5,
    priority_label == "stage_and_binary" ~ 4,
    priority_label == "stage_replicated_all5" ~ 3,
    priority_label == "stage_all5" ~ 2,
    priority_label == "binary_all5" ~ 2,
    TRUE ~ 1
  )
}

effect_direction <- function(estimate) {
  dplyr::case_when(
    estimate > 0 ~ "higher methylation in case/higher stage",
    estimate < 0 ~ "lower methylation in case/higher stage",
    TRUE ~ "no direction"
  )
}

candidate_summary <- robust_results |>
  dplyr::group_by(cpg) |>
  dplyr::summarise(
    robust_contrasts = dplyr::n_distinct(contrast),
    tier1_contrasts = sum(evidence_tier == "tier1_consistent_low_dataset_signal"),
    tier2_contrasts = sum(evidence_tier == "tier2_consistent_exploratory"),
    tier3_contrasts = sum(evidence_tier == "tier3_two_dataset_signal"),
    same_direction_contrasts = sum(direction_consistency == 1),
    max_n_datasets = max(n_datasets, na.rm = TRUE),
    all_datasets_seen = collapse_unique(datasets),
    best_fdr = min(fdr, na.rm = TRUE),
    best_neg_log10_fdr = -log10(pmax(best_fdr, .Machine$double.xmin)),
    max_abs_meta_effect = max(abs(random_effect), na.rm = TRUE),
    mean_i2 = mean(i2, na.rm = TRUE),
    min_i2 = min(i2, na.rm = TRUE),
    top_contrast = contrast[which.min(fdr)][[1]],
    top_random_effect = random_effect[which.min(fdr)][[1]],
    top_effect_direction = effect_direction(top_random_effect),
    robust_contrast_summary = paste(
      paste0(
        contrast,
        ": ",
        ifelse(random_effect > 0, "+", "-"),
        sprintf("%.3f", abs(random_effect)),
        ", FDR=",
        formatC(fdr, format = "e", digits = 2),
        ", I2=",
        sprintf("%.1f", i2)
      ),
      collapse = " | "
    ),
    evidence_outcomes = first_non_empty(evidence_outcomes),
    evidence_runs = first_non_empty(evidence_runs),
    priority_label = first_non_empty(priority_label),
    best_nonzero_folds = max(best_nonzero_folds, na.rm = TRUE),
    all5_any_run = any(all5_any_run, na.rm = TRUE),
    max_mean_abs_coefficient = max(max_mean_abs_coefficient, na.rm = TRUE),
    likely_dataset_sensitive = any(likely_dataset_sensitive, na.rm = TRUE),
    r2_disease_group = dplyr::first(r2_disease_group),
    r2_dataset = dplyr::first(r2_dataset),
    disease_r2_after_dataset = dplyr::first(disease_r2_after_dataset),
    annotation_source = first_non_empty(annotation_source),
    in_450k_annotation = any(in_450k_annotation, na.rm = TRUE),
    in_epic_annotation = any(in_epic_annotation, na.rm = TRUE),
    chr = first_non_empty(chr),
    pos = dplyr::first(pos),
    strand = first_non_empty(strand),
    UCSC_RefGene_Name = first_non_empty(UCSC_RefGene_Name),
    UCSC_RefGene_Group = first_non_empty(UCSC_RefGene_Group),
    Relation_to_Island = first_non_empty(Relation_to_Island),
    Regulatory_Feature_Group = first_non_empty(Regulatory_Feature_Group),
    DMR = first_non_empty(DMR),
    Enhancer = first_non_empty(Enhancer),
    DHS = first_non_empty(DHS),
    Probe_rs = first_non_empty(Probe_rs),
    Probe_maf = dplyr::first(Probe_maf),
    CpG_rs = first_non_empty(CpG_rs),
    CpG_maf = dplyr::first(CpG_maf),
    SBE_rs = first_non_empty(SBE_rs),
    SBE_maf = dplyr::first(SBE_maf),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    model_priority_score = model_priority_score(priority_label),
    dataset_warning = dplyr::case_when(
      likely_dataset_sensitive ~ "dataset_r2_gt_disease_r2",
      mean_i2 >= 50 ~ "high_heterogeneity",
      TRUE ~ "none"
    ),
    has_snp_annotation = (!is.na(Probe_rs) & Probe_rs != "") |
      (!is.na(CpG_rs) & CpG_rs != "") |
      (!is.na(SBE_rs) & SBE_rs != ""),
    candidate_priority_score =
      8 * robust_contrasts +
      5 * tier1_contrasts +
      3 * same_direction_contrasts +
      2 * model_priority_score +
      1.5 * pmin(max_n_datasets, 5) +
      0.7 * pmin(best_neg_log10_fdr, 50) +
      1.5 * as.integer(all5_any_run) -
      0.05 * pmin(mean_i2, 100) -
      4 * as.integer(likely_dataset_sensitive),
    candidate_selection_note = paste(
      "Selected/ranked from robust meta-analysis CpGs; not a final validated biomarker.",
      "Score favors multiple robust contrasts, Tier-1 evidence, same direction, more datasets, and Elastic Net stability."
    )
  ) |>
  dplyr::arrange(
    dplyr::desc(candidate_priority_score),
    dplyr::desc(robust_contrasts),
    best_fdr,
    dplyr::desc(max_abs_meta_effect),
    cpg
  )

if (nzchar(top_cpg_override)) {
  override_cpgs <- strsplit(top_cpg_override, ",", fixed = TRUE)[[1]]
  override_cpgs <- trimws(override_cpgs)
  override_cpgs <- override_cpgs[nzchar(override_cpgs)]

  top_candidates <- candidate_summary |>
    dplyr::mutate(override_rank = match(cpg, override_cpgs)) |>
    dplyr::filter(!is.na(override_rank)) |>
    dplyr::arrange(override_rank) |>
    dplyr::select(-override_rank)

  missing_override_cpgs <- setdiff(override_cpgs, top_candidates$cpg)

  if (length(missing_override_cpgs) > 0) {
    warning(
      "These override CpGs were not found in the candidate summary: ",
      paste(missing_override_cpgs, collapse = ", ")
    )
  }
} else {
  top_candidates <- candidate_summary |>
    dplyr::slice_head(n = top_n)
}

summary_path <- file.path(out_paths$results_features, "elastic_net_cpg_candidate_summary.csv")
top_file_name <- if (nzchar(top_cpg_override)) {
  paste0("elastic_net_cpg_candidate_", plot_suffix, "_plot_selection.csv")
} else {
  paste0("elastic_net_cpg_candidate_", plot_suffix, ".csv")
}
top_path <- file.path(out_paths$results_features, top_file_name)

readr::write_csv(candidate_summary, summary_path)
readr::write_csv(top_candidates, top_path)

beta <- readRDS(beta_path)
meta <- readRDS(meta_path)

if (!sample_id_col %in% colnames(meta)) {
  stop("Sample ID column missing from metadata: ", sample_id_col)
}

common_samples <- intersect(colnames(beta), as.character(meta[[sample_id_col]]))

meta <- meta |>
  dplyr::filter(.data[[sample_id_col]] %in% common_samples) |>
  dplyr::mutate(
    sample_id = as.character(.data[[sample_id_col]]),
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
    Array = factor(Array)
  ) |>
  dplyr::arrange(match(sample_id, common_samples))

top_cpgs <- intersect(top_candidates$cpg, rownames(beta))
beta_top <- beta[top_cpgs, meta$sample_id, drop = FALSE]

plot_df <- as.data.frame(t(beta_top)) |>
  tibble::rownames_to_column("sample_id") |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(top_cpgs),
    names_to = "cpg",
    values_to = "beta_value"
  ) |>
  dplyr::left_join(
    meta |>
      dplyr::select(sample_id, Dataset, Array, DiseaseGroup, Progression3, Age, Sex),
    by = "sample_id"
  ) |>
  dplyr::left_join(
    top_candidates |>
      dplyr::select(
        cpg,
        UCSC_RefGene_Name,
        top_contrast,
        top_effect_direction,
        robust_contrasts,
        candidate_priority_score,
        dataset_warning
      ),
    by = "cpg"
  ) |>
  dplyr::mutate(
    gene_short = dplyr::if_else(
      !is.na(UCSC_RefGene_Name) & UCSC_RefGene_Name != "",
      stringr::str_extract(UCSC_RefGene_Name, "^[^;]+"),
      "intergenic/no_gene"
    ),
    cpg_label = paste0(cpg, "\n", gene_short)
  )

plot_data_path <- file.path(
  out_paths$results_features,
  paste0("elastic_net_cpg_candidate_", plot_suffix, "_beta_values.csv")
)
readr::write_csv(plot_df, plot_data_path)

dataset_palette <- c(
  Ahrens = "#1B9E77",
  Horvath = "#D95F02",
  Johnson = "#7570B3",
  Murphy = "#E7298A",
  VanDijck = "#66A61E",
  ITEN = "#A6761D"
)

disease_group_palette <- c(
  Healthy = "#2C7BB6",
  Healthy_Obese = "#74ADD1",
  MASL = "#FDAE61",
  MASH = "#F46D43",
  Mild_Fibrosis = "#B2182B",
  Advanced_Fibrosis = "#5E3C99"
)

theme_pipeline <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size - 1, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1),
      legend.position = "bottom"
    )
}

p_top_group <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = DiseaseGroup, y = beta_value, color = Dataset)
) +
  ggplot2::geom_boxplot(
    ggplot2::aes(group = DiseaseGroup),
    outlier.shape = NA,
    color = "grey35",
    fill = "grey96",
    linewidth = 0.25
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.16, height = 0, seed = 123),
    size = 0.75,
    alpha = 0.65
  ) +
  ggplot2::facet_wrap(~ cpg_label, scales = "free_y", ncol = 5) +
  ggplot2::scale_color_manual(values = dataset_palette, drop = FALSE) +
  ggplot2::labs(
    x = NULL,
    y = "Beta value",
    color = "Dataset",
    title = "Top CpG candidates by DiseaseGroup",
    subtitle = "Candidates ranked from robust meta-analysis and Elastic Net stability; colors show dataset composition."
  ) +
  theme_pipeline()

top_effects <- dataset_effects |>
  dplyr::filter(cpg %in% top_cpgs) |>
  dplyr::left_join(
    top_candidates |>
      dplyr::select(cpg, UCSC_RefGene_Name),
    by = "cpg"
  ) |>
  dplyr::mutate(
    gene_short = dplyr::if_else(
      !is.na(UCSC_RefGene_Name) & UCSC_RefGene_Name != "",
      stringr::str_extract(UCSC_RefGene_Name, "^[^;]+"),
      "intergenic/no_gene"
    ),
    cpg_label = paste0(cpg, "\n", gene_short),
    ci_low = estimate - 1.96 * std_error,
    ci_high = estimate + 1.96 * std_error,
    contrast = factor(
      contrast,
      levels = c(
        "disease_vs_healthy",
        "fibrosis_vs_nonfibrosis",
        "advanced_vs_nonadvanced",
        "progression_score"
      )
    )
  )

p_effects <- ggplot2::ggplot(
  top_effects,
  ggplot2::aes(x = estimate, y = dataset, color = contrast)
) +
  ggplot2::geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    width = 0.15,
    linewidth = 0.25,
    alpha = 0.8
  ) +
  ggplot2::geom_point(size = 1.2, alpha = 0.9) +
  ggplot2::facet_wrap(~ cpg_label, scales = "free_x", ncol = 5) +
  ggplot2::labs(
    x = "Dataset-level effect estimate",
    y = NULL,
    color = "Contrast",
    title = "Dataset-level effects for top CpG candidates",
    subtitle = "Positive estimates mean higher methylation in disease/fibrosis/higher progression score."
  ) +
  theme_pipeline()

combined_group_path <- file.path(
  out_paths$results_plots,
  paste0("elastic_net_cpg_candidate_", plot_suffix, "_beta_by_disease_group.png")
)
combined_effect_path <- file.path(
  out_paths$results_plots,
  paste0("elastic_net_cpg_candidate_", plot_suffix, "_dataset_effects.png")
)

ggplot2::ggsave(combined_group_path, p_top_group, width = 15, height = 10, dpi = 300)
ggplot2::ggsave(combined_effect_path, p_effects, width = 15, height = 10, dpi = 300)

single_plot_dir <- file.path(out_paths$results_plots, paste0("candidate_", plot_suffix, "_cpgs"))
fs::dir_create(single_plot_dir)

make_single_plot <- function(cpg_id) {
  cpg_info <- top_candidates |>
    dplyr::filter(cpg == cpg_id) |>
    dplyr::slice(1)

  gene_short <- if (!is.na(cpg_info$UCSC_RefGene_Name) && cpg_info$UCSC_RefGene_Name != "") {
    stringr::str_extract(cpg_info$UCSC_RefGene_Name, "^[^;]+")
  } else {
    "intergenic/no_gene"
  }

  title <- paste(cpg_id, gene_short, sep = " / ")
  subtitle <- paste0(
    cpg_info$top_contrast,
    "; ",
    cpg_info$top_effect_direction,
    "; robust contrasts = ",
    cpg_info$robust_contrasts,
    "; warning = ",
    cpg_info$dataset_warning
  )

  cpg_plot_df <- plot_df |>
    dplyr::filter(cpg == cpg_id)

  cpg_effects <- top_effects |>
    dplyr::filter(cpg == cpg_id)

  p1 <- ggplot2::ggplot(
    cpg_plot_df,
    ggplot2::aes(x = DiseaseGroup, y = beta_value, color = Dataset)
  ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = DiseaseGroup),
      outlier.shape = NA,
      color = "grey35",
      fill = "grey96",
      linewidth = 0.25
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.16, height = 0, seed = 123),
      size = 0.9,
      alpha = 0.7
    ) +
    ggplot2::scale_color_manual(values = dataset_palette, drop = FALSE) +
    ggplot2::labs(x = NULL, y = "Beta value", color = "Dataset") +
    theme_pipeline(9)

  p2 <- ggplot2::ggplot(
    cpg_plot_df,
    ggplot2::aes(x = Dataset, y = beta_value, color = DiseaseGroup)
  ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = Dataset),
      outlier.shape = NA,
      color = "grey35",
      fill = "grey96",
      linewidth = 0.25
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.16, height = 0, seed = 123),
      size = 0.9,
      alpha = 0.7
    ) +
    ggplot2::scale_color_manual(values = disease_group_palette, drop = FALSE) +
    ggplot2::labs(x = NULL, y = "Beta value", color = "DiseaseGroup") +
    theme_pipeline(9)

  p3 <- ggplot2::ggplot(
    cpg_effects,
    ggplot2::aes(x = estimate, y = dataset, color = contrast)
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_low, xmax = ci_high),
      orientation = "y",
      width = 0.16,
      linewidth = 0.25,
      alpha = 0.85
    ) +
    ggplot2::geom_point(size = 1.3) +
    ggplot2::labs(
      x = "Dataset-level effect",
      y = NULL,
      color = "Contrast"
    ) +
    theme_pipeline(9)

  (p1 / p2 / p3) +
    patchwork::plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13),
        plot.subtitle = ggplot2::element_text(size = 10)
      )
    )
}

single_paths <- purrr::map_chr(top_cpgs, function(cpg_id) {
  p <- make_single_plot(cpg_id)
  path_out <- file.path(single_plot_dir, paste0(cpg_id, "_candidate_plot.png"))
  ggplot2::ggsave(path_out, p, width = 11, height = 12, dpi = 300)
  path_out
})

message("Candidate summary: ", summary_path)
message("Top CpGs: ", top_path)
message("Top CpG beta values: ", plot_data_path)
message("Combined DiseaseGroup plot: ", combined_group_path)
message("Combined dataset-effect plot: ", combined_effect_path)
message("Single CpG plots written to: ", single_plot_dir)

message("Top CpG candidates:")
print(
  top_candidates |>
    dplyr::select(
      cpg,
      UCSC_RefGene_Name,
      candidate_priority_score,
      robust_contrasts,
      tier1_contrasts,
      best_fdr,
      top_contrast,
      top_effect_direction,
      priority_label,
      dataset_warning
    )
)
