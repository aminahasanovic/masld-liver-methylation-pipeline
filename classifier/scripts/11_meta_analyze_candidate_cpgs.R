# 11_meta_analyze_candidate_cpgs.R
# Cross-study robustness analysis for prioritized Elastic Net CpG candidates.

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

candidate_path <- file.path(out_paths$results_features, "elastic_net_cpg_candidate_priority.csv")

if (!file.exists(candidate_path)) {
  stop("Missing candidate priority table. Run scripts/10_plot_candidate_cpgs.R first.")
}

beta_path <- check_file(config$inputs$beta_matrix_path)
meta_path <- check_file(config$inputs$metadata_path)
sample_id_col <- config$inputs$sample_id_col %||% "Sample_Name"
min_group_n <- as.integer(Sys.getenv("CLASSIFIER_META_MIN_GROUP_N", "3"))

if (is.na(min_group_n) || min_group_n < 2) {
  stop("CLASSIFIER_META_MIN_GROUP_N must be an integer >= 2.")
}

message("Running cross-study candidate CpG meta-analysis.")
message("Minimum per-group n for binary contrasts: ", min_group_n)

candidate_tbl <- readr::read_csv(candidate_path, show_col_types = FALSE)
candidate_cpgs <- unique(candidate_tbl$cpg)

beta <- readRDS(beta_path)
meta <- readRDS(meta_path)

candidate_cpgs <- intersect(candidate_cpgs, rownames(beta))

if (length(candidate_cpgs) == 0) {
  stop("No candidate CpGs are present in the configured beta matrix.")
}

if (!sample_id_col %in% colnames(meta)) {
  stop("Sample ID column missing from metadata: ", sample_id_col)
}

common_samples <- intersect(colnames(beta), as.character(meta[[sample_id_col]]))

meta <- meta |>
  dplyr::filter(.data[[sample_id_col]] %in% common_samples) |>
  dplyr::mutate(
    sample_id = as.character(.data[[sample_id_col]]),
    Age_num = suppressWarnings(as.numeric(as.character(Age))),
    Sex_clean = dplyr::if_else(is.na(Sex) | Sex == "", NA_character_, as.character(Sex)),
    disease_binary = dplyr::case_when(
      DiseaseGroup %in% c("Healthy", "Healthy_Obese") ~ 0,
      DiseaseGroup %in% c("MASL", "MASH", "Mild_Fibrosis", "Advanced_Fibrosis") ~ 1,
      TRUE ~ NA_real_
    ),
    fibrosis_binary = dplyr::case_when(
      DiseaseGroup %in% c("Healthy", "Healthy_Obese", "MASL", "MASH") ~ 0,
      DiseaseGroup %in% c("Mild_Fibrosis", "Advanced_Fibrosis") ~ 1,
      TRUE ~ NA_real_
    ),
    advanced_binary = dplyr::case_when(
      DiseaseGroup %in% c("Healthy", "Healthy_Obese", "MASL", "MASH", "Mild_Fibrosis") ~ 0,
      DiseaseGroup == "Advanced_Fibrosis" ~ 1,
      TRUE ~ NA_real_
    ),
    progression_score = dplyr::case_when(
      DiseaseGroup %in% c("Healthy", "Healthy_Obese") ~ 0,
      DiseaseGroup %in% c("MASL", "MASH") ~ 1,
      DiseaseGroup == "Mild_Fibrosis" ~ 2,
      DiseaseGroup == "Advanced_Fibrosis" ~ 3,
      TRUE ~ NA_real_
    )
  ) |>
  dplyr::arrange(match(sample_id, common_samples))

beta_sub <- beta[candidate_cpgs, meta$sample_id, drop = FALSE]

beta_df <- as.data.frame(t(beta_sub)) |>
  tibble::rownames_to_column("sample_id") |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(candidate_cpgs),
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
        Age_num,
        Sex_clean,
        disease_binary,
        fibrosis_binary,
        advanced_binary,
        progression_score
      ),
    by = "sample_id"
  )

make_formula <- function(df, predictor) {
  covariates <- character()

  age_ok <- sum(!is.na(df$Age_num)) >= 0.8 * nrow(df) &&
    stats::sd(df$Age_num, na.rm = TRUE) > 0 &&
    nrow(df) >= 20

  sex_ok <- sum(!is.na(df$Sex_clean)) >= 0.8 * nrow(df) &&
    dplyr::n_distinct(df$Sex_clean[!is.na(df$Sex_clean)]) >= 2 &&
    nrow(df) >= 20

  if (age_ok) {
    covariates <- c(covariates, "Age_num")
  }

  if (sex_ok) {
    covariates <- c(covariates, "Sex_clean")
  }

  rhs <- paste(c(predictor, covariates), collapse = " + ")
  stats::as.formula(paste("beta_value ~", rhs))
}

fit_one_dataset <- function(df, contrast_name, predictor, binary = TRUE) {
  df <- df |>
    dplyr::filter(!is.na(beta_value), !is.na(.data[[predictor]]))

  if (nrow(df) < 8) {
    return(tibble::tibble())
  }

  if (binary) {
    counts <- table(df[[predictor]])

    if (!all(c("0", "1") %in% names(counts))) {
      return(tibble::tibble())
    }

    if (any(counts[c("0", "1")] < min_group_n)) {
      return(tibble::tibble())
    }
  } else {
    score_counts <- table(df[[predictor]])

    if (length(score_counts) < 2) {
      return(tibble::tibble())
    }

    if (min(score_counts) < min_group_n) {
      return(tibble::tibble())
    }
  }

  fit <- tryCatch(
    stats::lm(make_formula(df, predictor), data = df),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble::tibble())
  }

  coef_tbl <- summary(fit)$coefficients

  if (!predictor %in% rownames(coef_tbl)) {
    return(tibble::tibble())
  }

  disease_counts <- table(df$DiseaseGroup)
  predictor_counts <- table(df[[predictor]])

  tibble::tibble(
    contrast = contrast_name,
    cpg = df$cpg[[1]],
    dataset = df$Dataset[[1]],
    array = paste(sort(unique(df$Array)), collapse = ";"),
    n = nrow(df),
    n_reference = if (binary) unname(predictor_counts["0"]) else NA_integer_,
    n_case = if (binary) unname(predictor_counts["1"]) else NA_integer_,
    predictor_min = min(df[[predictor]], na.rm = TRUE),
    predictor_max = max(df[[predictor]], na.rm = TRUE),
    estimate = unname(coef_tbl[predictor, "Estimate"]),
    std_error = unname(coef_tbl[predictor, "Std. Error"]),
    statistic = unname(coef_tbl[predictor, "t value"]),
    p_value = unname(coef_tbl[predictor, "Pr(>|t|)"]),
    covariates_used = paste(attr(stats::terms(fit), "term.labels")[-1], collapse = ";"),
    disease_group_counts = paste(
      paste(names(disease_counts), as.integer(disease_counts), sep = "="),
      collapse = ";"
    )
  )
}

contrast_specs <- tibble::tribble(
  ~contrast, ~predictor, ~binary, ~description,
  "disease_vs_healthy", "disease_binary", TRUE,
  "Disease groups (MASL/MASH/fibrosis) versus Healthy/Healthy_Obese.",
  "fibrosis_vs_nonfibrosis", "fibrosis_binary", TRUE,
  "Early/Advanced fibrosis versus Healthy/Healthy_Obese/MASL/MASH.",
  "advanced_vs_nonadvanced", "advanced_binary", TRUE,
  "Advanced_Fibrosis versus all non-advanced groups.",
  "progression_score", "progression_score", FALSE,
  "Linear trend: Healthy/Healthy_Obese=0, MASL/MASH=1, Mild_Fibrosis=2, Advanced_Fibrosis=3."
)

candidate_dataset_groups <- beta_df |>
  dplyr::group_by(cpg, Dataset) |>
  dplyr::group_split(.keep = TRUE)

dataset_effects <- purrr::map_dfr(seq_len(nrow(contrast_specs)), function(i) {
  spec <- contrast_specs[i, ]

  purrr::map_dfr(candidate_dataset_groups, function(.x) {
      fit_one_dataset(
        df = .x,
        contrast_name = spec$contrast,
        predictor = spec$predictor,
        binary = spec$binary
      )
    })
})

meta_analyze <- function(df) {
  df <- df |>
    dplyr::filter(
      is.finite(estimate),
      is.finite(std_error),
      std_error > 0
    )

  k <- nrow(df)

  if (k < 2) {
    return(tibble::tibble())
  }

  yi <- df$estimate
  vi <- df$std_error^2
  wi <- 1 / vi

  fixed_effect <- sum(wi * yi) / sum(wi)
  q_stat <- sum(wi * (yi - fixed_effect)^2)
  q_df <- k - 1
  c_val <- sum(wi) - (sum(wi^2) / sum(wi))
  tau2 <- if (c_val > 0) max(0, (q_stat - q_df) / c_val) else 0
  wi_random <- 1 / (vi + tau2)

  random_effect <- sum(wi_random * yi) / sum(wi_random)
  random_se <- sqrt(1 / sum(wi_random))
  z_value <- random_effect / random_se
  p_value <- 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)
  ci_low <- random_effect - stats::qnorm(0.975) * random_se
  ci_high <- random_effect + stats::qnorm(0.975) * random_se
  i2 <- if (q_stat > 0) max(0, (q_stat - q_df) / q_stat) * 100 else 0
  q_p_value <- stats::pchisq(q_stat, df = q_df, lower.tail = FALSE)

  n_positive <- sum(yi > 0)
  n_negative <- sum(yi < 0)
  direction_consistency <- max(n_positive, n_negative) / k

  tibble::tibble(
    n_datasets = k,
    datasets = paste(sort(unique(df$dataset)), collapse = ";"),
    total_n = sum(df$n),
    random_effect = random_effect,
    random_se = random_se,
    ci_low = ci_low,
    ci_high = ci_high,
    z_value = z_value,
    p_value = p_value,
    tau2 = tau2,
    i2 = i2,
    q_p_value = q_p_value,
    n_positive = n_positive,
    n_negative = n_negative,
    direction_consistency = direction_consistency
  )
}

meta_results <- dataset_effects |>
  dplyr::group_by(contrast, cpg) |>
  dplyr::group_modify(~ meta_analyze(.x)) |>
  dplyr::ungroup() |>
  dplyr::group_by(contrast) |>
  dplyr::mutate(fdr = stats::p.adjust(p_value, method = "BH")) |>
  dplyr::ungroup() |>
  dplyr::left_join(candidate_tbl, by = "cpg") |>
  dplyr::mutate(
    evidence_tier = dplyr::case_when(
      n_datasets >= 3 &
        fdr < 0.05 &
        direction_consistency == 1 &
        !likely_dataset_sensitive ~ "tier1_consistent_low_dataset_signal",
      n_datasets >= 3 &
        fdr < 0.10 &
        direction_consistency >= 0.75 &
        disease_r2_after_dataset >= 0.02 ~ "tier2_consistent_exploratory",
      n_datasets >= 2 &
        fdr < 0.10 &
        direction_consistency == 1 ~ "tier3_two_dataset_signal",
      TRUE ~ "not_robust"
    )
  ) |>
  dplyr::arrange(
    contrast,
    dplyr::desc(evidence_tier != "not_robust"),
    fdr,
    dplyr::desc(abs(random_effect)),
    cpg
  )

shortlist <- meta_results |>
  dplyr::filter(evidence_tier != "not_robust") |>
  dplyr::arrange(
    contrast,
    evidence_tier,
    fdr,
    dplyr::desc(abs(random_effect)),
    cpg
  )

strict_shortlist <- meta_results |>
  dplyr::filter(
    priority_label %in% c("stage_and_binary_all5", "stage_and_binary"),
    n_datasets >= 3,
    direction_consistency == 1,
    fdr < 0.05,
    i2 < 50,
    !likely_dataset_sensitive
  ) |>
  dplyr::mutate(
    strict_reason = paste(
      "stage/binary Elastic Net evidence",
      ">=3 datasets",
      "same effect direction",
      "FDR < 0.05",
      "I2 < 50",
      "dataset R2 <= disease-group R2",
      sep = "; "
    )
  ) |>
  dplyr::arrange(
    contrast,
    fdr,
    dplyr::desc(abs(random_effect)),
    cpg
  )

strict_unique <- strict_shortlist |>
  dplyr::group_by(cpg) |>
  dplyr::summarise(
    n_strict_contrasts = dplyr::n_distinct(contrast),
    strict_contrasts = paste(sort(unique(contrast)), collapse = ";"),
    best_fdr = min(fdr, na.rm = TRUE),
    max_abs_random_effect = max(abs(random_effect), na.rm = TRUE),
    directions = paste(
      paste(contrast, ifelse(random_effect > 0, "positive", "negative"), sep = "="),
      collapse = ";"
    ),
    datasets = paste(sort(unique(unlist(strsplit(datasets, ";", fixed = TRUE)))), collapse = ";"),
    UCSC_RefGene_Name = dplyr::first(UCSC_RefGene_Name),
    chr = dplyr::first(chr),
    pos = dplyr::first(pos),
    Relation_to_Island = dplyr::first(Relation_to_Island),
    priority_label = dplyr::first(priority_label),
    evidence_outcomes = dplyr::first(evidence_outcomes),
    best_nonzero_folds = dplyr::first(best_nonzero_folds),
    likely_dataset_sensitive = dplyr::first(likely_dataset_sensitive),
    r2_disease_group = dplyr::first(r2_disease_group),
    r2_dataset = dplyr::first(r2_dataset),
    disease_r2_after_dataset = dplyr::first(disease_r2_after_dataset),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(n_strict_contrasts),
    best_fdr,
    dplyr::desc(max_abs_random_effect),
    cpg
  )

dataset_effect_path <- file.path(
  out_paths$results_features,
  "elastic_net_candidate_cpg_meta_dataset_effects.csv"
)
meta_result_path <- file.path(
  out_paths$results_features,
  "elastic_net_candidate_cpg_meta_results.csv"
)
shortlist_path <- file.path(
  out_paths$results_features,
  "elastic_net_candidate_cpg_meta_shortlist.csv"
)
strict_shortlist_path <- file.path(
  out_paths$results_features,
  "elastic_net_candidate_cpg_meta_strict_shortlist.csv"
)
strict_unique_path <- file.path(
  out_paths$results_features,
  "elastic_net_candidate_cpg_meta_strict_unique.csv"
)

readr::write_csv(dataset_effects, dataset_effect_path)
readr::write_csv(meta_results, meta_result_path)
readr::write_csv(shortlist, shortlist_path)
readr::write_csv(strict_shortlist, strict_shortlist_path)
readr::write_csv(strict_unique, strict_unique_path)

plot_df <- meta_results |>
  dplyr::mutate(
    neg_log10_fdr = -log10(pmax(fdr, .Machine$double.xmin)),
    is_shortlisted = evidence_tier != "not_robust",
    label = dplyr::if_else(
      is_shortlisted & dplyr::row_number() <= 200,
      dplyr::if_else(
        !is.na(UCSC_RefGene_Name) & UCSC_RefGene_Name != "",
        paste0(cpg, " / ", stringr::str_extract(UCSC_RefGene_Name, "^[^;]+")),
        cpg
      ),
      NA_character_
    )
  )

p_volcano <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = random_effect,
    y = neg_log10_fdr,
    color = evidence_tier,
    shape = direction_consistency == 1
  )
) +
  ggplot2::geom_hline(yintercept = -log10(0.10), color = "grey55", linetype = "dashed") +
  ggplot2::geom_vline(xintercept = 0, color = "grey65") +
  ggplot2::geom_point(alpha = 0.8, size = 1.8) +
  ggplot2::facet_wrap(~ contrast, scales = "free") +
  ggplot2::scale_color_manual(
    values = c(
      tier1_consistent_low_dataset_signal = "#1B9E77",
      tier2_consistent_exploratory = "#7570B3",
      tier3_two_dataset_signal = "#E6AB02",
      not_robust = "grey70"
    ),
    drop = FALSE
  ) +
  ggplot2::labs(
    x = "Random-effects meta-analysis estimate",
    y = "-log10(FDR)",
    color = "Evidence tier",
    shape = "Same direction\nall datasets",
    title = "Cross-study robustness of Elastic Net CpG candidates",
    subtitle = "Effects are estimated within each dataset first, then combined across datasets."
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_volcano <- p_volcano +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = 2.2,
      max.overlaps = 35,
      na.rm = TRUE,
      show.legend = FALSE
    )
}

heatmap_source <- if (nrow(strict_shortlist) > 0) strict_shortlist else shortlist

top_heatmap_cpgs <- heatmap_source |>
  dplyr::arrange(fdr, dplyr::desc(abs(random_effect))) |>
  dplyr::slice_head(n = 40) |>
  dplyr::pull(cpg) |>
  unique()

heatmap_df <- dataset_effects |>
  dplyr::filter(cpg %in% top_heatmap_cpgs) |>
  dplyr::left_join(
    candidate_tbl |>
      dplyr::select(cpg, UCSC_RefGene_Name),
    by = "cpg"
  ) |>
  dplyr::mutate(
    cpg_label = dplyr::if_else(
      !is.na(UCSC_RefGene_Name) & UCSC_RefGene_Name != "",
      paste0(cpg, " / ", stringr::str_extract(UCSC_RefGene_Name, "^[^;]+")),
      cpg
    )
  )

p_heatmap <- ggplot2::ggplot(
  heatmap_df,
  ggplot2::aes(x = dataset, y = cpg_label, fill = estimate)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.2) +
  ggplot2::facet_wrap(~ contrast, scales = "free_y") +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "Effect",
    title = "Dataset-level effects for shortlisted CpGs",
    subtitle = "Consistent color direction across datasets is stronger evidence than model selection alone."
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

volcano_path <- file.path(out_paths$results_plots, "elastic_net_candidate_cpg_meta_volcano.png")
heatmap_path <- file.path(out_paths$results_plots, "elastic_net_candidate_cpg_meta_effect_heatmap.png")

ggplot2::ggsave(volcano_path, p_volcano, width = 12, height = 8, dpi = 300)
ggplot2::ggsave(heatmap_path, p_heatmap, width = 13, height = 10, dpi = 300)

message("Dataset-level effects: ", dataset_effect_path)
message("Meta-analysis results: ", meta_result_path)
message("Robust shortlist: ", shortlist_path)
message("Strict biomarker-oriented shortlist: ", strict_shortlist_path)
message("Strict unique CpG shortlist: ", strict_unique_path)
message("Meta volcano plot: ", volcano_path)
message("Dataset effect heatmap: ", heatmap_path)

message("Meta-analysis rows by contrast:")
print(table(meta_results$contrast))

message("Shortlist by contrast and tier:")
print(table(shortlist$contrast, shortlist$evidence_tier))

message("Strict shortlist by contrast:")
print(table(strict_shortlist$contrast))

message("Strict unique CpGs by number of robust contrasts:")
print(table(strict_unique$n_strict_contrasts))

message("Top shortlist entries:")
print(
  strict_unique |>
    dplyr::select(
      cpg,
      n_strict_contrasts,
      strict_contrasts,
      best_fdr,
      max_abs_random_effect,
      UCSC_RefGene_Name,
      priority_label,
      likely_dataset_sensitive
    ) |>
    utils::head(30)
)
