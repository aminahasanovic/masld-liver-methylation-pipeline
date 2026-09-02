# 25_kurokawa_transfer_analysis.R
# External cross-etiology transfer analysis on Kurokawa/GSE60753 liver DNAm.

if (!file.exists(Sys.getenv("CLASSIFIER_CONFIG", "config/classifier_config.yaml")) &&
    file.exists("classifier/config/classifier_config.yaml")) {
  Sys.setenv(CLASSIFIER_CONFIG = "classifier/config/classifier_config.yaml")
}

setup_path <- if (file.exists("scripts/00_setup_paths.R")) {
  "scripts/00_setup_paths.R"
} else {
  "classifier/scripts/00_setup_paths.R"
}

source(setup_path)

suppressPackageStartupMessages({
  library(caret)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
  library(readr)
  library(tibble)
  library(tidyr)
})

main_root <- normalizePath(file.path(config$project_root, ".."), mustWork = TRUE)
results_dir <- file.path(main_root, "results")
results_fig_dir <- file.path(results_dir, "figures")
results_tab_dir <- file.path(results_dir, "tables")
dir.create(results_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_tab_dir, recursive = TRUE, showWarnings = FALSE)

analysis_run_id <- Sys.getenv(
  "CLASSIFIER_KUROKAWA_RUN_ID",
  "disease_group_4_cpg_only_train_test"
)
model_name <- Sys.getenv("CLASSIFIER_KUROKAWA_MODEL", "elastic_net")

if (model_name != "elastic_net") {
  stop("This transfer analysis is currently implemented for the primary Elastic Net model.", call. = FALSE)
}

model_path <- file.path(out_paths$results_models, paste0(model_name, "_model_", analysis_run_id, ".rds"))
feature_path <- file.path(out_paths$results_features, paste0(model_name, "_features_", analysis_run_id, ".rds"))
input_path <- file.path(out_paths$data_processed, paste0("classifier_inputs_", analysis_run_id, ".rds"))

kuro_beta_path <- file.path(main_root, "data", "processed", "Kurokawa_GSE60753", "beta_60753.rds")
kuro_pheno_path <- file.path(main_root, "data", "processed", "Kurokawa_GSE60753", "pheno_60753_clean.rds")

fit <- readRDS(check_file(model_path))
feat <- readRDS(check_file(feature_path))
inp <- readRDS(check_file(input_path))
kuro_beta <- readRDS(check_file(kuro_beta_path))
kuro_pheno <- readRDS(check_file(kuro_pheno_path))

selected_cpgs <- as.character(feat$selected_cpgs)
class_levels <- as.character(fit$levels)

if (!all(c("Healthy", "Advanced_Fibrosis") %in% class_levels)) {
  stop("Expected Healthy and Advanced_Fibrosis probability columns in the model.", call. = FALSE)
}

if (is.null(rownames(kuro_beta)) || is.null(colnames(kuro_beta))) {
  stop("Kurokawa beta matrix must have CpG rownames and sample colnames.", call. = FALSE)
}

kuro_pheno <- kuro_pheno |>
  dplyr::filter(Sample_Name %in% colnames(kuro_beta)) |>
  dplyr::arrange(match(Sample_Name, colnames(kuro_beta)))

kuro_beta <- kuro_beta[, kuro_pheno$Sample_Name, drop = FALSE]
stopifnot(identical(colnames(kuro_beta), kuro_pheno$Sample_Name))

make_kuro_group <- function(group, disease_state) {
  dplyr::case_when(
    disease_state == "Normal" ~ "Normal/reference",
    group == "CirrEtOH" ~ "Alcohol cirrhosis",
    group == "CirrHCV" ~ "HCV cirrhosis",
    group == "CirrHBV" ~ "HBV cirrhosis",
    group %in% c("CC", "Cbil", "CG", "CI", "Hbil", "HM") ~ "Other cirrhosis",
    TRUE ~ NA_character_
  )
}

group_levels <- c(
  "Normal/reference",
  "Alcohol cirrhosis",
  "HCV cirrhosis",
  "HBV cirrhosis",
  "Other cirrhosis"
)

analysis_meta <- kuro_pheno |>
  dplyr::mutate(
    kurokawa_group = make_kuro_group(Group, DiseaseState),
    kurokawa_group = factor(kurokawa_group, levels = group_levels)
  ) |>
  dplyr::filter(!is.na(kurokawa_group)) |>
  dplyr::arrange(kurokawa_group, Sample_Name)

if (!any(analysis_meta$kurokawa_group == "Alcohol cirrhosis")) {
  stop("No Kurokawa alcohol cirrhosis samples available after filtering.", call. = FALSE)
}

# The external 450K dataset lacks some selected CpGs after QC filtering.
# Missing model predictors are set to the training median, making their
# contribution neutral relative to the training distribution.
x_train_all <- readRDS(check_file(inp$x_cpg_path))
if (!inherits(x_train_all, "Matrix")) {
  x_train_all <- Matrix::Matrix(x_train_all, sparse = TRUE)
}

train_ids <- inp$meta$sample_id[feat$train_idx]
train_ids <- intersect(train_ids, rownames(x_train_all))
selected_present_train <- intersect(selected_cpgs, colnames(x_train_all))

if (!setequal(selected_present_train, selected_cpgs)) {
  missing_train <- setdiff(selected_cpgs, selected_present_train)
  stop(
    "Selected model CpGs are missing from the training matrix: ",
    paste(head(missing_train, 10), collapse = ", "),
    call. = FALSE
  )
}

train_selected <- as.matrix(x_train_all[train_ids, selected_cpgs, drop = FALSE])
training_medians <- apply(train_selected, 2, stats::median, na.rm = TRUE)
training_medians[is.na(training_medians)] <- 0.5

sample_ids <- analysis_meta$Sample_Name
external_x <- matrix(
  rep(training_medians, each = length(sample_ids)),
  nrow = length(sample_ids),
  ncol = length(selected_cpgs),
  dimnames = list(sample_ids, selected_cpgs)
)

present_external_cpgs <- intersect(selected_cpgs, rownames(kuro_beta))
external_values <- t(kuro_beta[present_external_cpgs, sample_ids, drop = FALSE])
external_x[, present_external_cpgs] <- external_values

na_positions <- which(is.na(external_x), arr.ind = TRUE)
if (nrow(na_positions) > 0) {
  external_x[na_positions] <- training_medians[colnames(external_x)[na_positions[, "col"]]]
}

external_x <- Matrix::Matrix(external_x, sparse = TRUE)

pred_class <- predict(fit, newdata = external_x)
pred_prob <- predict(fit, newdata = external_x, type = "prob")

prob_tbl <- tibble::as_tibble(pred_prob)
for (nm in setdiff(class_levels, colnames(prob_tbl))) {
  prob_tbl[[nm]] <- NA_real_
}
prob_tbl <- prob_tbl[, class_levels, drop = FALSE]

predictions <- analysis_meta |>
  dplyr::transmute(
    sample_id = Sample_Name,
    sample_raw = Sample_raw,
    kurokawa_group = as.character(kurokawa_group),
    original_group = Group,
    disease_state = DiseaseState,
    etiology_group = Etiology,
    metadata_etiology = Etiology_age,
    age = Age,
    sex = Sex,
    predicted_class = as.character(pred_class)
  ) |>
  dplyr::bind_cols(prob_tbl) |>
  dplyr::mutate(
    prob_disease = 1 - .data[["Healthy"]],
    prob_advanced_fibrosis = .data[["Advanced_Fibrosis"]]
  )

auc_rank <- function(labels, scores) {
  labels <- as.logical(labels)
  ok <- !is.na(labels) & !is.na(scores)
  labels <- labels[ok]
  scores <- scores[ok]
  n_pos <- sum(labels)
  n_neg <- sum(!labels)
  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

normal_scores <- predictions |>
  dplyr::filter(kurokawa_group == "Normal/reference") |>
  dplyr::pull(prob_advanced_fibrosis)

group_auc <- lapply(setdiff(group_levels, "Normal/reference"), function(grp) {
  df <- predictions |>
    dplyr::filter(kurokawa_group %in% c("Normal/reference", grp))
  tibble::tibble(
    kurokawa_group = grp,
    auc_advanced_vs_normal = auc_rank(
      labels = df$kurokawa_group == grp,
      scores = df$prob_advanced_fibrosis
    ),
    wilcox_p_advanced_vs_normal = tryCatch(
      stats::wilcox.test(
        df$prob_advanced_fibrosis[df$kurokawa_group == grp],
        df$prob_advanced_fibrosis[df$kurokawa_group == "Normal/reference"],
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
  )
}) |>
  dplyr::bind_rows()

class_composition <- predictions |>
  dplyr::count(kurokawa_group, predicted_class, name = "n_predicted") |>
  tidyr::complete(
    kurokawa_group = group_levels,
    predicted_class = class_levels,
    fill = list(n_predicted = 0)
  ) |>
  dplyr::group_by(kurokawa_group) |>
  dplyr::mutate(
    n_group = sum(n_predicted),
    fraction_predicted = dplyr::if_else(n_group > 0, n_predicted / n_group, NA_real_)
  ) |>
  dplyr::ungroup()

prediction_summary <- predictions |>
  dplyr::group_by(kurokawa_group) |>
  dplyr::summarise(
    n_samples = dplyr::n(),
    pct_pred_healthy = mean(predicted_class == "Healthy") * 100,
    pct_pred_masl_mash = mean(predicted_class == "MASL_MASH") * 100,
    pct_pred_mild_fibrosis = mean(predicted_class == "Mild_Fibrosis") * 100,
    pct_pred_advanced_fibrosis = mean(predicted_class == "Advanced_Fibrosis") * 100,
    median_prob_advanced_fibrosis = stats::median(prob_advanced_fibrosis, na.rm = TRUE),
    q25_prob_advanced_fibrosis = stats::quantile(prob_advanced_fibrosis, 0.25, na.rm = TRUE),
    q75_prob_advanced_fibrosis = stats::quantile(prob_advanced_fibrosis, 0.75, na.rm = TRUE),
    median_prob_healthy = stats::median(.data[["Healthy"]], na.rm = TRUE),
    selected_cpgs_present = length(present_external_cpgs),
    selected_cpgs_total = length(selected_cpgs),
    selected_cpg_coverage = length(present_external_cpgs) / length(selected_cpgs),
    .groups = "drop"
  ) |>
  dplyr::left_join(group_auc, by = "kurokawa_group") |>
  dplyr::mutate(
    model = model_name,
    run_id = analysis_run_id,
    note = dplyr::case_when(
      kurokawa_group == "Normal/reference" ~ "Reference group; HCC and unclassified samples excluded.",
      TRUE ~ "Compared against Normal/reference with Advanced_Fibrosis probability."
    )
  ) |>
  dplyr::select(
    model, run_id, kurokawa_group, n_samples,
    pct_pred_healthy, pct_pred_masl_mash, pct_pred_mild_fibrosis,
    pct_pred_advanced_fibrosis, median_prob_advanced_fibrosis,
    q25_prob_advanced_fibrosis, q75_prob_advanced_fibrosis,
    median_prob_healthy, auc_advanced_vs_normal,
    wilcox_p_advanced_vs_normal, selected_cpgs_present,
    selected_cpgs_total, selected_cpg_coverage, note
  ) |>
  dplyr::arrange(match(kurokawa_group, group_levels))

predictions_path <- file.path(
  out_paths$results_metrics,
  paste0("kurokawa_transfer_predictions_", analysis_run_id, ".csv")
)
summary_path <- file.path(
  results_tab_dir,
  "T18_kurokawa_transfer_prediction_summary.csv"
)

readr::write_csv(predictions, predictions_path)
readr::write_csv(prediction_summary, summary_path)

class_palette <- c(
  Healthy = "#2E7D32",
  MASL_MASH = "#F2A93B",
  Mild_Fibrosis = "#D97A2B",
  Advanced_Fibrosis = "#B33A3A"
)

class_label <- c(
  Healthy = "Healthy",
  MASL_MASH = "MASL/MASH",
  Mild_Fibrosis = "Mild fibrosis",
  Advanced_Fibrosis = "Advanced fibrosis"
)

group_palette <- c(
  "Normal/reference" = "#4C8C6B",
  "Alcohol cirrhosis" = "#B85C38",
  "HCV cirrhosis" = "#5B6FA9",
  "HBV cirrhosis" = "#7A6A9E",
  "Other cirrhosis" = "#737373"
)

percent_label <- function(x) {
  paste0(round(100 * x), "%")
}

plot_theme <- theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8)
  )

prob_profile <- predictions |>
  dplyr::select(kurokawa_group, dplyr::all_of(class_levels)) |>
  tidyr::pivot_longer(cols = dplyr::all_of(class_levels), names_to = "class", values_to = "probability") |>
  dplyr::group_by(kurokawa_group, class) |>
  dplyr::summarise(
    median_probability = stats::median(probability, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    kurokawa_group = factor(kurokawa_group, levels = group_levels),
    class = factor(class, levels = rev(class_levels), labels = rev(unname(class_label[class_levels])))
  )

p_prob <- ggplot(prob_profile, aes(x = kurokawa_group, y = class, fill = median_probability)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = percent_label(median_probability)), size = 2.6) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1), labels = percent_label) +
  labs(
    title = "Median predicted probability",
    x = "Kurokawa group",
    y = "Predicted class",
    fill = "Median\nprobability"
  ) +
  plot_theme +
  theme(axis.text.y = element_text(size = 8.5))

p_comp <- class_composition |>
  dplyr::mutate(
    kurokawa_group = factor(kurokawa_group, levels = group_levels),
    predicted_class = factor(predicted_class, levels = class_levels, labels = unname(class_label[class_levels])),
    segment_label = dplyr::if_else(fraction_predicted > 0, percent_label(fraction_predicted), "")
  ) |>
  ggplot(aes(x = kurokawa_group, y = fraction_predicted, fill = predicted_class)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = segment_label),
    position = position_stack(vjust = 0.5),
    size = 2.4,
    color = "#1F1F1F"
  ) +
  scale_fill_manual(values = setNames(class_palette[class_levels], unname(class_label[class_levels])), drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    title = "Predicted class composition",
    x = "Kurokawa group",
    y = "Fraction of samples"
  ) +
  plot_theme +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

combined_plot <- p_prob + p_comp + plot_layout(widths = c(1.1, 1))

figure_stem <- "F25_kurokawa_transfer_classifier_predictions"
ggsave(
  filename = file.path(results_fig_dir, paste0(figure_stem, ".pdf")),
  plot = combined_plot,
  width = 10.2,
  height = 4.9,
  device = cairo_pdf
)
ggsave(
  filename = file.path(results_fig_dir, paste0(figure_stem, ".png")),
  plot = combined_plot,
  width = 10.2,
  height = 4.9,
  dpi = 600
)

update_manifest <- function(path, new_rows, key_col) {
  existing <- if (file.exists(path)) {
    readr::read_csv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  } else {
    tibble::tibble()
  }

  out <- dplyr::bind_rows(existing, new_rows) |>
    dplyr::distinct(.data[[key_col]], .keep_all = TRUE)

  readr::write_csv(out, path)
}

update_manifest(
  file.path(results_tab_dir, "figure_manifest.csv"),
  tibble::tibble(
    figure = "F25",
    file_stem = figure_stem,
    main_message = "The MASLD classifier assigns high Advanced_Fibrosis probability to Kurokawa cirrhosis samples, especially alcohol-related cirrhosis."
  ),
  key_col = "figure"
)

update_manifest(
  file.path(results_tab_dir, "table_manifest.csv"),
  tibble::tibble(
    table = "T18",
    file = "T18_kurokawa_transfer_prediction_summary.csv",
    main_message = "External Kurokawa transfer summary for the primary Elastic Net MASLD classifier."
  ),
  key_col = "table"
)

message("Kurokawa transfer analysis complete.")
message("Model: ", model_name, " / ", analysis_run_id)
message("Analysis samples: ", nrow(predictions))
message("Selected CpG coverage in Kurokawa: ", length(present_external_cpgs), "/", length(selected_cpgs))
message("Saved detailed predictions: ", predictions_path)
message("Saved result table: ", summary_path)
message("Saved result figure: ", file.path(results_fig_dir, paste0(figure_stem, ".pdf")))
print(prediction_summary, n = Inf)
