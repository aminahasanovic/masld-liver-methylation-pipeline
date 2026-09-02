# 13_add_biological_context_top_cpgs.R
# Join finalized manual biological curation to the fixed statistical CpG ranking.

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

candidate_summary_path <- file.path(
  out_paths$results_features,
  "elastic_net_cpg_candidate_summary.csv"
)
curation_path <- file.path(
  out_paths$results_features,
  "elastic_net_top30_biological_curation_final.csv"
)
dataset_effects_path <- file.path(
  out_paths$results_features,
  "elastic_net_candidate_cpg_meta_dataset_effects.csv"
)

review_path <- file.path(
  out_paths$results_features,
  "elastic_net_cpg_biological_review.csv"
)
interpretation_top15_path <- file.path(
  out_paths$results_features,
  "elastic_net_cpg_interpretation_top15.csv"
)
plot_data_path <- file.path(
  out_paths$results_features,
  "elastic_net_cpg_interpretation_top15_beta_values.csv"
)

interpretation_group_plot_path <- file.path(
  out_paths$results_plots,
  "elastic_net_cpg_interpretation_top15_beta_by_disease_group.png"
)
interpretation_effect_plot_path <- file.path(
  out_paths$results_plots,
  "elastic_net_cpg_interpretation_top15_dataset_effects.png"
)

# The curation table is a manual research input, not a pipeline product. The
# version used in the thesis is tracked in the repository, so the chain also
# completes on a fresh clone; a file placed in the results directory takes
# precedence and is used for a repeated review round.
tracked_curation_path <- file.path(classifier_root, "curation", basename(curation_path))
if (!file.exists(curation_path) && file.exists(tracked_curation_path)) {
  message("Using the curation table tracked in the repository: ", tracked_curation_path)
  curation_path <- tracked_curation_path
}

required_files <- c(candidate_summary_path, curation_path, dataset_effects_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input(s): ", paste(missing_files, collapse = ", "))
}

parse_bool <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  x <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("true", "t", "yes", "y", "1") ~ TRUE,
    x %in% c("false", "f", "no", "n", "0") ~ FALSE,
    TRUE ~ NA
  )
}

empty_to_na <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}

fmt_context <- function(original_annotation, current_gene_symbol, nearest_gene) {
  has_direct_gene <- !is.na(original_annotation) & nzchar(original_annotation)
  has_current_gene <- !is.na(current_gene_symbol) & nzchar(current_gene_symbol)
  has_nearest_gene <- !is.na(nearest_gene) & nzchar(nearest_gene)

  dplyr::case_when(
    has_direct_gene & has_current_gene ~ paste0(original_annotation, " / current: ", current_gene_symbol),
    has_direct_gene ~ original_annotation,
    has_nearest_gene ~ paste0("intergenic/no-gene; near ", nearest_gene),
    TRUE ~ "intergenic/no-gene"
  )
}

compare_numeric <- function(x, y, tolerance = 1e-12) {
  same_na <- is.na(x) & is.na(y)
  close <- !is.na(x) & !is.na(y) & abs(as.numeric(x) - as.numeric(y)) <= tolerance
  all(same_na | close)
}

compare_character <- function(x, y) {
  x <- empty_to_na(x)
  y <- empty_to_na(y)
  all((is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y))
}

message("Joining finalized manual biological curation to fixed statistical candidate ranking.")
message("Candidate summary: ", candidate_summary_path)
message("Finalized curation: ", curation_path)

candidate_summary <- readr::read_csv(candidate_summary_path, show_col_types = FALSE) |>
  dplyr::mutate(statistical_rank = dplyr::row_number())
curation <- readr::read_csv(curation_path, show_col_types = FALSE) |>
  dplyr::mutate(include_in_interpretation_panel = parse_bool(.data$include_in_interpretation_panel))
dataset_effects <- readr::read_csv(dataset_effects_path, show_col_types = FALSE)

required_curation_columns <- c(
  "statistical_rank",
  "cpg",
  "final_interpretation_priority",
  "include_in_interpretation_panel",
  "biological_evidence_class",
  "evidence_directness",
  "human_masld_evidence",
  "human_liver_fibrosis_evidence",
  "experimental_masld_evidence",
  "experimental_liver_fibrosis_evidence",
  "hepatic_metabolism_evidence",
  "literature_note",
  "evidence_source_1",
  "evidence_source_2",
  "evidence_source_3",
  "selection_reason"
)
missing_curation_columns <- setdiff(required_curation_columns, colnames(curation))
if (length(missing_curation_columns) > 0) {
  stop("Finalized curation file is missing column(s): ", paste(missing_curation_columns, collapse = ", "))
}

if (nrow(curation) != 30 || dplyr::n_distinct(curation$cpg) != 30) {
  stop("Finalized curation file must contain exactly 30 rows and 30 unique CpGs.")
}
if (!identical(as.integer(sort(curation$statistical_rank)), 1:30) ||
    anyDuplicated(curation$statistical_rank) > 0) {
  stop("Finalized curation file must contain statistical_rank exactly 1 through 30 without duplicates.")
}
if (sum(curation$include_in_interpretation_panel, na.rm = TRUE) != 15) {
  stop("Finalized curation file must contain exactly 15 TRUE include_in_interpretation_panel rows.")
}

selected_priorities <- curation |>
  dplyr::filter(.data$include_in_interpretation_panel) |>
  dplyr::pull(.data$final_interpretation_priority)

if (any(is.na(selected_priorities)) ||
    dplyr::n_distinct(selected_priorities) != 15 ||
    !identical(as.integer(sort(selected_priorities)), 1:15)) {
  stop("Selected interpretation candidates must have unique priorities exactly 1 through 15.")
}

top30_statistical <- candidate_summary |>
  dplyr::slice_head(n = 30) |>
  dplyr::select(
    statistical_rank,
    cpg,
    candidate_priority_score,
    robust_contrasts,
    tier1_contrasts,
    tier2_contrasts,
    tier3_contrasts,
    same_direction_contrasts,
    max_n_datasets,
    all_datasets_seen,
    best_fdr,
    best_neg_log10_fdr,
    max_abs_meta_effect,
    mean_i2,
    min_i2,
    strongest_contrast = top_contrast,
    top_random_effect,
    top_effect_direction,
    robust_contrast_summary,
    evidence_outcomes,
    evidence_runs,
    priority_label,
    best_nonzero_folds,
    all5_any_run,
    max_mean_abs_coefficient,
    likely_dataset_sensitive,
    r2_disease_group,
    r2_dataset,
    disease_r2_after_dataset,
    model_priority_score,
    dataset_warning,
    has_snp_annotation
  )

rank_check <- curation |>
  dplyr::select(cpg, curation_statistical_rank = statistical_rank) |>
  dplyr::left_join(
    candidate_summary |> dplyr::select(cpg, script12_statistical_rank = statistical_rank),
    by = "cpg"
  )

if (any(is.na(rank_check$script12_statistical_rank)) ||
    any(as.integer(rank_check$curation_statistical_rank) != rank_check$script12_statistical_rank)) {
  bad_rows <- rank_check |>
    dplyr::filter(
      is.na(.data$script12_statistical_rank) |
        as.integer(.data$curation_statistical_rank) != .data$script12_statistical_rank
    )
  print(bad_rows)
  stop("Finalized curation CpG/rank values do not match script-12 candidate summary.")
}

stat_check <- curation |>
  dplyr::inner_join(top30_statistical, by = "cpg", suffix = c("_curation", "_script12"))

numeric_stat_fields <- c(
  "candidate_priority_score",
  "robust_contrasts",
  "tier1_contrasts",
  "same_direction_contrasts",
  "max_n_datasets",
  "best_fdr",
  "best_neg_log10_fdr",
  "mean_i2",
  "max_abs_meta_effect"
)
character_stat_fields <- c("top_effect_direction", "strongest_contrast", "priority_label")

numeric_mismatch <- numeric_stat_fields[vapply(numeric_stat_fields, function(field) {
  !compare_numeric(stat_check[[paste0(field, "_curation")]], stat_check[[paste0(field, "_script12")]])
}, logical(1))]
character_mismatch <- character_stat_fields[vapply(character_stat_fields, function(field) {
  !compare_character(stat_check[[paste0(field, "_curation")]], stat_check[[paste0(field, "_script12")]])
}, logical(1))]

if (length(numeric_mismatch) > 0 || length(character_mismatch) > 0) {
  stop(
    "Finalized curation file differs from script-12 statistical values in: ",
    paste(c(numeric_mismatch, character_mismatch), collapse = ", ")
  )
}

curation_layer <- curation |>
  dplyr::transmute(
    cpg,
    evidence_tier,
    strict_subset_membership,
    direction_consistency,
    original_annotation,
    UCSC_RefGene_Name,
    current_gene_symbol,
    current_gene_aliases,
    gene_symbol_mapping_source,
    chromosome,
    genomic_position_hg19,
    MAPINFO,
    strand,
    UCSC_RefGene_Group,
    gene_region_annotation,
    UCSC_RefGene_Accession,
    Relation_to_UCSC_CpG_Island,
    CpG_island_name,
    Regulatory_Feature_Name,
    Regulatory_Feature_Group,
    Enhancer,
    DHS,
    DMR,
    HMM_Island,
    Phantom,
    nearby_regulatory_annotation,
    annotation_source,
    in_450k_annotation,
    in_epic_annotation,
    Probe_rs,
    Probe_maf,
    CpG_rs,
    CpG_maf,
    SBE_rs,
    SBE_maf,
    nearest_gene,
    nearest_gene_entrez,
    nearest_tss_transcript_id,
    nearest_tss_chr,
    nearest_tss_position,
    nearest_tss_strand,
    nearest_tss_distance,
    nearest_tss_orientation,
    biological_evidence_class,
    human_masld_evidence,
    human_liver_fibrosis_evidence,
    experimental_masld_evidence,
    experimental_liver_fibrosis_evidence,
    hepatic_metabolism_evidence,
    evidence_directness,
    evidence_source_1,
    evidence_source_2,
    evidence_source_3,
    literature_note,
    final_interpretation_priority,
    include_in_interpretation_panel,
    selection_reason
  )

reviewed <- top30_statistical |>
  dplyr::left_join(curation_layer, by = "cpg") |>
  dplyr::mutate(
    gene_or_genomic_context = fmt_context(
      .data$original_annotation,
      .data$current_gene_symbol,
      .data$nearest_gene
    ),
    interpretation_label = dplyr::if_else(
      !is.na(.data$final_interpretation_priority),
      paste0("Interpretation priority ", .data$final_interpretation_priority),
      NA_character_
    )
  ) |>
  dplyr::arrange(.data$statistical_rank)

if (nrow(reviewed) != 30 || dplyr::n_distinct(reviewed$cpg) != 30) {
  stop("Joined biological review must contain exactly 30 rows and 30 unique CpGs.")
}

interpretation_top15 <- reviewed |>
  dplyr::filter(.data$include_in_interpretation_panel) |>
  dplyr::arrange(.data$final_interpretation_priority, .data$statistical_rank, .data$cpg)

if (nrow(interpretation_top15) != 15 ||
    dplyr::n_distinct(interpretation_top15$cpg) != 15) {
  stop("Interpretation panel must contain exactly 15 unique CpGs.")
}

interpretation_top15_output <- interpretation_top15 |>
  dplyr::select(
    final_interpretation_priority,
    statistical_rank,
    cpg,
    original_annotation,
    current_gene_symbol,
    nearest_gene,
    nearest_tss_distance,
    gene_or_genomic_context,
    chromosome,
    genomic_position_hg19,
    Relation_to_UCSC_CpG_Island,
    gene_region_annotation,
    Regulatory_Feature_Group,
    Enhancer,
    DHS,
    DMR,
    evidence_tier,
    strict_subset_membership,
    strongest_contrast,
    best_fdr,
    mean_i2,
    max_n_datasets,
    direction_consistency,
    top_effect_direction,
    likely_dataset_sensitive,
    biological_evidence_class,
    evidence_directness,
    human_masld_evidence,
    human_liver_fibrosis_evidence,
    experimental_masld_evidence,
    experimental_liver_fibrosis_evidence,
    hepatic_metabolism_evidence,
    literature_note,
    evidence_source_1,
    evidence_source_2,
    evidence_source_3,
    selection_reason,
    dplyr::everything()
  )

readr::write_csv(reviewed, review_path)
readr::write_csv(interpretation_top15_output, interpretation_top15_path)

beta_path <- check_file(config$inputs$beta_matrix_path)
meta_path <- check_file(config$inputs$metadata_path)
sample_id_col <- config$inputs$sample_id_col %||% "Sample_Name"

beta <- readRDS(beta_path)
meta <- readRDS(meta_path)

if (!sample_id_col %in% colnames(meta)) {
  stop("Sample ID column missing from metadata: ", sample_id_col)
}

common_samples <- intersect(colnames(beta), as.character(meta[[sample_id_col]]))
top_cpgs <- interpretation_top15$cpg

missing_beta_cpgs <- setdiff(top_cpgs, rownames(beta))
if (length(missing_beta_cpgs) > 0) {
  stop("Interpretation panel CpGs missing from beta matrix: ", paste(missing_beta_cpgs, collapse = ", "))
}

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
    Progression3 = factor(Progression3)
  ) |>
  dplyr::arrange(match(sample_id, common_samples))

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
      dplyr::select(sample_id, Dataset, DiseaseGroup, Progression3),
    by = "sample_id"
  ) |>
  dplyr::left_join(
    interpretation_top15 |>
      dplyr::select(
        cpg,
        final_interpretation_priority,
        original_annotation,
        current_gene_symbol,
        nearest_gene,
        gene_or_genomic_context,
        strongest_contrast,
        top_effect_direction
      ),
    by = "cpg"
  ) |>
  dplyr::mutate(
    cpg = factor(cpg, levels = top_cpgs),
    cpg_label = paste0(
      final_interpretation_priority,
      ". ",
      as.character(cpg),
      "\n",
      stringr::str_trunc(gene_or_genomic_context, 35)
    )
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
    title = "Candidates prioritized for biological interpretation by disease group",
    subtitle = "Fifteen candidates selected from the statistically ranked CpG set for detailed biological interpretation."
  ) +
  theme_pipeline()

top_effects <- dataset_effects |>
  dplyr::filter(.data$cpg %in% top_cpgs) |>
  dplyr::left_join(
    interpretation_top15 |>
      dplyr::select(cpg, final_interpretation_priority, gene_or_genomic_context),
    by = "cpg"
  ) |>
  dplyr::mutate(
    cpg = factor(cpg, levels = top_cpgs),
    cpg_label = paste0(
      final_interpretation_priority,
      ". ",
      as.character(cpg),
      "\n",
      stringr::str_trunc(gene_or_genomic_context, 35)
    ),
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
    title = "Dataset-level effects for candidates prioritized for biological interpretation",
    subtitle = "Positive estimates mean higher methylation in disease/fibrosis/higher progression score."
  ) +
  theme_pipeline()

ggplot2::ggsave(interpretation_group_plot_path, p_top_group, width = 15, height = 10, dpi = 300)
ggplot2::ggsave(interpretation_effect_plot_path, p_effects, width = 15, height = 10, dpi = 300)

message("Full reviewed Top-30 table: ", review_path)
message("Interpretation Top-15 table: ", interpretation_top15_path)
message("Interpretation Top-15 beta values: ", plot_data_path)
message("Interpretation DiseaseGroup plot: ", interpretation_group_plot_path)
message("Interpretation dataset-effect plot: ", interpretation_effect_plot_path)

message("Candidates prioritized for biological interpretation:")
print(
  interpretation_top15_output |>
    dplyr::select(
      final_interpretation_priority,
      statistical_rank,
      cpg,
      gene_or_genomic_context,
      biological_evidence_class,
      strict_subset_membership,
      strongest_contrast,
      best_fdr,
      mean_i2,
      selection_reason
    )
)
