# 28_build_final_candidate_thesis_outputs.R
# Candidate-only thesis outputs from the finalized sparse-stability pipeline.

if (!file.exists("scripts/00_setup_paths.R") && file.exists("classifier/scripts/00_setup_paths.R")) {
  setwd("classifier")
}

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
})

message("Building candidate-only thesis outputs from finalized sparse-stability files.")

if (!identical(config$analysis_mode %||% "", "sparse_stability")) {
  stop("This thesis-output script must be run with the sparse-stability config.", call. = FALSE)
}

if (!grepl("selected_features_sparse_stability$", out_paths$results_features)) {
  stop("Configured results_features path is not the sparse-stability directory.", call. = FALSE)
}

if (!grepl("plots_sparse_stability$", out_paths$results_plots)) {
  stop("Configured results_plots path is not the sparse-stability directory.", call. = FALSE)
}

master_root <- normalizePath(file.path(project_root, ".."), mustWork = TRUE)

public_results_root <- file.path(master_root, "results")
public_fig_dir <- file.path(public_results_root, "figures")
public_tab_dir <- file.path(public_results_root, "tables")
public_sup_dir <- file.path(public_results_root, "supplementary")
public_sup_fig_dir <- file.path(public_sup_dir, "figures")
public_sup_tab_dir <- file.path(public_sup_dir, "tables")

thesis_out_root <- file.path(master_root, "outputs", "thesis_outputs")
thesis_out_fig_dir <- file.path(thesis_out_root, "figures")
thesis_out_tab_dir <- file.path(thesis_out_root, "tables")
pub_ready_sup_dir <- file.path(thesis_out_root, "supplementary")
pub_ready_sup_fig_dir <- file.path(pub_ready_sup_dir, "figures")
thesis_out_sup_tab_dir <- file.path(pub_ready_sup_dir, "tables")

invisible(lapply(
  c(
    public_fig_dir, public_tab_dir, public_sup_fig_dir, public_sup_tab_dir,
    thesis_out_fig_dir, thesis_out_tab_dir, pub_ready_sup_fig_dir, thesis_out_sup_tab_dir
  ),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

feature_path <- function(filename) file.path(out_paths$results_features, filename)
plot_path <- function(filename) file.path(out_paths$results_plots, filename)

source_files <- list(
  stable_unique = feature_path("elastic_net_cpg_annotation_priority_unique.csv"),
  candidate_priority = feature_path("elastic_net_cpg_candidate_priority.csv"),
  meta_results = feature_path("elastic_net_candidate_cpg_meta_results.csv"),
  meta_shortlist = feature_path("elastic_net_candidate_cpg_meta_shortlist.csv"),
  strict_unique = feature_path("elastic_net_candidate_cpg_meta_strict_unique.csv"),
  candidate_summary = feature_path("elastic_net_cpg_candidate_summary.csv"),
  biological_review = feature_path("elastic_net_cpg_biological_review.csv"),
  interpretation_top15 = feature_path("elastic_net_cpg_interpretation_top15.csv"),
  interpretation_beta = feature_path("elastic_net_cpg_interpretation_top15_beta_values.csv"),
  dataset_effects = feature_path("elastic_net_candidate_cpg_meta_dataset_effects.csv"),
  sparse_meta_heatmap_png = plot_path("elastic_net_candidate_cpg_meta_effect_heatmap.png")
)

missing_sources <- names(source_files)[!file.exists(unlist(source_files))]
if (length(missing_sources) > 0) {
  stop("Missing finalized sparse candidate sources: ", paste(missing_sources, collapse = ", "), call. = FALSE)
}

read_required_csv <- function(path_in) {
  readr::read_csv(path_in, show_col_types = FALSE)
}

stable_unique <- read_required_csv(source_files$stable_unique)
candidate_priority <- read_required_csv(source_files$candidate_priority)
meta_results <- read_required_csv(source_files$meta_results)
meta_shortlist <- read_required_csv(source_files$meta_shortlist)
strict_unique <- read_required_csv(source_files$strict_unique)
candidate_summary <- read_required_csv(source_files$candidate_summary)
biological_review <- read_required_csv(source_files$biological_review)
interpretation_top15 <- read_required_csv(source_files$interpretation_top15)
interpretation_beta <- read_required_csv(source_files$interpretation_beta)
dataset_effects <- read_required_csv(source_files$dataset_effects)

assert_unique_cpg_count <- function(tbl, expected, label) {
  if (!"cpg" %in% names(tbl)) {
    stop(label, " is missing required column cpg.", call. = FALSE)
  }
  observed <- dplyr::n_distinct(tbl$cpg)
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop(label, " has ", observed, " unique CpGs; expected ", expected, ".", call. = FALSE)
  }
  invisible(observed)
}

counts <- list(
  stable = assert_unique_cpg_count(stable_unique, 488, "Stable CpG table"),
  candidate_priority = assert_unique_cpg_count(candidate_priority, 488, "Candidate priority table"),
  cross_supported_summary = assert_unique_cpg_count(candidate_summary, 324, "Candidate summary table"),
  cross_supported_shortlist = assert_unique_cpg_count(meta_shortlist, 324, "Meta-analysis shortlist"),
  strict = assert_unique_cpg_count(strict_unique, 33, "Strict unique table"),
  biological_review = assert_unique_cpg_count(biological_review, 30, "Biological review table"),
  interpretation_panel = assert_unique_cpg_count(interpretation_top15, 15, "Interpretation panel table")
)

if (nrow(interpretation_top15) != 15 || anyDuplicated(interpretation_top15$cpg) > 0) {
  stop("Interpretation panel must contain exactly 15 rows and 15 unique CpGs.", call. = FALSE)
}

if (!setequal(interpretation_top15$cpg, interpretation_beta$cpg)) {
  stop("Interpretation beta table CpGs do not match the final Top-15 interpretation panel.", call. = FALSE)
}

if (!all(interpretation_top15$cpg %in% candidate_summary$cpg)) {
  stop("Final interpretation panel contains CpGs absent from the script-12 candidate summary.", call. = FALSE)
}

legacy_symbol_map <- c(
  SEPW1 = "SELENOW",
  KIAA0391 = "PRORP",
  LASS4 = "CERS4",
  SPG20 = "SPART"
)

pretty_gene_symbols <- function(x) {
  vapply(as.character(x), function(value) {
    if (is.na(value) || !nzchar(trimws(value))) {
      return(NA_character_)
    }

    genes <- unlist(strsplit(value, ";", fixed = TRUE), use.names = FALSE)
    genes <- trimws(genes)
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) == 0) {
      return(NA_character_)
    }

    genes <- ifelse(genes %in% names(legacy_symbol_map), unname(legacy_symbol_map[genes]), genes)
    paste(unique(genes), collapse = "/")
  }, character(1))
}

make_gene_context <- function(tbl) {
  current <- if ("current_gene_symbol" %in% names(tbl)) pretty_gene_symbols(tbl$current_gene_symbol) else rep(NA_character_, nrow(tbl))
  original <- if ("original_annotation" %in% names(tbl)) pretty_gene_symbols(tbl$original_annotation) else rep(NA_character_, nrow(tbl))
  ucsc <- if ("UCSC_RefGene_Name" %in% names(tbl)) pretty_gene_symbols(tbl$UCSC_RefGene_Name) else rep(NA_character_, nrow(tbl))

  dplyr::case_when(
    tbl$cpg == "cg10987682" ~ "intergenic; unannotated",
    tbl$cpg == "cg09870609" ~ "intergenic; unannotated",
    tbl$cpg == "cg08928408" ~ "intergenic; WWTR1-proximal",
    !is.na(current) & nzchar(current) ~ current,
    !is.na(original) & nzchar(original) ~ original,
    !is.na(ucsc) & nzchar(ucsc) ~ ucsc,
    TRUE ~ "intergenic/no-gene"
  )
}

strip_thesis_unannotated_context <- function(tbl) {
  if (!"cpg" %in% colnames(tbl)) {
    return(tbl)
  }

  unannotated <- tbl$cpg %in% c("cg10987682", "cg09870609")
  if (!any(unannotated)) {
    return(tbl)
  }

  nearest_cols <- intersect(
    c(
      "nearest_gene",
      "nearest_gene_entrez",
      "nearest_tss_transcript_id",
      "nearest_tss_chr",
      "nearest_tss_position",
      "nearest_tss_strand",
      "nearest_tss_distance",
      "nearest_tss_orientation"
    ),
    colnames(tbl)
  )
  for (col in nearest_cols) {
    tbl[[col]][unannotated] <- NA
  }

  set_if_present <- function(col, value) {
    if (col %in% colnames(tbl)) {
      tbl[[col]][unannotated] <<- value
    }
  }

  set_if_present("gene_or_genomic_context", "intergenic/no-gene")
  set_if_present(
    "evidence_directness",
    "Intergenic/no-gene candidate; no nearest-gene assignment is made in thesis-facing outputs."
  )
  set_if_present("human_masld_evidence", "No gene/locus-specific human MASLD evidence identified.")
  set_if_present("human_liver_fibrosis_evidence", "No gene/locus-specific human liver-fibrosis evidence identified.")
  set_if_present("experimental_masld_evidence", "No gene/locus-specific experimental MASLD evidence identified.")
  set_if_present("experimental_liver_fibrosis_evidence", "No gene/locus-specific experimental liver-fibrosis evidence identified.")
  set_if_present("hepatic_metabolism_evidence", "No convincing hepatic-metabolism evidence assigned.")
  set_if_present(
    "literature_note",
    "Biologically unresolved intergenic/no-gene candidate retained as a discovery-oriented statistical signal without assigning a nearest gene."
  )

  character_cols <- names(tbl)[vapply(tbl, is.character, logical(1))]
  for (col in character_cols) {
    tbl[[col]][unannotated] <- stringr::str_replace_all(
      tbl[[col]][unannotated],
      c("LINC02923" = "unannotated", "TMC3-AS1" = "unannotated")
    )
  }

  tbl
}

clean_contrast <- function(x) {
  dplyr::recode(
    as.character(x),
    disease_vs_healthy = "Disease vs healthy",
    fibrosis_vs_nonfibrosis = "Fibrosis vs non-fibrosis",
    advanced_vs_nonadvanced = "Advanced vs non-advanced",
    progression_score = "Ordinal progression score",
    .default = str_replace_all(as.character(x), "_", " ")
  )
}

display_dataset <- function(x) {
  dplyr::recode(as.character(x), Kim = "ITEN", .default = as.character(x))
}

display_dataset_text <- function(x) {
  gsub("\\bKim\\b", "ITEN", as.character(x), perl = TRUE)
}

format_yes_no <- function(x) {
  dplyr::if_else(as.logical(x), "yes", "no", missing = "no")
}

meta_key <- meta_results |>
  dplyr::select(
    cpg,
    strongest_contrast = contrast,
    strongest_contrast_i2 = i2,
    strongest_contrast_random_effect = random_effect,
    strongest_contrast_fdr = fdr,
    strongest_contrast_n_datasets = n_datasets,
    strongest_contrast_direction_consistency = direction_consistency
  )

if (anyDuplicated(meta_key[c("cpg", "strongest_contrast")]) > 0) {
  stop("Meta-results contain duplicate CpG/contrast rows.", call. = FALSE)
}

top15_joined <- interpretation_top15 |>
  dplyr::left_join(meta_key, by = c("cpg", "strongest_contrast"))

top15_aug <- top15_joined |>
  dplyr::mutate(
    gene_context_display = make_gene_context(top15_joined),
    contrast_display = clean_contrast(strongest_contrast),
    strict_subset_display = format_yes_no(strict_subset_membership),
    pooled_direction = dplyr::case_when(
      strongest_contrast_random_effect > 0 ~ "higher methylation in case/higher stage",
      strongest_contrast_random_effect < 0 ~ "lower methylation in case/higher stage",
      TRUE ~ top_effect_direction
    ),
    cpg_panel_label = paste0(
      final_interpretation_priority,
      ". ",
      cpg,
      "\n",
      stringr::str_wrap(gene_context_display, width = 28)
    )
  )

if (any(is.na(top15_aug$strongest_contrast_i2))) {
  stop("Could not retrieve strongest-contrast I2 for all final interpretation CpGs.", call. = FALSE)
}

if (any(abs(top15_aug$best_fdr - top15_aug$strongest_contrast_fdr) > 1e-12, na.rm = TRUE)) {
  stop("Top-15 best_fdr does not match the CpG + strongest_contrast FDR in script-11 meta-results.", call. = FALSE)
}

if (any(top15_aug$pooled_direction != top15_aug$top_effect_direction, na.rm = TRUE)) {
  stop("Top-15 pooled direction does not match the script-12/script-13 direction field.", call. = FALSE)
}

top15_table <- top15_aug |>
  dplyr::arrange(final_interpretation_priority, statistical_rank, cpg) |>
  dplyr::transmute(
    interpretation_priority = as.integer(final_interpretation_priority),
    statistical_rank = as.integer(statistical_rank),
    CpG = cpg,
    gene_genomic_context = gene_context_display,
    strongest_contrast = strongest_contrast,
    strongest_supported_contrast = contrast_display,
    direction = pooled_direction,
    FDR = best_fdr,
    n_datasets = as.integer(strongest_contrast_n_datasets),
    strongest_contrast_i2 = strongest_contrast_i2,
    strict_internal_support_subset = strict_subset_display,
    biological_evidence_class = biological_evidence_class
  )

strict_table <- strict_unique |>
  dplyr::mutate(
    gene_genomic_context = make_gene_context(strict_unique),
    datasets = display_dataset_text(datasets),
    likely_dataset_sensitive = format_yes_no(likely_dataset_sensitive)
  ) |>
  dplyr::transmute(
    CpG = cpg,
    gene_genomic_context,
    n_strict_contrasts,
    strict_contrasts,
    best_FDR = best_fdr,
    max_abs_random_effect,
    datasets,
    priority_label,
    best_nonzero_folds,
    likely_dataset_sensitive,
    r2_disease_group,
    r2_dataset,
    disease_r2_after_dataset
  ) |>
  dplyr::arrange(best_FDR, CpG)

supported_candidates_table <- candidate_summary |>
  dplyr::mutate(
    gene_genomic_context = make_gene_context(candidate_summary),
    strongest_supported_contrast = clean_contrast(top_contrast),
    strict_internal_support_subset = cpg %in% strict_unique$cpg
  ) |>
  dplyr::transmute(
    cpg,
    gene_genomic_context,
    statistical_rank = dplyr::row_number(),
    candidate_priority_score,
    robust_contrasts,
    tier1_contrasts,
    tier2_contrasts,
    tier3_contrasts,
    same_direction_contrasts,
    max_n_datasets,
    best_fdr,
    best_neg_log10_fdr,
    mean_i2,
    min_i2,
    strongest_contrast = top_contrast,
    strongest_supported_contrast,
    top_random_effect,
    top_effect_direction,
    max_abs_meta_effect,
    strict_internal_support_subset,
    likely_dataset_sensitive,
    priority_label,
    best_nonzero_folds,
    all5_any_run,
    r2_disease_group,
    r2_dataset,
    disease_r2_after_dataset,
    dataset_warning,
    UCSC_RefGene_Name,
    chr,
    pos,
    Relation_to_Island,
    UCSC_RefGene_Group,
    Regulatory_Feature_Group,
    Enhancer,
    DHS,
    DMR
  ) |>
  dplyr::arrange(statistical_rank, cpg)

biological_review_table <- biological_review |>
  strip_thesis_unannotated_context() |>
  dplyr::mutate(gene_genomic_context = make_gene_context(biological_review)) |>
  dplyr::arrange(statistical_rank, cpg)

selection_summary <- tibble::tibble(
  stage_id = c(
    "stable_elastic_net_selected_cpgs",
    "cross_study_supported_candidates",
    "strict_internal_support_subset",
    "manual_biological_review_top30",
    "interpretation_panel_top15"
  ),
  stage_label = c(
    "Stable Elastic-Net-selected CpGs",
    "Cross-study-supported Elastic-Net-derived candidates",
    "Strict internal-support subset",
    "Manual biological review of statistical Top 30",
    "Candidates prioritized for biological interpretation"
  ),
  unique_cpgs = as.integer(c(
    counts$stable,
    counts$cross_supported_summary,
    counts$strict,
    counts$biological_review,
    counts$interpretation_panel
  )),
  source_file = c(
    basename(source_files$stable_unique),
    basename(source_files$candidate_summary),
    basename(source_files$strict_unique),
    basename(source_files$biological_review),
    basename(source_files$interpretation_top15)
  ),
  filter_summary = c(
    "Sparse stability union after applying nonzero_folds >= 4 within each eligible standard LOSO analysis.",
    "Post-selection study-specific models and random-effects meta-analysis; not independent validation.",
    "Conservative internal-support subset from the candidate meta-analysis.",
    "Manual structured literature curation of the top 30 statistically ranked candidates.",
    "Manual interpretation panel selected only by include_in_interpretation_panel and ordered by final_interpretation_priority."
  )
)

write_main_table <- function(tbl, filename) {
  paths <- c(
    results = file.path(public_tab_dir, filename),
    thesis_outputs = file.path(thesis_out_tab_dir, filename)
  )
  readr::write_csv(tbl, paths[["results"]])
  readr::write_csv(tbl, paths[["thesis_outputs"]])
  paths
}

write_supp_table <- function(tbl, filename) {
  paths <- c(
    results = file.path(public_sup_tab_dir, filename),
    thesis_outputs = file.path(thesis_out_sup_tab_dir, filename)
  )
  readr::write_csv(tbl, paths[["results"]])
  readr::write_csv(tbl, paths[["thesis_outputs"]])
  paths
}

save_plot_set <- function(plot, stem, dirs, width, height, dpi = 600) {
  paths <- character()
  for (dir_out in dirs) {
    png_path <- file.path(dir_out, paste0(stem, ".png"))
    pdf_path <- file.path(dir_out, paste0(stem, ".pdf"))
    ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
    ggplot2::ggsave(pdf_path, plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
    paths <- c(paths, png_path, pdf_path)
  }
  paths
}

copy_existing <- function(src, dst) {
  if (!file.exists(src)) {
    warning("Missing source file: ", src, call. = FALSE)
    return(FALSE)
  }
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  file.copy(src, dst, overwrite = TRUE)
}

png_to_pdf <- function(png_path, pdf_path) {
  convert_bin <- Sys.which("convert")
  if (nzchar(convert_bin) && file.exists(png_path)) {
    status <- suppressWarnings(system2(convert_bin, c(png_path, pdf_path)))
    return(identical(status, 0L))
  }

  if (requireNamespace("png", quietly = TRUE) && file.exists(png_path)) {
    img <- png::readPNG(png_path)
    grDevices::cairo_pdf(pdf_path, width = 12, height = 8)
    grid::grid.newpage()
    grid::grid.raster(img)
    grDevices::dev.off()
    return(TRUE)
  }

  FALSE
}

update_manifest <- function(manifest_path, rows, key_col) {
  rows <- rows |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))

  existing <- if (file.exists(manifest_path)) {
    readr::read_csv(manifest_path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  } else {
    tibble::tibble()
  }

  if (nrow(existing) > 0) {
    missing_cols <- setdiff(names(rows), names(existing))
    for (col in missing_cols) {
      existing[[col]] <- NA_character_
    }
    existing <- existing |>
      dplyr::select(dplyr::all_of(names(rows))) |>
      dplyr::filter(!.data[[key_col]] %in% rows[[key_col]])
  }

  dplyr::bind_rows(existing, rows) |>
    dplyr::arrange(.data[[key_col]]) |>
    readr::write_csv(manifest_path)
}

table_outputs <- c(
  write_main_table(top15_table, "T08_interpretation_top15_cpg_candidates.csv"),
  write_main_table(strict_table, "T09_strict_internal_support_cpg_candidates.csv"),
  write_supp_table(selection_summary, "S04_cpg_filtering_selection_summary.csv"),
  write_supp_table(supported_candidates_table, "S08_cross_study_supported_elastic_net_candidates.csv"),
  write_supp_table(biological_review_table, "S09_top30_biological_review.csv")
)

theme_pub <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1, margin = margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = base_size, margin = margin(b = 6)),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey65", linewidth = 0.25),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_line(color = "grey88", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = margin(7, 7, 7, 7)
    )
}

nice_disease <- c(
  Healthy = "Healthy",
  Healthy_Obese = "Healthy obese",
  MASL = "MASL",
  MASH = "MASH",
  Mild_Fibrosis = "Mild fibrosis",
  Advanced_Fibrosis = "Advanced fibrosis"
)

pal_disease <- c(
  Healthy = "#4E9F6D",
  Healthy_Obese = "#91C77B",
  MASL = "#D8B447",
  MASH = "#E59C67",
  Mild_Fibrosis = "#D9822B",
  Advanced_Fibrosis = "#B64A4A"
)

class4_palette <- c(
  Healthy = pal_disease[["Healthy"]],
  `MASL/MASH` = "#F0CC61",
  `Mild fibrosis` = pal_disease[["Mild_Fibrosis"]],
  `Advanced fibrosis` = pal_disease[["Advanced_Fibrosis"]]
)

dataset_levels <- c("Ahrens", "Horvath", "Johnson", "ITEN", "Murphy", "VanDijck")
shape_dataset <- c(Ahrens = 15, Horvath = 16, Johnson = 17, ITEN = 8, Murphy = 3, VanDijck = 18)

beta_plot <- interpretation_beta |>
  dplyr::select(sample_id, cpg, beta_value, Dataset, DiseaseGroup, Progression3) |>
  dplyr::left_join(
    top15_aug |>
      dplyr::select(cpg, final_interpretation_priority, cpg_panel_label),
    by = "cpg"
  ) |>
  dplyr::mutate(
    DiseaseGroup = factor(as.character(DiseaseGroup), levels = names(nice_disease), labels = nice_disease),
    Dataset = factor(display_dataset(Dataset), levels = dataset_levels),
    cpg_panel_label = factor(cpg_panel_label, levels = top15_aug$cpg_panel_label[order(top15_aug$final_interpretation_priority)])
  )

p_beta_disease <- ggplot(beta_plot, aes(x = DiseaseGroup, y = beta_value, fill = DiseaseGroup)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "black", linewidth = 0.22, alpha = 0.78) +
  geom_point(
    aes(shape = Dataset),
    position = position_jitter(width = 0.13, height = 0, seed = 123),
    color = "grey25",
    alpha = 0.35,
    size = 0.55,
    na.rm = TRUE
  ) +
  facet_wrap(~ cpg_panel_label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = setNames(pal_disease, nice_disease[names(pal_disease)]), drop = FALSE) +
  scale_shape_manual(values = shape_dataset, drop = FALSE) +
  labs(
    title = "Candidates prioritized for biological interpretation by disease group",
    subtitle = "Fifteen candidates selected from the statistically ranked CpG set for detailed biological interpretation.",
    x = "Disease group",
    y = "DNA methylation beta value",
    fill = "Disease group",
    shape = "Study"
  ) +
  theme_pub(7.4) +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1),
    legend.key.size = grid::unit(0.34, "cm")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 1, byrow = TRUE))

beta_class4 <- beta_plot |>
  dplyr::mutate(
    Class4 = dplyr::case_when(
      as.character(DiseaseGroup) %in% c("Healthy", "Healthy obese") ~ "Healthy",
      as.character(DiseaseGroup) %in% c("MASL", "MASH") ~ "MASL/MASH",
      as.character(DiseaseGroup) == "Mild fibrosis" ~ "Mild fibrosis",
      as.character(DiseaseGroup) == "Advanced fibrosis" ~ "Advanced fibrosis",
      TRUE ~ NA_character_
    ),
    Class4 = factor(Class4, levels = names(class4_palette))
  ) |>
  dplyr::filter(!is.na(Class4))

p_beta_class4 <- ggplot(beta_class4, aes(x = Class4, y = beta_value, fill = Class4)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "black", linewidth = 0.22, alpha = 0.78) +
  geom_point(
    aes(shape = Dataset),
    position = position_jitter(width = 0.12, height = 0, seed = 123),
    color = "grey25",
    alpha = 0.33,
    size = 0.55,
    na.rm = TRUE
  ) +
  facet_wrap(~ cpg_panel_label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = class4_palette, drop = FALSE) +
  scale_shape_manual(values = shape_dataset, drop = FALSE) +
  labs(
    title = "Candidates prioritized for biological interpretation by 4-class disease group",
    subtitle = "Disease labels are collapsed to match the primary 4-class classifier outcome.",
    x = "4-class disease group",
    y = "DNA methylation beta value",
    fill = "4-class group",
    shape = "Study"
  ) +
  theme_pub(7.4) +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1),
    legend.key.size = grid::unit(0.34, "cm")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE), shape = guide_legend(nrow = 1, byrow = TRUE))

effects_plot <- dataset_effects |>
  dplyr::inner_join(
    top15_aug |>
      dplyr::select(cpg, strongest_contrast, cpg_panel_label),
    by = "cpg"
  ) |>
  dplyr::filter(contrast == strongest_contrast) |>
  dplyr::mutate(
    ci_low = estimate - 1.96 * std_error,
    ci_high = estimate + 1.96 * std_error,
    dataset = factor(display_dataset(dataset), levels = rev(dataset_levels)),
    cpg_panel_label = factor(cpg_panel_label, levels = top15_aug$cpg_panel_label[order(top15_aug$final_interpretation_priority)])
  )

if (nrow(effects_plot) == 0) {
  stop("No dataset-level effects matched the final Top-15 strongest contrasts.", call. = FALSE)
}

p_effects <- ggplot(effects_plot, aes(x = estimate, y = dataset)) +
  geom_vline(xintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.14, linewidth = 0.4, color = "grey35", orientation = "y") +
  geom_point(size = 1.6, color = "#0072B2") +
  facet_wrap(~ cpg_panel_label, scales = "free_x", ncol = 5) +
  labs(
    title = "Dataset-level effects for candidates prioritized for biological interpretation",
    x = "Within-study regression estimate (95% CI)",
    y = "Study"
  ) +
  theme_pub(7.4)

meta_plot <- meta_results |>
  dplyr::filter(cpg %in% candidate_summary$cpg) |>
  dplyr::mutate(
    neg_log10_fdr = -log10(pmax(fdr, .Machine$double.xmin)),
    interpretation_panel = cpg %in% top15_aug$cpg,
    contrast_display = factor(
      clean_contrast(contrast),
      levels = c("Disease vs healthy", "Fibrosis vs non-fibrosis", "Advanced vs non-advanced", "Ordinal progression score")
    )
  )

meta_labels <- meta_plot |>
  dplyr::filter(interpretation_panel, fdr < 0.01) |>
  dplyr::left_join(
    top15_aug |> dplyr::select(cpg, gene_context_display),
    by = "cpg"
  ) |>
  dplyr::group_by(contrast_display) |>
  dplyr::arrange(fdr, dplyr::desc(abs(random_effect)), .by_group = TRUE) |>
  dplyr::slice_head(n = 4) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    point_label = if_else(
      stringr::str_detect(gene_context_display, regex("^intergenic", ignore_case = TRUE)),
      cpg,
      paste0(gene_context_display, "\n", cpg)
    )
  )

p_meta <- ggplot(meta_plot, aes(x = random_effect, y = neg_log10_fdr)) +
  geom_vline(xintercept = 0, color = "grey65", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), color = "grey55", linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(color = interpretation_panel), alpha = 0.75, size = 1.25)

if (requireNamespace("ggrepel", quietly = TRUE) && nrow(meta_labels) > 0) {
  p_meta <- p_meta +
    ggrepel::geom_text_repel(
      data = meta_labels,
      aes(label = point_label),
      color = "grey10",
      size = 2.15,
      lineheight = 0.9,
      box.padding = 0.22,
      point.padding = 0.18,
      min.segment.length = 0,
      segment.color = "grey45",
      segment.size = 0.2,
      max.overlaps = Inf,
      seed = 123,
      show.legend = FALSE
    )
} else if (nrow(meta_labels) > 0) {
  p_meta <- p_meta +
    geom_text(
      data = meta_labels,
      aes(label = point_label),
      color = "grey10",
      size = 2.1,
      lineheight = 0.9,
      vjust = -0.55,
      check_overlap = TRUE,
      show.legend = FALSE
    )
}

p_meta <- p_meta +
  facet_wrap(~ contrast_display, scales = "free", nrow = 2) +
  scale_color_manual(
    values = c(`FALSE` = "grey72", `TRUE` = "#B64A4A"),
    labels = c("Other cross-study-supported candidates", "Interpretation panel")
  ) +
  labs(
    title = "Cross-study-supported Elastic-Net-derived candidates",
    subtitle = "Final interpretation-panel candidates are highlighted within the post-selection meta-analysis results.",
    x = "Random-effects estimate",
    y = expression(-log[10]("FDR")),
    color = NULL
  ) +
  theme_pub(8)

flow_nodes <- selection_summary |>
  dplyr::mutate(
    x = c(1, 2, 3, 4, 5),
    y = c(1, 1, 1, 1, 1),
    label = paste0(stage_label, "\n", unique_cpgs, " CpGs")
  )

flow_edges <- tibble::tibble(
  x = flow_nodes$x[-nrow(flow_nodes)] + 0.28,
  xend = flow_nodes$x[-1] - 0.28,
  y = 1,
  yend = 1
)

p_flow <- ggplot() +
  geom_segment(
    data = flow_edges,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed"),
    linewidth = 0.55,
    color = "grey35"
  ) +
  geom_label(
    data = flow_nodes,
    aes(x = x, y = y, label = stringr::str_wrap(label, width = 22), fill = stage_id),
    color = "black",
    fontface = "bold",
    size = 2.6,
    lineheight = 0.92,
    label.padding = grid::unit(0.30, "lines"),
    label.r = grid::unit(0.08, "lines"),
    linewidth = 0.25
  ) +
  annotate(
    "text",
    x = 3,
    y = 0.48,
    label = "Post-selection evidence from the same study collection; not independent validation.",
    size = 2.8,
    fontface = "italic",
    color = "grey25"
  ) +
  scale_fill_manual(
    values = c(
      stable_elastic_net_selected_cpgs = "#DCEBFA",
      cross_study_supported_candidates = "#E7E2F0",
      strict_internal_support_subset = "#E9F4D8",
      manual_biological_review_top30 = "#FFF0C7",
      interpretation_panel_top15 = "#F3D5DA"
    ),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0.52, 5.48), ylim = c(0.35, 1.55), clip = "off") +
  labs(
    title = "Candidate CpG evidence hierarchy",
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 8)))

hierarchy_nodes <- tibble::tibble(
  x = c(0.09, 0.33, 0.57, 0.81),
  y = c(0.64, 0.64, 0.64, 0.64),
  w = c(0.18, 0.20, 0.19, 0.18),
  h = c(0.18, 0.18, 0.18, 0.18),
  title = c(
    "Stable Elastic-Net-selected CpGs",
    "Cross-study-supported candidates",
    "Manual review Top 30",
    "Interpretation panel"
  ),
  number = c(counts$stable, counts$cross_supported_summary, counts$biological_review, counts$interpretation_panel),
  color = c("#2F5D8C", "#6A4C93", "#C79A2B", "#1B9AAA")
)

branch_node <- tibble::tibble(
  x = 0.43,
  y = 0.24,
  w = 0.22,
  h = 0.16,
  title = "Strict internal-support subset",
  number = counts$strict,
  color = "#C44536"
)

p_hierarchy <- ggplot() +
  geom_segment(aes(x = 0.27, xend = 0.33, y = 0.73, yend = 0.73), arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed"), color = "grey35", linewidth = 0.55) +
  geom_segment(aes(x = 0.53, xend = 0.57, y = 0.73, yend = 0.73), arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed"), color = "grey35", linewidth = 0.55) +
  geom_segment(aes(x = 0.76, xend = 0.81, y = 0.73, yend = 0.73), arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed"), color = "grey35", linewidth = 0.55) +
  geom_segment(aes(x = 0.43, xend = 0.51, y = 0.64, yend = 0.40), arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed"), color = "grey35", linewidth = 0.55) +
  geom_label(
    data = hierarchy_nodes,
    aes(x = x + w / 2, y = y + h / 2, label = paste0(number, "\n", stringr::str_wrap(title, 21)), fill = title),
    color = "white",
    fontface = "bold",
    size = 3.1,
    lineheight = 0.95,
    label.padding = grid::unit(0.34, "lines"),
    label.r = grid::unit(0.10, "lines"),
    linewidth = 0
  ) +
  geom_label(
    data = branch_node,
    aes(x = x + w / 2, y = y + h / 2, label = paste0(number, "\n", stringr::str_wrap(title, 24)), fill = title),
    color = "white",
    fontface = "bold",
    size = 3.0,
    lineheight = 0.95,
    label.padding = grid::unit(0.34, "lines"),
    label.r = grid::unit(0.10, "lines"),
    linewidth = 0
  ) +
  annotate(
    "text",
    x = 0.50,
    y = 0.10,
    label = "The strict subset and interpretation panel are complementary internal views of post-selection candidates.",
    color = "grey20",
    size = 3.0,
    fontface = "italic"
  ) +
  scale_fill_manual(
    values = c(
      setNames(hierarchy_nodes$color, hierarchy_nodes$title),
      setNames(branch_node$color, branch_node$title)
    ),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.95), clip = "off") +
  labs(
    title = "Candidate CpG evidence map",
    subtitle = "Final sparse-stability candidate outputs reported in the thesis tables and figures.",
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 3)),
    plot.subtitle = element_text(hjust = 0.5, color = "grey25", margin = margin(b = 8))
  )

plot_outputs <- c(
  save_plot_set(
    p_beta_disease,
    "F09_interpretation_top15_beta_by_disease_group",
    c(public_fig_dir, thesis_out_fig_dir),
    width = 12,
    height = 8.2
  ),
  save_plot_set(
    p_effects,
    "F10_interpretation_top15_dataset_effects",
    c(public_fig_dir, thesis_out_fig_dir),
    width = 12,
    height = 8.0
  ),
  save_plot_set(
    p_meta,
    "F11_candidate_cpg_meta_analysis",
    c(public_fig_dir, thesis_out_fig_dir),
    width = 9.4,
    height = 6.2
  ),
  save_plot_set(
    p_beta_class4,
    "S01_interpretation_top15_beta_by_4class_compact",
    c(public_sup_fig_dir, pub_ready_sup_fig_dir),
    width = 12,
    height = 8.1
  ),
  save_plot_set(
    p_flow,
    "S04_cpg_filtering_selection_flow",
    c(public_sup_fig_dir, pub_ready_sup_fig_dir),
    width = 12,
    height = 3.7
  ),
  save_plot_set(
    p_hierarchy,
    "F33_candidate_cpg_hierarchical_evidence_map",
    c(public_fig_dir, thesis_out_fig_dir),
    width = 12,
    height = 6.8
  )
)

s02_outputs <- character()
for (dir_out in c(public_sup_fig_dir, pub_ready_sup_fig_dir)) {
  png_dst <- file.path(dir_out, "S02_candidate_cpg_meta_effect_heatmap.png")
  pdf_dst <- file.path(dir_out, "S02_candidate_cpg_meta_effect_heatmap.pdf")
  copy_existing(source_files$sparse_meta_heatmap_png, png_dst)
  png_to_pdf(png_dst, pdf_dst)
  s02_outputs <- c(s02_outputs, png_dst, pdf_dst)
}

candidate_table_manifest <- tibble::tribble(
  ~table, ~file, ~main_message,
  "T08", "T08_interpretation_top15_cpg_candidates.csv", "Fifteen candidates prioritized for biological interpretation from the corrected sparse-stability candidate pipeline.",
  "T09", "T09_strict_internal_support_cpg_candidates.csv", "Strict internal-support subset among cross-study-supported Elastic-Net-derived candidates."
)

main_table_manifest <- candidate_table_manifest

candidate_figure_manifest <- tibble::tribble(
  ~figure, ~file_stem, ~main_message,
  "F09", "F09_interpretation_top15_beta_by_disease_group", "Candidates prioritized for biological interpretation show beta-value distributions across disease groups.",
  "F10", "F10_interpretation_top15_dataset_effects", "Dataset-level effects for candidates prioritized for biological interpretation.",
  "F11", "F11_candidate_cpg_meta_analysis", "Cross-study-supported Elastic-Net-derived candidates with the final interpretation panel highlighted.",
  "F33", "F33_candidate_cpg_hierarchical_evidence_map", "Candidate CpG evidence map for the corrected sparse-stability pipeline."
)

main_figure_manifest <- candidate_figure_manifest

supp_figure_manifest <- tibble::tribble(
  ~figure, ~file_stem, ~main_message, ~placement,
  "S01", "S01_interpretation_top15_beta_by_4class_compact", "Final interpretation-panel CpGs across the primary 4-class disease grouping.", "supplementary",
  "S02", "S02_candidate_cpg_meta_effect_heatmap", "Dataset-level candidate CpG effects from sparse-stability meta-analysis outputs.", "supplementary",
  "S04", "S04_cpg_filtering_selection_flow", "Candidate CpG workflow documenting sparse stability, post-selection support, strict internal support, manual review, and interpretation-panel selection.", "supplementary"
)

supp_table_manifest <- tibble::tribble(
  ~table, ~file, ~main_message,
  "S04", "S04_cpg_filtering_selection_summary.csv", "Candidate CpG workflow counts from finalized sparse-stability outputs.",
  "S08", "S08_cross_study_supported_elastic_net_candidates.csv", "Machine-readable statistical summary for cross-study-supported Elastic-Net-derived candidates.",
  "S09", "S09_top30_biological_review.csv", "Machine-readable manually curated biological review for the statistical Top 30 candidates."
)

update_manifest(file.path(public_tab_dir, "table_manifest.csv"), main_table_manifest, "table")
update_manifest(file.path(public_tab_dir, "figure_manifest.csv"), main_figure_manifest, "figure")
update_manifest(file.path(thesis_out_tab_dir, "thesis_table_manifest.csv"), main_table_manifest, "table")
update_manifest(file.path(thesis_out_tab_dir, "thesis_figure_manifest.csv"), main_figure_manifest, "figure")
update_manifest(file.path(public_sup_tab_dir, "supplementary_figure_manifest.csv"), supp_figure_manifest, "figure")
update_manifest(file.path(public_sup_tab_dir, "supplementary_table_manifest.csv"), supp_table_manifest, "table")
update_manifest(file.path(thesis_out_sup_tab_dir, "supplementary_figure_manifest.csv"), supp_figure_manifest, "figure")
update_manifest(file.path(thesis_out_sup_tab_dir, "supplementary_table_manifest.csv"), supp_table_manifest, "table")

readr::write_lines(
  c(
    "# Publication-ready supplementary candidate CpG outputs",
    "",
    "These files are derived from the corrected sparse-stability candidate pipeline.",
    "",
    "- S01_interpretation_top15_beta_by_4class_compact: final interpretation-panel CpGs across the primary 4-class disease grouping.",
    "- S02_candidate_cpg_meta_effect_heatmap: dataset-level meta-analysis effect heatmap from sparse-stability outputs.",
    "- S04_cpg_filtering_selection_flow: workflow counts for stable CpGs, cross-study-supported candidates, strict internal-support subset, manual review, and final interpretation panel.",
    "",
    "The candidate evidence is post-selection evidence from the same study collection, not independent validation."
  ),
  file.path(pub_ready_sup_dir, "README_supplementary.md")
)

all_output_paths <- c(table_outputs, plot_outputs, s02_outputs)
output_manifest <- tibble::tibble(
  output_path = all_output_paths,
  exists = file.exists(all_output_paths),
  size_bytes = ifelse(file.exists(all_output_paths), file.info(all_output_paths)$size, NA_real_)
)

readr::write_csv(
  output_manifest,
  file.path(thesis_out_tab_dir, "candidate_thesis_output_manifest.csv")
)

message("Candidate-only thesis outputs complete.")
message("Final counts: stable=", counts$stable,
        " cross-study-supported=", counts$cross_supported_summary,
        " strict=", counts$strict,
        " interpretation-panel=", counts$interpretation_panel)
message("Main Top-15 table: ", file.path(thesis_out_tab_dir, "T08_interpretation_top15_cpg_candidates.csv"))
message("Output manifest: ", file.path(thesis_out_tab_dir, "candidate_thesis_output_manifest.csv"))

print(top15_table, n = Inf)
