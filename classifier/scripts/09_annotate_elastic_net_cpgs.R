# 09_annotate_elastic_net_cpgs.R
# Annotate stable Elastic Net CpGs from LOSO coefficient outputs.

source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(minfi)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

min_nonzero_folds <- as.integer(Sys.getenv("CLASSIFIER_ANNOTATION_MIN_NONZERO_FOLDS", "4"))

if (is.na(min_nonzero_folds) || min_nonzero_folds < 1) {
  stop("CLASSIFIER_ANNOTATION_MIN_NONZERO_FOLDS must be a positive integer.")
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

parse_run_id <- function(path_in) {
  run_id <- basename(path_in)
  run_id <- sub("^elastic_net_loso_nonzero_feature_frequency_", "", run_id)
  sub("\\.csv$", "", run_id)
}

infer_outcome <- function(run_id) {
  sub("_cpg_only_loso$", "", run_id)
}

label_priority_set <- function(outcome) {
  dplyr::case_when(
    outcome == "disease_group_4" ~ "main_stage_model",
    outcome == "disease_group_3" ~ "stage_3_reference",
    outcome == "binary_healthy_vs_disease" ~ "binary_standard_reference",
    TRUE ~ "deprioritized_sensitivity"
  )
}

priority_rank <- function(priority_set) {
  dplyr::case_when(
    priority_set == "main_stage_model" ~ 1L,
    priority_set == "stage_3_reference" ~ 3L,
    priority_set == "binary_standard_reference" ~ 5L,
    TRUE ~ 6L
  )
}

load_annotation <- function() {
  suppressPackageStartupMessages({
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  })

  anno_450k_raw <- minfi::getAnnotation(
    IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19
  ) |>
    as.data.frame() |>
    tibble::rownames_to_column("cpg")

  anno_epic_raw <- NULL
  epic_cpgs <- character()

  if (requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
    suppressPackageStartupMessages({
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    })

    anno_epic_raw <- minfi::getAnnotation(
      IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19
    ) |>
      as.data.frame() |>
      tibble::rownames_to_column("cpg")

    epic_cpgs <- anno_epic_raw$cpg
  }

  anno_450k_cpgs <- anno_450k_raw$cpg

  anno_450k <- minfi::getAnnotation(
    IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19
  ) |>
    as.data.frame() |>
    tibble::rownames_to_column("cpg") |>
    dplyr::mutate(
      annotation_source = "450K_hg19",
      in_450k_annotation = .data$cpg %in% anno_450k_cpgs,
      in_epic_annotation = .data$cpg %in% epic_cpgs
    )

  if (is.null(anno_epic_raw)) {
    return(anno_450k)
  }

  anno_epic <- anno_epic_raw |>
    dplyr::mutate(
      annotation_source = "EPIC_hg19_fallback",
      in_450k_annotation = .data$cpg %in% anno_450k_cpgs,
      in_epic_annotation = .data$cpg %in% epic_cpgs
    )

  dplyr::bind_rows(
    anno_450k,
    dplyr::anti_join(anno_epic, anno_450k |> dplyr::select(cpg), by = "cpg")
  )
}

feature_files <- list.files(
  out_paths$results_features,
  pattern = "^elastic_net_loso_nonzero_feature_frequency_.*\\.csv$",
  full.names = TRUE
)

if (length(feature_files) == 0) {
  stop("No LOSO non-zero feature frequency files found in: ", out_paths$results_features)
}

feature_file_tbl <- tibble::tibble(
  path = feature_files,
  run_id = vapply(feature_files, parse_run_id, character(1))
)

eligible_run_ids <- config$stability_selection$eligible_run_ids %||% character()
eligible_run_ids <- as.character(unlist(eligible_run_ids, use.names = FALSE))
eligible_run_ids <- eligible_run_ids[!is.na(eligible_run_ids) & nzchar(eligible_run_ids)]

if (length(eligible_run_ids) > 0) {
  missing_run_ids <- setdiff(eligible_run_ids, feature_file_tbl$run_id)
  if (length(missing_run_ids) > 0) {
    stop(
      "Configured stability_selection eligible run ID(s) have no matching non-zero feature frequency file: ",
      paste(missing_run_ids, collapse = ", "),
      call. = FALSE
    )
  }

  feature_file_tbl <- feature_file_tbl |>
    dplyr::filter(run_id %in% eligible_run_ids) |>
    dplyr::mutate(run_order = match(run_id, eligible_run_ids)) |>
    dplyr::arrange(run_order)

  if (!all(feature_file_tbl$run_id %in% eligible_run_ids)) {
    stop("Internal whitelist error: matched feature files outside eligible run IDs.", call. = FALSE)
  }
}

feature_files <- feature_file_tbl$path

message("Annotating stable Elastic Net CpGs.")
message("Minimum non-zero LOSO folds: ", min_nonzero_folds)
message("Feature files selected: ", length(feature_files))
message("Contributing sparse-stability runs:")
purrr::walk2(
  feature_file_tbl$run_id,
  feature_file_tbl$path,
  ~ message("  ", .x, " -> ", .y)
)

feature_evidence <- purrr::map_dfr(feature_files, function(path_in) {
  run_id <- parse_run_id(path_in)
  outcome <- infer_outcome(run_id)
  priority_set <- label_priority_set(outcome)

  readr::read_csv(path_in, show_col_types = FALSE) |>
    dplyr::mutate(
      run_id = run_id,
      outcome = outcome,
      priority_set = priority_set,
      priority_rank = priority_rank(priority_set),
      .before = 1
    )
})

stable_evidence <- feature_evidence |>
  dplyr::filter(nonzero_folds >= min_nonzero_folds) |>
  dplyr::mutate(
    cpg = feature,
    coefficient_direction = dplyr::case_when(
      mean_coefficient > 0 ~ "positive",
      mean_coefficient < 0 ~ "negative",
      TRUE ~ "zero_or_mixed"
    )
  ) |>
  dplyr::arrange(
    priority_rank,
    dplyr::desc(nonzero_folds),
    dplyr::desc(mean_abs_coefficient),
    cpg
  )

if (nrow(stable_evidence) == 0) {
  stop("No CpGs passed nonzero_folds >= ", min_nonzero_folds, ".")
}

annotation <- load_annotation()

annotation_columns <- c(
  "cpg",
  "annotation_source",
  "in_450k_annotation",
  "in_epic_annotation",
  "chr",
  "pos",
  "strand",
  "Probe_rs",
  "Probe_maf",
  "CpG_rs",
  "CpG_maf",
  "SBE_rs",
  "SBE_maf",
  "Islands_Name",
  "Relation_to_Island",
  "UCSC_RefGene_Name",
  "UCSC_RefGene_Accession",
  "UCSC_RefGene_Group",
  "Phantom",
  "DMR",
  "Enhancer",
  "HMM_Island",
  "Regulatory_Feature_Name",
  "Regulatory_Feature_Group",
  "DHS"
)

annotation <- annotation |>
  dplyr::select(dplyr::any_of(annotation_columns))

annotated_long <- stable_evidence |>
  dplyr::left_join(annotation, by = "cpg") |>
  dplyr::relocate(cpg, .after = feature)

annotation_fields <- intersect(annotation_columns, colnames(annotated_long))
annotation_fields <- setdiff(annotation_fields, "cpg")

annotated_unique <- annotated_long |>
  dplyr::group_by(cpg) |>
  dplyr::summarise(
    evidence_rows = dplyr::n(),
    evidence_runs = collapse_unique(run_id),
    evidence_outcomes = collapse_unique(outcome),
    priority_sets = collapse_unique(priority_set),
    best_priority_rank = min(priority_rank, na.rm = TRUE),
    best_nonzero_folds = max(nonzero_folds, na.rm = TRUE),
    all5_runs = collapse_unique(run_id[nonzero_folds == 5]),
    max_mean_abs_coefficient = max(mean_abs_coefficient, na.rm = TRUE),
    mean_abs_coefficient_across_runs = mean(mean_abs_coefficient, na.rm = TRUE),
    coefficient_direction_summary = collapse_unique(
      paste(run_id, classes, coefficient_direction, sep = ":")
    ),
    dplyr::across(dplyr::all_of(annotation_fields), ~ dplyr::first(.x)),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    best_priority_rank,
    dplyr::desc(best_nonzero_folds),
    dplyr::desc(max_mean_abs_coefficient),
    cpg
  )

gene_summary <- annotated_unique |>
  dplyr::filter(!is.na(UCSC_RefGene_Name), UCSC_RefGene_Name != "") |>
  dplyr::select(
    cpg,
    best_nonzero_folds,
    max_mean_abs_coefficient,
    evidence_runs,
    evidence_outcomes,
    UCSC_RefGene_Name
  ) |>
  tidyr::separate_rows(UCSC_RefGene_Name, sep = ";") |>
  dplyr::filter(!is.na(UCSC_RefGene_Name), UCSC_RefGene_Name != "") |>
  dplyr::distinct(UCSC_RefGene_Name, cpg, .keep_all = TRUE) |>
  dplyr::group_by(gene = UCSC_RefGene_Name) |>
  dplyr::summarise(
    n_cpgs = dplyr::n_distinct(cpg),
    cpgs = paste(sort(unique(cpg)), collapse = ";"),
    best_nonzero_folds = max(best_nonzero_folds, na.rm = TRUE),
    max_mean_abs_coefficient = max(max_mean_abs_coefficient, na.rm = TRUE),
    evidence_runs = collapse_unique(evidence_runs),
    evidence_outcomes = collapse_unique(evidence_outcomes),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(n_cpgs),
    dplyr::desc(best_nonzero_folds),
    dplyr::desc(max_mean_abs_coefficient),
    gene
  )

long_path <- file.path(out_paths$results_features, "elastic_net_cpg_annotation_priority_long.csv")
unique_path <- file.path(out_paths$results_features, "elastic_net_cpg_annotation_priority_unique.csv")
gene_path <- file.path(out_paths$results_features, "elastic_net_cpg_annotation_gene_summary.csv")

readr::write_csv(annotated_long, long_path)
readr::write_csv(annotated_unique, unique_path)
readr::write_csv(gene_summary, gene_path)

message("Annotation complete.")
message("Stable evidence rows: ", nrow(annotated_long))
message("Unique CpGs: ", nrow(annotated_unique))
message("Genes with at least one annotated CpG: ", nrow(gene_summary))
message("Unannotated CpGs: ", sum(is.na(annotated_unique$chr)))
message("Saved long annotation: ", long_path)
message("Saved unique CpG annotation: ", unique_path)
message("Saved gene summary: ", gene_path)

print(
  annotated_unique |>
    dplyr::select(
      cpg,
      evidence_runs,
      best_nonzero_folds,
      max_mean_abs_coefficient,
      chr,
      pos,
      UCSC_RefGene_Name,
      Relation_to_Island,
      Regulatory_Feature_Group
    ) |>
    utils::head(20)
)
