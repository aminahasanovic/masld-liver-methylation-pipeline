suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tibble)
  library(tidyr)
})

source(file.path("data_preprocessing", "00_config.R"))
source(file.path(dir_helpers, "helpers_io.R"))
source(file.path(dir_helpers, "helpers_plotting.R"))

analysis_table_dir <- dir_out_epi_tables
supp_results_dir <- file.path(dir_results, "supplementary", "tables")
supp_thesis_out_dir <- file.path(project_root, "outputs", "thesis_outputs", "supplementary", "tables")
main_results_fig_dir <- file.path(dir_results, "figures")
main_thesis_out_fig_dir <- file.path(project_root, "outputs", "thesis_outputs", "figures")
main_results_table_dir <- file.path(dir_results, "tables")
main_thesis_out_table_dir <- file.path(project_root, "outputs", "thesis_outputs", "tables")
report_dir <- file.path(project_root, "outputs", "reports")

ensure_dirs(c(
  analysis_table_dir,
  supp_results_dir,
  supp_thesis_out_dir,
  dir_out_epi_plots,
  main_results_fig_dir,
  main_thesis_out_fig_dir,
  main_results_table_dir,
  main_thesis_out_table_dir,
  report_dir
))

correction_roles <- tibble(
  Correction = c("Raw betas", "ComBat betas"),
  matrix = c("Raw/pre-ComBat", "ComBat"),
  analysis_role = c("primary", "sensitivity")
)

methylcipher_version <- as.character(utils::packageVersion("methylCIPHER"))
ctsclocks_version <- as.character(utils::packageVersion("CTSclocks"))

read_required_csv <- function(path, label) {
  check_file_exists(path, label)
  readr::read_csv(path, show_col_types = FALSE)
}

write_table_all_locations <- function(x, filename, analysis_filename = filename) {
  write_csv_safe(x, file.path(analysis_table_dir, analysis_filename))
  write_csv_safe(x, file.path(supp_results_dir, filename))
  write_csv_safe(x, file.path(supp_thesis_out_dir, filename))
  invisible(filename)
}

update_supp_manifest <- function(manifest_path, rows) {
  rows <- rows |>
    mutate(across(everything(), as.character))

  existing <- if (file.exists(manifest_path)) {
    readr::read_csv(manifest_path, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
  } else {
    tibble()
  }

  if (nrow(existing) > 0) {
    missing_cols <- setdiff(names(rows), names(existing))
    for (col in missing_cols) {
      existing[[col]] <- NA_character_
    }
    existing <- existing |>
      select(all_of(names(rows))) |>
      filter(!.data$table %in% rows$table)
  }

  bind_rows(existing, rows) |>
    arrange(.data$table) |>
    readr::write_csv(manifest_path)
}

update_manifest_row <- function(manifest_path, key_col, row) {
  row <- tibble::as_tibble(row) |>
    mutate(across(everything(), as.character))

  existing <- if (file.exists(manifest_path)) {
    readr::read_csv(manifest_path, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
  } else {
    tibble()
  }

  if (nrow(existing) > 0) {
    missing_cols <- setdiff(names(row), names(existing))
    for (col in missing_cols) {
      existing[[col]] <- NA_character_
    }
    existing <- existing |>
      select(all_of(names(row))) |>
      filter(.data[[key_col]] != row[[key_col]][[1]])
  }

  bind_rows(existing, row) |>
    arrange(.data[[key_col]]) |>
    readr::write_csv(manifest_path)
}

save_main_figure_all_locations <- function(plot, file_stem, width, height, dpi = 450) {
  analysis_png <- file.path(dir_out_epi_plots, paste0(file_stem, ".png"))
  analysis_pdf <- file.path(dir_out_epi_plots, paste0(file_stem, ".pdf"))

  save_plot(plot, analysis_png, width = width, height = height, dpi = dpi)
  save_plot(plot, analysis_pdf, width = width, height = height, dpi = dpi, device = grDevices::cairo_pdf)

  targets <- c(main_results_fig_dir, main_thesis_out_fig_dir)
  for (target_dir in targets) {
    ensure_dir(target_dir)
    file.copy(analysis_png, file.path(target_dir, basename(analysis_png)), overwrite = TRUE)
    file.copy(analysis_pdf, file.path(target_dir, basename(analysis_pdf)), overwrite = TRUE)
  }

  invisible(c(analysis_png, analysis_pdf))
}

clock_table_path <- file.path(dir_out_epi_rds, "pheno_with_clocks_and_EAA_raw_combat.rds")
check_file_exists(clock_table_path, "Combined raw/ComBat EAA clock table")

clock_table <- readRDS(clock_table_path) |>
  as_tibble() |>
  mutate(
    Sample_Name = as.character(.data$Sample_Name),
    Dataset = as.character(.data$Dataset),
    DiseaseGroup = as.character(.data$DiseaseGroup),
    Progression3 = as.character(.data$Progression3),
    Sex = as.character(.data$Sex),
    Correction = as.character(.data$Correction)
  ) |>
  left_join(correction_roles, by = "Correction")

standardize_effects <- function(x) {
  x |>
    as_tibble() |>
    mutate(
      Correction = as.character(.data$Correction),
      StageVar = as.character(.data$StageVar),
      Clock = as.character(.data$Clock),
      Comparison = as.character(.data$Comparison)
    ) |>
    left_join(correction_roles, by = "Correction") |>
    mutate(
      reference_group = "Healthy",
      contrast = paste0(.data$Comparison, " vs ", .data$reference_group),
      beta = .data$estimate,
      ci_lower = .data$conf.low,
      ci_upper = .data$conf.high,
      p_value = .data[["p.value"]],
      std_error = .data[["std.error"]],
      model = paste0(.data$Clock, "_EAA_resid ~ ", .data$StageVar, " + Sex + Dataset"),
      covariates = "Sex + Dataset",
      residualization = "Global clock-age residuals within each matrix: lm(Clock ~ Age)"
    ) |>
    group_by(.data$Correction, .data$matrix, .data$analysis_role, .data$StageVar) |>
    mutate(
      q_value = stats::p.adjust(.data$p_value, method = "BH"),
      fdr_significant = .data$q_value < 0.05
    ) |>
    ungroup()
}

model_n_for <- function(correction_value, clock, stage_var) {
  outcome <- paste0(clock, "_EAA_resid")
  needed <- c(outcome, stage_var, "Sex", "Dataset")
  if (!all(needed %in% names(clock_table))) {
    return(NA_integer_)
  }
  d <- clock_table |>
    filter(.data$Correction == correction_value)
  sum(stats::complete.cases(d[, needed, drop = FALSE]))
}

effect_dg <- read_required_csv(
  file.path(analysis_table_dir, "effect_sizes_DiseaseGroup_vs_Healthy_RAW_vs_ComBat.csv"),
  "DiseaseGroup EAA effect-size table"
)

effect_prog3 <- read_required_csv(
  file.path(analysis_table_dir, "effect_sizes_Progression3_vs_Healthy_RAW_vs_ComBat.csv"),
  "Progression3 EAA effect-size table"
)

all_effects <- bind_rows(effect_dg, effect_prog3) |>
  standardize_effects() |>
  mutate(
    n_model_samples = mapply(model_n_for, .data$Correction, .data$Clock, .data$StageVar),
    Clock = factor(.data$Clock, levels = main_clock_cols),
    matrix = factor(.data$matrix, levels = c("Raw/pre-ComBat", "ComBat")),
    analysis_role = factor(.data$analysis_role, levels = c("primary", "sensitivity")),
    StageVar = factor(.data$StageVar, levels = c("DiseaseGroup", "Progression3"))
  ) |>
  arrange(.data$StageVar, .data$matrix, .data$Clock, .data$Comparison) |>
  mutate(
    Clock = as.character(.data$Clock),
    matrix = as.character(.data$matrix),
    analysis_role = as.character(.data$analysis_role),
    StageVar = as.character(.data$StageVar)
  ) |>
  select(
    analysis_role,
    matrix,
    Correction,
    StageVar,
    Clock,
    reference_group,
    Comparison,
    contrast,
    beta,
    ci_lower,
    ci_upper,
    p_value,
    q_value,
    fdr_significant,
    n_model_samples,
    std_error,
    statistic,
    term,
    model,
    covariates,
    residualization
  )

stopifnot(nrow(filter(all_effects, StageVar == "DiseaseGroup", matrix == "Raw/pre-ComBat")) == 60L)
stopifnot(nrow(filter(all_effects, StageVar == "DiseaseGroup", matrix == "ComBat")) == 60L)
stopifnot(nrow(filter(all_effects, StageVar == "Progression3", matrix == "Raw/pre-ComBat")) == 24L)
stopifnot(nrow(filter(all_effects, StageVar == "Progression3", matrix == "ComBat")) == 24L)

diseasegroup_effects <- all_effects |>
  filter(.data$StageVar == "DiseaseGroup")

multiple_testing_summary <- all_effects |>
  group_by(.data$StageVar, .data$analysis_role, .data$matrix, .data$Correction) |>
  summarise(
    n_tests = n(),
    n_clocks = n_distinct(.data$Clock),
    n_group_contrasts = n_distinct(.data$Comparison),
    min_nominal_p = min(.data$p_value, na.rm = TRUE),
    min_bh_fdr = min(.data$q_value, na.rm = TRUE),
    n_nominal_p_lt_0_05 = sum(.data$p_value < 0.05, na.rm = TRUE),
    n_bh_fdr_lt_0_05 = sum(.data$q_value < 0.05, na.rm = TRUE),
    adjustment = "Benjamini-Hochberg within each StageVar x matrix family",
    .groups = "drop"
  ) |>
  arrange(.data$StageVar, .data$analysis_role, .data$matrix)

coverage_files <- tibble(
  Correction = c("Raw betas", "ComBat betas"),
  path = file.path(
    analysis_table_dir,
    c("clock_probe_coverage_Raw_betas.csv", "clock_probe_coverage_ComBat_betas.csv")
  )
)

coverage <- bind_rows(lapply(seq_len(nrow(coverage_files)), function(i) {
  read_required_csv(coverage_files$path[[i]], basename(coverage_files$path[[i]])) |>
    mutate(Correction = coverage_files$Correction[[i]])
})) |>
  left_join(correction_roles, by = "Correction") |>
  mutate(
    ClockName = as.character(.data$ClockName),
    missing_probe_count = .data$TotalProbes - .data$ProbesAvailable
  )

pc_clock_names <- c("PCHannum", "PCHorvath1", "PCHorvath2", "PCPhenoAge", "PCGrimAge")

clock_map <- tibble(
  Clock = c(main_clock_cols, "GrimAgeV2", "HRSInChPhenoAge"),
  coverage_source = c(
    "Horvath1",
    "Horvath2",
    "Hannum",
    "PhenoAge",
    "GrimAge1",
    rep("PCClocks", length(pc_clock_names)),
    "HepClock",
    "LiverClock",
    "GrimAge2",
    "HRSInCHPhenoAge"
  )
) |>
  mutate(
    included_in_inference = .data$Clock %in% main_clock_cols,
    missing_cpg_strategy = case_when(
      .data$Clock %in% c("Horvath1", "Horvath2", "Hannum", "PhenoAge", "HRSInChPhenoAge") ~
        "available_cpg_weighted_sum_no_imputation",
      .data$Clock %in% c("GrimAgeV1", "GrimAgeV2") ~
        "zero_fill_absent_cpgs",
      .data$Clock %in% pc_clock_names ~
        "reference_mean_imputation",
      .data$Clock %in% c("HepClock", "LiverClock") ~
        "available_cpg_glmnet_prediction_no_imputation",
      TRUE ~ NA_character_
    ),
    clock_family = case_when(
      .data$Clock %in% pc_clock_names ~ "PC clock",
      .data$Clock %in% c("HepClock", "LiverClock") ~ "CTS tissue-specific clock",
      TRUE ~ "methylCIPHER clock"
    ),
    clock_function = case_when(
      .data$Clock == "Horvath1" ~ "methylCIPHER::calcHorvath1(imputation = FALSE)",
      .data$Clock == "Horvath2" ~ "methylCIPHER::calcHorvath2(imputation = FALSE)",
      .data$Clock == "Hannum" ~ "methylCIPHER::calcHannum(imputation = FALSE)",
      .data$Clock == "PhenoAge" ~ "methylCIPHER::calcPhenoAge(imputation = FALSE)",
      .data$Clock == "GrimAgeV1" ~ "methylCIPHER::calcGrimAgeV1()",
      .data$Clock == "GrimAgeV2" ~ "methylCIPHER::calcGrimAgeV2()",
      .data$Clock %in% pc_clock_names ~ "local calc_pcclocks_zenodo(); methylCIPHER::impute_DNAm(method = 'mean')",
      .data$Clock == "HepClock" ~ "CTSclocks::CTSclockAge(CTSclock = 'Hep')",
      .data$Clock == "LiverClock" ~ "CTSclocks::CTSclockAge(CTSclock = 'Liver')",
      .data$Clock == "HRSInChPhenoAge" ~ "methylCIPHER::calcHRSInChPhenoAge(imputation = FALSE)",
      TRUE ~ NA_character_
    ),
    absent_cpg_action = case_when(
      .data$missing_cpg_strategy == "available_cpg_weighted_sum_no_imputation" ~
        "CpGs absent from the common beta matrix are not imputed. The methylCIPHER function matches available CpGs only, applies the corresponding published coefficients, and keeps the original intercept/age transform.",
      .data$missing_cpg_strategy == "zero_fill_absent_cpgs" ~
        "CpGs absent from the common beta matrix are added as 0-valued beta columns via methylCIPHER::impute_DNAm(CpGs = zero_cpgs(...)) before GrimAge submodel scoring.",
      .data$missing_cpg_strategy == "reference_mean_imputation" ~
        "CpGs absent from the common beta matrix are added from PCClocks_data.qs2$imputeMissingCpGs and the matrix is then subset/reordered to the full PCClocks reference CpG vector.",
      .data$missing_cpg_strategy == "available_cpg_glmnet_prediction_no_imputation" ~
        "CTSclockAge matches clock CpGs to rownames(data.m), drops unmatched CpG rows and glmnet coefficients, updates the model dimension, and predicts on the available CpGs only.",
      TRUE ~ NA_character_
    ),
    present_na_action = case_when(
      .data$missing_cpg_strategy %in% c("reference_mean_imputation", "zero_fill_absent_cpgs") ~
        "methylCIPHER::impute_DNAm would mean-impute NA values in present CpGs; the pipeline input matrices contained no NAs after preprocessing.",
      .data$missing_cpg_strategy == "available_cpg_weighted_sum_no_imputation" ~
        "methylCIPHER rowSums(..., na.rm = TRUE) would omit NA contributions in present CpGs; the pipeline input matrices contained no NAs after preprocessing.",
      .data$missing_cpg_strategy == "available_cpg_glmnet_prediction_no_imputation" ~
        "CTSclockAge does not perform explicit NA imputation for Hep/Liver clocks; the pipeline input matrices contained no NAs after preprocessing.",
      TRUE ~ NA_character_
    ),
    score_implication = case_when(
      .data$missing_cpg_strategy == "available_cpg_weighted_sum_no_imputation" ~
        "Clock ages are available-CpG implementations of the published weighted formula, without refitting, reweighting, or reference-mean imputation of absent CpGs.",
      .data$missing_cpg_strategy == "zero_fill_absent_cpgs" ~
        "GrimAge submodels are scored after zero-filling CpGs not present in the common matrix.",
      .data$missing_cpg_strategy == "reference_mean_imputation" ~
        "PC-clock scores are computed on the full PCClocks reference CpG vector after reference-value imputation of absent CpGs.",
      .data$missing_cpg_strategy == "available_cpg_glmnet_prediction_no_imputation" ~
        "CTS-clock predictions use a reduced glmnet coefficient/input matrix over available CpGs, with the retained model intercept.",
      TRUE ~ NA_character_
    ),
    implementation_evidence = case_when(
      .data$missing_cpg_strategy == "available_cpg_weighted_sum_no_imputation" ~
        paste0("Local methylCIPHER ", methylcipher_version, " code path: present <- clock_CpGs %in% colnames(DNAm); betas <- DNAm[, na.omit(match(...))]; coefficients[present]; rowSums(..., na.rm = TRUE)."),
      .data$missing_cpg_strategy == "zero_fill_absent_cpgs" ~
        paste0("Local methylCIPHER ", methylcipher_version, " code path: calcGrimAgeV1/V2 calls impute_DNAm(method = 'mean', CpGs = zero_cpgs(CpGs_GrimAge1/2), subset = TRUE)."),
      .data$missing_cpg_strategy == "reference_mean_imputation" ~
        "Project code path: calc_pcclocks_zenodo calls methylCIPHER::impute_DNAm(method = 'mean', CpGs = RData$imputeMissingCpGs, subset = TRUE).",
      .data$missing_cpg_strategy == "available_cpg_glmnet_prediction_no_imputation" ~
        paste0("Local CTSclocks ", ctsclocks_version, " code path: idx <- match(ClockCpGs.v, rownames(data.m)); beta <- beta[!is.na(idx), ]; data.m <- data.m[na.omit(idx), ]; predict.glmnet(...)."),
      TRUE ~ NA_character_
    ),
    clock_handling = case_when(
      .data$Clock %in% pc_clock_names ~ "Calculated with local PCClocks_data.qs2 wrapper after methylCIPHER reference-mean imputation to the PCClocks reference CpG set.",
      .data$Clock %in% c("HepClock", "LiverClock") ~ "Calculated with CTSclocks::CTSclockAge on the common-CpG beta matrix; absent CpGs were dropped internally.",
      .data$Clock %in% c("GrimAgeV1", "GrimAgeV2") ~ "Calculated with methylCIPHER GrimAge function; absent GrimAge CpGs were added as 0-valued beta columns internally.",
      .data$Clock %in% c("GrimAgeV2", "HRSInChPhenoAge") ~ "Clock value calculated but residual EAA was not included in the curated inferential clock panel.",
      TRUE ~ "Calculated with methylCIPHER clock function using imputation = FALSE; absent CpGs were omitted from the weighted sum."
    ),
    exclusion_reason = case_when(
      .data$Clock == "GrimAgeV2" ~ "Calculated for audit context but excluded from curated EAA inference to keep the main clock panel fixed to the predefined 12 clocks.",
      .data$Clock == "HRSInChPhenoAge" ~ "Calculated for audit context but excluded from curated EAA inference because the installed package did not provide a resolved bibliographic source for citation.",
      TRUE ~ NA_character_
    )
  )

nonmissing_count <- function(correction_value, column_name) {
  if (!column_name %in% names(clock_table)) {
    return(0L)
  }
  sum(!is.na(clock_table[[column_name]][clock_table$Correction == correction_value]))
}

clock_audit <- tidyr::crossing(
  clock_map,
  correction_roles
) |>
  left_join(
    coverage |>
      select(
        Correction,
        ClockName,
        TotalProbes,
        ProbesAvailable,
        PercentProbesAvailable,
        missing_probe_count
      ),
    by = c("Correction", "coverage_source" = "ClockName")
  ) |>
  mutate(
    clock_value_nonmissing_n = mapply(nonmissing_count, .data$Correction, .data$Clock),
    eaa_resid_column = paste0(.data$Clock, "_EAA_resid"),
    eaa_resid_nonmissing_n = mapply(nonmissing_count, .data$Correction, .data$eaa_resid_column),
    matrix_probe_set_note = "Raw/pre-ComBat and ComBat matrices share the same common-CpG row set, so coverage counts should match across matrices."
  ) |>
  arrange(desc(.data$included_in_inference), .data$Clock, .data$analysis_role) |>
  select(
    analysis_role,
    matrix,
    Correction,
    Clock,
    coverage_source,
    clock_family,
    included_in_inference,
    clock_function,
    TotalProbes,
    ProbesAvailable,
    missing_probe_count,
    PercentProbesAvailable,
    clock_value_nonmissing_n,
    eaa_resid_nonmissing_n,
    missing_cpg_strategy,
    absent_cpg_action,
    present_na_action,
    score_implication,
    implementation_evidence,
    clock_handling,
    exclusion_reason,
    matrix_probe_set_note
  )

extract_johnson_term <- function(fit) {
  broom::tidy(fit, conf.int = TRUE) |>
    filter(.data$term == "FibrosisGroupAdvanced_Fibrosis") |>
    transmute(
      beta = .data$estimate,
      ci_lower = .data$conf.low,
      ci_upper = .data$conf.high,
      p_value = .data$p.value,
      std_error = .data$std.error,
      statistic = .data$statistic
    )
}

fit_johnson_residualization <- function(clock, correction_value, residualization_scope) {
  meta <- correction_roles |>
    filter(.data$Correction == correction_value)

  clock_col <- clock
  global_resid_col <- paste0(clock, "_EAA_resid")

  d <- clock_table |>
    filter(
      .data$Correction == correction_value,
      .data$Dataset == "Johnson",
      .data$DiseaseGroup %in% c("Healthy_Obese", "Advanced_Fibrosis")
    ) |>
    mutate(
      FibrosisGroup = factor(.data$DiseaseGroup, levels = c("Healthy_Obese", "Advanced_Fibrosis")),
      Sex = factor(.data$Sex)
    )

  if (identical(residualization_scope, "global_age_residuals")) {
    if (!global_resid_col %in% names(d)) {
      return(tibble())
    }
    d <- d |>
      mutate(rEAA = .data[[global_resid_col]])
    n_age_residualization <- sum(!is.na(clock_table[[global_resid_col]][clock_table$Correction == correction_value]))
  } else {
    if (!clock_col %in% names(d)) {
      return(tibble())
    }
    resid_data <- d |>
      drop_na(all_of(c(clock_col, "Age")))
    if (nrow(resid_data) < 10 || sd(resid_data[[clock_col]], na.rm = TRUE) == 0) {
      return(tibble())
    }
    age_fit <- lm(stats::reformulate("Age", response = clock_col), data = resid_data)
    d$rEAA <- NA_real_
    d$rEAA[match(resid_data$Sample_Name, d$Sample_Name)] <- resid(age_fit)
    n_age_residualization <- nrow(resid_data)
  }

  model_data <- d |>
    drop_na(all_of(c("rEAA", "FibrosisGroup", "Sex")))

  if (nrow(model_data) < 10 || n_distinct(model_data$FibrosisGroup) < 2) {
    return(tibble())
  }

  fit <- lm(rEAA ~ FibrosisGroup + Sex, data = model_data)
  extract_johnson_term(fit) |>
    mutate(
      Clock = clock,
      Correction = correction_value,
      matrix = meta$matrix,
      analysis_role = meta$analysis_role,
      residualization_scope = residualization_scope,
      comparison = "Advanced_Fibrosis vs Healthy_Obese",
      n_age_residualization = n_age_residualization,
      n_model_samples = nrow(model_data),
      model = "rEAA ~ FibrosisGroup + Sex",
      covariates = "Sex"
    )
}

johnson_residualization_audit <- bind_rows(lapply(correction_roles$Correction, function(correction_value) {
  bind_rows(lapply(main_clock_cols, function(clock) {
    bind_rows(
      fit_johnson_residualization(clock, correction_value, "global_age_residuals"),
      fit_johnson_residualization(clock, correction_value, "johnson_only_age_residuals")
    )
  }))
})) |>
  group_by(.data$matrix, .data$analysis_role, .data$residualization_scope) |>
  mutate(
    q_value = stats::p.adjust(.data$p_value, method = "BH"),
    fdr_significant = .data$q_value < 0.05
  ) |>
  ungroup() |>
  mutate(
    Clock = factor(.data$Clock, levels = main_clock_cols),
    matrix = factor(.data$matrix, levels = c("Raw/pre-ComBat", "ComBat")),
    analysis_role = factor(.data$analysis_role, levels = c("primary", "sensitivity")),
    residualization_scope = factor(
      .data$residualization_scope,
      levels = c("global_age_residuals", "johnson_only_age_residuals")
    )
  ) |>
  arrange(.data$matrix, .data$Clock, .data$residualization_scope) |>
  mutate(
    Clock = as.character(.data$Clock),
    matrix = as.character(.data$matrix),
    analysis_role = as.character(.data$analysis_role),
    residualization_scope = as.character(.data$residualization_scope)
  ) |>
  select(
    analysis_role,
    matrix,
    Correction,
    Clock,
    residualization_scope,
    comparison,
    beta,
    ci_lower,
    ci_upper,
    p_value,
    q_value,
    fdr_significant,
    n_age_residualization,
    n_model_samples,
    std_error,
    statistic,
    model,
    covariates
  )

stopifnot(nrow(johnson_residualization_audit) == 48L)

write_table_all_locations(
  all_effects,
  "T29_eaa_effects_all_contrasts_with_bh_fdr.csv",
  "eaa_effects_all_contrasts_with_bh_fdr_primary_sensitivity.csv"
)
write_table_all_locations(
  multiple_testing_summary,
  "T30_eaa_multiple_testing_summary.csv",
  "eaa_multiple_testing_summary_primary_sensitivity.csv"
)
write_table_all_locations(
  diseasegroup_effects,
  "T39_eaa_diseasegroup_primary_sensitivity_contrasts.csv",
  "eaa_diseasegroup_primary_sensitivity_contrasts.csv"
)
write_table_all_locations(
  clock_audit,
  "T40_eaa_clock_probe_handling_audit.csv",
  "eaa_clock_probe_handling_audit.csv"
)
write_table_all_locations(
  johnson_residualization_audit,
  "T41_johnson_eaa_residualization_audit.csv",
  "johnson_eaa_residualization_audit.csv"
)

make_f22_precombat_plot <- function(effects) {
  disease_order <- c("Healthy_Obese", "MASL", "MASH", "Mild_Fibrosis", "Advanced_Fibrosis")
  disease_labels <- c(
    Healthy_Obese = "Healthy obese",
    MASL = "MASL",
    MASH = "MASH",
    Mild_Fibrosis = "Mild fibrosis",
    Advanced_Fibrosis = "Advanced fibrosis"
  )
  disease_palette <- setNames(stage_colors[disease_order], disease_labels[disease_order])
  forest_clock_order <- c(
    "Horvath1", "PCHorvath1",
    "Horvath2", "PCHorvath2",
    "Hannum", "PCHannum",
    "PhenoAge", "PCPhenoAge",
    "GrimAgeV1", "PCGrimAge",
    "HepClock", "LiverClock"
  )
  forest_clock_labels <- c(
    Horvath1 = "Horvath1",
    PCHorvath1 = "PC-Horvath1",
    Horvath2 = "Horvath2",
    PCHorvath2 = "PC-Horvath2",
    Hannum = "Hannum",
    PCHannum = "PC-Hannum",
    PhenoAge = "PhenoAge",
    PCPhenoAge = "PC-PhenoAge",
    GrimAgeV1 = "GrimAgeV1",
    PCGrimAge = "PC-GrimAge",
    HepClock = "HepClock",
    LiverClock = "LiverClock"
  )

  df_sub <- effects |>
    filter(
      .data$StageVar == "DiseaseGroup",
      .data$analysis_role == "primary",
      .data$matrix == "Raw/pre-ComBat"
    ) |>
    mutate(
      Clock_chr = as.character(.data$Clock),
      Comparison_chr = as.character(.data$Comparison)
    )

  stopifnot(nrow(df_sub) == 60L)

  clock_levels <- forest_clock_order[forest_clock_order %in% unique(df_sub$Clock_chr)]
  forest_rows <- tidyr::expand_grid(
    Clock_chr = clock_levels,
    Comparison_chr = disease_order
  ) |>
    mutate(
      y = dplyr::n() - dplyr::row_number() + 1,
      Clock_label = dplyr::coalesce(unname(forest_clock_labels[Clock_chr]), Clock_chr),
      Comparison_label = unname(disease_labels[Comparison_chr])
    )

  plot_df <- df_sub |>
    select(Clock_chr, Comparison_chr, beta, ci_lower, ci_upper) |>
    right_join(forest_rows, by = c("Clock_chr", "Comparison_chr")) |>
    mutate(
      Comparison_label = factor(.data$Comparison_label, levels = unname(disease_labels[disease_order])),
      estimate_label = if_else(is.na(.data$beta), NA_character_, sprintf("%+.2f", .data$beta)),
      estimate_label_x = .data$beta + if_else(.data$beta >= 0, 0.25, -0.25),
      estimate_label_hjust = if_else(.data$beta >= 0, 0, 1)
    )

  clock_groups <- forest_rows |>
    group_by(.data$Clock_chr, .data$Clock_label) |>
    summarise(
      ymin = min(.data$y) - 0.42,
      ymax = max(.data$y) + 0.42,
      ymid = mean(.data$y),
      .groups = "drop"
    )

  x_vals <- c(plot_df$ci_lower, plot_df$ci_upper, plot_df$beta)
  x_vals <- x_vals[is.finite(x_vals)]
  if (length(x_vals) == 0) x_vals <- c(-1, 1)
  x_left <- floor(min(c(x_vals, 0), na.rm = TRUE) / 5) * 5
  x_right <- ceiling(max(c(x_vals, 0), na.rm = TRUE) / 5) * 5
  if (identical(x_left, x_right)) {
    x_left <- x_left - 1
    x_right <- x_right + 1
  }
  x_span <- x_right - x_left
  x_right <- x_right + max(0.8, x_span * 0.05)
  stage_label_x <- x_left - x_span * 0.012
  bracket_x <- x_left - x_span * 0.21
  bracket_tick_x <- x_left - x_span * 0.185
  clock_label_x <- x_left - x_span * 0.32

  forest_plot <- ggplot(plot_df, aes(x = .data$beta, y = .data$y, color = .data$Comparison_label)) +
    geom_hline(
      data = forest_rows,
      aes(yintercept = .data$y),
      inherit.aes = FALSE,
      color = "#E8EDF2",
      linewidth = 0.18
    ) +
    geom_vline(xintercept = 0, color = "#315C8A", linewidth = 0.42) +
    geom_errorbar(
      aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
      orientation = "y",
      width = 0.20,
      linewidth = 0.42,
      na.rm = TRUE
    ) +
    geom_point(size = 1.25, na.rm = TRUE) +
    geom_text(
      aes(
        x = .data$estimate_label_x,
        y = .data$y + 0.27,
        label = .data$estimate_label,
        hjust = .data$estimate_label_hjust
      ),
      color = "#111827",
      size = 1.55,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    geom_text(
      data = forest_rows,
      aes(y = .data$y, label = .data$Comparison_label),
      inherit.aes = FALSE,
      x = stage_label_x,
      hjust = 1,
      size = 1.62,
      color = "#40516A"
    ) +
    geom_segment(
      data = clock_groups,
      aes(y = .data$ymin, yend = .data$ymax),
      inherit.aes = FALSE,
      x = bracket_x,
      xend = bracket_x,
      color = "#5F6F82",
      linewidth = 0.28
    ) +
    geom_segment(
      data = clock_groups,
      aes(y = .data$ymin, yend = .data$ymin),
      inherit.aes = FALSE,
      x = bracket_x,
      xend = bracket_tick_x,
      color = "#5F6F82",
      linewidth = 0.28
    ) +
    geom_segment(
      data = clock_groups,
      aes(y = .data$ymax, yend = .data$ymax),
      inherit.aes = FALSE,
      x = bracket_x,
      xend = bracket_tick_x,
      color = "#5F6F82",
      linewidth = 0.28
    ) +
    geom_text(
      data = clock_groups,
      aes(y = .data$ymid, label = .data$Clock_label),
      inherit.aes = FALSE,
      x = clock_label_x,
      hjust = 1,
      size = 1.82,
      fontface = "bold",
      color = "#2F3B52"
    ) +
    scale_color_manual(values = disease_palette, drop = FALSE) +
    scale_y_continuous(breaks = NULL, expand = expansion(add = c(0.55, 0.55))) +
    scale_x_continuous(breaks = seq(x_left, floor(x_right / 5) * 5, 5), expand = expansion(mult = c(0, 0))) +
    coord_cartesian(xlim = c(x_left, x_right), clip = "off") +
    theme_pipeline(7.4) +
    labs(
      x = "Adjusted residual EAA effect vs Healthy (years)",
      y = NULL,
      color = "Disease group"
    ) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(face = "plain", hjust = 0.36, size = 7.8, color = "#17233C"),
      plot.title.position = "plot",
      plot.title = element_text(size = 9.8, face = "bold", hjust = 0.5, margin = margin(l = -80)),
      plot.subtitle = element_text(size = 7.4, hjust = 0.5, margin = margin(l = -80)),
      legend.title = element_text(size = 8.8, face = "bold"),
      legend.text = element_text(size = 8.1),
      legend.key.size = grid::unit(0.46, "cm"),
      legend.justification = "center",
      legend.box.just = "center",
      legend.box.margin = margin(t = 0, r = 190, b = 0, l = 0),
      legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
      plot.margin = margin(7, 8, 7, 170)
    ) +
    guides(color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(size = 2.6, linewidth = 0.8)
    ))

  title_plot <- ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.70,
      label = "DiseaseGroup effect on residual EAA - Primary pre-ComBat",
      fontface = "bold",
      size = 3.4,
      color = "#17233C"
    ) +
    annotate(
      "text",
      x = 0.5,
      y = 0.18,
      label = "Reference = Healthy; adjusted for Sex + Dataset",
      size = 2.7,
      color = "#4B5563"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(2, 0, 0, 0))

  patchwork::wrap_elements(full = title_plot) / forest_plot + patchwork::plot_layout(heights = c(0.065, 1))
}

f22_stem <- "F22_beta_CI_diseasegroup_precombat_all_clocks"
f22_plot <- make_f22_precombat_plot(diseasegroup_effects)
save_main_figure_all_locations(f22_plot, f22_stem, width = 7.2, height = 9.4)

f22_manifest_row <- tibble(
  figure = "F22",
  file_stem = f22_stem,
  main_message = "Adjusted DiseaseGroup effects on residual EAA in the primary Raw/pre-ComBat analysis."
)

update_manifest_row(file.path(main_results_table_dir, "figure_manifest.csv"), "figure", f22_manifest_row)
update_manifest_row(file.path(main_thesis_out_table_dir, "thesis_figure_manifest.csv"), "figure", f22_manifest_row)

manifest_rows <- tibble::tribble(
  ~table, ~file, ~main_message,
  "T29", "T29_eaa_effects_all_contrasts_with_bh_fdr.csv", "Complete EAA effect table with Raw/pre-ComBat primary and ComBat sensitivity group contrasts.",
  "T30", "T30_eaa_multiple_testing_summary.csv", "EAA multiple-testing summary by outcome family and matrix.",
  "T37", "T37_bmi_availability_audit_eaa.csv", "BMI availability in the final EAA-eligible cohort after technical replicate resolution.",
  "T38", "T38_johnson_bmi_eaa_sensitivity.csv", "Johnson-only BMI sensitivity for residual EAA advanced-fibrosis effects in Raw/pre-ComBat primary and ComBat sensitivity matrices.",
  "T39", "T39_eaa_diseasegroup_primary_sensitivity_contrasts.csv", "Focused 60-per-matrix DiseaseGroup EAA contrasts with beta, confidence interval, p-value, and q-value.",
  "T40", "T40_eaa_clock_probe_handling_audit.csv", "Clock CpG coverage and handling audit for the inferential and calculated-but-excluded EAA clocks.",
  "T41", "T41_johnson_eaa_residualization_audit.csv", "Johnson-only audit comparing global age residuals with Johnson-only age residuals."
)

update_supp_manifest(file.path(supp_results_dir, "supplementary_table_manifest.csv"), manifest_rows)
update_supp_manifest(file.path(supp_thesis_out_dir, "supplementary_table_manifest.csv"), manifest_rows)

decision_counts <- multiple_testing_summary |>
  select(StageVar, analysis_role, matrix, n_tests, n_bh_fdr_lt_0_05)

coverage_selected <- clock_audit |>
  filter(.data$included_in_inference, .data$analysis_role == "primary") |>
  transmute(
    Clock,
    coverage = sprintf("%.1f%%", .data$PercentProbesAvailable),
    missing_probe_count
  )

handling_summary <- clock_audit |>
  filter(.data$analysis_role == "primary") |>
  distinct(.data$Clock, .data$included_in_inference, .data$clock_family, .data$missing_cpg_strategy) |>
  count(.data$included_in_inference, .data$clock_family, .data$missing_cpg_strategy, name = "n_clocks") |>
  arrange(desc(.data$included_in_inference), .data$clock_family, .data$missing_cpg_strategy)

johnson_scope_summary <- johnson_residualization_audit |>
  group_by(.data$matrix, .data$analysis_role, .data$residualization_scope) |>
  summarise(
    n_tests = n(),
    n_fdr_lt_0_05 = sum(.data$q_value < 0.05, na.rm = TRUE),
    min_q = min(.data$q_value, na.rm = TRUE),
    .groups = "drop"
  )

bmi_path <- file.path(supp_thesis_out_dir, "T38_johnson_bmi_eaa_sensitivity.csv")
bmi_summary <- if (file.exists(bmi_path)) {
  readr::read_csv(bmi_path, show_col_types = FALSE) |>
    group_by(.data$matrix, .data$analysis_role) |>
    summarise(
      n_tests = n(),
      n_direction_concordant = sum(.data$direction_concordant, na.rm = TRUE),
      n_fdr_status_changed = sum(.data$FDR_status_change != "unchanged", na.rm = TRUE),
      max_abs_relative_beta_change_percent = max(abs(.data$relative_beta_change_percent), na.rm = TRUE),
      .groups = "drop"
    )
} else {
  tibble()
}

report_table <- function(x) {
  old_width <- getOption("width")
  on.exit(options(width = old_width), add = TRUE)
  options(width = 240)

  paste(
    capture.output(print(as.data.frame(x), row.names = FALSE, right = FALSE)),
    collapse = "\n"
  )
}

report_path <- file.path(report_dir, "eaa_audit_primary_raw_sensitivity_combat.md")
report_lines <- c(
  "# EAA Audit: Primary Raw/pre-ComBat, ComBat Sensitivity",
  "",
  "## Decision",
  "",
  "Recommendation: present Raw/pre-ComBat beta-derived EAA models with Dataset adjustment as the primary inferential analysis, and present ComBat-derived EAA models as sensitivity analysis.",
  "",
  "Rationale: the Raw/pre-ComBat matrix preserves the original within- and across-study methylation structure, while the linear EAA model explicitly adjusts for Dataset. ComBat remains useful as a sensitivity analysis because it tests whether the inference is robust after empirical batch harmonization.",
  "",
  "## Model",
  "",
  "- Residualization: global clock-age residuals are calculated within each matrix using `lm(Clock ~ Age)`.",
  "- Inference: `Clock_EAA_resid ~ DiseaseGroup + Sex + Dataset` and `Clock_EAA_resid ~ Progression3 + Sex + Dataset`.",
  "- Reference group: Healthy.",
  "- Multiple testing: Benjamini-Hochberg within each StageVar x matrix family.",
  "",
  "## Contrast Counts",
  "",
  "```",
  report_table(decision_counts),
  "```",
  "",
  "## Primary Clock Coverage",
  "",
  "```",
  report_table(coverage_selected),
  "```",
  "",
  "## Missing Clock-CpG Handling",
  "",
  "Only the 12 predefined clocks were included in EAA inference. GrimAgeV2 and HRSInChPhenoAge are retained in T40 only as audit-only calculated/excluded rows.",
  "",
  "- methylCIPHER Horvath1, Horvath2, Hannum, and PhenoAge were called with `imputation = FALSE`; absent CpGs were not imputed, but omitted from the weighted sum with the original intercept/transform retained.",
  "- methylCIPHER GrimAgeV1 filled absent GrimAge CpGs as 0-valued beta columns internally via `impute_DNAm(CpGs = zero_cpgs(...))` before submodel scoring.",
  "- PC clocks were calculated by the local `calc_pcclocks_zenodo()` wrapper; absent CpGs were filled from `PCClocks_data.qs2$imputeMissingCpGs` and the matrix was then ordered to the full PCClocks reference CpG vector.",
  "- CTSclocks HepClock and LiverClock used `CTSclockAge`; absent CpGs were dropped from both the glmnet coefficient vector and the input matrix before prediction, without explicit CpG imputation.",
  "- The input Raw/pre-ComBat and ComBat matrices had no remaining NA beta values, so these rules concern CpGs absent from the 450K/EPIC common-CpG matrix rather than sample-level missing beta values.",
  "",
  "```",
  report_table(handling_summary),
  "```",
  "",
  "## Johnson Residualization Check",
  "",
  "```",
  report_table(johnson_scope_summary),
  "```",
  "",
  "## BMI Sensitivity Check",
  "",
  "```",
  if (nrow(bmi_summary) > 0) report_table(bmi_summary) else "BMI sensitivity table was not present when this report was written.",
  "```",
  "",
  "## Main Output Tables",
  "",
  "- `outputs/thesis_outputs/supplementary/tables/T29_eaa_effects_all_contrasts_with_bh_fdr.csv`",
  "- `outputs/thesis_outputs/supplementary/tables/T39_eaa_diseasegroup_primary_sensitivity_contrasts.csv`",
  "- `outputs/thesis_outputs/supplementary/tables/T40_eaa_clock_probe_handling_audit.csv`",
  "- `outputs/thesis_outputs/supplementary/tables/T41_johnson_eaa_residualization_audit.csv`",
  "- `outputs/thesis_outputs/supplementary/tables/T37_bmi_availability_audit_eaa.csv`",
  "- `outputs/thesis_outputs/supplementary/tables/T38_johnson_bmi_eaa_sensitivity.csv`"
)

writeLines(report_lines, report_path)

message("EAA audit outputs written.")
message("Report: ", report_path)
