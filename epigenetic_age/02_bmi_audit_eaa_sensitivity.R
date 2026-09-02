suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

source(file.path("data_preprocessing", "00_config.R"))
source(file.path(dir_helpers, "helpers_io.R"))

eaa_pheno_file <- file.path(dir_out_epi_rds, "clock_cohort_pheno.rds")
eaa_clock_file <- file.path(dir_out_epi_rds, "pheno_with_clocks_and_EAA_raw_combat.rds")

check_file_exists(eaa_pheno_file, "EAA clock cohort phenotype table")
check_file_exists(eaa_clock_file, "EAA clock result table")

dataset_order <- c("Ahrens", "Horvath", "Johnson", "ITEN", "VanDijck", "Murphy")

pheno_files <- c(
  Ahrens = file.path(dir_processed, "Ahrens_GSE48325", "pheno_ahrens_clean.rds"),
  Horvath = file.path(dir_processed, "Horvath_GSE61258", "pheno_horvath_clean.rds"),
  Johnson = file.path(dir_processed, "Johnson_GSE180474", "pheno_johnson_clean.rds"),
  ITEN = file.path(dir_processed, "Kim_KSE102917", "pheno_kim_clean.rds"),
  VanDijck = file.path(dir_processed, "VanDijck_GSE294806", "pheno_vandijck_clean.rds"),
  Murphy = file.path(dir_processed, "Murphy_GSE49542", "pheno_murphy_clean.rds")
)

for (f in pheno_files) check_file_exists(f, basename(f))

ahrens_bmi_file <- file.path(
  dir_processed,
  "Ahrens_GSE48325",
  "ahrens_bmi_from_series_matrix.csv"
)
check_file_exists(ahrens_bmi_file, "Ahrens BMI extracted from GSE48325 series matrix")

bmi_source_variable <- c(
  Ahrens = "BMI (GSE48325 series_matrix !Sample_characteristics_ch1:bmi)",
  Horvath = "BMI (GEO characteristics_ch1:bmi)",
  Johnson = "BMI (GEO characteristics_ch1:bmi; retained after Individual_ID technical replicate resolution)",
  ITEN = "Not available in ITEN BioSample workbook, BioSample-IDAT map, KSE102917 IDF, or processed ITEN phenotype",
  VanDijck = "BMI (GEO characteristics_ch1:bmi)",
  Murphy = "Not available; Murphy excluded from EAA because age/sex are missing"
)

read_bmi_lookup <- function(dataset, path) {
  x <- readRDS(path) |>
    as_tibble() |>
    mutate(Sample_Name = as.character(.data$Sample_Name))

  if (identical(dataset, "Ahrens")) {
    ahrens_bmi <- readr::read_csv(ahrens_bmi_file, show_col_types = FALSE) |>
      transmute(
        Sample_Name = as.character(.data$Sample_Name),
        BMI = suppressWarnings(as.numeric(.data$BMI))
      )

    if (anyDuplicated(ahrens_bmi$Sample_Name) || anyNA(ahrens_bmi$BMI)) {
      stop("Ahrens BMI lookup has duplicated sample IDs or missing BMI values.", call. = FALSE)
    }

    x <- x |>
      select(-any_of("BMI")) |>
      left_join(ahrens_bmi, by = "Sample_Name")
  } else if ("BMI" %in% names(x)) {
    x <- x |>
      mutate(BMI = suppressWarnings(as.numeric(.data$BMI)))
  } else {
    x <- x |>
      mutate(BMI = NA_real_)
  }

  tibble(
    Sample_Name = x$Sample_Name,
    Dataset = dataset,
    BMI = x$BMI,
    BMI_source_variable = unname(bmi_source_variable[dataset])
  )
}

bmi_lookup <- bind_rows(Map(read_bmi_lookup, names(pheno_files), pheno_files))

clock_pheno <- readRDS(eaa_pheno_file) |>
  as_tibble() |>
  mutate(
    Sample_Name = as.character(.data$Sample_Name),
    Dataset = as.character(.data$Dataset)
  )

audit <- tibble(Dataset = dataset_order) |>
  left_join(
    clock_pheno |>
      count(Dataset, name = "final_EAA_eligible_n"),
    by = "Dataset"
  ) |>
  mutate(final_EAA_eligible_n = tidyr::replace_na(final_EAA_eligible_n, 0L)) |>
  left_join(
    clock_pheno |>
      left_join(bmi_lookup, by = c("Sample_Name", "Dataset")) |>
      group_by(Dataset) |>
      summarise(BMI_available_n = sum(!is.na(BMI)), .groups = "drop"),
    by = "Dataset"
  ) |>
  mutate(
    BMI_available_n = tidyr::replace_na(BMI_available_n, 0L),
    BMI_missing_n = final_EAA_eligible_n - BMI_available_n,
    BMI_source_variable = unname(bmi_source_variable[Dataset])
  ) |>
  select(Dataset, final_EAA_eligible_n, BMI_available_n, BMI_missing_n, BMI_source_variable)

correction_roles <- tibble::tribble(
  ~Correction, ~matrix, ~analysis_role,
  "Raw betas", "Raw/pre-ComBat", "primary",
  "ComBat betas", "ComBat", "sensitivity"
)

clock_data <- readRDS(eaa_clock_file) |>
  as_tibble() |>
  mutate(
    Sample_Name = as.character(.data$Sample_Name),
    Dataset = as.character(.data$Dataset),
    DiseaseGroup = as.character(.data$DiseaseGroup),
    Sex = as.character(.data$Sex),
    Correction = as.character(.data$Correction)
  ) |>
  left_join(bmi_lookup |> select(Sample_Name, Dataset, BMI), by = c("Sample_Name", "Dataset")) |>
  left_join(correction_roles, by = "Correction")

extract_group_term <- function(fit, model_label) {
  broom::tidy(fit, conf.int = TRUE) |>
    filter(.data$term == "FibrosisGroupAdvanced_Fibrosis") |>
    transmute(
      model = model_label,
      beta = .data$estimate,
      conf_low = .data$conf.low,
      conf_high = .data$conf.high,
      p_value = .data$p.value
    )
}

fit_johnson_clock <- function(clock, correction_value) {
  correction_meta <- correction_roles |>
    filter(.data$Correction == correction_value)

  johnson <- clock_data |>
    filter(
      .data$Correction == correction_value,
      .data$Dataset == "Johnson",
      .data$DiseaseGroup %in% c("Healthy_Obese", "Advanced_Fibrosis")
    ) |>
    mutate(
      FibrosisGroup = factor(
        .data$DiseaseGroup,
        levels = c("Healthy_Obese", "Advanced_Fibrosis")
      ),
      Sex = factor(.data$Sex)
    )

  bmi_group_summary <- johnson |>
    distinct(Sample_Name, FibrosisGroup, BMI) |>
    group_by(FibrosisGroup) |>
    summarise(
      n = sum(!is.na(BMI)),
      BMI_mean = mean(BMI, na.rm = TRUE),
      BMI_SD = sd(BMI, na.rm = TRUE),
      BMI_median = median(BMI, na.rm = TRUE),
      BMI_IQR = stats::IQR(BMI, na.rm = TRUE, type = 8),
      .groups = "drop"
    ) |>
    ungroup()

  outcome <- paste0(clock, "_EAA_resid")
  if (!outcome %in% names(johnson)) {
    return(tibble(
      Clock = clock,
      Correction = correction_value,
      matrix = correction_meta$matrix,
      analysis_role = correction_meta$analysis_role,
      n = 0L,
      note = paste("Missing EAA column:", outcome)
    ))
  }

  d <- johnson |>
    select(Sample_Name, FibrosisGroup, Sex, BMI, rEAA = all_of(outcome)) |>
    drop_na(FibrosisGroup, Sex, BMI, rEAA)

  if (nrow(d) < 10 || n_distinct(d$FibrosisGroup) < 2) {
    return(tibble(
      Clock = clock,
      n = nrow(d),
      Correction = correction_value,
      matrix = correction_meta$matrix,
      analysis_role = correction_meta$analysis_role,
      note = "Insufficient complete Johnson BMI cases"
    ))
  }

  fit_no_bmi <- lm(rEAA ~ FibrosisGroup + Sex, data = d)
  fit_bmi <- lm(rEAA ~ FibrosisGroup + Sex + BMI, data = d)

  no_bmi <- extract_group_term(fit_no_bmi, "no_BMI") |>
    rename_with(~ paste0(.x, "_no_BMI"), -model) |>
    select(-model)
  bmi <- extract_group_term(fit_bmi, "BMI_adjusted") |>
    rename_with(~ paste0(.x, "_BMI_adjusted"), -model) |>
    select(-model)

  bind_cols(
    tibble(
      Clock = clock,
      Correction = correction_value,
      matrix = correction_meta$matrix,
      analysis_role = correction_meta$analysis_role,
      comparison = "Advanced_Fibrosis vs Healthy_Obese",
      n = nrow(d)
    ),
    no_bmi,
    bmi,
    tibble(note = NA_character_)
  )
}

bmi_group_wide <- bind_rows(lapply(correction_roles$Correction, function(correction_value) {
  clock_data |>
    filter(
      .data$Correction == correction_value,
      .data$Dataset == "Johnson",
      .data$DiseaseGroup %in% c("Healthy_Obese", "Advanced_Fibrosis")
    ) |>
    mutate(
      FibrosisGroup = factor(
        .data$DiseaseGroup,
        levels = c("Healthy_Obese", "Advanced_Fibrosis")
      )
    ) |>
    distinct(Correction, matrix, analysis_role, Sample_Name, FibrosisGroup, BMI) |>
    group_by(Correction, matrix, analysis_role, FibrosisGroup) |>
    summarise(
      BMI_n = sum(!is.na(BMI)),
      BMI_mean = mean(BMI, na.rm = TRUE),
      BMI_SD = sd(BMI, na.rm = TRUE),
      BMI_median = median(BMI, na.rm = TRUE),
      BMI_IQR = stats::IQR(BMI, na.rm = TRUE, type = 8),
      .groups = "drop"
    )
})) |>
  tidyr::pivot_wider(
    names_from = FibrosisGroup,
    values_from = c(BMI_n, BMI_mean, BMI_SD, BMI_median, BMI_IQR),
    names_glue = "{FibrosisGroup}_{.value}"
  ) |>
  rename(
    Healthy_Obese_BMI_n = Healthy_Obese_BMI_n,
    Healthy_Obese_BMI_mean = Healthy_Obese_BMI_mean,
    Healthy_Obese_BMI_SD = Healthy_Obese_BMI_SD,
    Healthy_Obese_BMI_median = Healthy_Obese_BMI_median,
    Healthy_Obese_BMI_IQR = Healthy_Obese_BMI_IQR,
    Advanced_Fibrosis_BMI_n = Advanced_Fibrosis_BMI_n,
    Advanced_Fibrosis_BMI_mean = Advanced_Fibrosis_BMI_mean,
    Advanced_Fibrosis_BMI_SD = Advanced_Fibrosis_BMI_SD,
    Advanced_Fibrosis_BMI_median = Advanced_Fibrosis_BMI_median,
    Advanced_Fibrosis_BMI_IQR = Advanced_Fibrosis_BMI_IQR
  )

sensitivity <- bind_rows(lapply(correction_roles$Correction, function(correction_value) {
  bind_rows(lapply(main_clock_cols, fit_johnson_clock, correction_value = correction_value))
})) |>
  group_by(Correction, matrix, analysis_role) |>
  mutate(
    BH_no_BMI = p.adjust(.data$p_value_no_BMI, method = "BH"),
    BH_BMI_adjusted = p.adjust(.data$p_value_BMI_adjusted, method = "BH"),
    beta_change_BMI_minus_no_BMI = .data$beta_BMI_adjusted - .data$beta_no_BMI,
    relative_beta_change_percent = 100 * beta_change_BMI_minus_no_BMI / abs(.data$beta_no_BMI),
    direction_concordant = sign(.data$beta_no_BMI) == sign(.data$beta_BMI_adjusted),
    FDR_status_change = case_when(
      .data$BH_no_BMI < 0.05 & .data$BH_BMI_adjusted >= 0.05 ~ "lost_FDR_lt_0.05_after_BMI",
      .data$BH_no_BMI >= 0.05 & .data$BH_BMI_adjusted < 0.05 ~ "gained_FDR_lt_0.05_after_BMI",
      TRUE ~ "unchanged"
    )
  ) |>
  ungroup() |>
  left_join(bmi_group_wide, by = c("Correction", "matrix", "analysis_role")) |>
  select(
    Clock,
    Correction,
    matrix,
    analysis_role,
    comparison,
    n,
    beta_no_BMI,
    conf_low_no_BMI,
    conf_high_no_BMI,
    p_value_no_BMI,
    BH_no_BMI,
    beta_BMI_adjusted,
    conf_low_BMI_adjusted,
    conf_high_BMI_adjusted,
    p_value_BMI_adjusted,
    BH_BMI_adjusted,
    beta_change_BMI_minus_no_BMI,
    relative_beta_change_percent,
    direction_concordant,
    FDR_status_change,
    Healthy_Obese_BMI_n,
    Healthy_Obese_BMI_mean,
    Healthy_Obese_BMI_SD,
    Healthy_Obese_BMI_median,
    Healthy_Obese_BMI_IQR,
    Advanced_Fibrosis_BMI_n,
    Advanced_Fibrosis_BMI_mean,
    Advanced_Fibrosis_BMI_SD,
    Advanced_Fibrosis_BMI_median,
    Advanced_Fibrosis_BMI_IQR,
    note
  )

analysis_table_dir <- dir_out_epi_tables
supp_table_dir <- file.path(dir_results, "supplementary", "tables")
thesis_out_supp_table_dir <- file.path(project_root, "outputs", "thesis_outputs", "supplementary", "tables")

audit_file <- file.path(analysis_table_dir, "bmi_availability_audit_eaa.csv")
sensitivity_file <- file.path(analysis_table_dir, "johnson_bmi_eaa_sensitivity.csv")
supp_audit_file <- file.path(supp_table_dir, "T37_bmi_availability_audit_eaa.csv")
supp_sensitivity_file <- file.path(supp_table_dir, "T38_johnson_bmi_eaa_sensitivity.csv")
pub_supp_audit_file <- file.path(thesis_out_supp_table_dir, "T37_bmi_availability_audit_eaa.csv")
pub_supp_sensitivity_file <- file.path(thesis_out_supp_table_dir, "T38_johnson_bmi_eaa_sensitivity.csv")

write_csv_safe(audit, audit_file)
write_csv_safe(sensitivity, sensitivity_file)
write_csv_safe(audit, supp_audit_file)
write_csv_safe(sensitivity, supp_sensitivity_file)
write_csv_safe(audit, pub_supp_audit_file)
write_csv_safe(sensitivity, pub_supp_sensitivity_file)

message("BMI audit written to: ", audit_file)
message("Johnson BMI sensitivity written to: ", sensitivity_file)
