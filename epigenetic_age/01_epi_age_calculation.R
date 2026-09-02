suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(forcats)
  library(patchwork)
  library(methylCIPHER)
  library(CTSclocks)
  library(broom)
  library(GGally)
  library(rlang)
})

source(file.path("data_preprocessing", "00_config.R"))
source(file.path(dir_helpers, "helpers_io.R"))
source(file.path(dir_helpers, "helpers_plotting.R"))

ensure_dirs(all_output_dirs)

check_file_exists(file_beta_raw, "Raw imputed beta matrix")
check_file_exists(file_beta_combat, "ComBat imputed beta matrix")
check_file_exists(file_pheno_all, "Combined phenotype table")

beta_raw <- readRDS(file_beta_raw)
beta_cb <- readRDS(file_beta_combat)
ph <- readRDS(file_pheno_all)

check_object_columns(
  ph,
  c("Sample_Name", "Dataset", "Age", "Sex", "DiseaseGroup", "Progression3"),
  "ph"
)

stopifnot(!anyDuplicated(ph$Sample_Name))

for (nm in c("beta_raw", "beta_cb")) {
  b <- get(nm)
  miss <- setdiff(ph$Sample_Name, colnames(b))
  if (length(miss) > 0) {
    stop(
      nm, ": missing samples in beta colnames (first 10): ",
      paste(head(miss, 10), collapse = ", "),
      call. = FALSE
    )
  }
  b <- b[, ph$Sample_Name, drop = FALSE]
  stopifnot(identical(colnames(b), ph$Sample_Name))
  assign(nm, b)
}

if (!identical(rownames(beta_raw), rownames(beta_cb))) {
  common_cpgs <- intersect(rownames(beta_raw), rownames(beta_cb))
  beta_raw <- beta_raw[common_cpgs, , drop = FALSE]
  beta_cb <- beta_cb[common_cpgs, , drop = FALSE]
  stopifnot(identical(rownames(beta_raw), rownames(beta_cb)))
}

normalize_diseasegroup <- function(x) {
  x <- as.character(x)
  x[!(x %in% stage_levels)] <- NA_character_
  factor(x, levels = stage_levels)
}

normalize_progression3 <- function(x) {
  x <- as.character(x)
  x[!(x %in% prog3_levels)] <- NA_character_
  factor(x, levels = prog3_levels)
}

ph_age <- ph |>
  filter(!is.na(Sample_Name), !is.na(Age), !is.na(Sex)) |>
  mutate(
    Dataset = droplevels(as.factor(Dataset)),
    DiseaseGroup = normalize_diseasegroup(DiseaseGroup),
    Progression3 = normalize_progression3(Progression3)
  )

keep_samples <- ph_age$Sample_Name
beta_age_raw <- beta_raw[, keep_samples, drop = FALSE]
beta_age_cb <- beta_cb[, keep_samples, drop = FALSE]

stopifnot(identical(colnames(beta_age_raw), keep_samples))
stopifnot(identical(colnames(beta_age_cb), keep_samples))

save_rds_safe(ph_age, file.path(dir_out_epi_rds, "clock_cohort_pheno.rds"))
write_csv_safe(ph_age, file.path(dir_out_epi_tables, "clock_cohort_pheno.csv"))

write_csv_safe(dplyr::count(ph_age, Dataset, name = "n"), file.path(dir_out_epi_tables, "clock_cohort_dataset_counts.csv"))
write_csv_safe(dplyr::count(ph_age, DiseaseGroup, name = "n"), file.path(dir_out_epi_tables, "clock_cohort_diseasegroup_counts.csv"))

ph_plot <- ph_age |>
  mutate(Dataset = factor(as.character(Dataset), levels = sort(unique(as.character(Dataset)))))

counts_df <- ph_plot |>
  keep_diseasegroup(stage_levels) |>
  count(Dataset, DiseaseGroup, name = "n") |>
  group_by(Dataset) |>
  mutate(DiseaseGroup_rev = forcats::fct_rev(DiseaseGroup)) |>
  ungroup()

p_counts <- ggplot(counts_df, aes(x = DiseaseGroup_rev, y = n, fill = DiseaseGroup)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  geom_text(aes(label = n), size = 3.4) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = stage_colors, drop = FALSE) +
  labs(title = "Cohort composition (DiseaseGroup)", x = NULL, y = "N samples") +
  facet_wrap(~ Dataset, scales = "free_x") +
  theme_pipeline(13) +
  theme(legend.position = "none", plot.margin = margin(5.5, 30, 5.5, 5.5))

p_age <- ggplot(
  ph_plot |> keep_diseasegroup(stage_levels) |> filter(!is.na(Age)),
  aes(x = DiseaseGroup, y = Age, fill = DiseaseGroup)
) +
  geom_boxplot(width = 0.6, outlier.shape = NA, color = "black", linewidth = 0.25) +
  geom_jitter(width = 0.15, size = 1.1, alpha = 0.35) +
  scale_fill_manual(values = stage_colors, drop = FALSE) +
  labs(title = "Age distribution by DiseaseGroup", x = NULL, y = "Age") +
  facet_wrap(~ Dataset) +
  theme_pipeline(13) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

fig <- p_counts | p_age
save_plot(fig, file.path(dir_out_epi_plots, "cohort_overview_faceted.png"), width = 14, height = 5, dpi = 600)
save_plot(fig, file.path(dir_out_epi_plots, "cohort_overview_faceted.pdf"), width = 14, height = 5, device = cairo_pdf)

df_mix <- ph_age |>
  keep_diseasegroup(stage_levels) |>
  mutate(Dataset = as.character(Dataset)) |>
  count(DiseaseGroup, Dataset, name = "n") |>
  group_by(DiseaseGroup) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

p_mix <- ggplot(df_mix, aes(x = DiseaseGroup, y = Dataset, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.25) +
  theme_pipeline(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "QC: Dataset composition per DiseaseGroup",
    subtitle = "If a DiseaseGroup is mostly one Dataset, batch and biology are confounded",
    x = "DiseaseGroup", y = "Dataset", fill = "Proportion"
  )

save_plot(
  p_mix,
  file.path(dir_out_epi_plots, "QC_dataset_composition_per_DiseaseGroup.png"),
  width = 12, height = 5, dpi = 400
)

calc_eaa_resid <- function(df, clock_col, age_col = "Age", id_col = "Sample_Name", min_n = 10) {
  df2 <- df |> filter(!is.na(.data[[clock_col]]), !is.na(.data[[age_col]]))
  out_col <- paste0(clock_col, "_EAA_resid")
  df[[out_col]] <- NA_real_

  if (nrow(df2) < min_n) return(df)

  sd_clock <- sd(df2[[clock_col]], na.rm = TRUE)
  if (!is.finite(sd_clock) || sd_clock == 0) return(df)

  fit <- tryCatch(lm(reformulate(age_col, clock_col), data = df2), error = function(e) NULL)
  if (is.null(fit)) return(df)

  idx <- match(df2[[id_col]], df[[id_col]])
  df[[out_col]][idx] <- resid(fit)
  df
}

quick_num_summary <- function(x) {
  if (!is.numeric(x)) return("non-numeric")
  if (all(is.na(x))) return("ALL_NA")
  paste0(
    "n=", length(x),
    " | NA=", sum(is.na(x)),
    " | min=", suppressWarnings(min(x, na.rm = TRUE)),
    " | med=", suppressWarnings(median(x, na.rm = TRUE)),
    " | max=", suppressWarnings(max(x, na.rm = TRUE))
  )
}

check_clock_cols <- function(df, cols) {
  cols <- intersect(cols, colnames(df))
  if (!length(cols)) return(invisible(NULL))
  out <- sapply(cols, function(nm) {
    x <- df[[nm]]
    c(
      class = paste(class(x), collapse = "/"),
      summary = quick_num_summary(x)
    )
  })
  as.data.frame(t(out))
}

resolve_optional_clock_paths <- function() {
  default_dir <- Sys.getenv(
    "METHYLCIPHER_DATA_DIR",
    file.path(project_root, "data", "references", "methylCIPHER")
  )

  pc_file <- if (exists("file_pcclocks_rdata", inherits = TRUE)) {
    get("file_pcclocks_rdata", inherits = TRUE)
  } else if (nzchar(Sys.getenv("PCCLOCKS_DATA"))) {
    Sys.getenv("PCCLOCKS_DATA")
  } else {
    file.path(default_dir, "PCClocks_data.qs2")
  }

  systems_file <- if (exists("file_systemsage_rdata", inherits = TRUE)) {
    get("file_systemsage_rdata", inherits = TRUE)
  } else if (nzchar(Sys.getenv("SYSTEMSAGE_DATA"))) {
    Sys.getenv("SYSTEMSAGE_DATA")
  } else {
    file.path(default_dir, "SystemsAge_data.qs2")
  }

  list(pc_file = pc_file, systems_file = systems_file)
}

pcclocks_cache <- new.env(parent = emptyenv())

load_pcclocks_reference <- function(path) {
  expected_md5 <- "0f7fa3a89b6559e98478a9e986d28db0"
  expected_hash <- "bd7932093eabfcea4b8a6da3f6662413"

  path <- normalizePath(path, mustWork = TRUE)
  if (identical(pcclocks_cache$path, path)) {
    return(pcclocks_cache$obj)
  }

  observed_md5 <- unname(tools::md5sum(path))
  if (!identical(observed_md5, expected_md5)) {
    stop(
      "PCClocks_data.qs2 MD5 mismatch. Expected ", expected_md5,
      ", observed ", observed_md5,
      call. = FALSE
    )
  }

  obj <- qs2::qs_read(path, validate_checksum = TRUE)
  observed_hash <- rlang::hash(obj)

  if (!identical(observed_hash, expected_hash)) {
    warning(
      "PCClocks_data.qs2 object hash differs from the Zenodo reference hash: ",
      observed_hash,
      call. = FALSE
    )
  }

  required_names <- c(
    "imputeMissingCpGs",
    "CalcPCHorvath1",
    "CalcPCHorvath2",
    "CalcPCHannum",
    "CalcPCPhenoAge",
    "CalcPCDNAmTL",
    "CalcPCGrimAge"
  )
  missing_names <- setdiff(required_names, names(obj))
  if (length(missing_names) > 0) {
    stop(
      "PCClocks_data.qs2 is missing required components: ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }

  write_csv_safe(
    tibble::tibble(
      reference = "PCClocks_data.qs2",
      source = "Zenodo DOI 10.5281/zenodo.17162604 / 10.5281/zenodo.19455622",
      file = path,
      md5 = observed_md5,
      object_hash = observed_hash,
      methylCIPHER_version = as.character(utils::packageVersion("methylCIPHER")),
      note = "Zenodo MD5 and qs2 checksum verified; local wrapper follows current methylCIPHER PCClocks calculation because installed methylCIPHER 0.2.0 has an obsolete hard-coded object hash."
    ),
    file.path(dir_out_epi_tables, "pcclocks_reference_manifest.csv")
  )

  pcclocks_cache$path <- path
  pcclocks_cache$obj <- obj
  obj
}

calc_pcclocks_zenodo <- function(DNAm, pheno, ID = "Sample_ID", RData) {
  if (is.character(RData)) {
    RData <- load_pcclocks_reference(RData)
  }

  if (!isTRUE(all.equal(rownames(DNAm), pheno[[ID]]))) {
    samples <- intersect(rownames(DNAm), pheno[[ID]])
    if (length(samples) == 0) {
      stop("DNAm and pheno have no ID in common.", call. = FALSE)
    }
    DNAm <- DNAm[samples, , drop = FALSE]
    pheno <- pheno[match(samples, pheno[[ID]]), , drop = FALSE]
    stopifnot(identical(rownames(DNAm), pheno[[ID]]))
  }

  DNAm <- methylCIPHER::impute_DNAm(
    DNAm = DNAm,
    method = "mean",
    CpGs = RData$imputeMissingCpGs,
    subset = TRUE
  )
  DNAm <- DNAm[, names(RData$imputeMissingCpGs), drop = FALSE]

  message("Calculating PC Clocks now")
  pheno$PCHorvath1 <- as.numeric(methylCIPHER:::anti.trafo(sweep(DNAm, 2, RData$CalcPCHorvath1$center) %*% RData$CalcPCHorvath1$rotation %*% RData$CalcPCHorvath1$model + RData$CalcPCHorvath1$intercept))
  pheno$PCHorvath2 <- as.numeric(methylCIPHER:::anti.trafo(sweep(DNAm, 2, RData$CalcPCHorvath2$center) %*% RData$CalcPCHorvath2$rotation %*% RData$CalcPCHorvath2$model + RData$CalcPCHorvath2$intercept))
  pheno$PCHannum <- as.numeric(sweep(DNAm, 2, RData$CalcPCHannum$center) %*% RData$CalcPCHannum$rotation %*% RData$CalcPCHannum$model + RData$CalcPCHannum$intercept)
  pheno$PCPhenoAge <- as.numeric(sweep(DNAm, 2, RData$CalcPCPhenoAge$center) %*% RData$CalcPCPhenoAge$rotation %*% RData$CalcPCPhenoAge$model + RData$CalcPCPhenoAge$intercept)
  pheno$PCDNAmTL <- as.numeric(sweep(DNAm, 2, RData$CalcPCDNAmTL$center) %*% RData$CalcPCDNAmTL$rotation %*% RData$CalcPCDNAmTL$model + RData$CalcPCDNAmTL$intercept)

  DNAm_grim <- cbind(
    sweep(DNAm, 2, RData$CalcPCGrimAge$center) %*% RData$CalcPCGrimAge$rotation,
    Female = pheno$Female,
    Age = pheno$Age
  )
  pheno$PCPACKYRS <- as.numeric(DNAm_grim[, names(RData$CalcPCGrimAge$PCPACKYRS.model)] %*% RData$CalcPCGrimAge$PCPACKYRS.model + RData$CalcPCGrimAge$PCPACKYRS.intercept)
  pheno$PCADM <- as.numeric(DNAm_grim[, names(RData$CalcPCGrimAge$PCADM.model)] %*% RData$CalcPCGrimAge$PCADM.model + RData$CalcPCGrimAge$PCADM.intercept)
  pheno$PCB2M <- as.numeric(DNAm_grim[, names(RData$CalcPCGrimAge$PCB2M.model)] %*% RData$CalcPCGrimAge$PCB2M.model + RData$CalcPCGrimAge$PCB2M.intercept)
  pheno$PCCystatinC <- as.numeric(DNAm_grim[, names(RData$CalcPCGrimAge$PCCystatinC.model)] %*% RData$CalcPCGrimAge$PCCystatinC.model + RData$CalcPCGrimAge$PCCystatinC.intercept)
  pheno$PCGDF15 <- as.numeric(DNAm_grim[, names(RData$CalcPCGrimAge$PCGDF15.model)] %*% RData$CalcPCGrimAge$PCGDF15.model + RData$CalcPCGrimAge$PCGDF15.intercept)
  pheno$PCLeptin <- as.numeric(DNAm_grim[, names(RData$CalcPCGrimAge$PCLeptin.model)] %*% RData$CalcPCGrimAge$PCLeptin.model + RData$CalcPCGrimAge$PCLeptin.intercept)
  pheno$PCPAI1 <- as.numeric(DNAm_grim[, names(RData$CalcPCGrimAge$PCPAI1.model)] %*% RData$CalcPCGrimAge$PCPAI1.model + RData$CalcPCGrimAge$PCPAI1.intercept)
  pheno$PCTIMP1 <- as.numeric(DNAm_grim[, names(RData$CalcPCGrimAge$PCTIMP1.model)] %*% RData$CalcPCGrimAge$PCTIMP1.model + RData$CalcPCGrimAge$PCTIMP1.intercept)
  pheno$PCGrimAge <- as.numeric(as.matrix(subset(pheno, select = RData$CalcPCGrimAge$components)) %*% RData$CalcPCGrimAge$PCGrimAge.model + RData$CalcPCGrimAge$PCGrimAge.intercept)

  message("PC Clocks successfully calculated!")
  pheno
}

get_cts_clock_probes <- function(clock) {
  env <- new.env(parent = emptyenv())
  data(list = paste0(clock, "Clock"), package = "CTSclocks", envir = env)
  clock_obj <- env[[paste0(clock, "Clock.glm")]]
  rownames(clock_obj$beta)
}

get_cts_clock_coverage <- function(beta_age, clocks = c(HepClock = "Hep", LiverClock = "Liver")) {
  bind_rows(lapply(names(clocks), function(clock_name) {
    probes <- get_cts_clock_probes(clocks[[clock_name]])
    n_available <- sum(probes %in% rownames(beta_age))

    tibble(
      ClockName = clock_name,
      TotalProbes = length(probes),
      ProbesAvailable = n_available,
      PercentProbesAvailable = 100 * n_available / length(probes)
    )
  }))
}

run_cts_clock <- function(beta_age, sample_ids, clock, clock_name, correction_label) {
  stopifnot(identical(colnames(beta_age), sample_ids))

  pred <- tryCatch(
    CTSclocks::CTSclockAge(
      data.m = beta_age,
      CTSclock = clock,
      dataType = "bulk",
      CTF.m = NULL,
      tissue = "otherTissue"
    ),
    error = function(e) {
      warning(
        "[", correction_label, "] ", clock_name, " failed: ",
        conditionMessage(e),
        call. = FALSE
      )
      NULL
    }
  )

  if (is.null(pred)) return(NULL)
  pred <- as.numeric(pred)
  if (length(pred) != length(sample_ids)) {
    stop(
      "[", correction_label, "] ", clock_name, " returned ", length(pred),
      " predictions for ", length(sample_ids), " samples",
      call. = FALSE
    )
  }
  pred
}

cap_iqr <- function(x) {
  qs <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE, type = 8)
  iqr <- qs[2] - qs[1]
  lo <- qs[1] - 1.5 * iqr
  hi <- qs[2] + 1.5 * iqr
  pmin(pmax(x, lo), hi)
}

run_clocks_and_eaa <- function(beta_age, ph_age_df, correction_label) {
  tag <- tagify(correction_label)

  beta_age_t <- t(beta_age)
  if (sum(is.na(beta_age_t)) > 0) {
    stop("[", correction_label, "] NAs in transposed beta matrix", call. = FALSE)
  }

  coverage <- methylCIPHER::getClockProbes(beta_age_t)
  coverage_std <- coverage |>
    transmute(
      ClockName = Clock,
      TotalProbes = Total.Probes,
      ProbesAvailable = Present.Probes,
      PercentProbesAvailable = as.numeric(gsub("%", "", Percent.Present))
    )
  coverage_std <- bind_rows(coverage_std, get_cts_clock_coverage(beta_age))

  write_csv_safe(coverage_std, file.path(dir_out_epi_tables, paste0("clock_probe_coverage_", tag, ".csv")))

  ph_clock <- ph_age_df |>
    mutate(
      Sample_ID = Sample_Name,
      Sex = as.character(Sex),
      Female = case_when(
        Sex == "Female" ~ 1L,
        Sex == "Male" ~ 0L,
        TRUE ~ NA_integer_
      )
    )

  stopifnot(!anyDuplicated(ph_clock$Sample_ID))
  stopifnot(is.numeric(ph_clock$Age), !anyNA(ph_clock$Age))
  stopifnot(!anyNA(ph_clock$Female), all(ph_clock$Female %in% c(0L, 1L)))

  ph_clock <- as.data.frame(ph_clock)
  rownames(ph_clock) <- ph_clock$Sample_ID

  stopifnot(all(ph_clock$Sample_ID %in% rownames(beta_age_t)))
  beta_clock <- beta_age_t[ph_clock$Sample_ID, , drop = FALSE]
  stopifnot(identical(rownames(beta_clock), ph_clock$Sample_ID))

  beta_clock <- as.matrix(beta_clock)
  storage.mode(beta_clock) <- "double"
  stopifnot(sum(is.na(beta_clock)) == 0)

  run_clock <- function(label, FUN, post_cols = NULL) {
    out <- tryCatch(
      FUN(),
      error = function(e) {
        warning("[", correction_label, "] ", label, " failed: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (!is.null(out)) {
      out <- as.data.frame(out)
      if (!is.null(post_cols)) {
        check_clock_cols(out, post_cols)
      }
    }
    out
  }

  run_horvath1 <- TRUE
  run_hannum <- TRUE
  run_phenoage <- TRUE
  run_horvath2 <- TRUE
  run_grimage <- isTRUE(use_grimage)
  run_pcclocks <- isTRUE(use_pc_clocks)
  run_systemsage <- isTRUE(use_systems_age)
  run_hrspheno <- TRUE

  opt_paths <- resolve_optional_clock_paths()

  ph_clocks <- ph_clock

  if (run_horvath1) {
    tmp <- run_clock(
      "calcHorvath1",
      function() methylCIPHER::calcHorvath1(DNAm = beta_clock, pheno = ph_clocks, imputation = FALSE),
      post_cols = c("Horvath1")
    )
    if (!is.null(tmp)) ph_clocks <- tmp
  }

  if (run_hannum) {
    tmp <- run_clock(
      "calcHannum",
      function() methylCIPHER::calcHannum(DNAm = beta_clock, pheno = ph_clocks, imputation = FALSE),
      post_cols = c("Hannum")
    )
    if (!is.null(tmp)) ph_clocks <- tmp
  }

  if (run_phenoage) {
    tmp <- run_clock(
      "calcPhenoAge",
      function() methylCIPHER::calcPhenoAge(DNAm = beta_clock, pheno = ph_clocks, imputation = FALSE),
      post_cols = c("PhenoAge")
    )
    if (!is.null(tmp)) ph_clocks <- tmp
  }

  if (run_horvath2) {
    tmp <- run_clock(
      "calcHorvath2",
      function() methylCIPHER::calcHorvath2(DNAm = beta_clock, pheno = ph_clocks, imputation = FALSE),
      post_cols = c("Horvath2")
    )
    if (!is.null(tmp)) ph_clocks <- tmp
  }

  if (run_grimage) {
    tmp <- run_clock(
      "calcGrimAgeV1",
      function() methylCIPHER::calcGrimAgeV1(DNAm = beta_clock, pheno = ph_clocks, ID = "Sample_ID"),
      post_cols = c("GrimAgeV1")
    )
    if (!is.null(tmp)) ph_clocks <- tmp

    if ("calcGrimAgeV2" %in% getNamespaceExports("methylCIPHER")) {
      tmp <- run_clock(
        "calcGrimAgeV2",
        function() methylCIPHER::calcGrimAgeV2(DNAm = beta_clock, pheno = ph_clocks, ID = "Sample_ID"),
        post_cols = c("GrimAgeV2")
      )
      if (!is.null(tmp)) ph_clocks <- tmp
    }
  }

  if (run_pcclocks) {
    stopifnot(file.exists(opt_paths$pc_file))
    tmp <- run_clock(
      "calcPCClocks",
      function() calc_pcclocks_zenodo(
        DNAm = beta_clock, pheno = ph_clocks, ID = "Sample_ID", RData = opt_paths$pc_file
      ),
      post_cols = c("PCHannum", "PCPhenoAge", "PCHorvath1", "PCHorvath2", "PCGrimAge")
    )
    if (!is.null(tmp)) ph_clocks <- tmp
  }

  if (run_systemsage) {
    stopifnot(file.exists(opt_paths$systems_file))
    tmp <- run_clock(
      "calcSystemsAge",
      function() methylCIPHER::calcSystemsAge(
        DNAm = beta_clock, pheno = ph_clocks, ID = "Sample_ID", RData = opt_paths$systems_file
      )
    )
    if (!is.null(tmp)) ph_clocks <- tmp
  }

  if (run_hrspheno) {
    tmp <- run_clock(
      "calcHRSInChPhenoAge",
      function() methylCIPHER::calcHRSInChPhenoAge(DNAm = beta_clock, pheno = ph_clocks, imputation = FALSE),
      post_cols = c("HRSInChPhenoAge")
    )
    if (!is.null(tmp)) ph_clocks <- tmp
  }

  hep_pred <- run_cts_clock(beta_age, ph_clock$Sample_ID, "Hep", "HepClock", correction_label)
  if (!is.null(hep_pred)) ph_clocks$HepClock <- hep_pred

  liver_pred <- run_cts_clock(beta_age, ph_clock$Sample_ID, "Liver", "LiverClock", correction_label)
  if (!is.null(liver_pred)) ph_clocks$LiverClock <- liver_pred

  clock_cols_present <- intersect(main_clock_cols, colnames(ph_clocks))
  diff_ok_cols <- intersect(c("Horvath1", "Horvath2", "Hannum", "PhenoAge"), clock_cols_present)

  for (cc in clock_cols_present) {
    resid_col <- paste0(cc, "_EAA_resid")
    if (!resid_col %in% colnames(ph_clocks)) {
      ph_clocks <- calc_eaa_resid(ph_clocks, cc, min_n = 10)
    }
  }

  for (cc in diff_ok_cols) {
    diff_col <- paste0(cc, "_EAA_diff")
    if (!diff_col %in% colnames(ph_clocks)) {
      ph_clocks[[diff_col]] <- ifelse(
        is.na(ph_clocks[[cc]]) | is.na(ph_clocks$Age),
        NA_real_,
        ph_clocks[[cc]] - ph_clocks$Age
      )
    }
  }

  ph_clocks <- ph_clocks |> mutate(Correction = correction_label)

  save_rds_safe(ph_clocks, file.path(dir_out_epi_rds, paste0("pheno_with_clocks_and_EAA_", tag, ".rds")))
  write_tsv_safe(ph_clocks, file.path(dir_out_epi_tables, paste0("pheno_with_clocks_and_EAA_", tag, ".tsv")))

  ph_clocks
}

ph_clocks_raw <- run_clocks_and_eaa(beta_age_raw, ph_age, "Raw betas")
ph_clocks_cb <- run_clocks_and_eaa(beta_age_cb, ph_age, "ComBat betas")
ph_clocks_both <- bind_rows(ph_clocks_raw, ph_clocks_cb)

save_rds_safe(ph_clocks_both, file.path(dir_out_epi_rds, "pheno_with_clocks_and_EAA_raw_combat.rds"))
write_tsv_safe(ph_clocks_both, file.path(dir_out_epi_tables, "pheno_with_clocks_and_EAA_raw_combat.tsv"))

ph_clocks_both <- ph_clocks_both |>
  mutate(
    Correction = factor(as.character(Correction), levels = c("ComBat betas", "Raw betas")),
    DiseaseGroup = factor(as.character(DiseaseGroup), levels = stage_levels),
    Progression3 = factor(as.character(Progression3), levels = prog3_levels)
  )

clock_cols_plot <- intersect(main_clock_cols, colnames(ph_clocks_both))

df_scatter_dg <- ph_clocks_both |>
  select(Correction, Age, DiseaseGroup, any_of(clock_cols_plot)) |>
  pivot_longer(cols = any_of(clock_cols_plot), names_to = "Clock", values_to = "ClockValue") |>
  filter(!is.na(Age), !is.na(ClockValue), !is.na(DiseaseGroup)) |>
  mutate(Clock = factor(Clock, levels = clock_cols_plot))

make_scatter <- function(df_long, correction_value, color_var, palette, subtitle) {
  df_sub <- df_long |> filter(Correction == correction_value)

  labels <- ph_clocks_both |>
    dplyr::filter(Correction == correction_value) |>
    dplyr::distinct(Sample_Name, .data[[color_var]]) |>
    make_count_labels(!!rlang::sym(color_var))

  ggplot(df_sub, aes(x = Age, y = ClockValue, color = .data[[color_var]])) +
    geom_jitter(width = 0.25, height = 0, alpha = 0.6, size = 0.9) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.35) +
    facet_wrap(~ Clock, scales = "free_y", ncol = 3) +
    scale_color_manual(values = palette, labels = labels, drop = FALSE) +
    labs(
      title = paste0("Clock outputs vs chronological age - ", correction_value),
      subtitle = subtitle,
      x = "Chronological age (years)",
      y = "Clock output",
      color = color_var
    ) +
    theme_pipeline(10) +
    theme(legend.position = "bottom", strip.text = element_text(size = 8)) +
    guides(color = guide_legend(nrow = 2, override.aes = list(size = 2, alpha = 1)))
}

p_dg_cb <- make_scatter(df_scatter_dg, "ComBat betas", "DiseaseGroup", stage_colors, "Color: DiseaseGroup")
p_dg_raw <- make_scatter(df_scatter_dg, "Raw betas", "DiseaseGroup", stage_colors, "Color: DiseaseGroup")

save_plot(p_dg_cb, file.path(dir_out_epi_plots, "scatter_DiseaseGroup_ComBat.png"), width = 12, height = 8, dpi = 450)
save_plot(p_dg_raw, file.path(dir_out_epi_plots, "scatter_DiseaseGroup_Raw.png"), width = 12, height = 8, dpi = 450)

p_all_scatter <- (p_dg_cb / p_dg_raw) + plot_layout(heights = c(1, 1))
save_plot(p_all_scatter, file.path(dir_out_epi_plots, "scatter_DiseaseGroup_RAW_vs_ComBat.png"), width = 12, height = 16, dpi = 450)
save_plot(p_all_scatter, file.path(dir_out_epi_plots, "scatter_DiseaseGroup_RAW_vs_ComBat.pdf"), width = 12, height = 16)

eaa_resid_cols <- intersect(paste0(clock_cols_plot, "_EAA_resid"), colnames(ph_clocks_both))
eaa_resid_cols <- eaa_resid_cols[!sapply(eaa_resid_cols, function(nm) all(is.na(ph_clocks_both[[nm]])))]

df_eaa_dg_both <- ph_clocks_both |>
  select(Correction, Dataset, DiseaseGroup, any_of(eaa_resid_cols)) |>
  pivot_longer(cols = any_of(eaa_resid_cols), names_to = "Clock", values_to = "EAA_resid") |>
  mutate(
    Clock = gsub("_EAA_resid$", "", Clock),
    Clock = factor(Clock, levels = clock_cols_plot),
    DiseaseGroup = factor(as.character(DiseaseGroup), levels = stage_levels)
  ) |>
  filter(!is.na(EAA_resid), !is.na(DiseaseGroup)) |>
  group_by(Correction, Clock) |>
  mutate(EAA_resid_plot = cap_iqr(EAA_resid)) |>
  ungroup()

make_violin <- function(df_long, correction_value) {
  df_sub <- df_long |> filter(Correction == correction_value)

  labels <- ph_clocks_both |>
    dplyr::filter(Correction == correction_value) |>
    dplyr::distinct(Sample_Name, DiseaseGroup) |>
    make_count_labels(DiseaseGroup)

  ggplot(df_sub, aes(x = DiseaseGroup, y = EAA_resid_plot, fill = DiseaseGroup)) +
    geom_violin(alpha = 0.55, color = NA, trim = TRUE) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.12, alpha = 0.20, size = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~ Clock, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = stage_colors, labels = labels, drop = FALSE) +
    labs(
      title = paste0("EAA (residual) by DiseaseGroup - ", correction_value),
      subtitle = "Values capped at 1.5×IQR within each Clock/Correction",
      x = "DiseaseGroup",
      y = "EAA (residuals from lm(Clock ~ Age))"
    ) +
    theme_pipeline(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 9))
}

p_violin_cb <- make_violin(df_eaa_dg_both, "ComBat betas")
p_violin_raw <- make_violin(df_eaa_dg_both, "Raw betas")
p_violin_all <- p_violin_cb / p_violin_raw

save_plot(p_violin_cb, file.path(dir_out_epi_plots, "EAA_violin_DiseaseGroup_ComBat.png"), width = 12, height = 8, dpi = 450)
save_plot(p_violin_raw, file.path(dir_out_epi_plots, "EAA_violin_DiseaseGroup_Raw.png"), width = 12, height = 8, dpi = 450)
save_plot(p_violin_all, file.path(dir_out_epi_plots, "EAA_violin_DiseaseGroup_RAW_vs_ComBat.png"), width = 12, height = 16, dpi = 450)

fit_stage_effect <- function(df, outcome, stage_col, covars = c("Sex", "Dataset")) {
  df2 <- df |> tidyr::drop_na(all_of(c(outcome, stage_col, covars)))
  if (nrow(df2) < 10) return(NULL)

  df2[[stage_col]] <- factor(df2[[stage_col]])
  if ("Healthy" %in% levels(df2[[stage_col]])) {
    df2[[stage_col]] <- relevel(df2[[stage_col]], ref = "Healthy")
  }

  fml <- as.formula(paste(outcome, "~", stage_col, "+", paste(covars, collapse = " + ")))
  fit <- lm(fml, data = df2)

  broom::tidy(fit, conf.int = TRUE) |>
    filter(grepl(paste0("^", stage_col), term)) |>
    mutate(
      Clock = gsub("_EAA_resid$", "", outcome),
      StageVar = stage_col,
      Comparison = gsub(paste0("^", stage_col), "", term)
    )
}

eaa_resid_cols2 <- eaa_resid_cols[!sapply(eaa_resid_cols, function(nm) all(is.na(ph_clocks_both[[nm]])))]

res_dg <- bind_rows(lapply(levels(ph_clocks_both$Correction), function(cc) {
  df_sub <- filter(ph_clocks_both, Correction == cc)
  out <- bind_rows(lapply(eaa_resid_cols2, function(x) fit_stage_effect(df_sub, x, "DiseaseGroup")))
  if (is.null(out) || nrow(out) == 0) return(out)
  out |> mutate(Correction = cc)
}))

if (!is.null(res_dg) && nrow(res_dg) > 0) {
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
  forest_clock_order <- c(intersect(forest_clock_order, clock_cols_plot), setdiff(clock_cols_plot, forest_clock_order))
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

  res_dg <- res_dg |>
    dplyr::group_by(Correction, StageVar) |>
    dplyr::mutate(
      p_adj_bh = stats::p.adjust(p.value, method = "BH"),
      fdr_bh = p_adj_bh,
      fdr_significant = p_adj_bh < 0.05
    ) |>
    dplyr::ungroup() |>
    mutate(
      Comparison = factor(Comparison, levels = disease_order),
      Clock = factor(Clock, levels = clock_cols_plot)
    )

  make_beta_plot <- function(df, correction_value, title_prefix) {
    df_sub <- df |>
      filter(Correction == correction_value) |>
      mutate(
        Clock_chr = as.character(Clock),
        Comparison_chr = as.character(Comparison)
      )

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
      select(Clock_chr, Comparison_chr, estimate, conf.low, conf.high) |>
      right_join(forest_rows, by = c("Clock_chr", "Comparison_chr")) |>
      mutate(
        Comparison_label = factor(Comparison_label, levels = unname(disease_labels[disease_order])),
        estimate_label = if_else(is.na(estimate), NA_character_, sprintf("%+.2f", estimate)),
        estimate_label_x = estimate + if_else(estimate >= 0, 0.25, -0.25),
        estimate_label_hjust = if_else(estimate >= 0, 0, 1)
      )

    clock_groups <- forest_rows |>
      group_by(Clock_chr, Clock_label) |>
      summarise(
        ymin = min(y) - 0.42,
        ymax = max(y) + 0.42,
        ymid = mean(y),
        .groups = "drop"
      )

    x_vals <- c(plot_df$conf.low, plot_df$conf.high, plot_df$estimate)
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

    forest_plot <- ggplot(plot_df, aes(x = estimate, y = y, color = Comparison_label)) +
      geom_hline(
        data = forest_rows,
        aes(yintercept = y),
        inherit.aes = FALSE,
        color = "#E8EDF2",
        linewidth = 0.18
      ) +
      geom_vline(xintercept = 0, color = "#315C8A", linewidth = 0.42) +
      geom_errorbar(
        aes(xmin = conf.low, xmax = conf.high),
        orientation = "y",
        width = 0.20,
        linewidth = 0.42,
        na.rm = TRUE
      ) +
      geom_point(size = 1.25, na.rm = TRUE) +
      geom_text(
        aes(x = estimate_label_x, y = y + 0.27, label = estimate_label, hjust = estimate_label_hjust),
        color = "#111827",
        size = 1.55,
        show.legend = FALSE,
        na.rm = TRUE
      ) +
      geom_text(
        data = forest_rows,
        aes(y = y, label = Comparison_label),
        inherit.aes = FALSE,
        x = stage_label_x,
        hjust = 1,
        size = 1.62,
        color = "#40516A"
      ) +
      geom_segment(
        data = clock_groups,
        aes(y = ymin, yend = ymax),
        inherit.aes = FALSE,
        x = bracket_x,
        xend = bracket_x,
        color = "#5F6F82",
        linewidth = 0.28
      ) +
      geom_segment(
        data = clock_groups,
        aes(y = ymin, yend = ymin),
        inherit.aes = FALSE,
        x = bracket_x,
        xend = bracket_tick_x,
        color = "#5F6F82",
        linewidth = 0.28
      ) +
      geom_segment(
        data = clock_groups,
        aes(y = ymax, yend = ymax),
        inherit.aes = FALSE,
        x = bracket_x,
        xend = bracket_tick_x,
        color = "#5F6F82",
        linewidth = 0.28
      ) +
      geom_text(
        data = clock_groups,
        aes(y = ymid, label = Clock_label),
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
        legend.key.size = unit(0.46, "cm"),
        legend.justification = "center",
        legend.box.just = "center",
        legend.box.margin = margin(t = 0, r = 190, b = 0, l = 0),
        legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
        plot.margin = margin(7, 8, 7, 170)
      ) +
      guides(color = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5, override.aes = list(size = 2.6, linewidth = 0.8)))

    title_plot <- ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.70,
        label = paste0(title_prefix, " - ", correction_value),
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

  p_beta_dg_cb <- make_beta_plot(res_dg, "ComBat betas", "DiseaseGroup effect on residual EAA")
  p_beta_dg_raw <- make_beta_plot(res_dg, "Raw betas", "DiseaseGroup effect on residual EAA")
  p_beta_dg_all <- patchwork::wrap_elements(full = p_beta_dg_cb) /
    patchwork::wrap_elements(full = p_beta_dg_raw)

  save_plot(p_beta_dg_cb, file.path(dir_out_epi_plots, "beta_CI_DiseaseGroup_ComBat.png"), width = 7.2, height = 9.4, dpi = 450)
  save_plot(p_beta_dg_cb, file.path(dir_out_epi_plots, "beta_CI_DiseaseGroup_ComBat.pdf"), width = 7.2, height = 9.4, device = cairo_pdf)
  save_plot(p_beta_dg_raw, file.path(dir_out_epi_plots, "beta_CI_DiseaseGroup_Raw.png"), width = 7.2, height = 9.4, dpi = 450)
  save_plot(p_beta_dg_raw, file.path(dir_out_epi_plots, "beta_CI_DiseaseGroup_Raw.pdf"), width = 7.2, height = 9.4, device = cairo_pdf)
  save_plot(p_beta_dg_all, file.path(dir_out_epi_plots, "beta_CI_DiseaseGroup_RAW_vs_ComBat.png"), width = 7.2, height = 18.8, dpi = 450)

  write_csv_safe(res_dg, file.path(dir_out_epi_tables, "effect_sizes_DiseaseGroup_vs_Healthy_RAW_vs_ComBat.csv"))
}

# Same model on the ordinal three-level staging variable. Only the effect-size
# table is written here; the corresponding figures are produced downstream by
# epigenetic_age/03_eaa_audit_outputs.R, which also assembles the supplementary
# EAA tables (T29, T30) from both effect-size tables.
res_prog3 <- bind_rows(lapply(levels(ph_clocks_both$Correction), function(cc) {
  df_sub <- filter(ph_clocks_both, Correction == cc)
  out <- bind_rows(lapply(eaa_resid_cols2, function(x) fit_stage_effect(df_sub, x, "Progression3")))
  if (is.null(out) || nrow(out) == 0) return(out)
  out |> mutate(Correction = cc)
}))

if (!is.null(res_prog3) && nrow(res_prog3) > 0) {
  prog_order <- c("Mild", "Advanced")
  res_prog3 <- res_prog3 |>
    dplyr::group_by(Correction, StageVar) |>
    dplyr::mutate(
      p_adj_bh = stats::p.adjust(p.value, method = "BH"),
      fdr_bh = p_adj_bh,
      fdr_significant = p_adj_bh < 0.05
    ) |>
    dplyr::ungroup() |>
    mutate(
      Comparison = factor(Comparison, levels = prog_order),
      Clock = factor(Clock, levels = clock_cols_plot)
    )

  write_csv_safe(res_prog3, file.path(dir_out_epi_tables, "effect_sizes_Progression3_vs_Healthy_RAW_vs_ComBat.csv"))
}

panel_cor_tile <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)

  r <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs", method = "pearson"))
  if (!is.finite(r)) r <- NA_real_
  r_clamped <- if (is.na(r)) NA_real_ else max(min(r, 1), -1)
  lab <- if (is.na(r_clamped)) "NA" else sprintf("%.2f", r_clamped)

  ggplot(data.frame(r = r_clamped), aes(x = 0, y = 0, fill = r)) +
    geom_tile(width = 1, height = 1) +
    geom_text(aes(label = lab), fontface = "bold", size = 5) +
    scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, na.value = "grey90") +
    theme_void() +
    theme(legend.position = "none")
}

panel_points <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.25, size = 0.7) +
    theme_pipeline(11) +
    theme(panel.grid.minor = element_blank())
}

panel_diag_name <- function(data, mapping, ...) {
  var <- rlang::as_label(mapping$x)
  ggplot() +
    annotate("text", x = 0, y = 0, label = var, fontface = "bold", size = 5) +
    theme_void()
}

run_pairs <- function(df, vars, out_name) {
  dfX <- df |> select(all_of(vars))

  var_ok <- sapply(dfX, function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) return(FALSE)
    stats::sd(x) > 0
  })

  dfX <- dfX[, var_ok, drop = FALSE]
  if (ncol(dfX) < 2) return(invisible(NULL))

  p <- GGally::ggpairs(
    dfX,
    upper = list(continuous = panel_cor_tile),
    lower = list(continuous = panel_points),
    diag = list(continuous = panel_diag_name)
  ) +
    theme_pipeline(11) +
    theme(panel.grid.minor = element_blank())

  save_plot(p, file.path(dir_out_epi_plots, out_name), width = 13, height = 10)
  cor_matrix <- cor(dfX, use = "pairwise.complete.obs")
  write_csv_safe(as.data.frame(cor_matrix), file.path(dir_out_epi_tables, sub("\\.png$", ".csv", out_name)))
}

clock_cols_ok <- clock_cols_plot[!sapply(clock_cols_plot, function(nm) all(is.na(ph_clocks_both[[nm]])))]

for (cc in unique(ph_clocks_both$Correction)) {
  df_sub <- filter(ph_clocks_both, Correction == cc)

  run_pairs(
    df_sub,
    c(clock_cols_ok, "Age"),
    paste0("pairs_epigenetic_clocks_", tagify(cc), ".png")
  )

  eaa_cols_ok <- eaa_resid_cols
  eaa_cols_ok <- eaa_cols_ok[!sapply(eaa_cols_ok, function(nm) all(is.na(df_sub[[nm]])))]

  run_pairs(
    df_sub,
    c(eaa_cols_ok, "Age"),
    paste0("pairs_EAA_resid_", tagify(cc), ".png")
  )
}

cap_iqr <- function(x) {
  qs <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE, type = 8)
  iqr <- qs[2] - qs[1]
  lo <- qs[1] - 1.5 * iqr
  hi <- qs[2] + 1.5 * iqr
  pmin(pmax(x, lo), hi)
}
