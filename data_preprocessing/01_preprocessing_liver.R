suppressPackageStartupMessages({
  library(GEOquery)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(minfi)
  library(readxl)
  library(R.utils)
  library(ChAMP)
  library(ggplot2)
})

source(file.path("data_preprocessing", "00_config.R"))
source(file.path(dir_helpers, "helpers_io.R"))
source(file.path(dir_helpers, "helpers_plotting.R"))

ensure_dirs(all_output_dirs)
ensure_dir(dir_combined)

# =========================================================
# 01) Shared config + utilities
# =========================================================

detP_thresh <- 0.01

`%||%` <- function(a, b) if (!is.null(a)) a else b

normalize_sex <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("", "na", "n/a", "unknown", "not available")] <- NA_character_
  dplyr::case_when(
    x %in% c("m", "male", "1") ~ "Male",
    x %in% c("f", "female", "2") ~ "Female",
    TRUE ~ NA_character_
  )
}

extract_field_ci_vec <- function(row_vec, key) {
  x <- as.character(row_vec)
  hit <- x[grepl(paste0("^", key, "\\s*:"), x, ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("^", key, "\\s*:[[:space:]]*"), "", hit[1], ignore.case = TRUE)
}

to_chr <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.list(x)) x <- vapply(x, function(z) paste(z, collapse = ";"), character(1))
  trimws(as.character(x))
}

normalize_sample_id <- function(x) {
  x <- to_chr(x)
  x[x == "" | is.na(x)] <- NA_character_
  x <- sub("\\.0+$", "", x)
  x <- gsub("\\s+", "", x)
  toupper(x)
}

get_col_or_na <- function(df, nm) {
  if (nm %in% colnames(df)) df[[nm]] else rep(NA_character_, nrow(df))
}

signals_to_beta <- function(M_mat, U_mat, offset = 100) {
  M <- pmax(M_mat, 0)
  U <- pmax(U_mat, 0)
  M / (M + U + offset)
}

sanitize_id_component <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[x == "" | is.na(x)] <- NA_character_
  x
}

collapse_matrix_columns_by_group <- function(mat, group_ids, collapsed_names) {
  stopifnot(ncol(mat) == length(group_ids), length(group_ids) == length(collapsed_names))

  group_ids <- as.character(group_ids)
  group_levels <- unique(group_ids)
  out <- matrix(
    NA_real_,
    nrow = nrow(mat),
    ncol = length(group_levels),
    dimnames = list(rownames(mat), unname(collapsed_names[match(group_levels, group_ids)]))
  )

  for (i in seq_along(group_levels)) {
    idx <- which(group_ids == group_levels[[i]])
    if (length(idx) == 1) {
      out[, i] <- mat[, idx]
    } else {
      vals <- rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
      vals[is.nan(vals)] <- NA_real_
      out[, i] <- vals
    }
  }

  out
}

collapse_technical_replicates <- function(beta, pheno, group_ids, dataset_label, id_prefix, detP = NULL) {
  stopifnot(ncol(beta) == nrow(pheno))
  stopifnot(identical(colnames(beta), as.character(pheno$Sample_Name)))

  group_ids <- as.character(group_ids)
  group_ids[is.na(group_ids) | trimws(group_ids) == ""] <- as.character(pheno$Sample_Name)[is.na(group_ids) | trimws(group_ids) == ""]

  group_levels <- unique(group_ids)
  collapsed_names <- paste0(id_prefix, "_", sanitize_id_component(group_ids))

  beta_collapsed <- collapse_matrix_columns_by_group(beta, group_ids, collapsed_names)

  detP_collapsed <- NULL
  if (!is.null(detP)) {
    stopifnot(identical(colnames(detP), as.character(pheno$Sample_Name)))
    detP_collapsed <- collapse_matrix_columns_by_group(detP, group_ids, collapsed_names)
  }

  first_idx <- match(group_levels, group_ids)
  pheno_collapsed <- pheno[first_idx, , drop = FALSE]
  pheno_collapsed$Sample_Name <- unname(collapsed_names[first_idx])
  pheno_collapsed$Technical_Replicate_Group <- group_levels
  pheno_collapsed$Source_Profile_IDs <- vapply(
    group_levels,
    function(g) paste(as.character(pheno$Sample_Name)[group_ids == g], collapse = ";"),
    character(1)
  )
  pheno_collapsed$N_Technical_Profiles <- as.integer(table(factor(group_ids, levels = group_levels)))

  audit <- tibble::tibble(
    Dataset = dataset_label,
    Biological_Sample_ID = pheno_collapsed$Sample_Name,
    Technical_Replicate_Group = pheno_collapsed$Technical_Replicate_Group,
    Source_Profile_IDs = pheno_collapsed$Source_Profile_IDs,
    n_technical_profiles = pheno_collapsed$N_Technical_Profiles
  )

  message(
    dataset_label,
    ": collapsed ",
    length(group_ids),
    " array profiles to ",
    length(group_levels),
    " biological samples (",
    sum(audit$n_technical_profiles > 1),
    " replicate groups)."
  )

  list(beta = beta_collapsed, pheno = pheno_collapsed, detP = detP_collapsed, audit = audit)
}

# ----------------------------
# IDAT helpers
# ----------------------------
raw_tar_has_idats <- function(raw_tar) {
  stopifnot(length(raw_tar) == 1, file.exists(raw_tar))
  files <- untar(raw_tar, list = TRUE)
  any(grepl("\\.idat(\\.gz)?$", files, ignore.case = TRUE))
}

extract_raw_tar <- function(raw_tar, exdir) {
  stopifnot(length(raw_tar) == 1, file.exists(raw_tar))
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
  untar(raw_tar, exdir = exdir)

  gz <- list.files(exdir, pattern = "\\.idat\\.gz$", full.names = TRUE, recursive = TRUE)
  if (length(gz) > 0) {
    invisible(lapply(gz, function(f) R.utils::gunzip(f, overwrite = TRUE, remove = FALSE)))
  }
}

dir_has_idats <- function(idat_dir) {
  if (!dir.exists(idat_dir)) return(FALSE)
  any(grepl("\\.idat(\\.gz)?$", list.files(idat_dir, recursive = TRUE), ignore.case = TRUE))
}

read_idats_noob <- function(idat_dir) {
  rg <- minfi::read.metharray.exp(base = idat_dir, extended = TRUE)
  detP <- minfi::detectionP(rg)
  mset <- minfi::preprocessNoob(rg)
  beta <- minfi::getBeta(mset)

  gsm <- stringr::str_extract(colnames(beta), "GSM\\d+")
  if (!any(is.na(gsm))) {
    colnames(beta) <- gsm
    colnames(detP) <- gsm
  }

  list(beta = beta, detP = detP)
}

# ----------------------------
# ChAMP wrapper
# ----------------------------
champ_prefilter_beta <- function(beta, arraytype, detP = NULL, sample_names = NULL, do_bmiq = TRUE) {
  stopifnot(!is.null(beta), is.matrix(beta) || is.data.frame(beta))
  beta <- as.matrix(beta)

  if (is.null(sample_names)) sample_names <- colnames(beta)
  pd <- data.frame(Sample_Name = sample_names, row.names = sample_names, stringsAsFactors = FALSE)

  if (!is.null(detP)) {
    detP <- as.matrix(detP)
    common_probes <- intersect(rownames(beta), rownames(detP))
    beta <- beta[common_probes, , drop = FALSE]
    detP <- detP[common_probes, , drop = FALSE]
    stopifnot(all(colnames(beta) == colnames(detP)))
  }

  res <- ChAMP::champ.filter(
    beta = beta,
    pd = pd,
    detP = detP,
    autoimpute = FALSE,
    filterDetP = !is.null(detP),
    filterBeads = FALSE,
    fixOutlier = FALSE,
    filterNoCG = TRUE,
    filterSNPs = TRUE,
    filterMultiHit = TRUE,
    filterXY = TRUE,
    arraytype = arraytype
  )

  beta_f <- as.matrix(res$beta)

  if (do_bmiq) {
    beta_f <- ChAMP::champ.norm(
      beta = beta_f,
      arraytype = arraytype,
      method = "BMIQ"
    )
  }

  as.matrix(beta_f)
}

# ----------------------------
# Generic loader
# ----------------------------
load_beta_idat_or_signal <- function(
  gse_id,
  raw_root_dir,
  idat_dir = file.path(raw_root_dir, "idat"),
  raw_tar_pattern = paste0(gse_id, "_RAW\\.tar$"),
  arraytype,
  signal_reader = NULL,
  prefer_processed_beta = NULL
) {
  if (dir_has_idats(idat_dir)) {
    obj <- read_idats_noob(idat_dir)
    return(list(beta = obj$beta, detP = obj$detP, route = "idat"))
  }

  raw_tar <- list.files(raw_root_dir, pattern = raw_tar_pattern, full.names = TRUE, recursive = TRUE)
  if (length(raw_tar) == 1 && raw_tar_has_idats(raw_tar)) {
    extract_raw_tar(raw_tar, idat_dir)
    obj <- read_idats_noob(idat_dir)
    return(list(beta = obj$beta, detP = obj$detP, route = "idat"))
  }

  if (!is.null(prefer_processed_beta)) {
    beta <- prefer_processed_beta()
    return(list(beta = beta, detP = NULL, route = "processed_beta"))
  }

  if (is.null(signal_reader)) stop(gse_id, ": no IDATs and no fallback reader provided", call. = FALSE)
  sig <- signal_reader()
  list(beta = sig$beta, detP = sig$detP %||% NULL, route = "signal")
}

# =========================================================
# 02) Ahrens (GSE48325)
# =========================================================

ahrens_name <- "Ahrens_GSE48325"
ahrens_root <- file.path(dir_raw, ahrens_name)
ahrens_out <- file.path(dir_processed, ahrens_name)
ensure_dir(ahrens_out)

ahrens_signal_reader <- function() {
  f <- file.path(ahrens_root, "GSE48325", "GSE48325_signal_intensities.txt.gz")
  check_file_exists(f, "Ahrens signal file")

  df <- readr::read_tsv(f, show_col_types = FALSE)
  cn <- colnames(df)
  signalA_cols <- grep("Signal_A$", cn, value = TRUE)
  signalB_cols <- grep("Signal_B$", cn, value = TRUE)
  detP_cols <- grep("Detection Pval$", cn, value = TRUE)

  sample_nums_char <- sub("\\.Signal_A$", "", signalA_cols)

  U <- as.matrix(df[, signalA_cols])
  M <- as.matrix(df[, signalB_cols])
  rownames(U) <- df$`ID-REF`
  rownames(M) <- df$`ID-REF`

  beta <- signals_to_beta(M, U, offset = 100)
  colnames(beta) <- sample_nums_char

  detP <- NULL
  if (length(detP_cols) > 0) {
    detP <- as.matrix(df[, detP_cols])
    rownames(detP) <- df$`ID-REF`
    detP_names <- sub("\\.Detection Pval$", "", detP_cols)
    colnames(detP) <- detP_names
    if (all(sample_nums_char %in% detP_names)) {
      detP <- detP[, sample_nums_char, drop = FALSE]
    }
  }

  list(beta = beta, detP = detP)
}

ahrens_obj <- load_beta_idat_or_signal(
  gse_id = "GSE48325",
  raw_root_dir = ahrens_root,
  arraytype = "450K",
  signal_reader = ahrens_signal_reader
)

beta_ahrens_raw <- ahrens_obj$beta
detP_ahrens_raw <- ahrens_obj$detP

gse_ahrens <- GEOquery::getGEO("GSE48325", GSEMatrix = TRUE)[[1]]
ph <- Biobase::pData(gse_ahrens)
char_cols <- grep("^characteristics_ch1", colnames(ph), value = TRUE)

extract_field_ahr <- function(x, key) {
  x_chr <- as.character(x)
  x_low <- tolower(x_chr)
  key_low <- tolower(key)
  prefix <- paste0(key_low, ":")

  idx <- which(startsWith(x_low, prefix))
  if (!length(idx)) return(NA_character_)
  trimws(substr(x_chr[idx[1]], nchar(prefix) + 1, nchar(x_chr[idx[1]])))
}

disease_vec <- apply(ph[, char_cols, drop = FALSE], 1, extract_field_ahr, key = "group")
age_vec <- apply(ph[, char_cols, drop = FALSE], 1, extract_field_ahr, key = "age")
sex_code <- apply(ph[, char_cols, drop = FALSE], 1, extract_field_ahr, key = "sex (1 - male, 2 - female)")
sex_code <- suppressWarnings(as.integer(sex_code))
sex_vec <- dplyr::case_when(
  sex_code == 1 ~ "Male",
  sex_code == 2 ~ "Female",
  TRUE ~ NA_character_
)

bariatric_vec <- apply(
  ph[, char_cols, drop = FALSE], 1, extract_field_ahr,
  key = "bariatric surgery (1 - before surgery, 2 - after)"
)
bariatric_num <- suppressWarnings(as.integer(bariatric_vec))

sample_num_from_title <- ph$title |>
  stringr::str_extract("sample\\s*[0-9]+") |>
  stringr::str_extract("[0-9]+") |>
  as.integer()

pheno_ahrens_tmp <- tibble::tibble(
  Sample_Name = as.character(ph$geo_accession),
  Sample_Num = sample_num_from_title,
  DiseaseState = disease_vec,
  Age = suppressWarnings(as.numeric(age_vec)),
  Sex = normalize_sex(sex_vec),
  Bariatric_raw = bariatric_vec,
  Bariatric_num = bariatric_num
)

if (all(grepl("^GSM\\d+$", colnames(beta_ahrens_raw)))) {
  beta_ahrens <- beta_ahrens_raw
  detP_ahrens <- detP_ahrens_raw

  pheno_ahrens_clean <- pheno_ahrens_tmp |>
    dplyr::filter(Sample_Name %in% colnames(beta_ahrens)) |>
    dplyr::arrange(match(Sample_Name, colnames(beta_ahrens)))

  keep_idx <- is.na(pheno_ahrens_clean$Bariatric_num) | pheno_ahrens_clean$Bariatric_num != 2
  pheno_ahrens_clean <- pheno_ahrens_clean[keep_idx, , drop = FALSE]

  beta_ahrens <- beta_ahrens[, pheno_ahrens_clean$Sample_Name, drop = FALSE]
  if (!is.null(detP_ahrens)) {
    detP_ahrens <- detP_ahrens[, pheno_ahrens_clean$Sample_Name, drop = FALSE]
  }
} else {
  beta_sample_num <- suppressWarnings(as.integer(colnames(beta_ahrens_raw)))

  mapping <- tibble::tibble(
    beta_col_idx = seq_along(beta_sample_num),
    Sample_Num = beta_sample_num
  ) |>
    dplyr::left_join(pheno_ahrens_tmp, by = "Sample_Num")

  stopifnot(!any(is.na(mapping$Sample_Name)))

  mapping <- mapping |>
    dplyr::filter(is.na(Bariatric_num) | Bariatric_num != 2)

  beta_ahrens <- beta_ahrens_raw[, mapping$beta_col_idx, drop = FALSE]
  colnames(beta_ahrens) <- mapping$Sample_Name

  detP_ahrens <- NULL
  if (!is.null(detP_ahrens_raw)) {
    detP_ahrens <- detP_ahrens_raw[, as.character(mapping$Sample_Num), drop = FALSE]
    colnames(detP_ahrens) <- mapping$Sample_Name
  }

  pheno_ahrens_clean <- mapping |>
    dplyr::select(
      Sample_Name, Sample_Num, DiseaseState, Age, Sex, Bariatric_raw, Bariatric_num
    )
}

stopifnot(all(pheno_ahrens_clean$Sample_Name == colnames(beta_ahrens)))

beta_ahrens_f <- champ_prefilter_beta(beta_ahrens, arraytype = "450K", detP = detP_ahrens)

save_rds_safe(beta_ahrens_f, file.path(ahrens_out, "beta_ahrens.rds"))
save_rds_safe(pheno_ahrens_clean, file.path(ahrens_out, "pheno_ahrens_clean.rds"))
if (!is.null(detP_ahrens)) save_rds_safe(detP_ahrens, file.path(ahrens_out, "detP_ahrens.rds"))

# =========================================================
# 03) Horvath (GSE61258)
# =========================================================

horvath_name <- "Horvath_GSE61258"
horvath_root <- file.path(dir_raw, horvath_name)
horvath_out <- file.path(dir_processed, horvath_name)
ensure_dir(horvath_out)

horvath_signal_reader <- function() {
  f <- file.path(horvath_root, "GSE61258_datSignalAndBFinalLiver.csv.gz")
  check_file_exists(f, "Horvath signal file")
  df <- readr::read_csv(f, show_col_types = FALSE)

  cn <- colnames(df)
  signalA_cols <- grep("Signal_A$|SignalA$", cn, value = TRUE)
  signalB_cols <- grep("Signal_B$|SignalB$", cn, value = TRUE)

  sample_ids_char <- sub("\\.Signal_A$", "", signalA_cols)
  sample_nums <- suppressWarnings(as.integer(sub("^Sample", "", sample_ids_char)))

  U <- as.matrix(df[, signalA_cols])
  M <- as.matrix(df[, signalB_cols])
  rownames(U) <- df$TargeID
  rownames(M) <- df$TargeID

  beta <- signals_to_beta(M, U, offset = 100)
  colnames(beta) <- sample_nums

  list(beta = beta, detP = NULL)
}

horvath_obj <- load_beta_idat_or_signal(
  gse_id = "GSE61258",
  raw_root_dir = horvath_root,
  arraytype = "450K",
  signal_reader = horvath_signal_reader
)

beta_horvath_raw <- horvath_obj$beta

gse_horvath <- GEOquery::getGEO("GSE61258", GSEMatrix = TRUE)[[1]]
ph <- Biobase::pData(gse_horvath)
char_cols <- grep("^characteristics_ch1", colnames(ph), value = TRUE)

extract_field_h <- function(x, key) {
  hit <- x[grepl(paste0("^", key, ":"), x)]
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("^", key, ":[[:space:]]*"), "", hit[1])
}

if (all(grepl("^GSM\\d+$", colnames(beta_horvath_raw)))) {
  beta_horvath <- beta_horvath_raw
} else {
  stopifnot(ncol(beta_horvath_raw) == nrow(ph))
  beta_horvath <- beta_horvath_raw
  colnames(beta_horvath) <- as.character(ph$geo_accession)
}

keys <- c(
  "Sex", "dna.extraction.method", "age", "dnamage",
  "ageaccelerationresidualbasedontrainingandtestdata", "bmi",
  "diseasestatus", "nafldactivityscore", "liversteatosis",
  "liverinflammation", "fibrosis", "nas.ballooning",
  "dnamageadjustedagebmi", "dnamageadjustedagebmigender",
  "smokingstatus", "smokingpackyears", "subjectid"
)

char_list <- lapply(keys, function(k) apply(ph[, char_cols, drop = FALSE], 1, extract_field_h, key = k))
names(char_list) <- keys

pheno_horvath_clean <- tibble::tibble(
  Sample_Name = as.character(ph$geo_accession),
  subjectid = char_list$subjectid,
  Sex = normalize_sex(char_list$Sex),
  Age = suppressWarnings(as.numeric(char_list$age)),
  DNAmAge = suppressWarnings(as.numeric(char_list$dnamage)),
  AgeAccelResidual = suppressWarnings(as.numeric(char_list$ageaccelerationresidualbasedontrainingandtestdata)),
  BMI = suppressWarnings(as.numeric(char_list$bmi)),
  DiseaseState = char_list$diseasestatus,
  NAFLDActivityScore = suppressWarnings(as.numeric(char_list$nafldactivityscore)),
  Steatosis = suppressWarnings(as.numeric(char_list$liversteatosis)),
  Inflammation = suppressWarnings(as.numeric(char_list$liverinflammation)),
  Fibrosis = trimws(as.character(char_list$fibrosis)),
  NAS_Ballooning = suppressWarnings(as.numeric(char_list$`nas.ballooning`)),
  DNAmAgeAdj_BMI = suppressWarnings(as.numeric(char_list$dnamageadjustedagebmi)),
  DNAmAgeAdj_BMI_Sex = suppressWarnings(as.numeric(char_list$dnamageadjustedagebmigender)),
  SmokingStatus = char_list$smokingstatus,
  SmokingPackYears = suppressWarnings(as.numeric(char_list$smokingpackyears))
) |>
  dplyr::filter(Sample_Name %in% colnames(beta_horvath)) |>
  dplyr::arrange(match(Sample_Name, colnames(beta_horvath)))

beta_horvath <- beta_horvath[, pheno_horvath_clean$Sample_Name, drop = FALSE]
stopifnot(all(pheno_horvath_clean$Sample_Name == colnames(beta_horvath)))

beta_horvath_f <- champ_prefilter_beta(beta_horvath, arraytype = "450K", detP = NULL)

save_rds_safe(beta_horvath_f, file.path(horvath_out, "beta_horvath.rds"))
save_rds_safe(pheno_horvath_clean, file.path(horvath_out, "pheno_horvath_clean.rds"))

# =========================================================
# 04) Johnson (GSE180474)
# =========================================================

johnson_name <- "Johnson_GSE180474"
johnson_root <- file.path(dir_raw, johnson_name)
johnson_out <- file.path(dir_processed, johnson_name)
ensure_dir(johnson_out)

gse_johnson <- GEOquery::getGEO("GSE180474", GSEMatrix = TRUE)[[1]]
ph_j <- Biobase::pData(gse_johnson)
char_cols_j <- grep("^characteristics_ch1", colnames(ph_j), value = TRUE)

extract_field_j <- function(x, key) {
  hit <- x[grepl(paste0("^", key, ":"), x, ignore.case = FALSE)]
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("^", key, ":[[:space:]]*"), "", hit[1])
}

tissue_vec <- apply(ph_j[, char_cols_j, drop = FALSE], 1, extract_field_j, key = "tissue")
sex_vec <- apply(ph_j[, char_cols_j, drop = FALSE], 1, extract_field_j, key = "Sex")
disease_state_vec <- apply(ph_j[, char_cols_j, drop = FALSE], 1, extract_field_j, key = "disease state")
age_vec <- apply(ph_j[, char_cols_j, drop = FALSE], 1, extract_field_j, key = "age")
bmi_vec <- apply(ph_j[, char_cols_j, drop = FALSE], 1, extract_field_j, key = "bmi")
t2d_med_vec <- apply(ph_j[, char_cols_j, drop = FALSE], 1, extract_field_j, key = "t2d medication")
chipid_vec <- apply(ph_j[, char_cols_j, drop = FALSE], 1, extract_field_j, key = "chip id")
indiv_id_vec <- apply(ph_j[, char_cols_j, drop = FALSE], 1, extract_field_j, key = "individual id")

pheno_johnson_full <- tibble::tibble(
  Sample_Name = as.character(ph_j$geo_accession),
  Tissue = tissue_vec,
  Sex = normalize_sex(sex_vec),
  DiseaseState = disease_state_vec,
  Age = suppressWarnings(as.numeric(age_vec)),
  BMI = suppressWarnings(as.numeric(bmi_vec)),
  T2D_med = t2d_med_vec,
  Individual_ID = indiv_id_vec,
  chipID = trimws(chipid_vec)
)

pheno_johnson_liver <- pheno_johnson_full |>
  dplyr::filter(tolower(Tissue) == "liver") |>
  dplyr::select(Sample_Name, Tissue, DiseaseState, Sex, Age, BMI, T2D_med, Individual_ID, chipID)

johnson_signal_reader <- function() {
  f <- file.path(johnson_root, "GSE180474", "GSE180474_MatrixSignalIntensities.csv.gz")
  check_file_exists(f, "Johnson signal file")
  df <- readr::read_csv(f, show_col_types = FALSE)

  cn <- colnames(df)
  u_cols <- grep("Unmethylated Signal$", cn, value = TRUE)
  m_cols <- grep("Methylated Signal$", cn, value = TRUE)

  probe_ids <- df[[1]]
  U <- as.matrix(df[, u_cols])
  M <- as.matrix(df[, m_cols])
  rownames(U) <- probe_ids
  rownames(M) <- probe_ids

  beta <- signals_to_beta(M, U, offset = 100)
  sample_ids_raw <- sub(" Unmethylated Signal$", "", u_cols)
  colnames(beta) <- sample_ids_raw

  list(beta = beta, detP = NULL)
}

johnson_obj <- load_beta_idat_or_signal(
  gse_id = "GSE180474",
  raw_root_dir = johnson_root,
  arraytype = "EPIC",
  signal_reader = johnson_signal_reader,
  prefer_processed_beta = NULL
)

beta_johnson_raw <- johnson_obj$beta
detP_johnson_raw <- johnson_obj$detP

if (all(grepl("^GSM\\d+$", colnames(beta_johnson_raw)))) {
  beta_johnson <- beta_johnson_raw[, pheno_johnson_liver$Sample_Name, drop = FALSE]
  detP_johnson <- if (!is.null(detP_johnson_raw)) detP_johnson_raw[, pheno_johnson_liver$Sample_Name, drop = FALSE] else NULL
} else {
  stopifnot(ncol(beta_johnson_raw) == nrow(ph_j))
  colnames(beta_johnson_raw) <- ph_j$geo_accession
  beta_johnson <- beta_johnson_raw[, pheno_johnson_liver$Sample_Name, drop = FALSE]

  detP_johnson <- NULL
  if (!is.null(detP_johnson_raw)) {
    colnames(detP_johnson_raw) <- ph_j$geo_accession
    detP_johnson <- detP_johnson_raw[, pheno_johnson_liver$Sample_Name, drop = FALSE]
  }
}

stopifnot(all(pheno_johnson_liver$Sample_Name == colnames(beta_johnson)))

beta_johnson_f <- champ_prefilter_beta(beta_johnson, arraytype = "EPIC", detP = detP_johnson)

johnson_replicate_resolution <- collapse_technical_replicates(
  beta = beta_johnson_f,
  pheno = pheno_johnson_liver,
  group_ids = pheno_johnson_liver$Individual_ID,
  dataset_label = "Johnson",
  id_prefix = "Johnson_individual",
  detP = detP_johnson
)

beta_johnson_f <- johnson_replicate_resolution$beta
pheno_johnson_liver <- johnson_replicate_resolution$pheno
detP_johnson <- johnson_replicate_resolution$detP

save_rds_safe(beta_johnson_f, file.path(johnson_out, "beta_johnson.rds"))
save_rds_safe(pheno_johnson_liver, file.path(johnson_out, "pheno_johnson_clean.rds"))
if (!is.null(detP_johnson)) save_rds_safe(detP_johnson, file.path(johnson_out, "detP_johnson.rds"))

# =========================================================
# 05) Murphy (GSE49542)
# =========================================================

murphy_name <- "Murphy_GSE49542"
murphy_root <- file.path(dir_raw, murphy_name)
murphy_out <- file.path(dir_processed, murphy_name)
ensure_dir(murphy_out)

murphy_signal_reader <- function() {
  f <- file.path(murphy_root, "GSE49542", "GSE49542_unmethyl_methyl_signals.txt.gz")
  check_file_exists(f, "Murphy signal file")
  df <- readr::read_tsv(f, show_col_types = FALSE)

  cn <- colnames(df)
  m_cols <- grep("Methylated Signal$", cn, value = TRUE)
  u_cols <- grep("Unmethylated Signal$", cn, value = TRUE)
  detP_cols <- grep("Detection Pval", cn, value = TRUE)

  probe_ids <- df[[1]]
  U <- as.matrix(df[, u_cols])
  M <- as.matrix(df[, m_cols])
  rownames(U) <- probe_ids
  rownames(M) <- probe_ids

  beta <- signals_to_beta(M, U, offset = 100)
  sample_ids_raw <- sub(" Methylated Signal$", "", m_cols)
  colnames(beta) <- sample_ids_raw

  detP <- NULL
  if (length(detP_cols) > 0) {
    detP <- as.matrix(df[, detP_cols])
    rownames(detP) <- probe_ids
    detP_names <- sub(" Detection Pval$", "", detP_cols)
    colnames(detP) <- detP_names
    if (all(sample_ids_raw %in% detP_names)) detP <- detP[, sample_ids_raw, drop = FALSE]
  }

  list(beta = beta, detP = detP)
}

murphy_obj <- load_beta_idat_or_signal(
  gse_id = "GSE49542",
  raw_root_dir = murphy_root,
  arraytype = "450K",
  signal_reader = murphy_signal_reader
)

beta_murphy_raw <- murphy_obj$beta
detP_murphy_raw <- murphy_obj$detP

gse_murphy <- GEOquery::getGEO("GSE49542", GSEMatrix = TRUE)[[1]]
ph_m <- Biobase::pData(gse_murphy)

if (all(grepl("^GSM\\d+$", colnames(beta_murphy_raw)))) {
  beta_murphy <- beta_murphy_raw
  detP_murphy <- detP_murphy_raw

  pheno_murphy_clean <- tibble::tibble(
    Sample_Name = as.character(ph_m$geo_accession),
    Title = as.character(ph_m$title),
    Diagnosis = as.character(get_col_or_na(ph_m, "diagnosis:ch1")),
    DiseaseState = as.character(get_col_or_na(ph_m, "Stage:ch1")),
    Tissue = as.character(get_col_or_na(ph_m, "tissue:ch1")),
    Age = NA_real_,
    Sex = NA_character_
  ) |>
    dplyr::filter(Sample_Name %in% colnames(beta_murphy)) |>
    dplyr::arrange(match(Sample_Name, colnames(beta_murphy)))

  beta_murphy <- beta_murphy[, pheno_murphy_clean$Sample_Name, drop = FALSE]
} else {
  stopifnot(ncol(beta_murphy_raw) == nrow(ph_m))
  colnames(beta_murphy_raw) <- ph_m$geo_accession
  beta_murphy <- beta_murphy_raw

  detP_murphy <- NULL
  if (!is.null(detP_murphy_raw)) {
    stopifnot(ncol(detP_murphy_raw) == nrow(ph_m))
    colnames(detP_murphy_raw) <- ph_m$geo_accession
    detP_murphy <- detP_murphy_raw
  }

  pheno_murphy_clean <- tibble::tibble(
    Sample_Name = as.character(ph_m$geo_accession),
    Title = as.character(ph_m$title),
    Diagnosis = as.character(get_col_or_na(ph_m, "diagnosis:ch1")),
    DiseaseState = as.character(get_col_or_na(ph_m, "Stage:ch1")),
    Tissue = as.character(get_col_or_na(ph_m, "tissue:ch1")),
    Age = NA_real_,
    Sex = NA_character_
  )
}

stopifnot(all(pheno_murphy_clean$Sample_Name == colnames(beta_murphy)))

beta_murphy_f <- champ_prefilter_beta(beta_murphy, arraytype = "450K", detP = detP_murphy)

murphy_biopsy_id <- stringr::str_extract(pheno_murphy_clean$Title, "\\d+")
murphy_replicate_resolution <- collapse_technical_replicates(
  beta = beta_murphy_f,
  pheno = pheno_murphy_clean,
  group_ids = murphy_biopsy_id,
  dataset_label = "Murphy",
  id_prefix = "Murphy_biopsy",
  detP = detP_murphy
)

beta_murphy_f <- murphy_replicate_resolution$beta
pheno_murphy_clean <- murphy_replicate_resolution$pheno
detP_murphy <- murphy_replicate_resolution$detP

save_rds_safe(beta_murphy_f, file.path(murphy_out, "beta_murphy.rds"))
save_rds_safe(pheno_murphy_clean, file.path(murphy_out, "pheno_murphy_clean.rds"))
if (!is.null(detP_murphy)) save_rds_safe(detP_murphy, file.path(murphy_out, "detP_murphy.rds"))

technical_replicate_audit <- dplyr::bind_rows(
  johnson_replicate_resolution$audit,
  murphy_replicate_resolution$audit
)

write_csv_safe(
  technical_replicate_audit,
  file.path(dir_combined, "technical_replicate_resolution_audit.csv")
)

# =========================================================
# 06) Kurokawa (GSE60753)
# =========================================================

kuro_name <- "Kurokawa_GSE60753"
kuro_root <- file.path(dir_raw, kuro_name)
kuro_out <- file.path(dir_processed, kuro_name)
ensure_dir(kuro_out)

kuro_signal_reader <- function() {
  f <- file.path(kuro_root, "GSE60753", "GSE60753_Unmethylated_and_methylated_signal_intensities.txt.gz")
  check_file_exists(f, "Kurokawa signal file")

  df <- readr::read_tsv(f, show_col_types = FALSE)
  cn <- colnames(df)

  u_cols <- grep("-Unmethylated-Signal$", cn, value = TRUE)
  m_cols <- grep("-Methylated-Signal$", cn, value = TRUE)
  detP_cols <- grep("-Detection-Pval$", cn, value = TRUE)

  stopifnot(length(u_cols) > 0, length(m_cols) > 0, length(u_cols) == length(m_cols))

  sample_ids_raw <- trimws(sub("-Unmethylated-Signal$", "", u_cols))
  probe_ids <- df[["ID_REF"]]
  stopifnot(!is.null(probe_ids))

  U <- as.matrix(df[, u_cols])
  M <- as.matrix(df[, m_cols])
  rownames(U) <- probe_ids
  rownames(M) <- probe_ids

  beta <- signals_to_beta(M, U, offset = 100)
  colnames(beta) <- sample_ids_raw

  detP <- NULL
  if (length(detP_cols) > 0) {
    detP <- as.matrix(df[, detP_cols])
    rownames(detP) <- probe_ids
    colnames(detP) <- sample_ids_raw
  }

  list(beta = beta, detP = detP)
}

kuro_obj <- load_beta_idat_or_signal(
  gse_id = "GSE60753",
  raw_root_dir = kuro_root,
  arraytype = "450K",
  signal_reader = kuro_signal_reader
)

beta_60753_raw <- kuro_obj$beta
detP_60753_raw <- kuro_obj$detP

gse_60753 <- GEOquery::getGEO("GSE60753", GSEMatrix = TRUE)[[1]]
ph_k <- Biobase::pData(gse_60753)

char_cols_k <- grep("^characteristics_ch1", colnames(ph_k), value = TRUE)
stopifnot(length(char_cols_k) > 0)

extract_field_k <- function(x, key) {
  hit <- x[grepl(paste0("^", key, ":"), x, ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  trimws(sub(paste0("^", key, ":[[:space:]]*"), "", hit[1], ignore.case = TRUE))
}

Group <- apply(ph_k[, char_cols_k, drop = FALSE], 1, extract_field_k, key = "disease status")
Tissue <- apply(ph_k[, char_cols_k, drop = FALSE], 1, extract_field_k, key = "tissue")

is_gsm_cols <- all(grepl("^GSM\\d+$", colnames(beta_60753_raw)))

if (is_gsm_cols) {
  beta_60753 <- beta_60753_raw
  detP_60753 <- detP_60753_raw

  pheno_60753_clean <- tibble::tibble(
    Sample_Name = to_chr(ph_k$geo_accession),
    Sample_raw = normalize_sample_id(ph_k$description),
    Tissue = Tissue,
    Group = Group,
    Title = to_chr(ph_k$title)
  ) |>
    dplyr::filter(Sample_Name %in% colnames(beta_60753)) |>
    dplyr::arrange(match(Sample_Name, colnames(beta_60753)))

  beta_60753 <- beta_60753[, pheno_60753_clean$Sample_Name, drop = FALSE]
  if (!is.null(detP_60753)) detP_60753 <- detP_60753[, pheno_60753_clean$Sample_Name, drop = FALSE]
} else {
  stopifnot(ncol(beta_60753_raw) == nrow(ph_k))
  sample_ids_raw <- trimws(as.character(colnames(beta_60753_raw)))

  pheno_60753_clean <- tibble::tibble(
    Sample_Name = to_chr(ph_k$geo_accession),
    Sample_raw = sample_ids_raw,
    Tissue = Tissue,
    Group = Group,
    Title = to_chr(ph_k$title)
  )

  stopifnot(all(trimws(colnames(beta_60753_raw)) == pheno_60753_clean$Sample_raw))
  colnames(beta_60753_raw) <- pheno_60753_clean$Sample_Name
  beta_60753 <- beta_60753_raw

  detP_60753 <- NULL
  if (!is.null(detP_60753_raw)) {
    stopifnot(all(trimws(colnames(detP_60753_raw)) == pheno_60753_clean$Sample_raw))
    colnames(detP_60753_raw) <- pheno_60753_clean$Sample_Name
    detP_60753 <- detP_60753_raw
  }
}

check_file_exists(file_metadata_xlsx, "Kurokawa age/sex Excel file")
age_df_raw <- readxl::read_excel(file_metadata_xlsx, skip = 13)

req_cols <- c("Sample", "Group", "Etiology", "Simple Group", "Age", "Sex")
missing_cols <- setdiff(req_cols, colnames(age_df_raw))
if (length(missing_cols) > 0) {
  stop("Excel age table missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

age_df <- tibble::tibble(
  Sample_raw = normalize_sample_id(age_df_raw[["Sample"]]),
  Group_age = to_chr(age_df_raw[["Group"]]),
  Etiology_age = to_chr(age_df_raw[["Etiology"]]),
  SimpleGroup_age = to_chr(age_df_raw[["Simple Group"]]),
  Age = suppressWarnings(as.numeric(age_df_raw[["Age"]])),
  Sex_excel = to_chr(age_df_raw[["Sex"]])
) |>
  dplyr::filter(!is.na(Sample_raw), Sample_raw != "")

geo_keys <- unique(pheno_60753_clean$Sample_raw)

age_df2 <- age_df |>
  dplyr::mutate(
    Sample_raw_altN = dplyr::if_else(
      grepl("^[0-9]{3,}$", Sample_raw) &
        !(Sample_raw %in% geo_keys) &
        (paste0(Sample_raw, "N") %in% geo_keys),
      paste0(Sample_raw, "N"),
      NA_character_
    )
  )

age_long <- dplyr::bind_rows(
  age_df2 |> dplyr::mutate(Sample_raw_join = Sample_raw),
  age_df2 |> dplyr::filter(!is.na(Sample_raw_altN)) |> dplyr::mutate(Sample_raw_join = Sample_raw_altN)
) |>
  dplyr::select(-Sample_raw_altN) |>
  dplyr::group_by(Sample_raw_join) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

pheno_60753_clean <- pheno_60753_clean |>
  dplyr::left_join(age_long, by = c("Sample_raw" = "Sample_raw_join")) |>
  dplyr::mutate(Sex = normalize_sex(Sex_excel))

pheno_60753_clean <- tibble::as_tibble(pheno_60753_clean) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::where(is.list),
      ~ vapply(.x, function(z) paste(z, collapse = "; "), character(1))
    )
  )

stopifnot(ncol(beta_60753) == nrow(pheno_60753_clean))
stopifnot(setequal(colnames(beta_60753), pheno_60753_clean$Sample_Name))

idx <- match(colnames(beta_60753), pheno_60753_clean$Sample_Name)
stopifnot(!anyNA(idx))
pheno_60753_clean <- pheno_60753_clean[idx, , drop = FALSE]
stopifnot(identical(colnames(beta_60753), pheno_60753_clean$Sample_Name))

pheno_60753_clean <- pheno_60753_clean |>
  dplyr::mutate(
    DiseaseState = dplyr::case_when(
      Group %in% c("Normal", "normal_hepatocyte") ~ "Normal",
      Group %in% c("CC", "Cbil", "CG", "CI", "CirrEtOH", "CirrHBV", "CirrHCV", "Hbil", "HM") ~ "Cirrhosis",
      Group %in% c("HC", "HCCEtOH", "HCCHCV", "HCC cell line") ~ "HCC",
      Group == "unclassified" ~ "Unclassified",
      TRUE ~ NA_character_
    ),
    Etiology = dplyr::case_when(
      Group %in% c("CirrEtOH", "HCCEtOH") ~ "Alcohol",
      Group %in% c("CirrHCV", "HCCHCV") ~ "HCV",
      Group %in% c("CirrHBV") ~ "HBV",
      Group %in% c("CC", "Cbil", "CG", "CI", "Hbil", "HM") ~ "Other/cryptogenic",
      TRUE ~ NA_character_
    )
  )

beta_60753_f <- champ_prefilter_beta(beta_60753, arraytype = "450K", detP = detP_60753)

save_rds_safe(beta_60753_f, file.path(kuro_out, "beta_60753.rds"))
save_rds_safe(pheno_60753_clean, file.path(kuro_out, "pheno_60753_clean.rds"))
if (!is.null(detP_60753)) save_rds_safe(detP_60753, file.path(kuro_out, "detP_60753.rds"))
write_tsv_safe(pheno_60753_clean, file.path(kuro_out, "pheno_60753_clean.tsv"))

# =========================================================
# 07) VanDijck (GSE294806)
# =========================================================

vandijck_name <- "VanDijck_GSE294806"
vandijck_root <- file.path(dir_raw, vandijck_name)
vandijck_out <- file.path(dir_processed, vandijck_name)
ensure_dir(vandijck_out)

vandijck_obj <- load_beta_idat_or_signal(
  gse_id = "GSE294806",
  raw_root_dir = vandijck_root,
  arraytype = "EPIC",
  signal_reader = NULL
)

beta_vandijck <- vandijck_obj$beta
detP_vandijck <- vandijck_obj$detP

gsm_col <- stringr::str_extract(colnames(beta_vandijck), "GSM\\d+")
if (!any(is.na(gsm_col))) {
  colnames(beta_vandijck) <- gsm_col
  if (!is.null(detP_vandijck)) colnames(detP_vandijck) <- gsm_col
}

gse_v <- GEOquery::getGEO("GSE294806", GSEMatrix = TRUE)[[1]]
ph_v <- Biobase::pData(gse_v)
char_cols_v <- grep("^characteristics_ch1", colnames(ph_v), value = TRUE)

source_name_raw_vec <- as.character(ph_v$source_name_ch1)
source_name_fixed_vec <- source_name_raw_vec |>
  stringr::str_replace_all("MASH\\s*F0-F\\s*no\\s*MASLD", "MASH F0-F1") |>
  stringr::str_replace_all("MASH\\s*F0-Fno\\s*MASLD", "MASH F0-F1") |>
  stringr::str_squish()

pheno_vandijck_clean <- tibble::tibble(
  Sample_Name = as.character(ph_v$geo_accession),
  Title = as.character(ph_v$title),
  Dataset = "VanDijck",
  Source_Name_raw = source_name_raw_vec,
  Source_Name = source_name_fixed_vec,
  Age = suppressWarnings(as.numeric(unname(apply(ph_v[, char_cols_v, drop = FALSE], 1, extract_field_ci_vec, key = "age")))),
  BMI = suppressWarnings(as.numeric(unname(apply(ph_v[, char_cols_v, drop = FALSE], 1, extract_field_ci_vec, key = "bmi")))),
  Gender_raw = unname(apply(ph_v[, char_cols_v, drop = FALSE], 1, extract_field_ci_vec, key = "gender")),
  MASLD = suppressWarnings(as.numeric(unname(apply(ph_v[, char_cols_v, drop = FALSE], 1, extract_field_ci_vec, key = "masld")))),
  Fibrosis = unname(apply(ph_v[, char_cols_v, drop = FALSE], 1, extract_field_ci_vec, key = "fibrosis")),
  Steatosis = suppressWarnings(as.numeric(unname(apply(ph_v[, char_cols_v, drop = FALSE], 1, extract_field_ci_vec, key = "steatosis")))),
  Ballooning = suppressWarnings(as.numeric(unname(apply(ph_v[, char_cols_v, drop = FALSE], 1, extract_field_ci_vec, key = "ballooning")))),
  Inflammation = suppressWarnings(as.numeric(unname(apply(ph_v[, char_cols_v, drop = FALSE], 1, extract_field_ci_vec, key = "inflammation"))))
) |>
  dplyr::mutate(
    Sex = dplyr::case_when(
      Gender_raw == "0" ~ "Female",
      Gender_raw == "1" ~ "Male",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(Sample_Name %in% colnames(beta_vandijck)) |>
  dplyr::arrange(match(Sample_Name, colnames(beta_vandijck)))

beta_vandijck <- beta_vandijck[, pheno_vandijck_clean$Sample_Name, drop = FALSE]
if (!is.null(detP_vandijck)) detP_vandijck <- detP_vandijck[, pheno_vandijck_clean$Sample_Name, drop = FALSE]
stopifnot(all(pheno_vandijck_clean$Sample_Name == colnames(beta_vandijck)))

beta_vandijck_f <- champ_prefilter_beta(beta_vandijck, arraytype = "EPIC", detP = detP_vandijck)

save_rds_safe(beta_vandijck_f, file.path(vandijck_out, "beta_vandijck.rds"))
save_rds_safe(pheno_vandijck_clean, file.path(vandijck_out, "pheno_vandijck_clean.rds"))
if (!is.null(detP_vandijck)) save_rds_safe(detP_vandijck, file.path(vandijck_out, "detP_vandijck.rds"))

# =========================================================
# 08) ITEN / K-BDS KSE102917
# =========================================================

kim_name <- "Kim_KSE102917"
kim_root <- file.path(dir_raw, kim_name)
kim_idat_dir <- file.path(kim_root, "idat")
kim_out <- file.path(dir_processed, kim_name)
ensure_dir(kim_out)

kim_metadata_file <- file.path(kim_root, "BioSample_metadata-2.xlsx")
kim_idat_map_file <- file.path(kim_root, "BioSample_idat_map.csv")

kim_available <- file.exists(kim_metadata_file) &&
  file.exists(kim_idat_map_file) &&
  dir_has_idats(kim_idat_dir)

beta_kim_f <- NULL
detP_kim <- NULL
pheno_kim_clean <- NULL

if (kim_available) {
  message("Processing Kim_KSE102917 from KBDS IDATs.")

  kim_xlsx_raw <- readxl::read_excel(
    kim_metadata_file,
    sheet = "Human",
    col_names = FALSE,
    col_types = "text"
  )

  kim_xlsx_names <- as.character(unlist(kim_xlsx_raw[2, ]))
  kim_xlsx_names <- make.unique(kim_xlsx_names, sep = "_")

  kim_biosample <- as.data.frame(kim_xlsx_raw[-c(1:4), ], stringsAsFactors = FALSE)
  colnames(kim_biosample) <- kim_xlsx_names

  kim_map <- readr::read_csv(kim_idat_map_file, show_col_types = FALSE) |>
    dplyr::mutate(
      KBDS_Sample_Name = as.character(KBDS_Sample_Name),
      Basename = as.character(Basename)
    )

  kim_meta <- kim_biosample |>
    dplyr::transmute(
      KBDS_Sample_Name = as.character(.data[["Sample name"]]),
      Age = suppressWarnings(as.numeric(.data[["Age"]])),
      Sex = normalize_sex(.data[["Sex"]]),
      Isolate = stringr::str_squish(as.character(.data[["Isolate"]])),
      Tissue = stringr::str_squish(as.character(.data[["Tissue"]])),
      Biomaterial_provider = as.character(.data[["Biomaterial provider"]]),
      Collection_Date = as.character(.data[["Collection Date"]])
    ) |>
    dplyr::left_join(kim_map, by = "KBDS_Sample_Name") |>
    dplyr::mutate(
      Sample_Name = paste0("Kim_", KBDS_Sample_Name),
      Source_Name = Isolate,
      DiseaseState = Source_Name,
      Group = Isolate,
      Array = "EPIC"
    )

  if (nrow(kim_meta) != 106L) {
    stop("Kim_KSE102917 metadata should contain 106 samples; found ", nrow(kim_meta), ".", call. = FALSE)
  }

  if (any(is.na(kim_meta$Basename)) || any(is.na(kim_meta$Sample_Name))) {
    stop("Kim_KSE102917 metadata is missing Basename or sample IDs.", call. = FALSE)
  }

  kim_missing_idats <- setdiff(
    paste0(rep(kim_meta$Basename, each = 2), c("_Grn.idat", "_Red.idat")),
    list.files(kim_idat_dir, pattern = "\\.idat$", full.names = FALSE)
  )

  if (length(kim_missing_idats) > 0) {
    stop(
      "Kim_KSE102917 IDAT files missing from ",
      kim_idat_dir,
      " (first 10): ",
      paste(head(kim_missing_idats, 10), collapse = ", "),
      call. = FALSE
    )
  }

  kim_obj <- load_beta_idat_or_signal(
    gse_id = "KSE102917",
    raw_root_dir = kim_root,
    idat_dir = kim_idat_dir,
    arraytype = "EPIC",
    signal_reader = NULL
  )

  beta_kim <- kim_obj$beta
  detP_kim <- kim_obj$detP

  kim_beta_key <- stringr::str_extract(colnames(beta_kim), "[0-9]+_R[0-9]{2}C[0-9]{2}")
  if (any(is.na(kim_beta_key))) {
    stop("Could not infer KBDS Sentrix basenames from Kim_KSE102917 IDAT column names.", call. = FALSE)
  }

  kim_mapping <- tibble::tibble(
    beta_col_idx = seq_along(kim_beta_key),
    Basename = kim_beta_key
  ) |>
    dplyr::left_join(kim_meta, by = "Basename")

  if (any(is.na(kim_mapping$Sample_Name))) {
    missing_keys <- paste(kim_mapping$Basename[is.na(kim_mapping$Sample_Name)], collapse = ", ")
    stop("Kim_KSE102917 IDAT basenames not found in metadata: ", missing_keys, call. = FALSE)
  }

  beta_kim <- beta_kim[, kim_mapping$beta_col_idx, drop = FALSE]
  colnames(beta_kim) <- kim_mapping$Sample_Name

  if (!is.null(detP_kim)) {
    detP_kim <- detP_kim[, kim_mapping$beta_col_idx, drop = FALSE]
    colnames(detP_kim) <- kim_mapping$Sample_Name
  }

  pheno_kim_clean <- kim_mapping |>
    dplyr::transmute(
      Sample_Name,
      KBDS_Sample_Name,
      KBDS_Accession_ID,
      BioSample_accession_ID,
      Basename,
      Tissue,
      DiseaseState,
      Source_Name,
      Group,
      Isolate,
      Age,
      Sex,
      Array,
      Biomaterial_provider,
      Collection_Date
    )

  stopifnot(all(pheno_kim_clean$Sample_Name == colnames(beta_kim)))

  beta_kim_f <- champ_prefilter_beta(beta_kim, arraytype = "EPIC", detP = detP_kim)

  save_rds_safe(beta_kim_f, file.path(kim_out, "beta_kim.rds"))
  save_rds_safe(pheno_kim_clean, file.path(kim_out, "pheno_kim_clean.rds"))
  if (!is.null(detP_kim)) save_rds_safe(detP_kim, file.path(kim_out, "detP_kim.rds"))
  write_tsv_safe(pheno_kim_clean, file.path(kim_out, "pheno_kim_clean.tsv"))

  kim_counts <- pheno_kim_clean |>
    dplyr::count(Source_Name, name = "n")

  write_csv_safe(kim_counts, file.path(kim_out, "source_name_counts.csv"))

  message(
    "Kim_KSE102917 processed with Excel Isolate labels: ",
    paste(names(table(pheno_kim_clean$Isolate)), as.integer(table(pheno_kim_clean$Isolate)), sep = "=", collapse = ", ")
  )
} else {
  message(
    "Skipping Kim_KSE102917: place 212 KBDS IDAT files in ",
    kim_idat_dir,
    " and keep ",
    kim_metadata_file,
    " plus ",
    kim_idat_map_file,
    " to enable processing."
  )
}

# =========================================================
# 09) Combine datasets (notebook order preserved)
# =========================================================

beta_list <- list(
  Ahrens = beta_ahrens_f,
  Horvath = beta_horvath_f,
  Johnson = beta_johnson_f,
  Murphy = beta_murphy_f,
  VanDijck = beta_vandijck_f
)

if (!is.null(beta_kim_f)) {
  beta_list$ITEN <- beta_kim_f
}

pheno_ahrens_clean$Dataset <- "Ahrens"
pheno_horvath_clean$Dataset <- "Horvath"
pheno_johnson_liver$Dataset <- "Johnson"
pheno_murphy_clean$Dataset <- "Murphy"
pheno_vandijck_clean$Dataset <- "VanDijck"
if (!is.null(pheno_kim_clean)) pheno_kim_clean$Dataset <- "ITEN"

pheno_ahrens_clean$Array <- "450K"
pheno_horvath_clean$Array <- "450K"
pheno_johnson_liver$Array <- "EPIC"
pheno_murphy_clean$Array <- "450K"
pheno_vandijck_clean$Array <- "EPIC"
if (!is.null(pheno_kim_clean)) pheno_kim_clean$Array <- "EPIC"

coerce_for_bind <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  needed_cols <- c(
    "Sample_Name", "Sex", "Fibrosis", "Steatosis", "DiseaseState",
    "Group", "Title", "Diagnosis", "Source_Name", "Source_Name_raw", "Age", "Array"
  )

  for (nm in needed_cols) {
    if (!nm %in% names(df)) df[[nm]] <- NA
  }

  df$Sample_Name <- as.character(df$Sample_Name)
  df$Sex <- as.character(df$Sex)
  df$Fibrosis <- as.character(df$Fibrosis)
  df$Steatosis <- as.character(df$Steatosis)
  df$DiseaseState <- as.character(df$DiseaseState)
  df$Group <- as.character(df$Group)
  df$Title <- as.character(df$Title)
  df$Diagnosis <- as.character(df$Diagnosis)
  df$Source_Name <- as.character(df$Source_Name)
  df$Source_Name_raw <- as.character(df$Source_Name_raw)
  df$Age <- suppressWarnings(as.numeric(df$Age))

  df
}

pheno_ahrens_clean <- coerce_for_bind(pheno_ahrens_clean)
pheno_horvath_clean <- coerce_for_bind(pheno_horvath_clean)
pheno_johnson_liver <- coerce_for_bind(pheno_johnson_liver)
pheno_murphy_clean <- coerce_for_bind(pheno_murphy_clean)
pheno_vandijck_clean <- coerce_for_bind(pheno_vandijck_clean)
if (!is.null(pheno_kim_clean)) pheno_kim_clean <- coerce_for_bind(pheno_kim_clean)

pheno_all_wide <- dplyr::bind_rows(
  pheno_ahrens_clean,
  pheno_horvath_clean,
  pheno_johnson_liver,
  pheno_murphy_clean,
  pheno_vandijck_clean,
  pheno_kim_clean
) |>
  dplyr::mutate(
    Sample_Name = as.character(Sample_Name),
    Dataset = as.factor(Dataset),
    Sex = normalize_sex(Sex),
    Age = suppressWarnings(as.numeric(Age))
  ) |>
  dplyr::select(-dplyr::any_of(c("Sample_raw.y", "Sex_excel")))

is_unknown_like <- function(x) {
  x2 <- stringr::str_to_lower(stringr::str_trim(dplyr::coalesce(as.character(x), "")))
  x2 %in% c("", "na", "n/a", "nan", "unknown", "unclassified", "other", "undef", "undefined", "not classified") |
    stringr::str_detect(x2, "^unknown\\b|^unclassified\\b|^other\\b|^na\\b")
}

to_num <- function(x) {
  suppressWarnings(as.numeric(stringr::str_extract(as.character(x), "\\d+(?:\\.\\d+)?")))
}

pheno_all_wide <- pheno_all_wide |>
  dplyr::mutate(
    DiseaseState = dplyr::if_else(
      Dataset == "VanDijck" & (is.na(DiseaseState) | DiseaseState == ""),
      as.character(Source_Name),
      as.character(DiseaseState)
    ),
    Group = dplyr::if_else(
      Dataset == "VanDijck" & (is.na(Group) | Group == ""),
      as.character(Source_Name),
      as.character(Group)
    )
  ) |>
  dplyr::mutate(
    DiseaseStateUnified = dplyr::coalesce(as.character(DiseaseState), as.character(Group)) |>
      stringr::str_squish(),
    .dx_text = dplyr::coalesce(
      as.character(DiseaseState),
      as.character(DiseaseStateUnified),
      as.character(Group),
      as.character(Diagnosis),
      as.character(Title),
      ""
    ) |>
      stringr::str_squish(),
    .dx_low = stringr::str_to_lower(.dx_text),
    source_name = dplyr::coalesce(as.character(Source_Name), as.character(Source_Name_raw), "") |>
      stringr::str_replace_all("MASH\\s*F0-F\\s*no\\s*MASLD", "MASH F0-F1") |>
      stringr::str_replace_all("MASH\\s*F0-Fno\\s*MASLD", "MASH F0-F1") |>
      stringr::str_squish(),
    source_low = stringr::str_to_lower(source_name),
    VanDijck_Group4 = dplyr::if_else(Dataset == "VanDijck", source_name, NA_character_),
    is_unclassified = is_unknown_like(.dx_text) |
      stringr::str_detect(.dx_low, "unclassified|not classified|unknown"),
    is_HCC = stringr::str_detect(.dx_low, "hcc|hepatocellular|hepatoma|carcinoma"),
    is_PSC_PBC = stringr::str_detect(.dx_low, "\\bpsc\\b|\\bpbc\\b|primary sclerosing|primary biliary"),
    Fibrosis_num = to_num(Fibrosis),
    Steatosis_num = to_num(Steatosis),
    DiseaseGroup = dplyr::case_when(
      is_unclassified ~ NA_character_,
      is_HCC ~ NA_character_,
      is_PSC_PBC ~ NA_character_,

      Dataset == "Ahrens" & stringr::str_detect(.dx_low, "healthy\\s*obese|healthyobese") ~ "Healthy_Obese",
      Dataset == "Ahrens" & stringr::str_detect(.dx_low, "\\bnormal\\b|\\bcontrol\\b|\\bhealthy\\b") ~ "Healthy",
      Dataset == "Ahrens" & stringr::str_detect(.dx_low, "\\bnash\\b|steatohepatitis") ~ "MASH",
      Dataset == "Ahrens" & stringr::str_detect(.dx_low, "\\bnafld\\b|fatty liver|steatosis|\\bmasl\\b") ~ "MASL",

      Dataset == "Horvath" & !is.na(Fibrosis_num) & Fibrosis_num >= 3 ~ "Advanced_Fibrosis",
      Dataset == "Horvath" & !is.na(Fibrosis_num) & Fibrosis_num > 0 ~ "Mild_Fibrosis",
      Dataset == "Horvath" & stringr::str_detect(.dx_low, "\\bnash\\b|steatohepatitis") ~ "MASH",
      Dataset == "Horvath" & !is.na(Fibrosis_num) & Fibrosis_num == 0 &
        (!is.na(Steatosis_num) & Steatosis_num > 0) ~ "MASL",
      Dataset == "Horvath" & stringr::str_detect(.dx_low, "\\bnafld\\b|fatty liver|steatosis|\\bmasl\\b") ~ "MASL",
      Dataset == "Horvath" & stringr::str_detect(.dx_low, "healthy\\s*obese|healthyobese") ~ "Healthy_Obese",
      Dataset == "Horvath" & stringr::str_detect(.dx_low, "\\bnormal\\b|\\bcontrol\\b|\\bhealthy\\b") ~ "Healthy",

      Dataset == "Johnson" & stringr::str_detect(.dx_low, "grade\\s*0|grade0") ~ "Healthy_Obese",
      Dataset == "Johnson" & stringr::str_detect(.dx_low, "grade\\s*34|grade34|grade\\s*4|grade4|grade\\s*3|grade3") ~ "Advanced_Fibrosis",

      Dataset == "Murphy" & stringr::str_detect(.dx_low, "\\badvanced\\b") ~ "Advanced_Fibrosis",
      Dataset == "Murphy" & stringr::str_detect(.dx_low, "\\bmild\\b") ~ "Mild_Fibrosis",

      Dataset == "VanDijck" & stringr::str_detect(source_low, "no\\s*masld") ~ "Healthy",
      Dataset == "VanDijck" & stringr::str_detect(source_low, "^masl\\b") ~ "MASL",
      # VanDijck MASH F2-F4 is treated as advanced/significant fibrosis;
      # this label may include F2 as well as F3/F4.
      Dataset == "VanDijck" & stringr::str_detect(source_low, "mash\\s*f2\\s*-\\s*f4|mash f2-f4") ~ "Advanced_Fibrosis",
      Dataset == "VanDijck" & !is.na(Fibrosis_num) & Fibrosis_num >= 3 ~ "Advanced_Fibrosis",
      Dataset == "VanDijck" & !is.na(Fibrosis_num) & Fibrosis_num > 0 ~ "Mild_Fibrosis",
      Dataset == "VanDijck" & stringr::str_detect(source_low, "^mash\\b") ~ "MASH",
      Dataset == "VanDijck" & !is.na(Fibrosis_num) & Fibrosis_num == 0 &
        (!is.na(Steatosis_num) & Steatosis_num > 0) ~ "MASL",

      Dataset == "ITEN" & stringr::str_detect(source_low, "^control$") ~ "Healthy",
      Dataset == "ITEN" & stringr::str_detect(source_low, "^masl$") ~ "MASL",
      Dataset == "ITEN" & stringr::str_detect(source_low, "^mash$") ~ "MASH",

      stringr::str_detect(.dx_low, "healthy\\s*obese|healthyobese") ~ "Healthy_Obese",
      stringr::str_detect(.dx_low, "\\bnormal\\b|\\bcontrol\\b|\\bhealthy\\b") ~ "Healthy",
      stringr::str_detect(.dx_low, "\\bnash\\b|steatohepatitis") ~ "MASH",
      stringr::str_detect(.dx_low, "\\bnafld\\b|fatty liver|steatosis|\\bmasl\\b") ~ "MASL",
      TRUE ~ NA_character_
    ),
    DiseaseGroup = factor(DiseaseGroup, levels = stage_levels),
    Progression3 = dplyr::case_when(
      DiseaseGroup %in% c("Healthy", "Healthy_Obese") ~ "Healthy",
      DiseaseGroup %in% c("MASL", "MASH", "Mild_Fibrosis") ~ "Mild",
      DiseaseGroup %in% c("Advanced_Fibrosis") ~ "Advanced",
      TRUE ~ NA_character_
    ),
    Progression3 = factor(Progression3, levels = prog3_levels)
  ) |>
  dplyr::select(-.dx_text, -.dx_low, -source_low)

pheno_all <- pheno_all_wide |>
  dplyr::transmute(
    Sample_Name = as.character(Sample_Name),
    DiseaseState = as.character(DiseaseState),
    Age = suppressWarnings(as.numeric(Age)),
    Sex = normalize_sex(Sex),
    Dataset = as.factor(Dataset),
    Array = as.factor(Array),
    DiseaseGroup = DiseaseGroup,
    Progression3 = Progression3
  ) |>
  dplyr::filter(!is.na(Progression3)) |>
  dplyr::filter(!is.na(DiseaseGroup)) |>
  droplevels()

common_cpgs <- Reduce(intersect, lapply(beta_list, rownames))
common_cpgs <- sort(common_cpgs)

beta_all <- do.call(
  cbind,
  lapply(beta_list, function(b) b[common_cpgs, , drop = FALSE])
)

n_samples_beta_before_final_intersect <- ncol(beta_all)
n_samples_pheno_before_final_intersect <- nrow(pheno_all)

keep_samples <- intersect(colnames(beta_all), pheno_all$Sample_Name)
beta_all <- beta_all[, keep_samples, drop = FALSE]

pheno_all <- pheno_all |>
  dplyr::filter(Sample_Name %in% colnames(beta_all)) |>
  dplyr::arrange(match(Sample_Name, colnames(beta_all))) |>
  droplevels()

n_samples_final <- nrow(pheno_all)

stopifnot(identical(pheno_all$Sample_Name, colnames(beta_all)))
stopifnot(identical(rownames(beta_all), common_cpgs))

# Final grouping audit files for classifier outcome validation.
sample_grouping_audit <- pheno_all_wide |>
  dplyr::mutate(
    Include_in_pheno_all = Sample_Name %in% pheno_all$Sample_Name,
    Reason_or_rule = dplyr::case_when(
      is_unclassified ~ "Excluded: unknown/unclassified",
      is_HCC ~ "Excluded: HCC/hepatocellular carcinoma",
      is_PSC_PBC ~ "Excluded: PBC/PSC autoimmune liver disease",
      !Sample_Name %in% colnames(beta_all) ~ "Excluded: no matched beta sample after final intersection",
      is.na(DiseaseGroup) ~ "Excluded: no final DiseaseGroup rule matched",
      Dataset == "Horvath" & DiseaseGroup == "Advanced_Fibrosis" &
        !is.na(Fibrosis_num) & Fibrosis_num >= 3 ~ "Horvath fibrosis >= 3 -> Advanced_Fibrosis",
      Dataset == "Horvath" & DiseaseGroup == "Mild_Fibrosis" &
        !is.na(Fibrosis_num) & Fibrosis_num > 0 ~ "Horvath fibrosis > 0 and < 3 -> Mild_Fibrosis",
      Dataset == "Horvath" & DiseaseGroup == "MASH" &
        stringr::str_detect(stringr::str_to_lower(dplyr::coalesce(as.character(DiseaseState), "")), "\\bnash\\b|steatohepatitis") ~
        "Horvath NASH/steatohepatitis without fibrosis criterion -> MASH",
      Dataset == "VanDijck" & DiseaseGroup == "Advanced_Fibrosis" &
        stringr::str_detect(stringr::str_to_lower(dplyr::coalesce(as.character(source_name), "")), "mash\\s*f2\\s*-\\s*f4|mash f2-f4") ~
        "VanDijck MASH F2-F4 -> Advanced_Fibrosis (significant/advanced fibrosis label)",
      Dataset == "VanDijck" & DiseaseGroup == "Mild_Fibrosis" &
        !is.na(Fibrosis_num) & Fibrosis_num > 0 ~ "VanDijck fibrosis > 0 and not F2-F4 advanced label -> Mild_Fibrosis",
      Dataset == "VanDijck" & DiseaseGroup == "MASH" &
        stringr::str_detect(stringr::str_to_lower(dplyr::coalesce(as.character(source_name), "")), "^mash\\b") ~
        "VanDijck MASH label without fibrosis criterion -> MASH",
      Dataset == "ITEN" & DiseaseGroup == "Healthy" &
        stringr::str_detect(stringr::str_to_lower(dplyr::coalesce(as.character(source_name), "")), "^control$") ~
        "ITEN Isolate Control -> Healthy",
      Dataset == "ITEN" & DiseaseGroup == "MASL" &
        stringr::str_detect(stringr::str_to_lower(dplyr::coalesce(as.character(source_name), "")), "^masl$") ~
        "ITEN Isolate MASL -> MASL",
      Dataset == "ITEN" & DiseaseGroup == "MASH" &
        stringr::str_detect(stringr::str_to_lower(dplyr::coalesce(as.character(source_name), "")), "^mash$") ~
        "ITEN Isolate MASH -> MASH",
      DiseaseGroup == "Healthy_Obese" ~ "Healthy obese/control label -> Healthy_Obese",
      DiseaseGroup == "Healthy" ~ "Healthy/control/normal/no MASLD label -> Healthy",
      DiseaseGroup == "MASL" ~ "MASL/NAFLD/steatosis label without MASH/fibrosis criterion -> MASL",
      DiseaseGroup == "MASH" ~ "NASH/MASH/steatohepatitis label without fibrosis criterion -> MASH",
      DiseaseGroup == "Mild_Fibrosis" ~ "Fibrosis > 0 and not advanced -> Mild_Fibrosis",
      DiseaseGroup == "Advanced_Fibrosis" ~ "Advanced fibrosis rule -> Advanced_Fibrosis",
      TRUE ~ "Unspecified final grouping rule"
    ),
    DiseaseGroup_current = as.character(DiseaseGroup),
    Progression3_current = as.character(Progression3),
    NAFLDActivityScore = suppressWarnings(as.numeric(.data$NAFLDActivityScore)),
    Inflammation = suppressWarnings(as.numeric(.data$Inflammation)),
    Ballooning = dplyr::coalesce(
      suppressWarnings(as.numeric(.data$Ballooning)),
      suppressWarnings(as.numeric(.data$NAS_Ballooning)),
      NA_real_
    )
  ) |>
  dplyr::select(
    Sample_Name,
    Dataset,
    Array,
    DiseaseState,
    Group,
    Source_Name,
    Source_Name_raw,
    Diagnosis,
    Fibrosis,
    Fibrosis_num,
    Steatosis,
    Steatosis_num,
    NAFLDActivityScore,
    Inflammation,
    Ballooning,
    DiseaseGroup_current,
    Progression3_current,
    Include_in_pheno_all,
    Reason_or_rule
  ) |>
  dplyr::arrange(Dataset, Sample_Name)

write_csv_safe(
  sample_grouping_audit,
  file.path(dir_combined, "sample_grouping_audit_final_check.csv")
)

sample_grouping_counts <- dplyr::bind_rows(
  pheno_all |>
    dplyr::mutate(Dataset = as.character(Dataset), DiseaseGroup = as.character(DiseaseGroup)) |>
    dplyr::count(Dataset, DiseaseGroup, name = "n") |>
    dplyr::mutate(count_type = "Dataset x DiseaseGroup", Array = NA_character_, Progression3 = NA) |>
    dplyr::select(count_type, Dataset, Array, DiseaseGroup, Progression3, n),
  pheno_all |>
    dplyr::mutate(Dataset = as.character(Dataset), Progression3 = as.character(Progression3)) |>
    dplyr::count(Dataset, Progression3, name = "n") |>
    dplyr::mutate(count_type = "Dataset x Progression3", Array = NA_character_, DiseaseGroup = NA) |>
    dplyr::select(count_type, Dataset, Array, DiseaseGroup, Progression3, n),
  pheno_all |>
    dplyr::mutate(Array = as.character(Array), DiseaseGroup = as.character(DiseaseGroup)) |>
    dplyr::count(Array, DiseaseGroup, name = "n") |>
    dplyr::mutate(count_type = "Array x DiseaseGroup", Dataset = NA, Progression3 = NA) |>
    dplyr::select(count_type, Dataset, Array, DiseaseGroup, Progression3, n)
)

write_csv_safe(
  sample_grouping_counts,
  file.path(dir_combined, "sample_grouping_counts_final_check.csv")
)

# =========================================================
# 09) Imputation + ComBat
# =========================================================

pd_final <- data.frame(
  Sample_Name = pheno_all$Sample_Name,
  row.names = pheno_all$Sample_Name,
  stringsAsFactors = FALSE
)

stopifnot(identical(pd_final$Sample_Name, colnames(beta_all)))

na_before <- mean(is.na(beta_all))
imp <- ChAMP::champ.impute(beta = beta_all, pd = pd_final, method = "KNN")
beta_all_imp <- as.matrix(imp$beta)
na_after <- mean(is.na(beta_all_imp))

pd_combat <- pheno_all |>
  dplyr::mutate(
    Sample_Name = as.character(Sample_Name),
    Dataset = droplevels(as.factor(Dataset)),
    DiseaseGroup = droplevels(as.factor(DiseaseGroup)),
    Progression3 = droplevels(as.factor(Progression3))
  ) |>
  dplyr::arrange(match(Sample_Name, colnames(beta_all_imp)))

stopifnot(identical(pd_combat$Sample_Name, colnames(beta_all_imp)))

prog_by_dataset <- table(pd_combat$Dataset, pd_combat$Progression3)
prog_present_in_batches <- colSums(prog_by_dataset > 0)
can_protect_progression3 <- all(prog_present_in_batches >= 2)

combat_ok <- FALSE
beta_all_combat <- NULL
combat_mode <- "unknown"

if (can_protect_progression3) {
  combat_try <- try(
    ChAMP::champ.runCombat(
      beta = beta_all_imp,
      pd = as.data.frame(pd_combat),
      batchname = c("Dataset"),
      variablename = "Progression3",
      logitTrans = TRUE
    ),
    silent = TRUE
  )

  if (!inherits(combat_try, "try-error")) {
    beta_all_combat <- if (is.list(combat_try) && "beta" %in% names(combat_try)) {
      as.matrix(combat_try$beta)
    } else {
      as.matrix(combat_try)
    }
    combat_ok <- TRUE
    combat_mode <- "protected_Progression3"
  }
}

if (!combat_ok) {
  combat_try2 <- try(
    ChAMP::champ.runCombat(
      beta = beta_all_imp,
      pd = as.data.frame(pd_combat),
      batchname = c("Dataset"),
      variablename = NULL,
      logitTrans = TRUE
    ),
    silent = TRUE
  )

  if (inherits(combat_try2, "try-error")) {
    stop("ComBat failed both with and without protected covariate.", call. = FALSE)
  }

  beta_all_combat <- if (is.list(combat_try2) && "beta" %in% names(combat_try2)) {
    as.matrix(combat_try2$beta)
  } else {
    as.matrix(combat_try2)
  }
  combat_mode <- "fallback_unprotected"
}

stopifnot(identical(colnames(beta_all_combat), pd_combat$Sample_Name))

# =========================================================
# 10) Save canonical outputs + compatibility checks
# =========================================================

save_rds_safe(beta_all_imp, file_beta_raw)
save_rds_safe(beta_all_combat, file_beta_combat)
save_rds_safe(pheno_all, file_pheno_all)

# Keep legacy filtered reference (beta_all only)
save_rds_safe(beta_all, file_beta_filtered)

check_file_exists(file_beta_raw, "beta raw (imputed)")
check_file_exists(file_beta_combat, "beta combat")
check_file_exists(file_pheno_all, "combined pheno")

beta_raw_chk <- readRDS(file_beta_raw)
beta_cb_chk <- readRDS(file_beta_combat)
pheno_chk <- readRDS(file_pheno_all)

check_object_columns(
  pheno_chk,
  c("Sample_Name", "Dataset", "DiseaseGroup", "Progression3", "Age", "Sex"),
  "pheno_all_liver_filtered"
)

stopifnot(!anyDuplicated(pheno_chk$Sample_Name))
stopifnot(identical(colnames(beta_raw_chk), pheno_chk$Sample_Name))
stopifnot(identical(colnames(beta_cb_chk), pheno_chk$Sample_Name))

# =========================================================
# 11) Reporting outputs
# =========================================================

summary_dims <- tibble::tibble(
  object = c("beta_all_pre_imputation", "beta_all_imputed_raw", "beta_all_combat", "pheno_all"),
  n_rows = c(nrow(beta_all), nrow(beta_all_imp), nrow(beta_all_combat), nrow(pheno_all)),
  n_cols = c(ncol(beta_all), ncol(beta_all_imp), ncol(beta_all_combat), ncol(pheno_all))
)

write_csv_safe(summary_dims, file.path(dir_out_preprocessing, "combined_liver_object_dimensions.csv"))

counts_dataset <- pheno_all |> dplyr::count(Dataset, name = "n")
counts_disease <- pheno_all |> dplyr::count(DiseaseGroup, name = "n")
counts_prog3 <- pheno_all |> dplyr::count(Progression3, name = "n")

write_csv_safe(counts_dataset, file.path(dir_out_preprocessing, "pheno_dataset_counts.csv"))
write_csv_safe(counts_disease, file.path(dir_out_preprocessing, "pheno_diseasegroup_counts.csv"))
write_csv_safe(counts_prog3, file.path(dir_out_preprocessing, "pheno_progression3_counts.csv"))

xtab_dataset_disease <- pheno_all |>
  dplyr::count(Dataset, DiseaseGroup, name = "n")

xtab_dataset_prog3 <- pheno_all |>
  dplyr::count(Dataset, Progression3, name = "n")

write_csv_safe(xtab_dataset_disease, file.path(dir_out_preprocessing, "xtab_dataset_diseasegroup.csv"))
write_csv_safe(xtab_dataset_prog3, file.path(dir_out_preprocessing, "xtab_dataset_progression3.csv"))

preprocessing_manifest <- tibble::tibble(
  metric = c(
    "n_common_cpgs",
    "n_johnson_profiles_before_replicate_resolution",
    "n_johnson_biological_samples_after_replicate_resolution",
    "n_murphy_profiles_before_replicate_resolution",
    "n_murphy_biological_samples_after_replicate_resolution",
    "n_technical_profiles_removed_by_replicate_resolution",
    "n_samples_beta_before_final_intersect",
    "n_samples_pheno_before_final_intersect",
    "n_samples_final",
    "combat_mode",
    "legacy_beta_filtered_written"
  ),
  value = c(
    length(common_cpgs),
    sum(technical_replicate_audit$n_technical_profiles[technical_replicate_audit$Dataset == "Johnson"]),
    sum(technical_replicate_audit$Dataset == "Johnson"),
    sum(technical_replicate_audit$n_technical_profiles[technical_replicate_audit$Dataset == "Murphy"]),
    sum(technical_replicate_audit$Dataset == "Murphy"),
    sum(technical_replicate_audit$n_technical_profiles - 1),
    n_samples_beta_before_final_intersect,
    n_samples_pheno_before_final_intersect,
    n_samples_final,
    combat_mode,
    TRUE
  )
)

write_csv_safe(preprocessing_manifest, file.path(dir_out_preprocessing, "preprocessing_manifest.csv"))

combat_qc <- tibble::tibble(
  metric = c(
    "na_fraction_before_imputation",
    "na_fraction_after_imputation",
    "na_fraction_after_combat",
    "beta_combat_min",
    "beta_combat_max"
  ),
  value = c(
    na_before,
    na_after,
    mean(is.na(beta_all_combat)),
    min(beta_all_combat, na.rm = TRUE),
    max(beta_all_combat, na.rm = TRUE)
  )
)

write_csv_safe(combat_qc, file.path(dir_out_preprocessing, "imputation_combat_qc.csv"))

plot_df <- dplyr::bind_rows(
  data.frame(beta = as.numeric(beta_all_imp), Version = "Raw"),
  data.frame(beta = as.numeric(beta_all_combat), Version = "ComBat")
) |>
  dplyr::filter(!is.na(beta))

p_compare <- ggplot(plot_df, aes(x = beta)) +
  geom_histogram(bins = 100) +
  facet_wrap(~Version, scales = "free_y", ncol = 1) +
  labs(
    title = "Beta value distributions",
    x = "Beta value",
    y = "Count"
  ) +
  theme_pipeline(11)

save_plot(
  p_compare,
  file.path(dir_out_preprocessing, "beta_histogram_raw_vs_combat.png"),
  width = 8,
  height = 8,
  dpi = 300
)

save_plot(
  p_compare,
  file.path(dir_combined, "beta_histogram_raw_vs_combat.png"),
  width = 8,
  height = 8,
  dpi = 300
)
