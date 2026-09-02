find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    has_config <- file.exists(file.path(current, "data_preprocessing", "00_config.R"))
    has_classifier <- file.exists(file.path(current, "classifier", "config", "classifier_config.yaml"))
    if (has_config && has_classifier) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop(
        "Could not locate project root. Run pipeline commands from the repository root ",
        "or set PROJECT_ROOT.",
        call. = FALSE
      )
    }
    current <- parent
  }
}

project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) {
  project_root <- find_project_root()
}
project_root <- normalizePath(project_root, mustWork = TRUE)

# ---------------------------
# Core directories
# ---------------------------
dir_data      <- file.path(project_root, "data")
dir_raw       <- file.path(dir_data, "raw")
dir_metadata  <- file.path(dir_data, "metadata")
dir_processed <- file.path(dir_data, "processed")

dir_preprocessing <- file.path(project_root, "data_preprocessing")
dir_helpers   <- file.path(project_root, "R")
dir_notebooks <- file.path(project_root, "notebooks")
dir_archive   <- file.path(project_root, "archive")
dir_results   <- file.path(project_root, "results")
dir_outputs   <- file.path(project_root, "outputs")

# ---------------------------
# Main output directories
# ---------------------------
dir_out_preprocessing <- file.path(dir_outputs, "preprocessing")
dir_out_pca           <- file.path(dir_outputs, "pca")
dir_out_epi           <- file.path(dir_outputs, "epi_age")

dir_out_pca_figures   <- file.path(dir_out_pca, "figures")
dir_out_pca_tables    <- file.path(dir_out_pca, "tables")
dir_out_pca_rds       <- file.path(dir_out_pca, "rds")

dir_out_epi_plots     <- file.path(dir_out_epi, "plots")
dir_out_epi_tables    <- file.path(dir_out_epi, "tables")
dir_out_epi_rds       <- file.path(dir_out_epi, "rds")

all_output_dirs <- c(
  dir_out_preprocessing,
  dir_out_pca,
  dir_out_pca_figures,
  dir_out_pca_tables,
  dir_out_pca_rds,
  dir_out_epi,
  dir_out_epi_plots,
  dir_out_epi_tables,
  dir_out_epi_rds
)

# ---------------------------
# Input data
# ---------------------------
file_metadata_xlsx <- file.path(dir_metadata, "S_table1_27_12_2024.xlsx")

dir_combined <- file.path(dir_processed, "Combined_liver")

file_beta_raw <- file.path(
  dir_combined,
  "beta_all_liver_commonCpGs_filtered_imputed_raw.rds"
)

file_beta_combat <- file.path(
  dir_combined,
  "beta_all_liver_commonCpGs_filtered_imputed_combat.rds"
)

file_beta_filtered <- file.path(
  dir_combined,
  "beta_all_liver_commonCpGs_filtered.rds"
)

file_beta_all <- file.path(
  dir_combined,
  "beta_all_liver_commonCpGs.rds"
)

file_pheno_all <- file.path(
  dir_combined,
  "pheno_all_liver_filtered.rds"
)

file_pca_scores <- file.path(
  dir_combined,
  "pca_scores_by_dataset_top50k.rds"
)

file_pca_df <- file.path(
  dir_combined,
  "pca_df_raw_vs_combat_top50k.rds"
)

# ---------------------------
# Shared analysis options
# ---------------------------
pca_top_n <- 50000L

use_grimage <- TRUE
use_systems_age <- FALSE
use_pc_clocks <- TRUE

# ---------------------------
# Group definitions
# ---------------------------
stage_colors <- c(
  "Healthy"           = "#5B9E6B",
  "Healthy_Obese"     = "#A8D08D",
  "MASL"              = "#F3E55C",
  "MASH"              = "#F3B35E",
  "Mild_Fibrosis"    = "#E7873C",
  "Advanced_Fibrosis" = "#B85C5C"
)

stage_levels <- names(stage_colors)

prog3_colors <- c(
  "Healthy"  = "#5B9E6B",
  "Mild"     = "#E7873C",
  "Advanced" = "#B85C5C"
)

prog3_levels <- names(prog3_colors)

shape_dataset_map <- c(
  "Ahrens"   = 15,
  "Horvath"  = 16,
  "Johnson"  = 17,
  "Murphy"   = 7,
  "VanDijck" = 18,
  "ITEN"      = 8
)

dataset_colors <- c(
  "Ahrens"   = "#1B9E77",
  "Horvath"  = "#D95F02",
  "Johnson"  = "#7570B3",
  "Murphy"   = "#E7298A",
  "VanDijck" = "#66A61E",
  "ITEN"     = "#A6761D"
)

dataset_levels <- names(dataset_colors)

main_clock_cols <- c(
  "Horvath1",
  "Horvath2",
  "Hannum",
  "PhenoAge",
  "GrimAgeV1",
  "PCHannum",
  "PCHorvath1",
  "PCHorvath2",
  "PCPhenoAge",
  "PCGrimAge",
  "HepClock",
  "LiverClock"
)

array_colors <- c("450K" = "#1F77B4", "EPIC" = "#FF7F0E")
array_shapes <- c("450K" = 16, "EPIC" = 17)
