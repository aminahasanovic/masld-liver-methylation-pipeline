#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Install the R packages required by this pipeline.
#
# The thesis results were produced with R 4.5.0 and Bioconductor 3.22; the
# versions of every package in that environment are listed in docs/SOFTWARE.md.
# This script installs whatever is missing from the current library and reports
# the resulting versions, so a fresh environment can be compared against the
# documented one.
#
#   Rscript workflow/install_dependencies.R
#
# Two epigenetic-clock packages are not on CRAN or Bioconductor and are
# installed from GitHub at the exact commits used for the thesis.
# ---------------------------------------------------------------------------

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "BiocManager", "remotes",
  "broom", "caret", "dplyr", "forcats", "fs", "GGally", "ggplot2",
  "ggrepel", "glmnet", "glue", "janitor", "matrixStats", "pROC", "patchwork",
  "png", "purrr", "qs2", "R.utils", "ranger", "readr", "readxl", "rlang",
  "Rtsne", "scales", "stringr", "tibble", "tidyr", "uwot", "viridis",
  "xgboost", "yaml", "yardstick"
)

bioc_pkgs <- c(
  "AnnotationDbi", "Biobase", "ChAMP", "GEOquery",
  "GenomicFeatures", "GenomicRanges", "HiBED", "impute", "IRanges", "limma",
  "minfi", "org.Hs.eg.db", "S4Vectors", "sva", "wateRmelon",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "TxDb.Hsapiens.UCSC.hg19.knownGene"
)

# Clock packages, pinned to the commits used for the thesis results.
github_pkgs <- c(
  "danbelsky/DunedinPoAm38",
  "MorganLevineLab/prcPhenoAge",
  "HigginsChenLab/methylCIPHER@73f2d7dc7701c4266bd35f91b3b61921996fa9a0",
  "HGT-UwU/CTSclocks@7c242cf66cee48041f14ebf512cf31cb07a4a7b4"
)

installed <- function(p) requireNamespace(p, quietly = TRUE)

missing_cran <- cran_pkgs[!vapply(cran_pkgs, installed, logical(1))]
if (length(missing_cran)) {
  message("Installing from CRAN: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
} else {
  message("All CRAN packages present")
}

if (!installed("BiocManager")) {
  stop("BiocManager could not be installed; cannot continue.", call. = FALSE)
}

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, installed, logical(1))]
if (length(missing_bioc)) {
  message("Installing from Bioconductor: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
} else {
  message("All Bioconductor packages present")
}

github_name <- function(spec) sub("@.*$", "", basename(spec))
missing_github <- github_pkgs[!vapply(vapply(github_pkgs, github_name, character(1)),
                                     installed, logical(1))]
if (length(missing_github)) {
  message("Installing from GitHub: ", paste(missing_github, collapse = ", "))
  remotes::install_github(missing_github, upgrade = "never")
} else {
  message("All GitHub packages present")
}

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
all_pkgs <- sort(unique(c(cran_pkgs, bioc_pkgs,
                          vapply(github_pkgs, github_name, character(1)))))
versions <- vapply(
  all_pkgs,
  function(p) tryCatch(as.character(utils::packageVersion(p)),
                       error = function(e) NA_character_),
  character(1)
)

cat("\n", R.version.string, "\n", sep = "")
cat("Bioconductor ",
    tryCatch(as.character(BiocManager::version()), error = function(e) "unknown"),
    "\n\n", sep = "")
for (p in all_pkgs) {
  cat(sprintf("%-50s %s\n", p, ifelse(is.na(versions[[p]]), "MISSING", versions[[p]])))
}

still_missing <- all_pkgs[is.na(versions)]
if (length(still_missing)) {
  cat("\n")
  stop("Not installed: ", paste(still_missing, collapse = ", "),
       call. = FALSE)
}
cat("\nAll required packages are installed.\n")
