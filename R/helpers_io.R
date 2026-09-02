ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

ensure_dirs <- function(paths) {
  invisible(lapply(paths, ensure_dir))
}

save_rds_safe <- function(object, path) {
  ensure_dir(dirname(path))
  saveRDS(object, path)
  invisible(path)
}

write_csv_safe <- function(x, path, row.names = FALSE) {
  ensure_dir(dirname(path))
  utils::write.csv(x, path, row.names = row.names)
  invisible(path)
}

write_tsv_safe <- function(x, path) {
  ensure_dir(dirname(path))
  readr::write_tsv(x, path)
  invisible(path)
}

check_file_exists <- function(path, label = NULL) {
  if (!file.exists(path)) {
    if (is.null(label)) label <- path
    stop("Missing required file: ", label, "\nPath: ", path, call. = FALSE)
  }
  invisible(path)
}

check_object_columns <- function(df, required_cols, object_name = deparse(substitute(df))) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}