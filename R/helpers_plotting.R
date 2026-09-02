# Reminder: helpers_plotting.R uses ensure_dir() from helpers_io.R --> source helpers_io.R before helpers_plotting.R

base_theme_classic <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "white", color = "black"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

base_theme_bw <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "white", color = "black"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

save_plot <- function(plot, path, width = 10, height = 7, dpi = 300, device = NULL) {
  ensure_dir(dirname(path))
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    device = device
  )
  invisible(path)
}

keep_diseasegroup <- function(df, stage_levels) {
  df %>%
    dplyr::filter(!is.na(DiseaseGroup), as.character(DiseaseGroup) %in% stage_levels) %>%
    dplyr::mutate(
      DiseaseGroup = factor(as.character(DiseaseGroup), levels = stage_levels)
    )
}

keep_progression3 <- function(df, prog3_levels) {
  df %>%
    dplyr::filter(!is.na(Progression3), as.character(Progression3) %in% prog3_levels) %>%
    dplyr::mutate(
      Progression3 = factor(as.character(Progression3), levels = prog3_levels)
    )
}

tagify <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

theme_pipeline <- function(base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey85"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

make_count_labels <- function(df, var) {
  var <- rlang::ensym(var)
  df |>
    dplyr::filter(!is.na(!!var)) |>
    dplyr::count(!!var, name = "n") |>
    dplyr::mutate(label = paste0(!!var, " (n=", n, ")")) |>
    dplyr::select(!!var, label) |>
    tibble::deframe()
}
