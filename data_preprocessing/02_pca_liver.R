suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(viridis)
  library(readr)
})

source(file.path("data_preprocessing", "00_config.R"))
source(file.path(dir_helpers, "helpers_io.R"))
source(file.path(dir_helpers, "helpers_plotting.R"))

ensure_dirs(all_output_dirs)

message("=== PCA script started ===")

check_file_exists(file_beta_raw, "Raw imputed beta matrix")
check_file_exists(file_beta_combat, "ComBat imputed beta matrix")
check_file_exists(file_pheno_all, "Combined phenotype table")

beta_raw <- readRDS(file_beta_raw)
beta_cb <- readRDS(file_beta_combat)
pheno_all <- readRDS(file_pheno_all)

check_object_columns(
  pheno_all,
  c("Sample_Name", "Dataset", "DiseaseGroup", "Progression3", "Age", "Array"),
  "pheno_all"
)

stopifnot(ncol(beta_raw) == nrow(pheno_all))
stopifnot(ncol(beta_cb) == nrow(pheno_all))
stopifnot(all(colnames(beta_raw) == pheno_all$Sample_Name))
stopifnot(all(colnames(beta_cb) == pheno_all$Sample_Name))

pheno_masld <- pheno_all %>%
  mutate(
    Dataset = droplevels(as.factor(Dataset)),
    Array = droplevels(as.factor(Array))
  )

run_global_pca <- function(beta_all_mat, pheno_pca_df, top_n = 50000L) {
  beta_pca_input <- beta_all_mat[, pheno_pca_df$Sample_Name, drop = FALSE]
  stopifnot(all(colnames(beta_pca_input) == pheno_pca_df$Sample_Name))

  keep_probes <- which(rowSums(is.na(beta_pca_input)) == 0)
  beta_complete <- beta_pca_input[keep_probes, , drop = FALSE]

  row_var <- function(mat) {
    m <- rowMeans(mat, na.rm = TRUE)
    rowMeans((mat - m)^2, na.rm = TRUE)
  }

  var_all <- row_var(beta_complete)
  top_idx <- order(var_all, decreasing = TRUE)[seq_len(min(top_n, length(var_all)))]
  beta_pca <- beta_complete[top_idx, , drop = FALSE]

  pca_res <- prcomp(t(beta_pca), center = TRUE, scale. = FALSE)
  pc_scores <- as.data.frame(pca_res$x[, 1:3, drop = FALSE])
  pc_scores$Sample_Name <- rownames(pc_scores)

  pca_df <- pheno_pca_df %>%
    left_join(pc_scores, by = "Sample_Name")

  stopifnot(nrow(pca_df) == ncol(beta_pca))

  list(
    pca = pca_res,
    pca_df = pca_df,
    beta_pca = beta_pca
  )
}

res_raw <- run_global_pca(beta_raw, pheno_masld, top_n = pca_top_n)
res_cb  <- run_global_pca(beta_cb,  pheno_masld, top_n = pca_top_n)

pca_df_raw <- res_raw$pca_df %>%
  mutate(Correction = "Raw")

pca_df_cb <- res_cb$pca_df %>%
  mutate(Correction = "ComBat")

pca_df_both <- bind_rows(pca_df_raw, pca_df_cb) %>%
  mutate(
    Correction = factor(Correction, levels = c("Raw", "ComBat")),
    Dataset = factor(as.character(Dataset), levels = names(shape_dataset_map)),
    Array = factor(as.character(Array), levels = names(array_shapes))
  )

save_rds_safe(
  pca_df_both,
  file_pca_df
)

write_csv_safe(
  pca_df_both,
  file.path(dir_combined, "pca_df_raw_vs_combat_top50k.csv")
)

plot_scree <- function(pca_obj, title = "") {
  ve <- summary(pca_obj)$importance["Proportion of Variance", ]
  ve_df <- tibble(
    PC = factor(paste0("PC", seq_along(ve)), levels = paste0("PC", seq_along(ve))),
    VarExpl = as.numeric(ve)
  ) %>%
    slice_head(n = 10) %>%
    mutate(Percent = 100 * VarExpl)

  ggplot(ve_df, aes(x = PC, y = Percent)) +
    geom_col() +
    theme_pipeline(14) +
    labs(title = title, x = NULL, y = "Variance explained (%)")
}

p_scree_raw <- plot_scree(res_raw$pca, "Explained variance (RAW, top PCs)")
p_scree_cb  <- plot_scree(res_cb$pca,  "Explained variance (ComBat, top PCs)")

save_plot(
  p_scree_raw,
  file.path(dir_out_pca_figures, "scree_raw_top10.png"),
  width = 8,
  height = 5,
  dpi = 300
)

save_plot(
  p_scree_cb,
  file.path(dir_out_pca_figures, "scree_combat_top10.png"),
  width = 8,
  height = 5,
  dpi = 300
)

# ---- legend labels from unique samples (not Raw+ComBat duplicated) ----
labels_dg <- pheno_all |>
  dplyr::distinct(Sample_Name, DiseaseGroup) |>
  make_count_labels(DiseaseGroup)

labels_p3 <- pheno_all |>
  dplyr::distinct(Sample_Name, Progression3) |>
  make_count_labels(Progression3)

labels_array <- pheno_all |>
  dplyr::distinct(Sample_Name, Array) |>
  make_count_labels(Array)

# ------------------------------------------------------------------
# Global plots
# ------------------------------------------------------------------

scores_all_plot <- pca_df_both %>%
  keep_diseasegroup(stage_levels) %>%
  filter(!is.na(PC1), !is.na(PC2))

p_pc12_dg <- ggplot(
  scores_all_plot,
  aes(PC1, PC2, color = DiseaseGroup, shape = Dataset)
) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = stage_colors, labels = labels_dg, drop = FALSE) +
  scale_shape_manual(values = shape_dataset_map, drop = FALSE) +
  facet_wrap(~Correction) +
  theme_pipeline(12) +
  labs(
    title = "PCA (PC1 vs PC2): DiseaseGroup - Raw vs ComBat",
    x = "PC1",
    y = "PC2"
  )

save_plot(
  p_pc12_dg,
  file.path(dir_out_pca_figures, "pca_pc1_pc2_diseasegroup_raw_vs_combat.png"),
  width = 12,
  height = 7,
  dpi = 300
)

plot_df_prog3 <- pca_df_both %>%
  keep_progression3(prog3_levels) %>%
  filter(!is.na(PC1), !is.na(PC2))

p_pc12_prog3 <- ggplot(
  plot_df_prog3,
  aes(PC1, PC2, color = Progression3, shape = Dataset)
) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = prog3_colors, labels = labels_p3, drop = FALSE) +
  scale_shape_manual(values = shape_dataset_map, drop = FALSE) +
  facet_wrap(~Correction) +
  theme_pipeline(12) +
  labs(
    title = "PCA (PC1 vs PC2): Progression3 - Raw vs ComBat",
    subtitle = "All datasets in one plot; shapes = Dataset",
    x = "PC1",
    y = "PC2"
  )

save_plot(
  p_pc12_prog3,
  file.path(dir_out_pca_figures, "pca_pc1_pc2_progression3_raw_vs_combat.png"),
  width = 12,
  height = 7,
  dpi = 300
)

plot_df_prog3_23 <- pca_df_both %>%
  keep_progression3(prog3_levels) %>%
  filter(!is.na(PC2), !is.na(PC3))

p_pc23_prog3 <- ggplot(
  plot_df_prog3_23,
  aes(PC2, PC3, color = Progression3, shape = Dataset)
) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = prog3_colors, labels = labels_p3, drop = FALSE) +
  scale_shape_manual(values = shape_dataset_map, drop = FALSE) +
  facet_wrap(~Correction) +
  theme_pipeline(12) +
  labs(
    title = "PCA (PC2 vs PC3): Progression3 - Raw vs ComBat",
    subtitle = "All datasets together; shapes = Dataset",
    x = "PC2",
    y = "PC3"
  )

save_plot(
  p_pc23_prog3,
  file.path(dir_out_pca_figures, "pca_pc2_pc3_progression3_raw_vs_combat.png"),
  width = 12,
  height = 7,
  dpi = 300
)

plot_df_age <- pca_df_both %>%
  filter(!is.na(PC1), !is.na(PC2), !is.na(Age))

p_pc12_age <- ggplot(
  plot_df_age,
  aes(PC1, PC2, color = Age, shape = Dataset)
) +
  geom_point(size = 2, alpha = 0.8) +
  scale_shape_manual(values = shape_dataset_map, drop = FALSE) +
  scale_color_viridis_c(option = "plasma", na.value = "grey80") +
  facet_wrap(~Correction) +
  theme_pipeline(12) +
  labs(
    title = "PCA colored by Age - Raw vs ComBat",
    subtitle = "All datasets together; shapes = Dataset",
    x = "PC1",
    y = "PC2",
    color = "Age"
  )

save_plot(
  p_pc12_age,
  file.path(dir_out_pca_figures, "pca_pc1_pc2_age_raw_vs_combat.png"),
  width = 12,
  height = 7,
  dpi = 300
)

# ---- Array type plot ----
p_pc12_array <- ggplot(
  scores_all_plot,
  aes(PC1, PC2, color = Array, shape = Array)
) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = array_colors, labels = labels_array, drop = FALSE) +
  scale_shape_manual(values = array_shapes, labels = labels_array, drop = FALSE) +
  facet_wrap(~Correction) +
  theme_pipeline(12) +
  labs(
    title = "PCA (PC1 vs PC2): Array type",
    x = "PC1",
    y = "PC2"
  )

save_plot(
  p_pc12_array,
  file.path(dir_out_pca_figures, "pca_pc1_pc2_array.png"),
  width = 10,
  height = 6,
  dpi = 450
)

# ------------------------------------------------------------------
# Within-dataset PCA (pre-ComBat/raw-imputed)
# ------------------------------------------------------------------

beta_for_within <- beta_raw
stopifnot(identical(colnames(beta_for_within), pheno_masld$Sample_Name))
stopifnot(ncol(beta_for_within) == nrow(pheno_masld))

datasets <- sort(unique(as.character(pheno_masld$Dataset)))

dir_thesis_out_figures <- file.path(project_root, "outputs", "thesis_outputs", "figures")
dir_thesis_out_tables <- file.path(project_root, "outputs", "thesis_outputs", "tables")
dir_results_figures <- file.path(project_root, "results", "figures")
dir_results_tables <- file.path(project_root, "results", "tables")

ensure_dirs(c(
  dir_thesis_out_figures,
  dir_thesis_out_tables,
  dir_results_figures,
  dir_results_tables
))

save_plot_pair <- function(plot, path_stem, width, height, dpi = 300) {
  save_plot(
    plot,
    paste0(path_stem, ".png"),
    width = width,
    height = height,
    dpi = dpi
  )
  save_plot(
    plot,
    paste0(path_stem, ".pdf"),
    width = width,
    height = height,
    dpi = dpi,
    device = grDevices::cairo_pdf
  )
}

update_manifest_row <- function(manifest_path, key_col, row) {
  row <- tibble::as_tibble(row) %>%
    mutate(across(everything(), as.character))

  existing <- if (file.exists(manifest_path)) {
    readr::read_csv(
      manifest_path,
      col_types = readr::cols(.default = "c"),
      show_col_types = FALSE
    )
  } else {
    tibble::tibble()
  }

  if (nrow(existing) > 0) {
    missing_cols <- setdiff(names(row), names(existing))
    for (col in missing_cols) {
      existing[[col]] <- NA_character_
    }
    existing <- existing %>%
      select(all_of(names(row))) %>%
      filter(.data[[key_col]] != row[[key_col]][[1]])
  }

  bind_rows(existing, row) %>%
    arrange(.data[[key_col]]) %>%
    readr::write_csv(manifest_path)
}

run_pca_one_dataset <- function(dataset_name, beta_mat, pheno_df, top_n = 50000L) {
  samples_ds <- pheno_df %>%
    filter(Dataset == dataset_name) %>%
    pull(Sample_Name)

  beta_ds <- beta_mat[, samples_ds, drop = FALSE]

  message("\n=== Dataset: ", dataset_name, " ===")
  message("Samples: ", ncol(beta_ds), " | CpGs: ", nrow(beta_ds))

  keep_probes <- rowSums(is.na(beta_ds)) == 0
  beta_ds <- beta_ds[keep_probes, , drop = FALSE]

  message("CpGs after removing NAs: ", nrow(beta_ds))

  probe_var <- apply(beta_ds, 1, var)
  probe_var <- sort(probe_var, decreasing = TRUE)
  top_n_eff <- min(top_n, length(probe_var))
  top_probes <- names(probe_var)[seq_len(top_n_eff)]

  message("Top variable CpGs used: ", top_n_eff)

  x <- t(beta_ds[top_probes, , drop = FALSE])
  pca <- prcomp(x, center = TRUE, scale. = FALSE)
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

  message(
    "Var explained PC1/PC2/PC3: ",
    round(var_expl[1], 3), " / ",
    round(var_expl[2], 3), " / ",
    round(var_expl[3], 3)
  )

  scores <- as.data.frame(pca$x[, 1:3, drop = FALSE])
  scores$Sample_Name <- rownames(scores)

  list(
    pca = pca,
    scores = scores,
    var_expl = var_expl
  )
}

pca_by_dataset_raw <- setNames(
  lapply(
    datasets,
    run_pca_one_dataset,
    beta_mat = beta_for_within,
    pheno_df = pheno_masld,
    top_n = pca_top_n
  ),
  datasets
)

scores_all_raw <- bind_rows(lapply(names(pca_by_dataset_raw), function(ds) {
  out <- pca_by_dataset_raw[[ds]]$scores
  out$Dataset <- ds
  out
})) %>%
  left_join(
    pheno_masld %>%
      select(
        Sample_Name,
        Dataset,
        DiseaseGroup,
        Progression3,
        Age,
        Sex
      ),
    by = c("Sample_Name", "Dataset")
  ) %>%
  mutate(
    Correction = "Raw"
  )

save_rds_safe(
  pca_by_dataset_raw,
  file.path(dir_out_pca_rds, "pca_by_dataset_top50k_raw.rds")
)

save_rds_safe(
  scores_all_raw,
  file.path(dir_out_pca_rds, "pca_scores_by_dataset_top50k_raw.rds")
)

var_expl_df <- bind_rows(lapply(names(pca_by_dataset_raw), function(ds) {
  ve <- pca_by_dataset_raw[[ds]]$var_expl
  tibble(
    Dataset = ds,
    PC1 = ve[1],
    PC2 = ve[2],
    PC3 = ve[3]
  )
}))

write_csv_safe(
  var_expl_df,
  file.path(dir_out_pca_tables, "within_dataset_variance_explained_raw.csv")
)

var_expl_pub <- var_expl_df %>%
  transmute(
    Dataset,
    `PC1 variance (%)` = round(100 * PC1, 1),
    `PC2 variance (%)` = round(100 * PC2, 1),
    `PC3 variance (%)` = round(100 * PC3, 1)
  )

readr::write_csv(
  var_expl_pub,
  file.path(dir_results_tables, "T03_within_dataset_pca_variance_explained.csv")
)

readr::write_csv(
  var_expl_pub,
  file.path(dir_thesis_out_tables, "T03_within_dataset_pca_variance_explained.csv")
)

scores_all_plot_ds <- scores_all_raw %>%
  keep_diseasegroup(stage_levels) %>%
  filter(!is.na(PC1), !is.na(PC2)) %>%
  mutate(
    Dataset = factor(as.character(Dataset), levels = names(shape_dataset_map))
  )

labels_dg_ds <- make_count_labels(scores_all_plot_ds, DiseaseGroup)

facet_annotations <- scores_all_plot_ds %>%
  count(Dataset, name = "n") %>%
  left_join(
    var_expl_df %>%
      mutate(Dataset = factor(as.character(Dataset), levels = names(shape_dataset_map))),
    by = "Dataset"
  ) %>%
  mutate(
    label = sprintf("n=%d\nPC1 %.1f%%; PC2 %.1f%%", n, 100 * PC1, 100 * PC2)
  )

p_ds_pc12_dg <- ggplot(
  scores_all_plot_ds,
  aes(PC1, PC2, color = DiseaseGroup)
) +
  geom_point(alpha = 0.8, size = 2) +
  geom_text(
    data = facet_annotations,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.15,
    size = 2.8,
    color = "grey20",
    lineheight = 0.95
  ) +
  scale_color_manual(values = stage_colors, labels = labels_dg_ds, drop = FALSE) +
  facet_wrap(~Dataset, scales = "free") +
  theme_pipeline(12) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))
  ) +
  theme(
    legend.box = "vertical",
    legend.box.just = "center"
  ) +
  labs(
    title = "Within-study PCA (PC1 vs PC2) - Raw/pre-ComBat",
    subtitle = "Top 50,000 variable CpGs selected within each study; color = DiseaseGroup",
    x = "PC1",
    y = "PC2",
    color = "DiseaseGroup"
  )

save_plot_pair(
  p_ds_pc12_dg,
  file.path(dir_out_pca_figures, "within_dataset_pca_pc1_pc2_diseasegroup_raw"),
  width = 12,
  height = 8,
  dpi = 300
)

save_plot_pair(
  p_ds_pc12_dg,
  file.path(dir_thesis_out_figures, "F16_within_dataset_pca_pc1_pc2_diseasegroup_raw"),
  width = 12,
  height = 8,
  dpi = 300
)

save_plot_pair(
  p_ds_pc12_dg,
  file.path(dir_results_figures, "F16_within_dataset_pca_pc1_pc2_diseasegroup_raw"),
  width = 12,
  height = 8,
  dpi = 300
)

scores_all_prog3 <- scores_all_raw %>%
  keep_progression3(prog3_levels) %>%
  filter(!is.na(PC1), !is.na(PC2)) %>%
  mutate(
    Dataset = factor(as.character(Dataset), levels = names(shape_dataset_map))
  )

labels_p3_ds <- make_count_labels(scores_all_prog3, Progression3)

p_ds_pc12_prog3 <- ggplot(scores_all_prog3, aes(PC1, PC2, color = Progression3)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = prog3_colors, labels = labels_p3_ds, drop = FALSE) +
  facet_wrap(~Dataset, scales = "free") +
  theme_pipeline(12) +
  labs(
    title = "Within-study PCA (PC1 vs PC2) - Raw/pre-ComBat",
    subtitle = "Color = Progression3",
    x = "PC1",
    y = "PC2"
  )

save_plot_pair(
  p_ds_pc12_prog3,
  file.path(dir_out_pca_figures, "within_dataset_pca_pc1_pc2_progression3_raw"),
  width = 12,
  height = 8,
  dpi = 300
)

scores_all_prog3_23 <- scores_all_raw %>%
  keep_progression3(prog3_levels) %>%
  filter(!is.na(PC2), !is.na(PC3)) %>%
  mutate(
    Dataset = factor(as.character(Dataset), levels = names(shape_dataset_map))
  )

labels_p3_ds_23 <- make_count_labels(scores_all_prog3_23, Progression3)

p_ds_pc23_prog3 <- ggplot(scores_all_prog3_23, aes(PC2, PC3, color = Progression3)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = prog3_colors, labels = labels_p3_ds_23, drop = FALSE) +
  facet_wrap(~Dataset, scales = "free") +
  theme_pipeline(12) +
  labs(
    title = "Within-study PCA (PC2 vs PC3) - Raw/pre-ComBat",
    subtitle = "Color = Progression3",
    x = "PC2",
    y = "PC3"
  )

save_plot_pair(
  p_ds_pc23_prog3,
  file.path(dir_out_pca_figures, "within_dataset_pca_pc2_pc3_progression3_raw"),
  width = 12,
  height = 8,
  dpi = 300
)

update_manifest_row(
  file.path(dir_results_tables, "figure_manifest.csv"),
  "figure",
  tibble(
    figure = "F16",
    file_stem = "F16_within_dataset_pca_pc1_pc2_diseasegroup_raw",
    main_message = "Within-study PCA on the raw/pre-ComBat matrix; color = DiseaseGroup."
  )
)

update_manifest_row(
  file.path(dir_thesis_out_tables, "thesis_figure_manifest.csv"),
  "figure",
  tibble(
    figure = "F16",
    file_stem = "F16_within_dataset_pca_pc1_pc2_diseasegroup_raw",
    main_message = "Within-study PCA on the raw/pre-ComBat matrix; color = DiseaseGroup."
  )
)

update_manifest_row(
  file.path(dir_results_tables, "table_manifest.csv"),
  "table",
  tibble(
    table = "T03",
    file = "T03_within_dataset_pca_variance_explained.csv",
    main_message = "PC1-PC3 variance explained for raw/pre-ComBat within-study PCAs."
  )
)

update_manifest_row(
  file.path(dir_thesis_out_tables, "thesis_table_manifest.csv"),
  "table",
  tibble(
    table = "T03",
    file = "T03_within_dataset_pca_variance_explained.csv",
    main_message = "PC1-PC3 variance explained for raw/pre-ComBat within-study PCAs."
  )
)

# ------
# Array type PCA plot
# ------

scores_all_plot <- pca_df_both %>%
  keep_diseasegroup(stage_levels) %>%
  filter(!is.na(PC1), !is.na(PC2))

labels_array <- make_count_labels(scores_all_plot, Array)

p_pc12_array <- ggplot(
  scores_all_plot,
  aes(PC1, PC2, color = Array, shape = Array)
) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = array_colors, labels = labels_array, drop = FALSE) +
  scale_shape_manual(values = array_shapes, labels = labels_array, drop = FALSE) +
  facet_wrap(~Correction) +
  theme_pipeline(12) +
  labs(
    title = "PCA (PC1 vs PC2): Array type",
    x = "PC1",
    y = "PC2"
  )

save_plot(
  p_pc12_array,
  file.path(dir_out_pca_figures, "pca_pc1_pc2_array.png"),
  width = 10,
  height = 6,
  dpi = 450
)

message("=== PCA script finished ===")

make_count_labels <- function(df, var) {
  var <- rlang::ensym(var)
  df |>
    dplyr::filter(!is.na(!!var)) |>
    dplyr::count(!!var, name = "n") |>
    dplyr::mutate(label = paste0(!!var, " (n=", n, ")")) |>
    tibble::deframe()
}
