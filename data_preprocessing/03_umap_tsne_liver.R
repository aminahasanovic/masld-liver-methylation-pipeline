# ---------------------------------------------------------------------------
# Exploratory nonlinear embeddings (UMAP and t-SNE) of the ComBat matrix
#
# Supplementary companion to 02_pca_liver.R. Both embeddings are computed from
# the first 30 principal components of the 50,000 most variable complete CpGs,
# so they describe the same feature space as the global PCA.
#
# Produces the supplementary UMAP/t-SNE figure reported in the thesis, plus the
# embedding coordinates and the parameter settings as CSV.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(readr)
})

source(file.path("data_preprocessing", "00_config.R"))
source(file.path(dir_helpers, "helpers_io.R"))
source(file.path(dir_helpers, "helpers_plotting.R"))

ensure_dirs(all_output_dirs)

message("=== UMAP / t-SNE script started ===")

for (pkg in c("uwot", "Rtsne")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for the exploratory embeddings.", call. = FALSE)
  }
}

check_file_exists(file_beta_combat, "ComBat imputed beta matrix")
check_file_exists(file_pheno_all, "Combined phenotype table")

beta_cb <- readRDS(file_beta_combat)
pheno_all <- readRDS(file_pheno_all)

check_object_columns(
  pheno_all,
  c("Sample_Name", "Dataset", "DiseaseGroup", "Progression3", "Age", "Array"),
  "pheno_all"
)

# The study labels used in the figures come from `dataset_levels` in the config.
# A phenotype table written before the cohort was renamed still carries the old
# label, which would silently drop those samples from every figure, so the
# mismatch is resolved here and any remaining unknown label stops the script.
if (any(as.character(pheno_all$Dataset) == "Kim")) {
  message(
    "Phenotype table uses the previous label 'Kim' for the ITEN cohort; ",
    "relabelling ", sum(as.character(pheno_all$Dataset) == "Kim"), " samples."
  )
  pheno_all$Dataset <- dplyr::recode(as.character(pheno_all$Dataset), Kim = "ITEN")
}

unknown_datasets <- setdiff(unique(as.character(pheno_all$Dataset)), dataset_levels)
if (length(unknown_datasets) > 0) {
  stop(
    "Phenotype table contains study labels that are not defined in the config: ",
    paste(unknown_datasets, collapse = ", "),
    call. = FALSE
  )
}
pheno_all$Dataset <- factor(as.character(pheno_all$Dataset), levels = dataset_levels)

stopifnot(ncol(beta_cb) == nrow(pheno_all))
stopifnot(all(colnames(beta_cb) == pheno_all$Sample_Name))

# ------------------------------------------------------------------
# Embedding input: PCs of the most variable complete CpGs
# ------------------------------------------------------------------

embedding_seed <- 20260611L
n_pcs <- 30L

build_embedding_input <- function(beta_mat, pheno_df, top_n = pca_top_n, n_pcs = 30L) {
  beta_in <- beta_mat[, pheno_df$Sample_Name, drop = FALSE]
  stopifnot(all(colnames(beta_in) == pheno_df$Sample_Name))

  keep_probes <- which(rowSums(is.na(beta_in)) == 0)
  beta_complete <- beta_in[keep_probes, , drop = FALSE]

  row_var <- function(mat) {
    m <- rowMeans(mat, na.rm = TRUE)
    rowMeans((mat - m)^2, na.rm = TRUE)
  }

  var_all <- row_var(beta_complete)
  top_idx <- order(var_all, decreasing = TRUE)[seq_len(min(top_n, length(var_all)))]

  max_rank <- min(50L, ncol(beta_complete) - 1L)
  pca_res <- prcomp(
    t(beta_complete[top_idx, , drop = FALSE]),
    center = TRUE,
    scale. = FALSE,
    rank. = max_rank
  )
  scores <- pca_res$x[, seq_len(min(n_pcs, ncol(pca_res$x))), drop = FALSE]

  rm(beta_in, beta_complete)
  gc()

  list(scores = scores, n_cpgs = length(top_idx), n_pcs = ncol(scores))
}

emb_input <- build_embedding_input(beta_cb, pheno_all, top_n = pca_top_n, n_pcs = n_pcs)
message(
  "Embedding input: ", emb_input$n_cpgs, " CpGs -> ",
  emb_input$n_pcs, " PCs for ", nrow(emb_input$scores), " samples."
)

# ------------------------------------------------------------------
# UMAP and t-SNE
# ------------------------------------------------------------------

set.seed(embedding_seed)
umap_xy <- uwot::umap(
  emb_input$scores,
  n_neighbors = 15L,
  min_dist = 0.1,
  metric = "euclidean",
  scale = FALSE,
  n_threads = 1L,
  ret_model = FALSE,
  verbose = FALSE
)

tsne_perplexity <- min(30L, floor((nrow(emb_input$scores) - 1L) / 3L) - 1L)
set.seed(embedding_seed)
tsne_xy <- Rtsne::Rtsne(
  emb_input$scores,
  dims = 2L,
  perplexity = tsne_perplexity,
  theta = 0.5,
  pca = FALSE,
  check_duplicates = FALSE,
  verbose = FALSE,
  max_iter = 1000L
)$Y

embedding_df <- bind_rows(
  tibble(
    Method = "UMAP",
    Sample_Name = rownames(emb_input$scores),
    Dim1 = umap_xy[, 1],
    Dim2 = umap_xy[, 2]
  ),
  tibble(
    Method = "t-SNE",
    Sample_Name = rownames(emb_input$scores),
    Dim1 = tsne_xy[, 1],
    Dim2 = tsne_xy[, 2]
  )
) %>%
  left_join(
    pheno_all %>%
      select(Sample_Name, Dataset, Array, DiseaseGroup, Progression3, Age, Sex),
    by = "Sample_Name"
  ) %>%
  mutate(Method = factor(Method, levels = c("UMAP", "t-SNE")))

embedding_params <- tibble::tribble(
  ~parameter, ~value,
  "input_matrix", "ComBat-corrected beta matrix",
  "feature_selection", paste0("Top ", format(emb_input$n_cpgs, big.mark = ","), " most variable complete CpGs"),
  "embedding_input", paste0("First ", emb_input$n_pcs, " principal components"),
  "random_seed", as.character(embedding_seed),
  "umap_package", "uwot",
  "umap_n_neighbors", "15",
  "umap_min_dist", "0.1",
  "umap_metric", "euclidean",
  "tsne_package", "Rtsne",
  "tsne_perplexity", as.character(tsne_perplexity),
  "tsne_theta", "0.5",
  "tsne_max_iter", "1000"
)

write_csv(embedding_df, file.path(dir_out_pca_tables, "umap_tsne_coordinates_combat_pca30.csv"))
write_csv(embedding_params, file.path(dir_out_pca_tables, "umap_tsne_parameters_combat_pca30.csv"))

# ------------------------------------------------------------------
# Figure: DiseaseGroup colour, study shape
# ------------------------------------------------------------------

labels_dg <- pheno_all %>%
  dplyr::distinct(Sample_Name, DiseaseGroup) %>%
  make_count_labels(DiseaseGroup)

plot_df <- embedding_df %>%
  keep_diseasegroup(stage_levels) %>%
  filter(!is.na(Dim1), !is.na(Dim2)) %>%
  mutate(Dataset = droplevels(as.factor(Dataset)))

p_embed_dg <- ggplot(plot_df, aes(Dim1, Dim2, color = DiseaseGroup, shape = Dataset)) +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(values = stage_colors, labels = labels_dg, drop = FALSE) +
  scale_shape_manual(values = shape_dataset_map, drop = FALSE) +
  facet_wrap(~Method, scales = "free") +
  theme_pipeline(12) +
  labs(
    title = "UMAP and t-SNE of the ComBat matrix (30 PCs of the top 50,000 variable CpGs)",
    x = "Embedding dimension 1",
    y = "Embedding dimension 2"
  )

save_plot(
  p_embed_dg,
  file.path(dir_out_pca_figures, "umap_tsne_diseasegroup_combat_pca30.png"),
  width = 12,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------------
# Figure: Progression3 colour, study shape
# ------------------------------------------------------------------

labels_p3 <- pheno_all %>%
  dplyr::distinct(Sample_Name, Progression3) %>%
  make_count_labels(Progression3)

plot_df_prog3 <- embedding_df %>%
  keep_progression3(prog3_levels) %>%
  filter(!is.na(Dim1), !is.na(Dim2)) %>%
  mutate(Dataset = droplevels(as.factor(Dataset)))

p_embed_prog3 <- ggplot(plot_df_prog3, aes(Dim1, Dim2, color = Progression3, shape = Dataset)) +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(values = prog3_colors, labels = labels_p3, drop = FALSE) +
  scale_shape_manual(values = shape_dataset_map, drop = FALSE) +
  facet_wrap(~Method, scales = "free") +
  theme_pipeline(12) +
  labs(
    title = "UMAP and t-SNE of the ComBat matrix: collapsed progression stage",
    x = "Embedding dimension 1",
    y = "Embedding dimension 2"
  )

save_plot(
  p_embed_prog3,
  file.path(dir_out_pca_figures, "umap_tsne_progression3_combat_pca30.png"),
  width = 12,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------------
# Figure: study of origin
# ------------------------------------------------------------------

labels_ds <- pheno_all %>%
  dplyr::distinct(Sample_Name, Dataset) %>%
  make_count_labels(Dataset)

plot_df_dataset <- embedding_df %>%
  filter(!is.na(Dim1), !is.na(Dim2), !is.na(Dataset)) %>%
  mutate(Dataset = factor(as.character(Dataset), levels = dataset_levels))

p_embed_dataset <- ggplot(plot_df_dataset, aes(Dim1, Dim2, color = Dataset, shape = Dataset)) +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(values = dataset_colors, labels = labels_ds, drop = FALSE) +
  scale_shape_manual(values = shape_dataset_map, labels = labels_ds, drop = FALSE) +
  facet_wrap(~Method, scales = "free") +
  theme_pipeline(12) +
  labs(
    title = "UMAP and t-SNE of the ComBat matrix: study of origin",
    x = "Embedding dimension 1",
    y = "Embedding dimension 2"
  )

save_plot(
  p_embed_dataset,
  file.path(dir_out_pca_figures, "umap_tsne_dataset_combat_pca30.png"),
  width = 12,
  height = 7,
  dpi = 300
)

message("=== UMAP / t-SNE script finished ===")
