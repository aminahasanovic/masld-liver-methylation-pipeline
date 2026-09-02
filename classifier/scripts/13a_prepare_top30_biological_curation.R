source("scripts/00_setup_paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

candidate_summary_path <- file.path(
  out_paths$results_features,
  "elastic_net_cpg_candidate_summary.csv"
)
meta_results_path <- file.path(
  out_paths$results_features,
  "elastic_net_candidate_cpg_meta_results.csv"
)
strict_unique_path <- file.path(
  out_paths$results_features,
  "elastic_net_candidate_cpg_meta_strict_unique.csv"
)
output_path <- file.path(
  out_paths$results_features,
  "elastic_net_top30_biological_curation_input.csv"
)
script13_path <- file.path("scripts", "13_add_biological_context_top_cpgs.R")

required_files <- c(candidate_summary_path, meta_results_path, strict_unique_path, script13_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input(s): ", paste(missing_files, collapse = ", "))
}

message("Preparing Top-30 biological curation input table.")
message("Candidate summary: ", candidate_summary_path)

candidate_summary <- readr::read_csv(candidate_summary_path, show_col_types = FALSE)
meta_results <- readr::read_csv(meta_results_path, show_col_types = FALSE)
strict_unique <- readr::read_csv(strict_unique_path, show_col_types = FALSE)

if (nrow(candidate_summary) < 30) {
  stop("Candidate summary contains fewer than 30 ranked candidates.")
}
if (anyDuplicated(candidate_summary$cpg) > 0) {
  stop("Candidate summary contains duplicated CpGs.")
}

top30 <- candidate_summary |>
  dplyr::mutate(statistical_rank = dplyr::row_number()) |>
  dplyr::slice_head(n = 30)

if (dplyr::n_distinct(top30$cpg) != 30) {
  stop("Top-30 candidate selection did not produce 30 unique CpGs.")
}

first_non_empty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

collapse_unique <- function(x, sep = ";") {
  x <- as.character(x)
  x <- unlist(strsplit(x, sep, fixed = TRUE), use.names = FALSE)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(sort(unique(x)), collapse = sep)
}

empty_to_na <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}

load_illumina_annotation <- function() {
  if (!requireNamespace("minfi", quietly = TRUE) ||
      !requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
    warning("Illumina annotation packages are unavailable; using existing pipeline annotation only.")
    return(tibble::tibble(cpg = character()))
  }

  suppressPackageStartupMessages({
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  })

  anno_450k_raw <- minfi::getAnnotation(
    IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19
  ) |>
    as.data.frame() |>
    tibble::rownames_to_column("cpg")

  anno_epic_raw <- NULL
  epic_cpgs <- character()

  if (requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
    suppressPackageStartupMessages({
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    })

    anno_epic_raw <- minfi::getAnnotation(
      IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19
    ) |>
      as.data.frame() |>
      tibble::rownames_to_column("cpg")

    epic_cpgs <- anno_epic_raw$cpg
  }

  anno_450k_cpgs <- anno_450k_raw$cpg
  anno_450k <- anno_450k_raw |>
    dplyr::mutate(
      annotation_source_full = "450K_hg19",
      in_450k_annotation_full = .data$cpg %in% anno_450k_cpgs,
      in_epic_annotation_full = .data$cpg %in% epic_cpgs
    )

  if (is.null(anno_epic_raw)) {
    return(anno_450k)
  }

  anno_epic <- anno_epic_raw |>
    dplyr::mutate(
      annotation_source_full = "EPIC_hg19_fallback",
      in_450k_annotation_full = .data$cpg %in% anno_450k_cpgs,
      in_epic_annotation_full = .data$cpg %in% epic_cpgs
    )

  dplyr::bind_rows(
    anno_450k,
    dplyr::anti_join(anno_epic, anno_450k |> dplyr::select(cpg), by = "cpg")
  )
}

safe_select_org <- function(keys_in, keytype, columns) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE) ||
      length(keys_in) == 0) {
    return(tibble::tibble())
  }

  available_keys <- tryCatch(
    AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = keytype),
    error = function(e) character()
  )
  keys_in <- unique(keys_in[keys_in %in% available_keys])
  if (length(keys_in) == 0) {
    return(tibble::tibble())
  }

  suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = keys_in,
      keytype = keytype,
      columns = columns
    )
  ) |>
    tibble::as_tibble()
}

split_gene_symbols <- function(x) {
  x <- as.character(x)
  x <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  x <- trimws(x)
  sort(unique(x[!is.na(x) & nzchar(x)]))
}

map_gene_symbols <- function(gene_values) {
  original_genes <- sort(unique(unlist(lapply(gene_values, split_gene_symbols), use.names = FALSE)))
  if (length(original_genes) == 0 ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(tibble::tibble(
      original_gene_symbol = original_genes,
      current_gene_symbol = NA_character_,
      current_gene_aliases = NA_character_,
      gene_symbol_mapping_source = "org.Hs.eg.db_unavailable"
    ))
  }

  symbol_hits <- safe_select_org(
    original_genes,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ALIAS", "GENENAME", "ENTREZID")
  )
  alias_hits <- safe_select_org(
    original_genes,
    keytype = "ALIAS",
    columns = c("SYMBOL", "ALIAS", "GENENAME", "ENTREZID")
  )

  dplyr::bind_rows(lapply(original_genes, function(gene_symbol) {
    symbol_match <- symbol_hits |>
      dplyr::filter(.data$SYMBOL == gene_symbol, !is.na(.data$ENTREZID))

    if (nrow(symbol_match) > 0) {
      return(tibble::tibble(
        original_gene_symbol = gene_symbol,
        current_gene_symbol = gene_symbol,
        current_gene_aliases = collapse_unique(symbol_match$ALIAS),
        gene_symbol_mapping_source = "org.Hs.eg.db:SYMBOL"
      ))
    }

    alias_match <- alias_hits |>
      dplyr::filter(.data$ALIAS == gene_symbol, !is.na(.data$SYMBOL), !is.na(.data$ENTREZID))

    if (nrow(alias_match) == 0) {
      return(tibble::tibble(
        original_gene_symbol = gene_symbol,
        current_gene_symbol = NA_character_,
        current_gene_aliases = NA_character_,
        gene_symbol_mapping_source = "not_resolved_locally"
      ))
    }

    candidate_symbols <- sort(unique(alias_match$SYMBOL))
    current_symbol <- if (gene_symbol %in% candidate_symbols) {
      gene_symbol
    } else {
      paste(candidate_symbols, collapse = ";")
    }

    tibble::tibble(
      original_gene_symbol = gene_symbol,
      current_gene_symbol = current_symbol,
      current_gene_aliases = collapse_unique(alias_match$ALIAS),
      gene_symbol_mapping_source = "org.Hs.eg.db:ALIAS"
    )
  }))
}

gene_mapping <- map_gene_symbols(top30$UCSC_RefGene_Name)

if (nrow(gene_mapping) > 0) {
  current_symbols <- sort(unique(unlist(
    strsplit(gene_mapping$current_gene_symbol, ";", fixed = TRUE),
    use.names = FALSE
  )))
  current_symbols <- current_symbols[!is.na(current_symbols) & nzchar(current_symbols)]

  current_symbol_aliases <- safe_select_org(
    current_symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ALIAS", "ENTREZID")
  ) |>
    dplyr::group_by(.data$SYMBOL) |>
    dplyr::summarise(
      current_symbol_aliases = collapse_unique(.data$ALIAS),
      .groups = "drop"
    )

  if (nrow(current_symbol_aliases) > 0) {
    gene_mapping <- gene_mapping |>
      dplyr::rowwise() |>
      dplyr::mutate(
        current_gene_aliases = {
          symbols <- split_gene_symbols(.data$current_gene_symbol)
          resolved_aliases <- current_symbol_aliases |>
            dplyr::filter(.data$SYMBOL %in% symbols) |>
            dplyr::pull(.data$current_symbol_aliases)
          collapse_unique(c(.data$current_gene_aliases, resolved_aliases))
        }
      ) |>
      dplyr::ungroup()
  }
}

summarise_gene_mapping <- function(gene_annotation, field) {
  genes <- split_gene_symbols(gene_annotation)
  if (length(genes) == 0 || nrow(gene_mapping) == 0) {
    return(NA_character_)
  }

  values <- gene_mapping |>
    dplyr::filter(.data$original_gene_symbol %in% genes) |>
    dplyr::pull(.data[[field]])

  collapse_unique(values)
}

parse_old_biological_context <- function(path_in) {
  lines <- readLines(path_in, warn = FALSE)
  start <- grep("^\\s*biological_context\\s*<-\\s*tibble::tribble\\(", lines)
  if (length(start) != 1) {
    warning("Could not locate biological_context table in script 13.")
    return(tibble::tibble())
  }

  balance <- 0L
  end <- start
  for (i in seq(from = start, to = length(lines))) {
    balance <- balance +
      stringr::str_count(lines[[i]], "\\(") -
      stringr::str_count(lines[[i]], "\\)")
    if (balance == 0L && i > start) {
      end <- i
      break
    }
  }

  expr <- paste(lines[start:end], collapse = "\n")
  env <- new.env(parent = globalenv())
  env$tibble <- asNamespace("tibble")
  eval(parse(text = expr), envir = env)
  env$biological_context |>
    tibble::as_tibble()
}

old_biological_context <- parse_old_biological_context(script13_path)

# The gene-level context table was hardcoded in an earlier version of script 13,
# which now reads the manually curated review file instead. Where that table is
# no longer present, the "previously reviewed" columns of the curation sheet stay
# empty and are filled during the manual review.
if (!"gene_key" %in% names(old_biological_context)) {
  message(
    "No hardcoded gene-level context table in script 13; ",
    "the previously-reviewed columns of the curation sheet stay empty."
  )
  old_biological_context <- tibble::tibble(gene_key = character(0))
}

old_specific_context <- old_biological_context |>
  dplyr::filter(.data$gene_key != "intergenic_or_no_gene")

lookup_old_context <- function(gene_annotation, column) {
  genes <- split_gene_symbols(gene_annotation)
  if (length(genes) == 0 || nrow(old_specific_context) == 0) {
    return(NA_character_)
  }
  hits <- old_specific_context |>
    dplyr::filter(.data$gene_key %in% genes)
  if (nrow(hits) == 0) {
    return(NA_character_)
  }
  values <- hits[[column]]
  if (is.numeric(values)) {
    return(paste(sort(unique(values)), collapse = ";"))
  }
  collapse_unique(values)
}

find_nearest_tss <- function(cpg_tbl) {
  no_gene_tbl <- cpg_tbl |>
    dplyr::filter(is.na(.data$original_annotation) | .data$original_annotation == "") |>
    dplyr::filter(!is.na(.data$chromosome), !is.na(.data$MAPINFO))

  empty_result <- tibble::tibble(
    cpg = cpg_tbl$cpg,
    nearest_gene = NA_character_,
    nearest_gene_entrez = NA_character_,
    nearest_tss_transcript_id = NA_character_,
    nearest_tss_chr = NA_character_,
    nearest_tss_position = NA_integer_,
    nearest_tss_strand = NA_character_,
    nearest_tss_distance = NA_integer_,
    nearest_tss_orientation = NA_character_
  )

  if (nrow(no_gene_tbl) == 0) {
    return(empty_result)
  }

  needed <- c("TxDb.Hsapiens.UCSC.hg19.knownGene", "GenomicFeatures", "GenomicRanges", "IRanges", "S4Vectors")
  if (!all(vapply(needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)))) {
    warning("TxDb/GenomicRanges resources are unavailable; nearest TSS fields left as NA.")
    return(empty_result)
  }

  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  tx <- GenomicFeatures::transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
  gene_id <- as.character(S4Vectors::mcols(tx)$gene_id)
  tx <- tx[!is.na(gene_id) & nzchar(gene_id)]
  gene_id <- as.character(S4Vectors::mcols(tx)$gene_id)
  tss_pos <- ifelse(as.character(GenomicRanges::strand(tx)) == "-", GenomicRanges::end(tx), GenomicRanges::start(tx))

  tss <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(tx),
    ranges = IRanges::IRanges(start = tss_pos, end = tss_pos),
    strand = GenomicRanges::strand(tx),
    tx_id = as.character(S4Vectors::mcols(tx)$tx_id),
    tx_name = as.character(S4Vectors::mcols(tx)$tx_name),
    gene_id = gene_id
  )

  query <- GenomicRanges::GRanges(
    seqnames = no_gene_tbl$chromosome,
    ranges = IRanges::IRanges(start = as.integer(no_gene_tbl$MAPINFO), end = as.integer(no_gene_tbl$MAPINFO))
  )

  hits <- GenomicRanges::nearest(query, tss, ignore.strand = TRUE)
  valid <- !is.na(hits)
  if (!any(valid)) {
    return(empty_result)
  }

  nearest_tss <- tss[hits[valid]]
  query_valid <- query[valid]
  distance_to_tss <- GenomicRanges::distance(query_valid, nearest_tss)
  nearest_gene_ids <- as.character(S4Vectors::mcols(nearest_tss)$gene_id)
  nearest_symbols <- safe_select_org(
    nearest_gene_ids,
    keytype = "ENTREZID",
    columns = c("SYMBOL", "ENTREZID")
  ) |>
    dplyr::distinct(.data$ENTREZID, .data$SYMBOL)

  nearest <- tibble::tibble(
    cpg = no_gene_tbl$cpg[valid],
    MAPINFO = no_gene_tbl$MAPINFO[valid],
    nearest_gene_entrez = nearest_gene_ids,
    nearest_tss_transcript_id = as.character(S4Vectors::mcols(nearest_tss)$tx_name),
    nearest_tss_chr = as.character(GenomicRanges::seqnames(nearest_tss)),
    nearest_tss_position = GenomicRanges::start(nearest_tss),
    nearest_tss_strand = as.character(GenomicRanges::strand(nearest_tss)),
    nearest_tss_distance = as.integer(distance_to_tss)
  ) |>
    dplyr::left_join(nearest_symbols, by = c("nearest_gene_entrez" = "ENTREZID")) |>
    dplyr::mutate(
      nearest_gene = .data$SYMBOL,
      nearest_tss_orientation = dplyr::case_when(
        .data$MAPINFO == .data$nearest_tss_position ~ "at_tss",
        .data$nearest_tss_strand == "+" & .data$MAPINFO < .data$nearest_tss_position ~ "upstream",
        .data$nearest_tss_strand == "+" & .data$MAPINFO > .data$nearest_tss_position ~ "downstream",
        .data$nearest_tss_strand == "-" & .data$MAPINFO > .data$nearest_tss_position ~ "upstream",
        .data$nearest_tss_strand == "-" & .data$MAPINFO < .data$nearest_tss_position ~ "downstream",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::select(
      cpg,
      nearest_gene,
      nearest_gene_entrez,
      nearest_tss_transcript_id,
      nearest_tss_chr,
      nearest_tss_position,
      nearest_tss_strand,
      nearest_tss_distance,
      nearest_tss_orientation
    )

  empty_result |>
    dplyr::select(-dplyr::any_of(names(nearest)[-1])) |>
    dplyr::left_join(nearest, by = "cpg")
}

illumina_annotation <- load_illumina_annotation() |>
  dplyr::rename_with(~ paste0(.x, "_illumina"), -cpg)

top_meta <- meta_results |>
  dplyr::semi_join(top30 |> dplyr::select(cpg), by = "cpg") |>
  dplyr::inner_join(
    top30 |> dplyr::select(cpg, top_contrast),
    by = c("cpg", "contrast" = "top_contrast")
  ) |>
  dplyr::select(
    cpg,
    evidence_tier,
    strongest_contrast = contrast,
    strongest_random_effect = random_effect,
    strongest_fdr = fdr,
    strongest_i2 = i2,
    strongest_n_datasets = n_datasets,
    direction_consistency
  )

curation_base <- top30 |>
  dplyr::left_join(top_meta, by = "cpg") |>
  dplyr::left_join(illumina_annotation, by = "cpg") |>
  dplyr::mutate(
    original_annotation = dplyr::coalesce(
      empty_to_na(.data$UCSC_RefGene_Name_illumina),
      empty_to_na(.data$UCSC_RefGene_Name)
    ),
    chromosome = dplyr::coalesce(empty_to_na(.data$chr_illumina), empty_to_na(.data$chr)),
    genomic_position_hg19 = dplyr::coalesce(.data$pos_illumina, .data$pos),
    MAPINFO = .data$genomic_position_hg19,
    strand = dplyr::coalesce(empty_to_na(.data$strand_illumina), empty_to_na(.data$strand)),
    UCSC_RefGene_Name = .data$original_annotation,
    UCSC_RefGene_Group = dplyr::coalesce(
      empty_to_na(.data$UCSC_RefGene_Group_illumina),
      empty_to_na(.data$UCSC_RefGene_Group)
    ),
    gene_region_annotation = .data$UCSC_RefGene_Group,
    UCSC_RefGene_Accession = empty_to_na(.data$UCSC_RefGene_Accession_illumina),
    Relation_to_UCSC_CpG_Island = dplyr::coalesce(
      empty_to_na(.data$Relation_to_Island_illumina),
      empty_to_na(.data$Relation_to_Island)
    ),
    CpG_island_name = empty_to_na(.data$Islands_Name_illumina),
    Regulatory_Feature_Name = empty_to_na(.data$Regulatory_Feature_Name_illumina),
    Regulatory_Feature_Group = dplyr::coalesce(
      empty_to_na(.data$Regulatory_Feature_Group_illumina),
      empty_to_na(.data$Regulatory_Feature_Group)
    ),
    DMR = dplyr::coalesce(empty_to_na(.data$DMR_illumina), empty_to_na(.data$DMR)),
    Enhancer = dplyr::coalesce(empty_to_na(.data$Enhancer_illumina), empty_to_na(.data$Enhancer)),
    HMM_Island = empty_to_na(.data$HMM_Island_illumina),
    Phantom = empty_to_na(.data$Phantom_illumina),
    DHS = dplyr::coalesce(empty_to_na(.data$DHS_illumina), empty_to_na(.data$DHS)),
    annotation_source = dplyr::coalesce(
      empty_to_na(.data$annotation_source_full_illumina),
      empty_to_na(.data$annotation_source)
    ),
    in_450k_annotation = dplyr::coalesce(.data$in_450k_annotation_full_illumina, .data$in_450k_annotation),
    in_epic_annotation = dplyr::coalesce(.data$in_epic_annotation_full_illumina, .data$in_epic_annotation),
    current_gene_symbol = vapply(.data$original_annotation, summarise_gene_mapping, character(1), field = "current_gene_symbol"),
    current_gene_aliases = vapply(.data$original_annotation, summarise_gene_mapping, character(1), field = "current_gene_aliases"),
    gene_symbol_mapping_source = vapply(.data$original_annotation, summarise_gene_mapping, character(1), field = "gene_symbol_mapping_source"),
    strict_subset_membership = .data$cpg %in% strict_unique$cpg,
    old_biological_support_tier = vapply(.data$original_annotation, lookup_old_context, character(1), column = "biological_support_tier"),
    old_biological_support_score = vapply(.data$original_annotation, lookup_old_context, character(1), column = "biological_support_score"),
    old_source_label = vapply(.data$original_annotation, lookup_old_context, character(1), column = "source_label"),
    previously_reviewed_in_old_script13 = !is.na(.data$old_biological_support_tier),
    nearby_regulatory_annotation = paste(
      "Regulatory_Feature_Group=", dplyr::coalesce(.data$Regulatory_Feature_Group, "NA"),
      "; Regulatory_Feature_Name=", dplyr::coalesce(.data$Regulatory_Feature_Name, "NA"),
      "; Enhancer=", dplyr::coalesce(.data$Enhancer, "NA"),
      "; DHS=", dplyr::coalesce(.data$DHS, "NA"),
      "; DMR=", dplyr::coalesce(.data$DMR, "NA"),
      sep = ""
    )
  )

nearest_tss <- find_nearest_tss(curation_base)

manual_fields <- tibble::tibble(
  biological_evidence_class = NA_character_,
  human_masld_evidence = NA_character_,
  human_liver_fibrosis_evidence = NA_character_,
  experimental_masld_evidence = NA_character_,
  experimental_liver_fibrosis_evidence = NA_character_,
  hepatic_metabolism_evidence = NA_character_,
  evidence_directness = NA_character_,
  evidence_source_1 = NA_character_,
  evidence_source_2 = NA_character_,
  evidence_source_3 = NA_character_,
  literature_note = NA_character_,
  final_interpretation_priority = NA_character_,
  include_in_interpretation_panel = NA_character_,
  selection_reason = NA_character_
)

curation <- curation_base |>
  dplyr::left_join(nearest_tss, by = "cpg") |>
  dplyr::bind_cols(manual_fields[rep(1, nrow(curation_base)), ]) |>
  dplyr::transmute(
    statistical_rank,
    cpg,
    original_annotation,
    UCSC_RefGene_Name,
    current_gene_symbol,
    current_gene_aliases,
    gene_symbol_mapping_source,
    candidate_priority_score,
    evidence_tier,
    strict_subset_membership,
    robust_contrasts,
    tier1_contrasts,
    same_direction_contrasts,
    max_n_datasets,
    best_fdr,
    best_neg_log10_fdr,
    mean_i2,
    max_abs_meta_effect,
    top_effect_direction,
    strongest_contrast,
    direction_consistency,
    likely_dataset_sensitive,
    priority_label,
    all5_any_run,
    chromosome,
    genomic_position_hg19,
    MAPINFO,
    strand,
    UCSC_RefGene_Group,
    gene_region_annotation,
    UCSC_RefGene_Accession,
    Relation_to_UCSC_CpG_Island,
    CpG_island_name,
    Regulatory_Feature_Name,
    Regulatory_Feature_Group,
    Enhancer,
    DHS,
    DMR,
    HMM_Island,
    Phantom,
    nearby_regulatory_annotation,
    annotation_source,
    in_450k_annotation,
    in_epic_annotation,
    Probe_rs = dplyr::coalesce(empty_to_na(.data$Probe_rs_illumina), empty_to_na(.data$Probe_rs)),
    Probe_maf = dplyr::coalesce(.data$Probe_maf_illumina, .data$Probe_maf),
    CpG_rs = dplyr::coalesce(empty_to_na(.data$CpG_rs_illumina), empty_to_na(.data$CpG_rs)),
    CpG_maf = dplyr::coalesce(.data$CpG_maf_illumina, .data$CpG_maf),
    SBE_rs = dplyr::coalesce(empty_to_na(.data$SBE_rs_illumina), empty_to_na(.data$SBE_rs)),
    SBE_maf = dplyr::coalesce(.data$SBE_maf_illumina, .data$SBE_maf),
    nearest_gene,
    nearest_gene_entrez,
    nearest_tss_transcript_id,
    nearest_tss_chr,
    nearest_tss_position,
    nearest_tss_strand,
    nearest_tss_distance,
    nearest_tss_orientation,
    previously_reviewed_in_old_script13,
    old_biological_support_tier,
    old_biological_support_score,
    old_source_label,
    biological_evidence_class,
    human_masld_evidence,
    human_liver_fibrosis_evidence,
    experimental_masld_evidence,
    experimental_liver_fibrosis_evidence,
    hepatic_metabolism_evidence,
    evidence_directness,
    evidence_source_1,
    evidence_source_2,
    evidence_source_3,
    literature_note,
    final_interpretation_priority,
    include_in_interpretation_panel,
    selection_reason
  ) |>
  dplyr::arrange(.data$statistical_rank)

if (nrow(curation) != 30 || dplyr::n_distinct(curation$cpg) != 30) {
  stop("Curation output must contain exactly 30 rows and 30 unique CpGs.")
}

readr::write_csv(curation, output_path)

message("Curation input table: ", output_path)
message("Rows: ", nrow(curation))
message("Unique CpGs: ", dplyr::n_distinct(curation$cpg))
message("Direct UCSC gene annotation: ", sum(!is.na(curation$UCSC_RefGene_Name) & curation$UCSC_RefGene_Name != ""))
message("Intergenic/no-gene CpGs: ", sum(is.na(curation$UCSC_RefGene_Name) | curation$UCSC_RefGene_Name == ""))
message("No-gene CpGs with nearest TSS resolved: ", sum((is.na(curation$UCSC_RefGene_Name) | curation$UCSC_RefGene_Name == "") & !is.na(curation$nearest_gene)))
message("Previously represented in old script-13 specific biological table: ", sum(curation$previously_reviewed_in_old_script13))
