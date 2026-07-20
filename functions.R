upper_mad_cutoff <- function(x, nmads = 3) {
  as.numeric(stats::median(x, na.rm = TRUE) + nmads * stats::mad(x, na.rm = TRUE))
}

lower_mad_cutoff <- function(x, nmads = 3) {
  as.numeric(stats::median(x, na.rm = TRUE) - nmads * stats::mad(x, na.rm = TRUE))
}

save_plot <- function(plot, filename, figure_dir, width = 8, height = 6) {
  ggplot2::ggsave(file.path(figure_dir, filename), plot = plot,
                  width = width, height = height, dpi = 300)
}

load_sparse_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path,
         "\nPlace the PanglaoDB .RData file in data/raw/.")
  }
  loaded_names <- load(path)
  if (length(loaded_names) != 1L) {
    stop("Expected exactly one object in the .RData file; found: ",
         paste(loaded_names, collapse = ", "))
  }
  matrix_object <- get(loaded_names[[1]])
  if (!inherits(matrix_object, "Matrix")) {
    warning("Loaded object is not a Matrix sparse matrix.")
  }
  matrix_object
}

clean_gene_names <- function(count_matrix) {
  gene_names <- sub("_[^_]*$", "", rownames(count_matrix))
  rownames(count_matrix) <- make.names(gene_names, unique = TRUE)
  count_matrix
}

add_qc_metrics <- function(object) {
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = "^mt\\.")
  object[["percent.rbp"]] <- Seurat::PercentageFeatureSet(object, pattern = "^Rp[ls]")
  object
}

calculate_qc_thresholds <- function(object, nmads = 3, min_features = 200,
                                    max_mito_percent = 20) {
  feature_low <- max(min_features, floor(lower_mad_cutoff(object$nFeature_RNA, nmads)))
  feature_high <- ceiling(upper_mad_cutoff(object$nFeature_RNA, nmads))
  count_high <- ceiling(upper_mad_cutoff(object$nCount_RNA, nmads))
  mito_high <- min(max_mito_percent,
                   max(5, upper_mad_cutoff(object$percent.mt, nmads)))
  data.frame(
    parameter = c("nFeature_RNA_min", "nFeature_RNA_max",
                  "nCount_RNA_max", "percent.mt_max"),
    value = c(feature_low, feature_high, count_high, mito_high)
  )
}

filter_cells <- function(object, thresholds) {
  value <- setNames(thresholds$value, thresholds$parameter)
  subset(
    object,
    subset = nFeature_RNA >= value[["nFeature_RNA_min"]] &
      nFeature_RNA <= value[["nFeature_RNA_max"]] &
      nCount_RNA <= value[["nCount_RNA_max"]] &
      percent.mt <= value[["percent.mt_max"]]
  )
}

score_cell_cycle_if_available <- function(object, path) {
  if (!file.exists(path)) {
    object$S.Score <- NA_real_
    object$G2M.Score <- NA_real_
    object$Phase <- "not_scored"
    message("Cell-cycle file not found; scoring skipped: ", path)
    return(object)
  }
  genes <- readRDS(path)
  required <- c("s.genes", "g2m.genes")
  if (!all(required %in% names(genes))) {
    stop("Cell-cycle file must contain s.genes and g2m.genes.")
  }
  Seurat::CellCycleScoring(
    object,
    s.features = intersect(genes$s.genes, rownames(object)),
    g2m.features = intersect(genes$g2m.genes, rownames(object)),
    set.ident = FALSE
  )
}

run_clustering_solution <- function(object, n_pcs, resolution, seed = 1234) {
  dims <- seq_len(n_pcs)
  object <- Seurat::FindNeighbors(object, dims = dims, k.param = 20, verbose = FALSE)
  object <- Seurat::FindClusters(object, resolution = resolution,
                                 random.seed = seed, verbose = FALSE)
  object <- Seurat::RunUMAP(object, dims = dims, seed.use = seed, verbose = FALSE)
  object <- Seurat::RunTSNE(object, dims = dims, seed.use = seed, verbose = FALSE)
  object
}

summarize_clusters <- function(object, solution_name) {
  out <- as.data.frame(table(object$seurat_clusters), stringsAsFactors = FALSE)
  names(out) <- c("cluster", "n_cells")
  out$solution <- solution_name
  out[, c("solution", "cluster", "n_cells")]
}

get_pca_loadings <- function(object, n = 20) {
  loadings <- Seurat::Loadings(object, reduction = "pca")
  components <- intersect(c("PC_1", "PC_2"), colnames(loadings))
  dplyr::bind_rows(lapply(components, function(component) {
    values <- loadings[, component]
    positive <- sort(values, decreasing = TRUE)[seq_len(min(n, length(values)))]
    negative <- sort(values, decreasing = FALSE)[seq_len(min(n, length(values)))]
    dplyr::bind_rows(
      data.frame(PC = component, direction = "positive",
                 gene = names(positive), loading = unname(positive)),
      data.frame(PC = component, direction = "negative",
                 gene = names(negative), loading = unname(negative))
    )
  }))
}

clean_marker_table <- function(markers) {
  markers |>
    dplyr::filter(
      !grepl("^(Rpl|Rps|Mrpl|Mrps|mt\\.|Mt-)", gene),
      !grepl("^RP[0-9]+\\.", gene),
      !grepl("^AC[0-9]+\\.", gene),
      !grepl("^Gm[0-9]+", gene),
      !grepl("Rik$", gene)
    ) |>
    dplyr::mutate(
      pct.diff = pct.1 - pct.2,
      marker_score = avg_log2FC * pct.diff
    ) |>
    dplyr::arrange(cluster, dplyr::desc(marker_score), dplyr::desc(avg_log2FC))
}

known_marker_panel <- function() {
  list(
    Interneurons = c("Gad1", "Gad2", "Slc32a1", "Dlx1", "Dlx2", "Pvalb", "Sst", "Vip", "Reln", "Npy", "Cck"),
    Astrocytes = c("Aldoc", "Aqp4", "Gja1", "Slc1a2", "Slc1a3", "Gfap", "S100b", "Fabp7", "Clu", "Mt1"),
    Oligodendrocytes = c("Mbp", "Plp1", "Mog", "Mag", "Mobp", "Ermn", "Cldn11", "Opalin", "Mal"),
    Microglia = c("P2ry12", "Tmem119", "Cx3cr1", "Csf1r", "Hexb", "Aif1", "C1qa", "C1qb", "C1qc", "Tyrobp", "Lyz2", "Ms4a7"),
    Endothelial_cells = c("Cldn5", "Pecam1", "Kdr", "Flt1", "Slco1a4", "Igfbp7", "Rbp7"),
    Oligodendrocyte_progenitor_cells = c("Pdgfra", "Cspg4", "Sox10", "Olig1", "Olig2", "Enpp6", "Bcan", "Nkx2.2"),
    Neurons = c("Snap25", "Syt1", "Syp", "Rbfox3", "Tubb3", "Map2", "Slc17a7", "Slc17a6", "Camk2a", "Camk2b", "Ncdn", "Nrgn", "Pcp4", "Adcy1", "Cplx2", "Olfm1"),
    Cajal_Retzius_cells = c("Reln", "Trp73", "Lhx5", "Calb2", "Ebf2", "Cux2", "Cxcl14"),
    Smooth_muscle_cells = c("Acta2", "Tagln", "Myl9", "Myh11", "Cnn1", "Rgs5", "Pdgfrb", "Vtn"),
    Ependymal_cells = c("Foxj1", "Tmem212", "Pifo", "Dnah5", "Cfap43", "Hydin", "Ttr")
  )
}

build_annotation_table <- function(top_markers, all_clusters) {
  marker_table <- utils::stack(known_marker_panel())
  names(marker_table) <- c("gene", "candidate_cell_type")
  matches <- top_markers |>
    dplyr::inner_join(marker_table, by = "gene") |>
    dplyr::group_by(cluster, candidate_cell_type) |>
    dplyr::summarise(n_marker_genes = dplyr::n(),
                     matched_genes = paste(unique(gene), collapse = ", "),
                     .groups = "drop") |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(n_marker_genes, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  top_compact <- top_markers |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(genes_for_annotation = paste(gene, collapse = ", "), .groups = "drop")
  data.frame(cluster = as.character(all_clusters)) |>
    dplyr::left_join(dplyr::mutate(top_compact, cluster = as.character(cluster)), by = "cluster") |>
    dplyr::left_join(dplyr::mutate(matches, cluster = as.character(cluster)), by = "cluster") |>
    dplyr::mutate(
      final_cell_type = dplyr::coalesce(as.character(candidate_cell_type), "Unknown / uncertain"),
      annotation_confidence = dplyr::if_else(is.na(candidate_cell_type), "Low", "Medium")
    )
}
