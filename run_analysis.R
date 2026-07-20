source(file.path("config", "config.R"))
source(file.path("R", "functions.R"))

required <- c("Seurat", "dplyr", "ggplot2", "ggrepel", "patchwork", "Matrix")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nRun: Rscript install_packages.R")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

set.seed(config$seed)
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$figure_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Import and object creation
sm <- load_sparse_matrix(config$data_file)
sm <- clean_gene_names(sm)
pbmc <- CreateSeuratObject(
  counts = sm,
  project = config$project_name,
  min.cells = config$min_cells,
  min.features = config$min_features
)
pbmc <- add_qc_metrics(pbmc)

# 2. Quality control
qc_before <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"),
                     ncol = 4, pt.size = 0)
save_plot(qc_before, "01_qc_violin_before_filtering.png", config$figure_dir, 12, 5)

scatter_mt <- FeatureScatter(pbmc, "nCount_RNA", "percent.mt")
scatter_features <- FeatureScatter(pbmc, "nCount_RNA", "nFeature_RNA")
scatter_ribo <- FeatureScatter(pbmc, "nCount_RNA", "percent.rbp")
save_plot(scatter_mt, "02_qc_counts_vs_mito.png", config$figure_dir)
save_plot(scatter_features, "03_qc_counts_vs_features.png", config$figure_dir)
save_plot(scatter_ribo, "04_qc_counts_vs_ribosomal.png", config$figure_dir)

thresholds <- calculate_qc_thresholds(pbmc, config$qc_nmads,
                                      config$min_features,
                                      config$max_mito_percent)
write.csv(thresholds, file.path(config$output_dir, "qc_thresholds.csv"), row.names = FALSE)

n_before <- ncol(pbmc)
pbmc <- filter_cells(pbmc, thresholds)
qc_summary <- data.frame(
  cells_before = n_before,
  cells_after = ncol(pbmc),
  cells_removed = n_before - ncol(pbmc),
  percent_removed = 100 * (n_before - ncol(pbmc)) / n_before
)
write.csv(qc_summary, file.path(config$output_dir, "qc_cell_summary.csv"), row.names = FALSE)

qc_after <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"),
                    ncol = 4, pt.size = 0)
save_plot(qc_after, "05_qc_violin_after_filtering.png", config$figure_dir, 12, 5)

# 3. Normalization and optional cell-cycle scoring
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
pbmc <- score_cell_cycle_if_available(pbmc, config$cell_cycle_file)
mean_expression <- Matrix::rowMeans(LayerData(pbmc, assay = "RNA", layer = "data"))
mean_expression <- sort(mean_expression, decreasing = TRUE)
write.csv(data.frame(gene = names(mean_expression), mean_expression = unname(mean_expression)),
          file.path(config$output_dir, "mean_gene_expression.csv"), row.names = FALSE)

# 4. Variable genes, scaling and PCA
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",
                             nfeatures = config$variable_features, verbose = FALSE)
top10_variable <- head(VariableFeatures(pbmc), 10)
variable_plot <- LabelPoints(VariableFeaturePlot(pbmc), points = top10_variable, repel = TRUE)
save_plot(variable_plot, "06_variable_features.png", config$figure_dir)
write.csv(data.frame(gene = top10_variable),
          file.path(config$output_dir, "top10_variable_genes.csv"), row.names = FALSE)

pbmc <- ScaleData(pbmc, features = rownames(pbmc), verbose = FALSE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc),
               npcs = config$pca_components, verbose = FALSE)
write.csv(get_pca_loadings(pbmc, n = 20),
          file.path(config$output_dir, "pca_top_loadings_pc1_pc2.csv"), row.names = FALSE)
save_plot(DimPlot(pbmc, reduction = "pca"), "07_pca_cells.png", config$figure_dir)
save_plot(ElbowPlot(pbmc, ndims = config$pca_components), "08_pca_elbow_plot.png", config$figure_dir)

variance <- pbmc[["pca"]]@stdev^2
variance <- variance / sum(variance)
pc_75_percent <- which(cumsum(variance) >= 0.75)[1]
pc_recommendation <- min(20, max(5, pc_75_percent))
write.csv(data.frame(pc_75_percent = pc_75_percent,
                     recommended_pc_count = pc_recommendation),
          file.path(config$output_dir, "pca_dimension_recommendation.csv"), row.names = FALSE)

# 5. Compare clustering solutions
solutions <- list()
cluster_sizes <- list()
for (i in seq_len(nrow(config$clustering_solutions))) {
  row <- config$clustering_solutions[i, ]
  message("Running ", row$name, ": ", row$n_pcs, " PCs, resolution ", row$resolution)
  solution <- run_clustering_solution(pbmc, row$n_pcs, row$resolution, config$seed)
  solutions[[row$name]] <- solution
  cluster_sizes[[row$name]] <- summarize_clusters(solution, row$name)

  umap_plot <- DimPlot(solution, reduction = "umap", label = TRUE, repel = TRUE) +
    NoLegend() + ggtitle(paste0(row$name, ": UMAP"))
  tsne_plot <- DimPlot(solution, reduction = "tsne", label = TRUE, repel = TRUE) +
    NoLegend() + ggtitle(paste0(row$name, ": t-SNE"))
  save_plot(umap_plot, paste0("clustering_", row$name, "_umap.png"), config$figure_dir)
  save_plot(tsne_plot, paste0("clustering_", row$name, "_tsne.png"), config$figure_dir)
}
write.csv(bind_rows(cluster_sizes), file.path(config$output_dir, "cluster_sizes_all_solutions.csv"), row.names = FALSE)

if (!config$final_solution %in% names(solutions)) {
  stop("Final solution not found in clustering solutions: ", config$final_solution)
}
pbmc <- solutions[[config$final_solution]]
Idents(pbmc) <- "seurat_clusters"

# 6. Bias checks for the final solution
for (feature in c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rbp")) {
  p <- VlnPlot(pbmc, features = feature, group.by = "seurat_clusters", pt.size = 0)
  save_plot(p, paste0("final_clusters_", feature, ".png"), config$figure_dir, 10, 5)
}
if ("Phase" %in% colnames(pbmc@meta.data) && any(pbmc$Phase != "not_scored")) {
  phase_table <- pbmc@meta.data |>
    count(seurat_clusters, Phase) |>
    group_by(seurat_clusters) |>
    mutate(percent = 100 * n / sum(n)) |>
    ungroup()
  phase_plot <- ggplot(phase_table, aes(seurat_clusters, percent, fill = Phase)) +
    geom_col() + labs(x = "Cluster", y = "Cells (%)", title = "Cell-cycle phase by cluster") +
    theme_classic()
  save_plot(phase_plot, "final_clusters_cell_cycle.png", config$figure_dir, 10, 5)
}

# 7. Marker discovery
markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = config$marker_min_pct,
  min.diff.pct = config$marker_min_diff_pct,
  logfc.threshold = config$marker_logfc,
  test.use = "wilcox",
  return.thresh = config$marker_adjusted_p
)
markers <- clean_marker_table(markers)
write.csv(markers, file.path(config$output_dir, "all_cluster_markers_cleaned.csv"), row.names = FALSE)

top10 <- markers |>
  filter(p_val_adj < config$marker_adjusted_p, avg_log2FC > 0,
         pct.1 >= config$marker_min_pct, pct.diff >= config$marker_min_diff_pct) |>
  group_by(cluster) |>
  slice_max(marker_score, n = 10, with_ties = FALSE) |>
  ungroup()
top5 <- top10 |> group_by(cluster) |> slice_max(marker_score, n = 5, with_ties = FALSE) |> ungroup()
top3 <- top10 |> group_by(cluster) |> slice_max(marker_score, n = 3, with_ties = FALSE) |> ungroup()
write.csv(top10, file.path(config$output_dir, "top10_markers_per_cluster.csv"), row.names = FALSE)
write.csv(top5, file.path(config$output_dir, "top5_markers_per_cluster.csv"), row.names = FALSE)

heatmap_cells <- WhichCells(pbmc, downsample = config$heatmap_cells_per_cluster)
pbmc_heatmap <- subset(pbmc, cells = heatmap_cells)
heatmap_plot <- DoHeatmap(pbmc_heatmap, features = unique(top3$gene), raster = TRUE) +
  NoLegend() + ggtitle("Top 3 markers per cluster")
save_plot(heatmap_plot, "09_marker_heatmap_top3.png", config$figure_dir, 12, 14)

# 8. Annotation support from known mouse CNS markers
annotation <- build_annotation_table(top10, sort(unique(as.character(pbmc$seurat_clusters))))
write.csv(annotation, file.path(config$output_dir, "cluster_annotation_to_review.csv"), row.names = FALSE)

cluster_to_type <- setNames(annotation$final_cell_type, annotation$cluster)
manual_type <- unname(cluster_to_type[as.character(pbmc$seurat_clusters)])
manual_type[is.na(manual_type)] <- "Unknown / uncertain"
names(manual_type) <- Cells(pbmc)
pbmc <- AddMetaData(pbmc, manual_type, col.name = "manual_celltype")

annotation_umap <- DimPlot(pbmc, reduction = "umap", group.by = "manual_celltype",
                           label = TRUE, repel = TRUE) + NoLegend() +
  ggtitle("UMAP - marker-supported annotation")
annotation_tsne <- DimPlot(pbmc, reduction = "tsne", group.by = "manual_celltype",
                           label = TRUE, repel = TRUE) + NoLegend() +
  ggtitle("t-SNE - marker-supported annotation")
save_plot(annotation_umap, "10_umap_marker_supported_annotation.png", config$figure_dir, 11, 8)
save_plot(annotation_tsne, "11_tsne_marker_supported_annotation.png", config$figure_dir, 11, 8)

# 9. Save the final object and reproducibility information
saveRDS(pbmc, file.path(config$output_dir, "final_seurat_object.rds"))
writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "sessionInfo.txt"))
message("Analysis completed. Results: ", config$output_dir, "; figures: ", config$figure_dir)
