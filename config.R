config <- list(
  data_file = file.path("data", "raw", "SRA667466_SRS3059954.sparse.RData"),
  cell_cycle_file = file.path("data", "raw", "cc.genes.mouse.RDS"),
  project_name = "SRA667466_SRS3059954",
  seed = 1234,
  min_cells = 3,
  min_features = 200,
  qc_nmads = 3,
  max_mito_percent = 20,
  variable_features = 2000,
  pca_components = 50,
  clustering_solutions = data.frame(
    name = c("pc10_res05", "pc20_res05", "pc20_res08", "pc25_res05", "pc30_res05"),
    n_pcs = c(10, 20, 20, 25, 30),
    resolution = c(0.5, 0.5, 0.8, 0.5, 0.5),
    stringsAsFactors = FALSE
  ),
  final_solution = "pc25_res05",
  marker_min_pct = 0.30,
  marker_min_diff_pct = 0.15,
  marker_logfc = 0.50,
  marker_adjusted_p = 0.05,
  heatmap_cells_per_cluster = 100,
  output_dir = "results",
  figure_dir = "figures"
)
