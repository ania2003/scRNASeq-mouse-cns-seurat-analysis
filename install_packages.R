cran_packages <- c("dplyr", "ggplot2", "ggrepel", "patchwork", "remotes")
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)

if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")

message("Core packages installed. Azimuth is optional and is not required by the main pipeline.")
