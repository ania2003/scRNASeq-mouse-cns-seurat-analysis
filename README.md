# Mouse CNS single-cell RNA-seq analysis with Seurat

This repository contains a cleaned and reproducible Seurat workflow for the PanglaoDB dataset `SRA667466_SRS3059954`. The analysis covers quality control, MAD-based cell filtering, log normalization, optional cell-cycle scoring, highly variable gene selection, PCA, comparison of clustering solutions, marker discovery, and marker-supported cell-type annotation.

## Main improvements over the original script

- Replaced absolute Windows paths with relative project paths.
- Moved package installation outside the analysis script.
- Removed duplicated cell-cycle, PCA, clustering, marker, and annotation blocks.
- Fixed inconsistent clustering dimensions and an accidental use of the wrong object in the 30-PC solution.
- Replaced repeated code with reusable functions.
- Added systematic output names and automatic saving of plots, tables, the final Seurat object, and `sessionInfo()`.
- Preserved the complete original script in `legacy/`.

## Repository structure

```text
.
├── analysis/run_analysis.R       # Complete cleaned pipeline
├── R/functions.R                 # Reusable analysis functions
├── config/config.R               # Parameters and clustering solutions
├── data/raw/                     # Local input files (not tracked by Git)
├── figures/                      # Generated plots
├── results/                      # Generated tables and Seurat object
├── legacy/                       # Original script
├── install_packages.R
├── .gitignore
├── LICENSE
└── README.md
```

## Required input

Place this file in `data/raw/`:

```text
SRA667466_SRS3059954.sparse.RData
```

Cell-cycle scoring is optional. To enable it, also place:

```text
cc.genes.mouse.RDS
```

The cell-cycle file must contain two named elements: `s.genes` and `g2m.genes`.

## Run the analysis

Open R or a terminal in the repository root and run:

```r
source("install_packages.R")
source("analysis/run_analysis.R")
```

or:

```bash
Rscript install_packages.R
Rscript analysis/run_analysis.R
```

## Clustering comparison

The default configuration evaluates:

- 10 PCs, resolution 0.5
- 20 PCs, resolution 0.5
- 20 PCs, resolution 0.8
- 25 PCs, resolution 0.5
- 30 PCs, resolution 0.5

The final solution is set to **25 PCs and resolution 0.5**, matching the selection in the original project. Change `final_solution` in `config/config.R` to select another solution.

## Important interpretation note

The generated annotations are based on overlap between cluster markers and a curated mouse CNS marker panel. They are intended as annotation support and must be reviewed using marker specificity, expression plots, biological context, and potential technical biases.

## Data availability

Large source data files are intentionally excluded from Git. See `data/README.md` for placement instructions.
