#!/usr/bin/env Rscript

cran_packages <- c(
  "circlize", "clustree", "cowplot", "data.table", "future", "ggplot2",
  "ggpubr", "ggrepel", "harmony", "msigdbr", "patchwork", "pheatmap",
  "remotes", "scales", "Seurat", "tidyverse", "viridis"
)

bioc_packages <- c(
  "celldex", "clusterProfiler", "ComplexHeatmap", "DESeq2", "DropletUtils",
  "enrichplot", "fgsea", "org.Mm.eg.db", "scater", "SingleR"
)

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

message(
  "CellChat, monocle3, SeuratWrappers, and nichenetr may require installation ",
  "from their current upstream repositories. See environment/README.md."
)
