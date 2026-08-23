#!/usr/bin/env Rscript


rm(list = ls())
gc()


source("scripts/_common.R")
base_dir <- project_root()
setwd(base_dir)
cat("Working directory:", getwd(), "\n")


dirs <- make_output_dirs(base_dir)

output_dir <- file.path(dirs$figures, "Figure6_panels")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


condition_colors <- c(
  "Ctrl" = "#6B9CD3",
  "NO"   = "#8B5AA3"
)

celltype_colors <- c(
  "PTC-1"          = "#2E5A87",
  "PTC-2"          = "#4A7BB0",
  "PTC-3"          = "#E74C3C",
  "Macrophage"     = "#8B5AA3",
  "NK cell"        = "#1ABC9C",
  "Medullary cell" = "#7FBAC4",
  "DCTC"           = "#3498DB",
  "VEC"            = "#2ECC71",
  "Neutrophil"     = "#F39C12",
  "CDPC"           = "#5DADE2",
  "Fibroblast"     = "#E67E22",
  "B cell"         = "#16A085"
)

feature_colors <- c("lightgrey", "#2E5A87", "#B2182B")
pseudotime_colors <- c("lightgrey", "#FEE08B", "#D73027")


theme_publication <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size + 2, hjust = 0.5, color = "grey40"),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      legend.title = element_text(size = base_size + 1, face = "bold"),
      legend.text = element_text(size = base_size),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 2, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}


log_message <- function(msg, log_file = NULL) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  message(paste(timestamp, msg))
}


suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggplot2)
  library(monocle3)
  library(SeuratWrappers)
  library(viridis)
})

log_message("Starting Trajectory Analysis")


log_message("Loading data...")

seurat_obj <- readRDS(file.path(dirs$processed, "seurat_integrated_annotated.rds"))
message(sprintf("Loaded: %d cells, %d genes", ncol(seurat_obj), nrow(seurat_obj)))
message("Available reductions: ", paste(names(seurat_obj@reductions), collapse = ", "))


trajectory_celltypes <- c("PTC-1", "PTC-2", "PTC-3", "DCTC", "CDPC", "Medullary cell")
trajectory_celltypes <- trajectory_celltypes[trajectory_celltypes %in% unique(seurat_obj$cell_type)]

message("Trajectory cell types: ", paste(trajectory_celltypes, collapse = ", "))

seurat_subset <- subset(seurat_obj, subset = cell_type %in% trajectory_celltypes)
message(sprintf("Subset: %d cells", ncol(seurat_subset)))

trajectory_colors <- celltype_colors[names(celltype_colors) %in% trajectory_celltypes]


log_message("Running Monocle3 trajectory analysis...")


cds <- as.cell_data_set(seurat_subset)
cds@colData$cell_type <- seurat_subset$cell_type
cds@colData$sample <- seurat_subset$sample


if ("umap" %in% names(seurat_subset@reductions)) {
  reducedDims(cds)[["UMAP"]] <- Embeddings(seurat_subset, "umap")
  message("✓ Using Harmony-integrated UMAP")
} else {
  stop("ERROR: 'umap' reduction not found!")
}


cds <- cluster_cells(cds, reduction_method = "UMAP", verbose = FALSE)
cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)


ptc1_cells <- colnames(cds)[cds@colData$cell_type == "PTC-1"]
if (length(ptc1_cells) > 0) {
  cds <- order_cells(cds, root_cells = head(ptc1_cells, 10))
  message("✓ Root set to PTC-1 cells")
}


seurat_subset$pseudotime <- pseudotime(cds)[colnames(seurat_subset)]


log_message("Creating Trajectory by Condition Plot...")

p_traj_condition <- plot_cells(cds,
                               color_cells_by = "cell_type",
                               label_cell_groups = FALSE,
                               label_leaves = FALSE,
                               label_branch_points = FALSE,
                               label_roots = FALSE,
                               label_principal_points = FALSE,
                               graph_label_size = 0,
                               cell_size = 0.8) +
  scale_color_manual(values = trajectory_colors) +
  facet_wrap(~sample, ncol = 2) +
  labs(title = "Trajectory Analysis by Condition",
       x = "UMAP 1", y = "UMAP 2",
       color = "Cell Type") +
  theme_publication(base_size = 14) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "Fig6A_trajectory_by_condition.pdf"),
       p_traj_condition, width = 14, height = 7, device = cairo_pdf)
ggsave(file.path(output_dir, "Fig6A_trajectory_by_condition.png"),
       p_traj_condition, width = 14, height = 7, dpi = 300)

ggsave(file.path(dirs$trajectory, "Fig6A_trajectory_by_condition.pdf"),
       p_traj_condition, width = 14, height = 7, device = cairo_pdf)

log_message("✓ Trajectory by Condition saved!")


log_message("Creating Pseudotime UMAP...")

p_pseudotime <- FeaturePlot(seurat_subset,
                            features = "pseudotime",
                            reduction = "umap",
                            order = TRUE,
                            pt.size = 0.5) +
  scale_color_gradientn(colors = pseudotime_colors, na.value = "lightgrey") +
  labs(title = "Pseudotime Trajectory",
       x = "UMAP 1", y = "UMAP 2",
       color = "Pseudotime") +
  theme_publication(base_size = 14)

ggsave(file.path(output_dir, "Fig6B_pseudotime_umap.pdf"),
       p_pseudotime, width = 9, height = 7, device = cairo_pdf)
ggsave(file.path(output_dir, "Fig6B_pseudotime_umap.png"),
       p_pseudotime, width = 9, height = 7, dpi = 300)

log_message("✓ Pseudotime UMAP saved!")


log_message("Calculating Fate Scores...")

fate_genes <- list(
  EMT = c("Vim", "Snai1", "Snai2", "Twist1", "Zeb1", "Zeb2",
          "Fn1", "Col1a1", "Col1a2", "Acta2"),
  Regeneration = c("Pcna", "Mki67", "Ccnd1", "Ccne1", "Cdk4",
                   "Slc34a1", "Lrp2", "Aqp1"),
  Apoptosis = c("Bax", "Bad", "Casp3", "Casp9", "Bcl2", "Bcl2l1")
)

for (fate in names(fate_genes)) {
  genes_present <- fate_genes[[fate]][fate_genes[[fate]] %in% rownames(seurat_subset)]
  if (length(genes_present) >= 3) {
    seurat_subset <- AddModuleScore(seurat_subset,
                                    features = list(genes_present),
                                    name = paste0(fate, "_score"),
                                    search = TRUE)
    message(sprintf("  Added %s score (%d genes)", fate, length(genes_present)))
  }
}


log_message("Creating Fate Score UMAPs...")

fate_score_cols <- grep("_score1$", colnames(seurat_subset@meta.data), value = TRUE)

if (length(fate_score_cols) > 0) {
  p_fate_list <- list()

  for (col in fate_score_cols) {
    fate_name <- gsub("_score1$", "", col)

    p <- FeaturePlot(seurat_subset,
                     features = col,
                     reduction = "umap",
                     order = TRUE,
                     pt.size = 0.5) +
      scale_color_gradientn(colors = pseudotime_colors) +
      labs(title = sprintf("%s Score", fate_name),
           color = "Score",
           x = "UMAP 1", y = "UMAP 2") +
      theme_publication(base_size = 14)

    p_fate_list[[fate_name]] <- p
  }

  p_fate_combined <- wrap_plots(p_fate_list, ncol = 3)

  ggsave(file.path(output_dir, "Fig6_fate_scores_umap.pdf"),
         p_fate_combined, width = 16, height = 6, device = cairo_pdf)
  ggsave(file.path(output_dir, "Fig6_fate_scores_umap.png"),
         p_fate_combined, width = 16, height = 6, dpi = 300)

  log_message("✓ Fate score UMAPs saved!")
}


log_message("Creating Pseudotime Density Plot...")

pseudotime_df <- data.frame(
  pseudotime = seurat_subset$pseudotime,
  sample = seurat_subset$sample,
  cell_type = seurat_subset$cell_type
) %>% filter(!is.na(pseudotime))

p_density <- ggplot(pseudotime_df, aes(x = pseudotime, fill = sample)) +
  geom_density(alpha = 0.6, linewidth = 0.8) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Pseudotime Distribution by Condition",
       x = "Pseudotime", y = "Density", fill = "Condition") +
  theme_publication(base_size = 14) +
  theme(legend.position = c(0.85, 0.85))

ggsave(file.path(output_dir, "Fig6_pseudotime_density.pdf"),
       p_density, width = 8, height = 6, device = cairo_pdf)
ggsave(file.path(output_dir, "Fig6_pseudotime_density.png"),
       p_density, width = 8, height = 6, dpi = 300)

log_message("✓ Pseudotime density saved!")


log_message("Creating Pseudotime by Cell Type Plot...")

p_pseudotime_celltype <- ggplot(pseudotime_df, aes(x = cell_type, y = pseudotime, fill = cell_type)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = trajectory_colors) +
  facet_wrap(~sample, ncol = 2) +
  labs(title = "Pseudotime by Cell Type and Condition",
       x = NULL, y = "Pseudotime") +
  theme_publication(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(file.path(output_dir, "Fig6_pseudotime_by_celltype.pdf"),
       p_pseudotime_celltype, width = 12, height = 6, device = cairo_pdf)
ggsave(file.path(output_dir, "Fig6_pseudotime_by_celltype.png"),
       p_pseudotime_celltype, width = 12, height = 6, dpi = 300)

log_message("✓ Pseudotime by cell type saved!")


log_message("Saving data...")

saveRDS(seurat_subset, file.path(dirs$processed, "seurat_trajectory_subset.rds"))
saveRDS(cds, file.path(dirs$processed, "monocle3_cds.rds"))


seurat_obj$pseudotime <- NA
seurat_obj$pseudotime[colnames(seurat_subset)] <- seurat_subset$pseudotime


for (col in fate_score_cols) {
  seurat_obj[[col]] <- NA
  seurat_obj[[col]][colnames(seurat_subset)] <- seurat_subset[[col]]
}

saveRDS(seurat_obj, file.path(dirs$processed, "seurat_with_trajectory.rds"))

log_message("✓ Data saved!")


message("\n========================================")
message("TRAJECTORY ANALYSIS COMPLETE!")
message("========================================")
message("\nOutput files:")
message("  Fig6A_trajectory_by_condition.pdf/png")
message("  Fig6B_pseudotime_umap.pdf/png")
message("  Fig6_fate_scores_umap.pdf/png")
message("  Fig6_pseudotime_density.pdf/png")
message("  Fig6_pseudotime_by_celltype.pdf/png")
message("\nData files:")
message("  seurat_trajectory_subset.rds")
message("  monocle3_cds.rds")
message("  seurat_with_trajectory.rds")
message("========================================")
message("All plots use Harmony-integrated UMAP!")
message("========================================")

sessionInfo()
