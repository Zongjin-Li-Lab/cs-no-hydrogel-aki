#!/usr/bin/env Rscript


rm(list = ls())
gc()


source("scripts/_common.R")
base_dir <- project_root()
setwd(base_dir)


dirs <- make_output_dirs(base_dir)


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


regulation_colors <- c(
  "Up"   = "#B2182B",
  "Down" = "#2E5A87",
  "NS"   = "grey70"
)


heatmap_colors <- colorRampPalette(c("#2E5A87", "white", "#B2182B"))(100)


feature_colors <- c("lightgrey", "#2E5A87", "#B2182B")


theme_publication <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "grey40"),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 2, color = "black"),
      axis.line = element_line(linewidth = 0.6, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.key.size = unit(0.6, "cm"),
      legend.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}


log_message <- function(msg, log_file = NULL) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  full_msg <- paste(timestamp, msg)
  message(full_msg)
  if (!is.null(log_file)) {
    cat(full_msg, "\n", file = log_file, append = TRUE)
  }
}


suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggplot2)
  library(harmony)
  library(SingleR)
  library(celldex)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(future)
})


n_workers <- min(4L, future::availableCores())
if (n_workers > 1L) future::plan("multisession", workers = n_workers) else future::plan("sequential")
options(future.globals.maxSize = 8000 * 1024^2)


log_file <- file.path(dirs$clustering, "integration_clustering_log.txt")
log_message("Starting Integration and Clustering Analysis", log_file)


log_message("Loading QC-filtered data...", log_file)

seurat_list <- readRDS(file.path(dirs$processed, "seurat_list_qc_filtered.rds"))

log_message(sprintf("Ctrl: %d cells", ncol(seurat_list$Ctrl)), log_file)
log_message(sprintf("NO: %d cells", ncol(seurat_list$NO)), log_file)


log_message("Merging and normalizing data...", log_file)


seurat_merged <- merge(seurat_list[[1]], y = seurat_list[[2]],
                       add.cell.ids = names(seurat_list),
                       project = "CS_NO_AKI")

log_message(sprintf("Merged object: %d cells, %d genes",
                    ncol(seurat_merged), nrow(seurat_merged)), log_file)


seurat_merged <- NormalizeData(seurat_merged, verbose = FALSE)
seurat_merged <- FindVariableFeatures(seurat_merged,
                                      selection.method = "vst",
                                      nfeatures = 3000,
                                      verbose = FALSE)
seurat_merged <- ScaleData(seurat_merged, verbose = FALSE)


log_message("Running PCA...", log_file)

seurat_merged <- RunPCA(seurat_merged, npcs = 50, verbose = FALSE)


p_elbow <- ElbowPlot(seurat_merged, ndims = 50) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(title = "PCA Elbow Plot", subtitle = "Red line indicates selected PCs (30)") +
  theme_publication(base_size = 14)

ggsave(file.path(dirs$clustering, "S01_pca_elbow_plot.pdf"),
       p_elbow, width = 8, height = 6, device = cairo_pdf)


seurat_merged <- RunUMAP(seurat_merged, dims = 1:30, verbose = FALSE,
                         reduction.name = "umap_uncorrected")

p_before_integration <- DimPlot(seurat_merged, reduction = "umap_uncorrected",
                                group.by = "sample", pt.size = 0.3) +
  scale_color_manual(values = condition_colors) +
  labs(title = "UMAP Before Integration", x = "UMAP 1", y = "UMAP 2", color = "Condition") +
  theme_publication(base_size = 14) +
  theme(legend.position = "right")

ggsave(file.path(dirs$clustering, "S02_umap_before_integration.pdf"),
       p_before_integration, width = 9, height = 7, device = cairo_pdf)


log_message("Running Harmony integration...", log_file)

seurat_merged <- RunHarmony(seurat_merged,
                            group.by.vars = "sample",
                            reduction.use = "pca",
                            dims.use = 1:30,
                            plot_convergence = FALSE,
                            verbose = FALSE)


seurat_merged <- RunUMAP(seurat_merged,
                         reduction = "harmony",
                         dims = 1:30,
                         verbose = FALSE)


p_after_integration <- DimPlot(seurat_merged, reduction = "umap",
                               group.by = "sample", pt.size = 0.3) +
  scale_color_manual(values = condition_colors) +
  labs(title = "UMAP After Harmony Integration", x = "UMAP 1", y = "UMAP 2", color = "Condition") +
  theme_publication(base_size = 14) +
  theme(legend.position = "right")

ggsave(file.path(dirs$clustering, "S03_umap_after_integration.pdf"),
       p_after_integration, width = 9, height = 7, device = cairo_pdf)


p_integration_comparison <- p_before_integration + p_after_integration +
  plot_annotation(
    title = "Integration Quality Assessment",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

ggsave(file.path(dirs$clustering, "S04_integration_comparison.pdf"),
       p_integration_comparison, width = 16, height = 7, device = cairo_pdf)
ggsave(file.path(dirs$clustering, "S04_integration_comparison.png"),
       p_integration_comparison, width = 16, height = 7, dpi = 300)


log_message("Performing clustering...", log_file)


seurat_merged <- FindNeighbors(seurat_merged,
                               reduction = "harmony",
                               dims = 1:30,
                               verbose = FALSE)


resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)

for (res in resolutions) {
  seurat_merged <- FindClusters(seurat_merged, resolution = res, verbose = FALSE)
  col_name <- paste0("RNA_snn_res.", res)
  log_message(sprintf("Resolution %.1f: %d clusters", res,
                      length(unique(seurat_merged@meta.data[[col_name]]))), log_file)
}


if (requireNamespace("clustree", quietly = TRUE)) {
  library(clustree)
  p_clustree <- clustree(seurat_merged, prefix = "RNA_snn_res.") +
    labs(title = "Cluster Resolution Analysis") +
    theme_publication(base_size = 12)
  ggsave(file.path(dirs$clustering, "S05_clustree_resolution.pdf"),
         p_clustree, width = 12, height = 10, device = cairo_pdf)
}


seurat_merged <- FindClusters(seurat_merged, resolution = 0.2, verbose = FALSE)


seurat_merged <- JoinLayers(seurat_merged)

n_clusters <- length(unique(seurat_merged$seurat_clusters))
log_message(sprintf("Selected resolution: 0.2 (%d initial clusters)", n_clusters), log_file)


log_message("Defining cell type markers...", log_file)

kidney_markers <- list(
  "PTC" = c("Slc34a1", "Lrp2", "Cubn", "Slc5a2", "Aqp1"),
  "PTC_injured" = c("Havcr1", "Lcn2", "Cd44", "Sox9", "Krt8", "Cxcl10", "Ccl5"),
  "DCTC" = c("Slc12a3", "Pvalb", "Calb1"),
  "CDPC" = c("Aqp2", "Aqp3", "Fxyd4"),
  "Medullary" = c("Slc12a1", "Umod"),
  "VEC" = c("Pecam1", "Cdh5", "Kdr", "Emcn"),
  "Fibroblast" = c("Col1a1", "Dcn", "Pdgfra"),
  "Macrophage" = c("Adgre1", "Cd68", "Csf1r"),
  "M1_Macrophage" = c("Cd86", "Nos2", "Il1b", "Tnf"),
  "M2_Macrophage" = c("Mrc1", "Arg1", "Cd163"),
  "B_cell" = c("Cd19", "Ms4a1", "Cd79a"),
  "Neutrophil" = c("S100a8", "S100a9", "Ly6g"),
  "NK_cell" = c("Klrb1c", "Ncr1", "Nkg7")
)

save(kidney_markers, file = file.path(dirs$clustering, "kidney_cell_type_markers.RData"))


log_message("Creating marker expression visualizations...", log_file)


merge_markers <- c(

  "Slc34a1", "Lrp2", "Aqp1", "Cubn",

  "Havcr1", "Lcn2", "Sox9", "Vcam1",

  "C3", "Ttc36", "Akr1c21", "Inmt", "Acy3", "Ghr", "Ccn1",

  "Slc12a3", "Pvalb",

  "Aqp2", "Aqp3",

  "Slc12a1", "Umod",

  "Pecam1", "Kdr", "Cdh5", "Emcn",

  "Col1a1", "Pdgfra", "Dcn",

  "Adgre1", "Cd68", "Csf1r",

  "Cd86", "Mrc1", "Arg1",

  "Klrb1c", "Ncr1", "Nkg7",

  "Cd79a", "Ms4a1",

  "S100a8", "S100a9"
)

merge_markers <- merge_markers[merge_markers %in% rownames(seurat_merged)]


p_dotplot <- DotPlot(seurat_merged, features = merge_markers) +
  coord_flip() +
  labs(title = "Cluster Marker Expression for Annotation") +
  theme_publication(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 10))

ggsave(file.path(dirs$clustering, "S06_dotplot_markers.pdf"),
       p_dotplot, width = 14, height = 16, device = cairo_pdf)
ggsave(file.path(dirs$clustering, "S06_dotplot_markers.png"),
       p_dotplot, width = 14, height = 16, dpi = 300)


avg_expr <- AverageExpression(seurat_merged, features = merge_markers,
                              group.by = "seurat_clusters")$RNA
avg_expr_scaled <- t(scale(t(as.matrix(avg_expr))))

pdf(file.path(dirs$clustering, "S07_cluster_heatmap.pdf"), width = 12, height = 18)
pheatmap(avg_expr_scaled,
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         color = colorRampPalette(c("#2E5A87", "white", "#B2182B"))(100),
         main = "Average Expression by Cluster",
         fontsize = 12,
         fontsize_row = 10,
         fontsize_col = 12)
dev.off()


log_message("Running automated cell type annotation with SingleR...", log_file)

tryCatch({
  ref_immgen <- celldex::ImmGenData()
  ref_mouse <- celldex::MouseRNAseqData()

  sce <- as.SingleCellExperiment(seurat_merged)

  singler_immgen <- SingleR(test = sce, ref = ref_immgen, labels = ref_immgen$label.main)
  seurat_merged$singler_immgen <- singler_immgen$labels

  singler_mouse <- SingleR(test = sce, ref = ref_mouse, labels = ref_mouse$label.main)
  seurat_merged$singler_mouse <- singler_mouse$labels

  log_message("SingleR annotation completed", log_file)
}, error = function(e) {
  log_message(sprintf("SingleR error: %s - proceeding with manual annotation", e$message), log_file)
})


log_message("Performing manual cell type annotation (12 cell types)...", log_file)


saveRDS(seurat_merged, file.path(dirs$processed, "seurat_integrated_clustered_original.rds"))
log_message("Original clustered object saved", log_file)


cat("\nInitial cluster distribution:\n")
print(table(seurat_merged$seurat_clusters))


seurat_annotated <- seurat_merged


if ("13" %in% levels(seurat_annotated$seurat_clusters)) {
  cells_to_keep <- colnames(seurat_annotated)[seurat_annotated$seurat_clusters != "13"]
  seurat_annotated <- subset(seurat_annotated, cells = cells_to_keep)
  seurat_annotated$seurat_clusters <- droplevels(seurat_annotated$seurat_clusters)
  log_message("Removed cluster 13 (low quality)", log_file)
}


cluster_to_celltype <- c(
  "0"  = "PTC-1",
  "1"  = "PTC-1",
  "2"  = "PTC-2",
  "3"  = "Macrophage",
  "4"  = "NK cell",
  "5"  = "Medullary cell",
  "6"  = "DCTC",
  "7"  = "VEC",
  "8"  = "Neutrophil",
  "9"  = "CDPC",
  "10" = "Fibroblast",
  "11" = "B cell",
  "12" = "Macrophage",
  "14" = "PTC-2",
  "15" = "PTC-3"
)


clusters <- as.character(seurat_annotated$seurat_clusters)
cell_types <- cluster_to_celltype[clusters]
names(cell_types) <- colnames(seurat_annotated)
seurat_annotated <- AddMetaData(seurat_annotated, metadata = cell_types, col.name = "cell_type")


cat("\nCell type distribution (12 cell types):\n")
print(table(seurat_annotated$cell_type))


Idents(seurat_annotated) <- "cell_type"


log_message("Finding markers for annotated cell types...", log_file)

all_markers <- FindAllMarkers(seurat_annotated,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.5,
                              verbose = FALSE)

write.csv(all_markers, file.path(dirs$clustering, "all_celltype_markers.csv"), row.names = FALSE)

top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(top_markers, file.path(dirs$clustering, "top10_markers_per_celltype.csv"), row.names = FALSE)


log_message("Creating publication-quality main figures...", log_file)


p_umap_celltype <- DimPlot(seurat_annotated,
                           reduction = "umap",
                           group.by = "cell_type",
                           label = TRUE, label.size = 5,
                           repel = TRUE, pt.size = 0.3) +
  scale_color_manual(values = celltype_colors) +
  labs(title = "Cell Type Annotations (12 Cell Types)",
       x = "UMAP 1", y = "UMAP 2") +
  theme_publication(base_size = 14) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12))

ggsave(file.path(dirs$clustering, "Fig01_umap_celltypes.pdf"),
       p_umap_celltype, width = 11, height = 9, device = cairo_pdf)
ggsave(file.path(dirs$clustering, "Fig01_umap_celltypes.png"),
       p_umap_celltype, width = 11, height = 9, dpi = 300)


p_umap_split <- DimPlot(seurat_annotated,
                        reduction = "umap",
                        group.by = "cell_type",
                        split.by = "sample",
                        label = TRUE, label.size = 4,
                        repel = TRUE, pt.size = 0.3) +
  scale_color_manual(values = celltype_colors) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_publication(base_size = 14) +
  theme(legend.title = element_blank(),
        strip.text = element_text(size = 16, face = "bold"))

ggsave(file.path(dirs$clustering, "Fig02_umap_split_condition.pdf"),
       p_umap_split, width = 16, height = 7, device = cairo_pdf)
ggsave(file.path(dirs$clustering, "Fig02_umap_split_condition.png"),
       p_umap_split, width = 16, height = 7, dpi = 300)


proportion_data <- seurat_annotated@meta.data %>%
  group_by(sample, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = n / sum(n) * 100) %>%
  ungroup()

p_proportion <- ggplot(proportion_data,
                       aes(x = cell_type, y = proportion, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Cell Type Proportions: Control vs NO Treatment",
       x = NULL, y = "Proportion (%)", fill = "Condition") +
  theme_publication(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.position = "top")

ggsave(file.path(dirs$clustering, "Fig03_celltype_proportions.pdf"),
       p_proportion, width = 12, height = 7, device = cairo_pdf)
ggsave(file.path(dirs$clustering, "Fig03_celltype_proportions.png"),
       p_proportion, width = 12, height = 7, dpi = 300)

write.csv(proportion_data, file.path(dirs$clustering, "celltype_proportions.csv"), row.names = FALSE)


key_features <- c("Slc34a1", "Havcr1", "Ghr", "Adgre1", "Cd86", "Mrc1")
key_features <- key_features[key_features %in% rownames(seurat_annotated)]

p_feature <- FeaturePlot(seurat_annotated,
                         reduction = "umap",
                         features = key_features,
                         cols = feature_colors,
                         order = TRUE, pt.size = 0.3, ncol = 3) &
  theme_publication(base_size = 12) &
  labs(x = "UMAP 1", y = "UMAP 2")

ggsave(file.path(dirs$clustering, "Fig04_feature_plots.pdf"),
       p_feature, width = 14, height = 10, device = cairo_pdf)
ggsave(file.path(dirs$clustering, "Fig04_feature_plots.png"),
       p_feature, width = 14, height = 10, dpi = 300)


top5_markers <- top_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

p_heatmap <- DoHeatmap(seurat_annotated, features = top5_markers,
                       group.by = "cell_type",
                       group.colors = celltype_colors,
                       size = 4, angle = 45) +
  scale_fill_gradientn(colors = c("#2E5A87", "white", "#B2182B")) +
  labs(title = "Top Marker Genes by Cell Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(file.path(dirs$clustering, "Fig05_marker_heatmap.pdf"),
       p_heatmap, width = 16, height = 14, device = cairo_pdf)
ggsave(file.path(dirs$clustering, "Fig05_marker_heatmap.png"),
       p_heatmap, width = 16, height = 14, dpi = 300)


log_message("Creating supplementary figures...", log_file)


p_umap_condition <- DimPlot(seurat_annotated,
                            reduction = "umap",
                            group.by = "sample",
                            pt.size = 0.3) +
  scale_color_manual(values = condition_colors) +
  labs(title = "UMAP by Condition", x = "UMAP 1", y = "UMAP 2", color = "Condition") +
  theme_publication(base_size = 14) +
  theme(legend.position = "right")

ggsave(file.path(dirs$clustering, "S08_umap_condition.pdf"),
       p_umap_condition, width = 9, height = 7, device = cairo_pdf)


p_proportion_stacked <- ggplot(proportion_data,
                               aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = celltype_colors) +
  labs(title = "Cell Type Composition by Condition",
       x = NULL, y = "Proportion (%)", fill = "Cell Type") +
  theme_publication(base_size = 14) +
  theme(legend.position = "right")

ggsave(file.path(dirs$clustering, "S09_proportion_stacked.pdf"),
       p_proportion_stacked, width = 9, height = 7, device = cairo_pdf)


fibroblast_markers <- c("Col1a1", "Pdgfra", "Dcn")
fibroblast_markers <- fibroblast_markers[fibroblast_markers %in% rownames(seurat_annotated)]

if (length(fibroblast_markers) >= 2) {
  p_fibroblast <- FeaturePlot(seurat_annotated,
                              reduction = "umap",
                              features = fibroblast_markers,
                              ncol = 3, order = TRUE, pt.size = 0.3) &
    scale_color_gradientn(colors = feature_colors) &
    theme_publication(base_size = 12) &
    labs(x = "UMAP 1", y = "UMAP 2")

  ggsave(file.path(dirs$clustering, "S10_fibroblast_evidence.pdf"),
         p_fibroblast, width = 14, height = 5, device = cairo_pdf)
}


log_message("Saving data and generating report...", log_file)


saveRDS(seurat_annotated, file.path(dirs$processed, "seurat_integrated_annotated.rds"))


celltype_info <- data.frame(
  cell_type = names(table(seurat_annotated$cell_type)),
  n_cells = as.vector(table(seurat_annotated$cell_type)),
  color = celltype_colors[names(table(seurat_annotated$cell_type))]
)
write.csv(celltype_info, file.path(dirs$clustering, "celltype_info.csv"), row.names = FALSE)


write.csv(seurat_annotated@meta.data, file.path(dirs$clustering, "cell_metadata.csv"))


sink(file.path(dirs$clustering, "Integration_Clustering_Report.txt"))
cat("===============================================================================\n")
cat("Integration and Clustering Analysis Report\n")
cat("CS-NO Hydrogel scRNA-seq Analysis\n")
cat("===============================================================================\n\n")

cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("1. INTEGRATION PARAMETERS:\n")
cat("--------------------------\n")
cat("   Method: Harmony\n")
cat("   PCs used: 30\n")
cat("   Batch variable: sample\n\n")

cat("2. CLUSTERING PARAMETERS:\n")
cat("--------------------------\n")
cat("   Resolution: 0.2\n")
cat("   Final cell types: 12\n\n")

cat("3. CELL TYPE ANNOTATIONS (12 cell types):\n")
cat("-----------------------------------------\n")
cat("   PTC-1: Healthy proximal tubule cells (clusters 0, 1)\n")
cat("   PTC-2: Proximal tubule cells (clusters 2, 14)\n")
cat("   PTC-3: Injured/failed-repair proximal tubule cells (cluster 15)\n")
cat("   Macrophage: Kidney macrophages (clusters 3, 12)\n")
cat("   NK cell: Natural killer cells (cluster 4)\n")
cat("   Medullary cell: Medullary cells (cluster 5)\n")
cat("   DCTC: Distal convoluted tubule cells (cluster 6)\n")
cat("   VEC: Vascular endothelial cells (cluster 7)\n")
cat("   Neutrophil: Neutrophils (cluster 8)\n")
cat("   CDPC: Collecting duct principal cells (cluster 9)\n")
cat("   Fibroblast: Fibroblasts (cluster 10)\n")
cat("   B cell: B lymphocytes (cluster 11)\n\n")

cat("4. CELL TYPE COUNTS:\n")
cat("--------------------\n")
print(table(seurat_annotated$cell_type))

cat("\n5. CELL COUNTS BY CONDITION:\n")
cat("-----------------------------\n")
print(table(seurat_annotated$cell_type, seurat_annotated$sample))

cat("\n6. OUTPUT FILES:\n")
cat("-----------------\n")
cat("   Main Figures:\n")
cat("   - Fig01_umap_celltypes.pdf/png\n")
cat("   - Fig02_umap_split_condition.pdf/png\n")
cat("   - Fig03_celltype_proportions.pdf/png\n")
cat("   - Fig04_feature_plots.pdf/png\n")
cat("   - Fig05_marker_heatmap.pdf/png\n\n")
cat("   Supplementary Figures:\n")
cat("   - S01-S10 (see clustering directory)\n\n")
cat("   Data Files:\n")
cat("   - seurat_integrated_annotated.rds\n")
cat("   - all_celltype_markers.csv\n")
cat("   - top10_markers_per_celltype.csv\n")
cat("   - celltype_proportions.csv\n")

cat("\n===============================================================================\n")
cat("Analysis Complete\n")
cat("Next step: Run 03_differential_expression.R\n")
cat("===============================================================================\n")
sink()

log_message("Integration and Clustering Analysis completed successfully!", log_file)


message("\n===============================================================================")
message("INTEGRATION AND CLUSTERING COMPLETE")
message("===============================================================================")
message("\nCell Type Summary (12 cell types):")
print(table(seurat_annotated$cell_type))
message("\nBy Condition:")
print(table(seurat_annotated$cell_type, seurat_annotated$sample))
message("\nOutput saved to:")
message(sprintf("  %s", dirs$clustering))
message("\nNext step: Run 03_differential_expression.R")
message("===============================================================================")

sessionInfo()
