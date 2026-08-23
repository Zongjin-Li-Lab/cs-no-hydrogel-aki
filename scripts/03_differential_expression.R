#!/usr/bin/env Rscript


rm(list = ls())
gc()


source("scripts/_common.R")
base_dir <- project_root()
setwd(base_dir)
cat("Working directory:", getwd(), "\n")


dirs <- make_output_dirs(base_dir)


condition_colors <- c(
  "Ctrl" = "#6B9CD3",
  "NO"   = "#8B5AA3"
)

regulation_colors <- c(
  "Up"   = "#B2182B",
  "Down" = "#2E5A87",
  "NS"   = "grey70"
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

heatmap_colors <- colorRampPalette(c("#2E5A87", "white", "#B2182B"))(100)
feature_colors <- c("lightgrey", "#2E5A87", "#B2182B")


theme_publication <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size + 2, hjust = 0.5, color = "grey40"),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.ticks = element_line(linewidth = 0.6, color = "black"),
      legend.title = element_text(size = base_size + 1, face = "bold"),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(0.6, "cm"),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 2, face = "bold"),
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
  library(DESeq2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(msigdbr)
  library(fgsea)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(future)
})

n_workers <- min(4L, future::availableCores())
if (n_workers > 1L) future::plan("multisession", workers = n_workers) else future::plan("sequential")
options(future.globals.maxSize = 8000 * 1024^2)

log_file <- file.path(dirs$de, "differential_expression_log.txt")
log_message("Starting Differential Expression Analysis", log_file)


log_message("Loading integrated data...", log_file)

seurat_obj <- readRDS(file.path(dirs$processed, "seurat_integrated_annotated.rds"))
message(sprintf("Loaded: %d cells, %d genes", ncol(seurat_obj), nrow(seurat_obj)))
message(sprintf("Cell types: %s", paste(unique(seurat_obj$cell_type), collapse = ", ")))

Idents(seurat_obj) <- "cell_type"


log_message("Running DE analysis for each cell type...", log_file)

run_de_by_celltype <- function(seurat_obj, celltype) {
  cells <- WhichCells(seurat_obj, idents = celltype)

  if (length(cells) < 20) {
    message(sprintf("  %s: Too few cells (%d), skipping", celltype, length(cells)))
    return(NULL)
  }

  subset_obj <- subset(seurat_obj, cells = cells)

  if (length(unique(subset_obj$sample)) < 2) {
    message(sprintf("  %s: Only one condition, skipping", celltype))
    return(NULL)
  }


  n_ctrl <- sum(subset_obj$sample == "Ctrl")
  n_no <- sum(subset_obj$sample == "NO")

  if (n_ctrl < 10 || n_no < 10) {
    message(sprintf("  %s: Insufficient cells (Ctrl=%d, NO=%d), skipping", celltype, n_ctrl, n_no))
    return(NULL)
  }

  tryCatch({
    Idents(subset_obj) <- "sample"

    markers <- FindMarkers(subset_obj,
                           ident.1 = "NO",
                           ident.2 = "Ctrl",
                           test.use = "wilcox",
                           min.pct = 0.1,
                           logfc.threshold = 0.1,
                           verbose = FALSE)

    markers$gene <- rownames(markers)
    markers$cell_type <- celltype

    markers <- markers %>%
      mutate(
        significance = case_when(
          p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Up",
          p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down",
          TRUE ~ "NS"
        )
      )

    message(sprintf("  %s: %d genes (Up=%d, Down=%d)",
                    celltype, nrow(markers),
                    sum(markers$significance == "Up"),
                    sum(markers$significance == "Down")))

    return(markers)

  }, error = function(e) {
    message(sprintf("  %s: Error - %s", celltype, e$message))
    return(NULL)
  })
}


celltypes <- unique(seurat_obj$cell_type)
de_results_list <- lapply(celltypes, function(ct) run_de_by_celltype(seurat_obj, ct))
names(de_results_list) <- celltypes


de_combined <- bind_rows(de_results_list[!sapply(de_results_list, is.null)])


write.csv(de_combined, file.path(dirs$de, "FindMarkers_NO_vs_Ctrl_all_celltypes.csv"),
          row.names = FALSE)

log_message(sprintf("DE complete: %d genes across %d cell types",
                    nrow(de_combined), length(unique(de_combined$cell_type))), log_file)


log_message("Creating volcano plots...", log_file)

create_volcano <- function(de_data, title, top_n = 15) {

  de_data <- de_data %>%
    mutate(
      neg_log10_pval = -log10(p_val_adj + 1e-300),
      label = ifelse(significance != "NS" &
                       (rank(-abs(avg_log2FC)) <= top_n | rank(p_val_adj) <= top_n),
                     gene, NA)
    )

  p <- ggplot(de_data, aes(x = avg_log2FC, y = neg_log10_pval, color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(aes(label = label), size = 4, max.overlaps = 20,
                    box.padding = 0.5, segment.color = "grey50") +
    scale_color_manual(values = regulation_colors) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.6) +
    labs(title = title,
         x = "Log2 Fold Change (NO vs Ctrl)",
         y = "-Log10 Adjusted P-value",
         color = "Regulation") +
    theme_publication(base_size = 14) +
    theme(legend.position = "right")

  return(p)
}


p_volcano_all <- create_volcano(de_combined, "Differential Expression: NO vs Ctrl (All Cell Types)")

ggsave(file.path(dirs$de, "Volcano_all_celltypes.pdf"),
       p_volcano_all, width = 12, height = 10, device = cairo_pdf)
ggsave(file.path(dirs$de, "Volcano_all_celltypes.png"),
       p_volcano_all, width = 12, height = 10, dpi = 300)


key_celltypes <- c("PTC-1", "PTC-2", "PTC-3", "Macrophage")

for (ct in key_celltypes) {
  if (ct %in% unique(de_combined$cell_type)) {
    ct_data <- de_combined %>% filter(cell_type == ct)
    p <- create_volcano(ct_data, sprintf("DE: %s (NO vs Ctrl)", ct))

    ggsave(file.path(dirs$de, sprintf("Volcano_%s.pdf", gsub("-", "", ct))),
           p, width = 10, height = 8, device = cairo_pdf)
    ggsave(file.path(dirs$de, sprintf("Volcano_%s.png", gsub("-", "", ct))),
           p, width = 10, height = 8, dpi = 300)
  }
}

log_message("Volcano plots saved", log_file)


log_message("Running GSEA analysis...", log_file)


ranks <- de_combined %>%
  filter(!is.na(avg_log2FC) & !is.na(p_val)) %>%
  group_by(gene) %>%
  summarise(avg_log2FC = mean(avg_log2FC), .groups = "drop") %>%
  arrange(desc(avg_log2FC))

gene_list <- setNames(ranks$avg_log2FC, ranks$gene)


hallmark <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol)

reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)


run_fgsea_analysis <- function(gene_list, pathways_df, name) {
  pathways <- split(pathways_df$gene_symbol, pathways_df$gs_name)

  fgsea_res <- fgsea(pathways = pathways,
                     stats = gene_list,
                     minSize = 15,
                     maxSize = 500,
                     nPermSimple = 10000)

  fgsea_res <- fgsea_res %>%
    arrange(pval) %>%
    mutate(database = name)

  return(fgsea_res)
}

gsea_hallmark <- run_fgsea_analysis(gene_list, hallmark, "Hallmark")
gsea_kegg <- run_fgsea_analysis(gene_list, kegg, "KEGG")
gsea_reactome <- run_fgsea_analysis(gene_list, reactome, "Reactome")


gsea_combined <- bind_rows(gsea_hallmark, gsea_kegg, gsea_reactome)
write.csv(gsea_combined %>% dplyr::select(-leadingEdge),
          file.path(dirs$de, "GSEA_results_combined.csv"), row.names = FALSE)

log_message(sprintf("GSEA complete: %d significant pathways",
                    sum(gsea_combined$padj < 0.05)), log_file)


log_message("Creating GSEA visualizations...", log_file)


plot_gsea_bar <- function(gsea_res, title, top_n = 20) {

  top_pathways <- gsea_res %>%
    filter(padj < 0.05) %>%
    arrange(NES) %>%
    head(top_n)

  if (nrow(top_pathways) == 0) {
    message(sprintf("No significant pathways for %s", title))
    return(NULL)
  }

  top_pathways$pathway <- gsub("^HALLMARK_|^KEGG_|^REACTOME_", "", top_pathways$pathway)
  top_pathways$pathway <- gsub("_", " ", top_pathways$pathway)
  top_pathways$pathway <- factor(top_pathways$pathway, levels = top_pathways$pathway)

  p <- ggplot(top_pathways, aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2E5A87"),
                      labels = c("TRUE" = "Up in NO", "FALSE" = "Down in NO")) +
    labs(title = title,
         x = "Normalized Enrichment Score",
         y = NULL,
         fill = "Direction") +
    theme_publication(base_size = 12) +
    theme(axis.text.y = element_text(size = 10))

  return(p)
}

p_hallmark <- plot_gsea_bar(gsea_hallmark, "GSEA: Hallmark Pathways")
p_kegg <- plot_gsea_bar(gsea_kegg, "GSEA: KEGG Pathways")

if (!is.null(p_hallmark)) {
  ggsave(file.path(dirs$de, "GSEA_Hallmark_bar.pdf"),
         p_hallmark, width = 12, height = 8, device = cairo_pdf)
}

if (!is.null(p_kegg)) {
  ggsave(file.path(dirs$de, "GSEA_KEGG_bar.pdf"),
         p_kegg, width = 12, height = 10, device = cairo_pdf)
}


log_message("Creating key genes heatmap...", log_file)


key_genes <- list(
  Injury = c("Havcr1", "Lcn2", "Sox9", "Krt8", "Krt18", "Cd44"),
  MAPK = c("Mapk1", "Mapk3", "Mapk14", "Jun", "Fos", "Dusp1"),
  Inflammation = c("Il1b", "Tnf", "Ccl2", "Cxcl1", "Cxcl10", "Il6"),
  M1_markers = c("Cd86", "Nos2", "Il1b", "Tnf"),
  M2_markers = c("Mrc1", "Arg1", "Cd163", "Retnla"),
  Angiogenesis = c("Vegfa", "Vegfb", "Flt1", "Kdr", "Angpt1", "Angpt2"),
  Repair = c("Pcna", "Mki67", "Ccnd1", "Sox9")
)


all_key_genes <- unique(unlist(key_genes))
genes_present <- all_key_genes[all_key_genes %in% rownames(seurat_obj)]

if (length(genes_present) >= 10) {

  avg_expr <- AverageExpression(seurat_obj,
                                features = genes_present,
                                group.by = c("cell_type", "sample"))$RNA


  avg_expr_scaled <- t(scale(t(as.matrix(avg_expr))))


  col_fun <- colorRamp2(c(-2, 0, 2), c("#2E5A87", "white", "#B2182B"))

  ht <- Heatmap(avg_expr_scaled,
                name = "Z-score",
                col = col_fun,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize = 10),
                column_names_gp = gpar(fontsize = 10),
                column_title = "Key Genes Expression",
                column_title_gp = gpar(fontsize = 14, fontface = "bold"))

  pdf(file.path(dirs$de, "Key_genes_heatmap.pdf"), width = 14, height = 12)
  draw(ht)
  dev.off()

  log_message("Key genes heatmap saved", log_file)
}


log_message("Creating DE summary...", log_file)

de_summary <- de_combined %>%
  group_by(cell_type, significance) %>%
  summarise(n_genes = n(), .groups = "drop") %>%
  pivot_wider(names_from = significance, values_from = n_genes, values_fill = 0)

write.csv(de_summary, file.path(dirs$de, "DE_summary_by_celltype.csv"), row.names = FALSE)


log_message("Saving results...", log_file)

de_results <- list(
  cell_level = de_combined,
  gsea_hallmark = gsea_hallmark,
  gsea_kegg = gsea_kegg,
  gsea_reactome = gsea_reactome,
  summary = de_summary
)

saveRDS(de_results, file.path(dirs$processed, "differential_expression_results.rds"))


sink(file.path(dirs$de, "Differential_Expression_Report.txt"))
cat("===============================================================================\n")
cat("Differential Expression Analysis Report\n")
cat("CS-NO Hydrogel scRNA-seq Analysis\n")
cat("===============================================================================\n\n")

cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Random seed: 20260614\n\n")

cat("1. DE SUMMARY (NO vs Ctrl):\n")
cat("---------------------------\n")
cat(sprintf("   Total DE genes: %d\n", nrow(de_combined)))
cat(sprintf("   Upregulated in NO: %d\n", sum(de_combined$significance == "Up")))
cat(sprintf("   Downregulated in NO: %d\n", sum(de_combined$significance == "Down")))

cat("\n2. DE BY CELL TYPE:\n")
cat("-------------------\n")
print(de_summary)

cat("\n3. GSEA SUMMARY:\n")
cat("----------------\n")
cat(sprintf("   Hallmark significant: %d\n", sum(gsea_hallmark$padj < 0.05)))
cat(sprintf("   KEGG significant: %d\n", sum(gsea_kegg$padj < 0.05)))
cat(sprintf("   Reactome significant: %d\n", sum(gsea_reactome$padj < 0.05)))

cat("\n4. OUTPUT FILES:\n")
cat("----------------\n")
cat("   - FindMarkers_NO_vs_Ctrl_all_celltypes.csv\n")
cat("   - GSEA_results_combined.csv\n")
cat("   - Volcano plots (PDF/PNG)\n")
cat("   - GSEA bar plots (PDF)\n")
cat("   - Key_genes_heatmap.pdf\n")

cat("\n===============================================================================\n")
cat("Analysis Complete\n")
cat("Next step: Run 04_trajectory_analysis.R\n")
cat("===============================================================================\n")
sink()

log_message("Differential Expression Analysis completed!", log_file)

sessionInfo()
