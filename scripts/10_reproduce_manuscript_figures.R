#!/usr/bin/env Rscript


rm(list = ls()); gc()


source("scripts/_common.R")
base_dir <- project_root()
PROC_DIR <- file.path(base_dir, "data", "processed")
OUT_DIR  <- file.path(base_dir, "figures", "manuscript")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


condition_colors <- c("Control" = "#6B9CD3", "NO_treated" = "#8B5AA3")
celltype_colors <- c(
  "PTC-1" = "#2E5A87", "PTC-2" = "#4A7BB0", "PTC-3" = "#E74C3C",
  "Macrophage" = "#8B5AA3", "NK cell" = "#1ABC9C", "Medullary cell" = "#7FBAC4",
  "DCTC" = "#3498DB", "VEC" = "#2ECC71", "Neutrophil" = "#F39C12",
  "CDPC" = "#5DADE2", "Fibroblast" = "#E67E22", "B cell" = "#16A085"
)
ptc_colors <- c("PTC-1" = "#2E5A87", "PTC-2" = "#4A7BB0", "PTC-3" = "#E74C3C")
regulation_colors <- c("Up" = "#B2182B", "Down" = "#2E5A87", "NS" = "grey70")
feature_colors <- c("lightgrey", "#2E5A87", "#B2182B")


theme_pub <- function(base_size = 18) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title    = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      axis.title    = element_text(size = base_size, face = "bold"),
      axis.text     = element_text(size = base_size - 2, color = "black"),
      axis.line     = element_line(linewidth = 0.8, color = "black"),
      axis.ticks    = element_line(linewidth = 0.6, color = "black"),
      legend.title  = element_text(size = base_size - 2, face = "bold"),
      legend.text   = element_text(size = base_size - 3),
      strip.text    = element_text(size = base_size, face = "bold"),
      plot.margin   = margin(10, 15, 10, 15),
      panel.grid    = element_blank()
    )
}


suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(patchwork)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
})


CELLTYPE_COL  <- "cell_type"
CONDITION_COL <- "condition"


cat("Loading Seurat object...\n")
seu <- readRDS(file.path(PROC_DIR, "seurat_with_kirita_scores.rds"))
de_results <- readRDS(file.path(PROC_DIR, "differential_expression_results.rds"))


cat("Conditions:", unique(seu@meta.data[[CONDITION_COL]]), "\n")


Idents(seu) <- CELLTYPE_COL


library(dplyr)
ct <- table(seu$cell_type, seu$sample)
results <- data.frame()
for (i in 1:nrow(ct)) {
  mat <- matrix(c(ct[i,"Ctrl"], ct[i,"NO"],
                  sum(ct[,"Ctrl"]) - ct[i,"Ctrl"],
                  sum(ct[,"NO"]) - ct[i,"NO"]), nrow=2)
  ft <- fisher.test(mat)
  results <- rbind(results, data.frame(
    cell_type = rownames(ct)[i],
    Ctrl = ct[i,"Ctrl"], NO = ct[i,"NO"],
    pvalue = ft$p.value))
}
results$FDR <- p.adjust(results$pvalue, method = "BH")
print(results)


cat("\n========== FIGURE 4 ==========\n")


pdf(file.path(OUT_DIR, "Fig4A_UMAP_split.pdf"), width = 16, height = 7)
print(
  DimPlot(seu, group.by = CELLTYPE_COL, split.by = CONDITION_COL,reduction = "umap",
          cols = celltype_colors, label = FALSE, repel = TRUE,
          label.size = 10, pt.size = 0.4) +
    theme_pub(20) +
    theme(legend.text = element_text(size = 20),
          legend.key.size = unit(0.6, "cm")) +
    ggtitle(NULL)
)
dev.off()


comp_df <- seu@meta.data %>%
  group_by(!!sym(CONDITION_COL), !!sym(CELLTYPE_COL)) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(!!sym(CONDITION_COL)) %>%
  mutate(pct = n / sum(n) * 100)

pdf(file.path(OUT_DIR, "Fig4B_composition_bar.pdf"), width = 8, height = 6)
print(
  ggplot(comp_df, aes(x = !!sym(CONDITION_COL), y = pct, fill = !!sym(CELLTYPE_COL))) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = celltype_colors) +
    labs(x = NULL, y = "Proportion (%)", fill = "Cell Type") +
    theme_pub(18) +
    theme(axis.text.x = element_text(size = 20),
          legend.text = element_text(size = 20))
)
dev.off()


marker_genes <- c("Aqp1", "Ccn1", "Slc34a1", "Inmt", "Akr1c21", "Ttc36",
                   "Akr1b8", "C1qa", "C1qb", "Flt1", "Pgam2",
                   "Csf1r", "Adgre1", "Pecam1", "Cdh5",
                   "Plvap", "Mmp9", "Col3a1", "Dcn",
                   "S100a9", "S100a8", "Ppp1r1b", "Nkg7",
                   "Cd79a", "Borkl")
marker_genes <- marker_genes[marker_genes %in% rownames(seu)]

pdf(file.path(OUT_DIR, "Fig4C_marker_violins.pdf"), width = 18, height = 14)

p <- VlnPlot(seu, features = marker_genes,
             pt.size = 0, stack = TRUE, flip = TRUE, fill.by = "ident") +
  scale_fill_manual(values = celltype_colors)

p <- p & theme(axis.title.x = element_text(size = 20, face = "bold"),
               axis.title.y = element_text(size = 20, face = "bold"),
               axis.text.x = element_text(size = 22, angle = 45, hjust = 1, color = "black"),
               axis.text.y = element_text(size = 20, color = "black"),
               strip.text.y.right = element_text(size = 20, face = "italic", angle = 0),
               legend.position = "none")

print(p)
dev.off()


ct <- table(seu@meta.data[[CELLTYPE_COL]], seu@meta.data[[CONDITION_COL]])
totals <- colSums(ct)
CTRL <- colnames(ct)[grep("Ctrl|Control", colnames(ct))]
NO   <- colnames(ct)[grep("NO", colnames(ct))]

results <- do.call(rbind, lapply(rownames(ct), function(cell_type) {
  a <- ct[cell_type, CTRL]; b <- ct[cell_type, NO]
  mat <- matrix(c(a, totals[CTRL] - a, b, totals[NO] - b), nrow = 2)
  ft <- fisher.test(mat)
  data.frame(cell_type = cell_type, fold_change = (b/totals[NO])/(a/totals[CTRL]),
             p_value = ft$p.value)
}))
results$FDR <- p.adjust(results$p_value, method = "BH")
results$sig <- ifelse(results$FDR < 0.001, "***",
               ifelse(results$FDR < 0.01, "**",
               ifelse(results$FDR < 0.05, "*", "ns")))
results$fill <- ifelse(results$cell_type %in% names(ptc_colors),
                       ptc_colors[results$cell_type],
                ifelse(results$FDR < 0.05, "#8B5AA3", "grey60"))
results <- results %>% arrange(desc(fold_change))
results$cell_type <- factor(results$cell_type, levels = results$cell_type)

pdf(file.path(OUT_DIR, "Fig4D_log2FC_composition.pdf"), width = 8, height = 7)
print(
  ggplot(results, aes(x = cell_type, y = log2(fold_change))) +
    geom_col(fill = results$fill, width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_text(aes(label = sig),
              vjust = ifelse(results$fold_change > 1, -0.5, 1.5),
              size = 6, fontface = "bold") +
    labs(x = NULL, y = expression(bold(log[2]~"(Fold Change: NO / Control)"))) +
    theme_pub(18) +
    theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1))
)
dev.off()


ptc <- subset(seu, idents = c("PTC-1", "PTC-2", "PTC-3"))
score_col <- grep("Kirita_FR_PTC", colnames(ptc@meta.data), value = TRUE)[1]


pdf(file.path(OUT_DIR, "Fig4E_Kirita_violin.pdf"), width = 7, height = 5)
print(
  VlnPlot(ptc, features = score_col, group.by = CELLTYPE_COL,
          split.by = CONDITION_COL, cols = condition_colors, pt.size = 0) +
    labs(x = NULL, y = "Kirita Failed-Repair\nPTC Score") +
    theme_pub(18) +
    theme(legend.position = "right",
          axis.text.x = element_text(size = 18),
          axis.title.y = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16))
)
dev.off()


pdf(file.path(OUT_DIR, "Fig4F_Mrc1_UMAP.pdf"), width = 12, height = 5.5)
print(
  FeaturePlot(seu, features = "Mrc1", split.by = CONDITION_COL,
              reduction = "umap", cols = c("lightgrey", "#B2182B"),
              pt.size = 0.3, order = TRUE) &
    labs(x = "umap_1", y = "umap_2") &
    theme_pub(18) &
    theme(plot.title = element_text(size = 22, face = "bold.italic"),
          axis.title = element_text(size = 18, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.key.height = unit(1.5, "cm"))
)
dev.off()


pdf(file.path(OUT_DIR, "Fig4G_Cd86_UMAP.pdf"), width = 12, height = 5.5)
print(
  FeaturePlot(seu, features = "Cd86", split.by = CONDITION_COL,
              reduction = "umap", cols = c("lightgrey", "#B2182B"),
              pt.size = 0.3, order = TRUE) &
    labs(x = "umap_1", y = "umap_2") &
    theme_pub(18) &
    theme(plot.title = element_text(size = 22, face = "bold.italic"),
          axis.title = element_text(size = 18, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.key.height = unit(1.5, "cm"))
)
dev.off()


library(CellChat)
library(grid)

pdf(file.path(OUT_DIR, "Fig5C_Ctrl_CCL.pdf"), width = 10, height = 9)
par(mar = c(1, 1, 4, 1), cex = 1.5)
netVisual_aggregate(cellchat_ctrl, signaling = "CCL",
                    layout = "chord", lab.cex = 2)
grid.text("Ctrl - CCL", x = 0.5, y = 0.97,
          gp = gpar(fontsize = 24, fontface = "bold"))
dev.off()

pdf(file.path(OUT_DIR, "Fig5C_NO_CCL.pdf"), width = 10, height = 9)
par(mar = c(1, 1, 4, 1), cex = 1.5)
netVisual_aggregate(cellchat_no, signaling = "CCL",
                    layout = "chord", lab.cex = 2)
grid.text("NO - CCL", x = 0.5, y = 0.97,
          gp = gpar(fontsize = 24, fontface = "bold"))
dev.off()


cat("\n========== FIGURE 5 ==========\n")


str(de_results, max.level = 1)


de_df <- bind_rows(de_results, .id = "cluster")


volcano_data <- de_df %>%
  filter(!is.na(avg_log2FC) & !is.na(p_val_adj)) %>%
  mutate(
    neg_log10 = -log10(p_val_adj + 1e-300),
    sig = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down",
      TRUE ~ "NS"
    )
  )

highlight_genes <- c("Cxcl10", "Ccl5", "Ccl2", "S100a9", "Cxcl1",
                      "Ccl3", "Ccl4", "Cxcl2")
top_sig <- volcano_data %>% filter(sig != "NS") %>% arrange(p_val_adj) %>% head(15)
label_genes <- unique(c(highlight_genes, top_sig$gene))
volcano_data$label <- ifelse(volcano_data$gene %in% label_genes, volcano_data$gene, NA)

pdf(file.path(OUT_DIR, "Fig5D_volcano.pdf"), width = 10, height = 8)
print(
  ggplot(volcano_data, aes(x = avg_log2FC, y = neg_log10, color = sig)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_text_repel(aes(label = label), size = 5.5, max.overlaps = 30,
                    box.padding = 0.6, segment.color = "grey40",
                    fontface = "italic", show.legend = FALSE) +
    scale_color_manual(values = regulation_colors, name = "Regulation") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    labs(x = expression(bold(log[2]~"Fold Change")),
         y = expression(bold(-log[10]~"Adjusted P-value")),
         title = "Inflammatory Cytokine,\nNO vs Control in all PTC clusters") +
    theme_pub(18) +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = "grey80"))
)
dev.off()


for (gene in c("Cxcl10", "Ccl5")) {
  pdf(file.path(OUT_DIR, paste0("Fig5E_", gene, "_UMAP.pdf")), width = 12, height = 5)
  print(
    FeaturePlot(seu, features = gene, split.by = CONDITION_COL,reduction = "umap",
                cols = c("lightgrey", "#B2182B"), pt.size = 0.3, order = TRUE) &
      theme_pub(18) &
      theme(plot.title = element_text(size = 22, face = "bold.italic"),
            axis.title = element_text(size = 18, face = "bold"),
            legend.title = element_blank(),
            legend.text = element_text(size = 14),
            legend.key.height = unit(1.5, "cm"))
  )
  dev.off()
}


sig_genes_up <- de_df %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5) %>%
  pull(gene) %>% unique()

if (length(sig_genes_up) > 10) {
  gene_entrez <- bitr(sig_genes_up, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)
  if (nrow(gene_entrez) > 10) {
    kegg_res <- enrichKEGG(gene = gene_entrez$ENTREZID, organism = "mmu",
                            pvalueCutoff = 0.05)
    if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
      pdf(file.path(OUT_DIR, "Fig5F_KEGG.pdf"), width = 10, height = 8)
      print(
        dotplot(kegg_res, showCategory = 12) +
          labs(title = "PTC Up-expressed Genes Pathway Enrichment,\nCtrl vs NO in all PTC clusters") +
          theme_pub(16) +
          theme(axis.text.y = element_text(size = 14))
      )
      dev.off()
    }
  }
}


mapk_genes <- c("Fos", "Jun", "Junb", "Egr1", "Dusp1", "Dusp6")
mapk_present <- mapk_genes[mapk_genes %in% rownames(seu)]

pdf(file.path(OUT_DIR, "Fig5G_MAPK_dotplot.pdf"), width = 10, height = 8)
print(
  DotPlot(ptc, features = mapk_present, group.by = CELLTYPE_COL,
          split.by = CONDITION_COL, cols = c("#6B9CD3", "#8B5AA3")) +
    coord_flip() +
    labs(title = "MAPK Pathway Genes", x = NULL, y = NULL) +
    theme_pub(18) +
    theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 16, face = "italic"),
          plot.margin = margin(10, 15, 10, 40))
)
dev.off()


cat("\n========== FIGURE 6 ==========\n")


library(monocle3)
library(ggplot2)
library(pheatmap)


cds <- readRDS(file.path(PROC_DIR, "monocle3_cds.rds"))


pdf(file.path(OUT_DIR, "Fig6A_trajectory.pdf"), width = 16, height = 7)
print(
  plot_cells(cds,
             color_cells_by = "cell_type",
             label_cell_groups = FALSE,
             label_leaves = FALSE,
             label_branch_points = FALSE,
             label_roots = FALSE,
             graph_label_size = 0,
             cell_size = 0.4,
             trajectory_graph_color = "black",
             trajectory_graph_segment_size = 1.2) +
    scale_color_manual(values = celltype_colors) +
    facet_wrap(~sample, ncol = 2) +
    labs(x = "umap_1", y = "umap_2", color = "Cell Type") +
    theme_pub(18) +
    theme(strip.text = element_text(size = 22, face = "bold"),
          legend.text = element_text(size = 18),
          legend.key.size = unit(0.8, "cm"))
)
dev.off()


pt_data <- data.frame(
  pseudotime = pseudotime(cds),
  condition = colData(cds)$orig.ident,
  cell_type = colData(cds)$cell_type
)

genes_to_plot <- c("Aqp1", "Angpt1", "Col1a1", "Cxcl10",
                   "Ghr", "Havcr1", "Lcn2", "Lrp2",
                   "Slc34a1", "Ttc36", "Vim", "Vegfa")
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(cds)]


expr_mat <- as.matrix(normalized_counts(cds)[genes_to_plot, , drop = FALSE])
pt_data <- cbind(pt_data, t(expr_mat))
pt_data <- pt_data[!is.infinite(pt_data$pseudotime), ]

library(tidyr)
pt_long <- pt_data %>%
  pivot_longer(cols = all_of(genes_to_plot), names_to = "gene", values_to = "expression")

pdf(file.path(OUT_DIR, "Fig6B_pseudotime_curves.pdf"), width = 16, height = 10)
print(
  ggplot(pt_long, aes(x = pseudotime, y = expression, color = condition)) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, alpha = 0.2) +
    scale_color_manual(values = c("Ctrl" = "#6B9CD3", "NO" = "#8B5AA3"),
                       labels = c("Ctrl", "NO")) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) +
    labs(x = "Pseudotime", y = "Expression", color = "Condition",
         title = "Gene Expression Along Pseudotime") +
    theme_pub(16) +
    theme(strip.text = element_text(size = 18, face = "bold.italic"),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 16),
          legend.text = element_text(size = 20),
          legend.position = "bottom",
          legend.direction = "horizontal")
)
dev.off()


library(ggplot2)
library(dplyr)
library(tidyr)


pathway_sets <- list(
  "Inflammation/Chemokine" = c("Cxcl10", "Ccl5", "Ccl2", "Tnf", "Il1b", "Cxcl1"),
  "MAPK Signaling" = c("Fos", "Jun", "Junb", "Egr1", "Dusp1", "Dusp6"),
  "Angiogenesis" = c("Vegfa", "Vegfc", "Pecam1", "Kdr", "Flt1", "Esm1"),
  "Fibrosis/EMT" = c("Col1a1", "Fn1", "Vim", "Acta2", "Tgfb1", "Col3a1"),
  "Healthy PTC" = c("Slc34a1", "Lrp2", "Aqp1", "Cubn"),
  "Apoptosis" = c("Bax", "Bad", "Casp3", "Bcl2")
)


for (pw in names(pathway_sets)) {
  genes <- pathway_sets[[pw]]
  genes <- genes[genes %in% rownames(cds)]
  if (length(genes) >= 3) {
    expr <- as.matrix(normalized_counts(cds)[genes, ])
    score <- colMeans(expr)
    colData(cds)[[paste0("pw_", gsub("[/ ]", "_", pw))]] <- score
  }
}


pt_pw <- data.frame(
  pseudotime = pseudotime(cds),
  condition = colData(cds)$orig.ident,
  colData(cds)[, grep("^pw_", colnames(colData(cds)))]
)
pt_pw <- pt_pw[!is.infinite(pt_pw$pseudotime), ]


pw_cols <- grep("^pw_", colnames(pt_pw), value = TRUE)
pt_long <- pt_pw %>%
  pivot_longer(cols = all_of(pw_cols), names_to = "pathway", values_to = "score") %>%
  mutate(pathway = gsub("^pw_", "", pathway),
         pathway = gsub("_", " ", pathway))


pt_long$pathway <- factor(pt_long$pathway, levels = c(
  "Healthy PTC", "Angiogenesis", "MAPK Signaling",
  "Inflammation Chemokine", "Fibrosis EMT", "Apoptosis"
))


pdf(file.path(OUT_DIR, "Fig6CD_pathway_pseudotime.pdf"), width = 14, height = 10)
print(
  ggplot(pt_long, aes(x = pseudotime, y = score, color = condition)) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.5, alpha = 0.2) +
    scale_color_manual(values = c("Ctrl" = "#6B9CD3", "NO" = "#8B5AA3"),
                       labels = c("Control", "NO")) +
    facet_wrap(~pathway, scales = "free_y", ncol = 3) +
    labs(x = "Pseudotime", y = "Pathway Score",
         color = "Condition") +
    theme_pub(18) +
    theme(strip.text = element_text(size = 16, face = "bold"),
          strip.background = element_rect(fill = "grey95", color = NA),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 14),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 16))
)
dev.off()


pt_df <- data.frame(
  pseudotime = pseudotime(cds),
  condition = colData(cds)$orig.ident,
  cell_type = colData(cds)$cell_type
)
pt_df <- pt_df[!is.infinite(pt_df$pseudotime), ]


pt_df$stage <- cut(pt_df$pseudotime, breaks = 5,
                   labels = c("Early", "Early-Mid", "Mid", "Mid-Late", "Late"))


prop_df <- pt_df %>%
  group_by(condition, stage, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition, stage) %>%
  mutate(pct = n / sum(n) * 100)

pdf(file.path(OUT_DIR, "Fig6D_celltype_pseudotime_bar.pdf"), width = 12, height = 6)
print(
  ggplot(prop_df, aes(x = stage, y = pct, fill = cell_type)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = celltype_colors) +
    facet_wrap(~condition) +
    labs(x = "Pseudotime Stage", y = "Proportion (%)", fill = "Cell Type") +
    theme_pub(18) +
    theme(axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
          legend.text = element_text(size = 13),
          strip.text = element_text(size = 18, face = "bold"),
          legend.position = "right")
)
dev.off()


library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)


seurat_traj <- readRDS(file.path(PROC_DIR, "seurat_with_trajectory.rds"))


traj_types <- c("PTC-1", "PTC-2", "PTC-3", "DCTC", "CDPC", "Medullary cell")
traj_types <- traj_types[traj_types %in% unique(seurat_traj$cell_type)]
seurat_traj <- subset(seurat_traj, cell_type %in% traj_types)


if (exists("graph_test_res")) {
  top_genes <- graph_test_res %>%
    as.data.frame() %>%
    filter(q_value < 0.05 & morans_I > 0.1) %>%
    arrange(desc(morans_I)) %>%
    head(100) %>%
    rownames()
  cat("Using", length(top_genes), "graph_test genes\n")
} else {
  seurat_traj <- FindVariableFeatures(seurat_traj, nfeatures = 500)
  top_genes <- head(VariableFeatures(seurat_traj), 100)
  cat("Using", length(top_genes), "VariableFeature genes\n")
}


cluster_colors <- c("1" = "#E67E22", "2" = "#5DADE2", "3" = "#8B5AA3", "4" = "#1ABC9C")
bg_colors     <- c("1" = "#FFF3E0", "2" = "#E8F5E9", "3" = "#FCE4EC", "4" = "#E0F7FA")
text_colors   <- c("1" = "#E65100", "2" = "#2E7D32", "3" = "#AD1457", "4" = "#00695C")


graph_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
saveRDS(graph_test_res, file.path(OUT_DIR, "graph_test_results.rds"))
graph_test_res <- readRDS(file.path(OUT_DIR, "graph_test_results.rds"))


for (cond in c("Ctrl", "NO")) {
  cat(paste0("\n=== ", cond, " ===\n"))


  cds_sub <- cds[, colData(cds)$orig.ident == cond]


  cds_sub <- cluster_cells(cds_sub, reduction_method = "UMAP")
  cds_sub <- learn_graph(cds_sub, use_partition = FALSE)


  ptc1_cells <- colnames(cds_sub)[colData(cds_sub)$cell_type == "PTC-1"]
  cell_closest <- cds_sub@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  ptc1_vertices <- as.numeric(cell_closest[ptc1_cells, 1])
  root_node <- names(which.min(table(ptc1_vertices)))
  root_node <- igraph::V(principal_graph(cds_sub)[["UMAP"]])$name[as.numeric(root_node)]
  cds_sub <- order_cells(cds_sub, root_pr_nodes = root_node)


  gt_res <- graph_test(cds_sub, neighbor_graph = "principal_graph", cores = 4)
  saveRDS(gt_res, file.path(OUT_DIR, paste0("graph_test_", cond, ".rds")))

  top_genes_cond <- gt_res %>%
    as.data.frame() %>%
    filter(q_value < 0.05 & morans_I > 0.1) %>%
    arrange(desc(morans_I)) %>%
    head(100) %>%
    rownames()

  cat(paste0(cond, ": ", length(top_genes_cond), " genes\n"))

  create_pseudotime_heatmap(seurat_traj, cond, top_genes_cond,
                            file.path(OUT_DIR, paste0("Fig6", ifelse(cond=="Ctrl","C","D"),
                                                      "_heatmap_", cond, "_final.pdf")))
}


cat("Fig6A, 6B, 6C, 6D saved to:", OUT_DIR, "\n")


library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)


seu <- readRDS(file.path(PROC_DIR, "seurat_with_kirita_scores.rds"))


ptc2 <- subset(seu, cell_type == "PTC-2")


Idents(ptc2) <- "sample"
de_ptc2 <- FindMarkers(ptc2, ident.1 = "Ctrl", ident.2 = "NO",
                       logfc.threshold = 0.25, min.pct = 0.1)


up_in_ctrl <- de_ptc2 %>%
  filter(avg_log2FC > 0.25 & p_val_adj < 0.05) %>%
  rownames()


ids <- bitr(up_in_ctrl, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)


kegg_res <- enrichKEGG(gene = ids$ENTREZID,
                       organism = "mmu",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1)


print(kegg_res@result %>% filter(grepl("TGF|HIF|Hypoxia|MAPK", Description)))


pdf("Fig_S8F_NR-PTC_KEGG.pdf", width = 8, height = 6)
dotplot(kegg_res, showCategory = 15,
        title = "KEGG: Pathways enriched in Control NR-PTCs vs NO") +
  theme(axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))
dev.off()


cat("\n=== All significant KEGG terms ===\n")
print(kegg_res@result %>% filter(p.adjust < 0.05) %>% select(Description, p.adjust, Count))


ct <- table(seu$cell_type, seu$sample)
props <- prop.table(ct, margin = 2) * 100

cat("\n=== S9A: Cell type proportions ===\n")
print(round(props, 2))
cat("\nMacrophage Ctrl:", round(props["Macrophage","Ctrl"], 2), "%\n")
cat("Macrophage NO:", round(props["Macrophage","NO"], 2), "%\n")
cat("Direction: Macrophage proportion is",
    ifelse(props["Macrophage","NO"] > props["Macrophage","Ctrl"], "HIGHER", "LOWER"),
    "in NO\n")


results_or <- data.frame()
for (i in 1:nrow(ct)) {
  a <- ct[i, "NO"]
  b <- ct[i, "Ctrl"]
  c <- sum(ct[, "NO"]) - a
  d <- sum(ct[, "Ctrl"]) - b
  or <- (a/c) / (b/d)
  log2or <- log2(or)
  results_or <- rbind(results_or, data.frame(
    cell_type = rownames(ct)[i],
    NO = a, Ctrl = b,
    odds_ratio = round(or, 3),
    log2OR = round(log2or, 3)
  ))
}

cat("\n=== S9B: Log2 Odds Ratio (positive = enriched in NO) ===\n")
print(results_or)


results_or$cell_type <- factor(results_or$cell_type,
                               levels = results_or$cell_type[order(results_or$log2OR)])
results_or$significant <- ifelse(results_or$cell_type %in%
                                   c("PTC-2","DCTC","CDPC","Fibroblast","Macrophage","Neutrophil","B cell",
                                     "Medullary cell","PTC-1"), "Yes", "No")

pdf("Fig_S9B_forest_plot.pdf", width = 8, height = 5)
ggplot(results_or, aes(x = log2OR, y = cell_type, color = significant)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Yes" = "red", "No" = "grey50")) +
  labs(x = "log2(Odds Ratio: NO/Ctrl)", y = "",
       title = "Cell type enrichment in NO vs Control") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")
dev.off()


ptcs <- subset(seu, cell_type %in% c("PTC-1", "PTC-2", "PTC-3"))
Idents(ptcs) <- "sample"
de_all <- FindMarkers(ptcs, ident.1 = "NO", ident.2 = "Ctrl",
                      logfc.threshold = 0, min.pct = 0.1)


de_all$gene <- rownames(de_all)
gene_list <- de_all$avg_log2FC
names(gene_list) <- de_all$gene
gene_list <- sort(gene_list, decreasing = TRUE)

cat("\n=== S10 sign check ===\n")
cat("Positive log2FC = upregulated in NO\n")
cat("Negative log2FC = upregulated in Ctrl (downregulated in NO)\n")
cat("\nTop 10 UP in NO:\n")
print(head(gene_list, 10))
cat("\nTop 10 DOWN in NO (UP in Ctrl):\n")
print(tail(gene_list, 10))


cat("\n=== Key gene directions ===\n")
check_genes <- c("Vegfa", "Esm1", "Pecam1",
                 "Cxcl10", "Ccl5", "Fos", "Jun",
                 "Vim", "Fn1", "Col1a1",
                 "Casp3", "Bax")
for (g in check_genes) {
  if (g %in% names(gene_list)) {
    cat(sprintf("  %s: log2FC = %.3f (%s in NO)\n", g, gene_list[g],
                ifelse(gene_list[g] > 0, "UP", "DOWN")))
  }
}


ptc2 <- subset(seu, cell_type == "PTC-2")
Idents(ptc2) <- "sample"
de_ptc2 <- FindMarkers(ptc2, ident.1 = "NO", ident.2 = "Ctrl",
                       features = c("Fos","Jun","Junb","Egr1","Dusp1",
                                    "Vegfa","Esm1","Cxcl10","Ccl5","Bax"),
                       logfc.threshold = 0, min.pct = 0)
print(de_ptc2)


ptc3 <- subset(seu, cell_type == "PTC-3")
Idents(ptc3) <- "sample"
de_ptc3 <- FindMarkers(ptc3, ident.1 = "NO", ident.2 = "Ctrl",
                       features = c("Fos","Jun","Junb","Egr1","Dusp1",
                                    "Vegfa","Esm1","Cxcl10","Ccl5"),
                       logfc.threshold = 0, min.pct = 0)
print(de_ptc3)


cat("\n")
cat("═══════════════════════════════════════════════════\n")
cat("  ALL PANELS SAVED TO:", OUT_DIR, "\n")
cat("═══════════════════════════════════════════════════\n")
cat("  Fig4A_UMAP_split.pdf\n")
cat("  Fig4B_composition_bar.pdf\n")
cat("  Fig4C_marker_violins.pdf\n")
cat("  Fig4D_log2FC_composition.pdf\n")
cat("  Fig4E_Kirita_violin.pdf\n")
cat("  Fig4F_Mrc1_UMAP.pdf\n")
cat("  Fig4G_Cd86_UMAP.pdf\n")
cat("  Fig5D_volcano.pdf\n")
cat("  Fig5E_Cxcl10_UMAP.pdf\n")
cat("  Fig5E_Ccl5_UMAP.pdf\n")
cat("  Fig5F_KEGG.pdf\n")
cat("  Fig5G_MAPK_dotplot.pdf\n")
cat("  Fig6A_trajectory_UMAP.pdf\n")
cat("═══════════════════════════════════════════════════\n")
cat("\n  Key changes vs old plots:\n")
cat("  - base_size: 14 -> 18\n")
cat("  - gene labels: 4pt -> 5.5pt (volcano)\n")
cat("  - axis text: 12pt -> 16pt\n")
cat("  - legend text: 10pt -> 14pt\n")
cat("  - UMAP labels: 3pt -> 6pt\n")
cat("  - pathway names: 11pt -> 14pt\n")
cat("═══════════════════════════════════════════════════\n")
