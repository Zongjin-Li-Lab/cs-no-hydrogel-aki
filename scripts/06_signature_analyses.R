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

phase_colors <- c(
  "G1"  = "#2E5A87",
  "S"   = "#8B5AA3",
  "G2M" = "#E74C3C"
)

polarization_colors <- c(
  "M1-like" = "#B2182B",
  "M2-like" = "#2166AC",
  "Intermediate" = "#999999"
)

feature_colors <- c("lightgrey", "#2E5A87", "#B2182B")
gradient_injury <- c("lightgrey", "#FDAE6B", "#E6550D", "#8B0000")
gradient_repair <- c("lightgrey", "#A1D99B", "#31A354", "#006D2C")


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
      strip.text = element_text(size = base_size + 2, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}


log_message <- function(msg, log_file = NULL) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  message(paste(timestamp, msg))
  if (!is.null(log_file)) {
    cat(paste(timestamp, msg, "\n"), file = log_file, append = TRUE)
  }
}


suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
  library(ggpubr)
})

log_file <- file.path(dirs$advanced, "advanced_analysis_log.txt")
log_message("Starting Advanced Analyses", log_file)


log_message("Loading data...", log_file)

seurat_obj <- readRDS(file.path(dirs$processed, "seurat_with_trajectory.rds"))
message(sprintf("Loaded: %d cells, %d genes", ncol(seurat_obj), nrow(seurat_obj)))


log_message("Performing cell cycle scoring...", log_file)


s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes


s.genes.mouse <- stringr::str_to_title(tolower(s.genes))
g2m.genes.mouse <- stringr::str_to_title(tolower(g2m.genes))


s.genes.present <- s.genes.mouse[s.genes.mouse %in% rownames(seurat_obj)]
g2m.genes.present <- g2m.genes.mouse[g2m.genes.mouse %in% rownames(seurat_obj)]

message(sprintf("S-phase genes: %d/%d found", length(s.genes.present), length(s.genes.mouse)))
message(sprintf("G2M-phase genes: %d/%d found", length(g2m.genes.present), length(g2m.genes.mouse)))


seurat_obj <- CellCycleScoring(seurat_obj,
                               s.features = s.genes.present,
                               g2m.features = g2m.genes.present,
                               set.ident = FALSE)


p_cc_umap <- DimPlot(seurat_obj,
                     reduction = "umap",
                     group.by = "Phase",
                     pt.size = 0.3) +
  scale_color_manual(values = phase_colors) +
  labs(title = "Cell Cycle Phase Distribution",
       x = "UMAP 1", y = "UMAP 2") +
  theme_publication(base_size = 14)

ggsave(file.path(dirs$advanced, "Fig01_cell_cycle_umap.pdf"),
       p_cc_umap, width = 10, height = 8, device = cairo_pdf)
ggsave(file.path(dirs$advanced, "Fig01_cell_cycle_umap.png"),
       p_cc_umap, width = 10, height = 8, dpi = 300)


cc_summary <- seurat_obj@meta.data %>%
  group_by(sample, cell_type, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample, cell_type) %>%
  mutate(prop = n / sum(n) * 100)

p_cc_bar <- ggplot(cc_summary, aes(x = cell_type, y = prop, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = phase_colors) +
  facet_wrap(~sample) +
  labs(title = "Cell Cycle Phase by Cell Type",
       x = NULL, y = "Proportion (%)", fill = "Phase") +
  theme_publication(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(dirs$advanced, "Fig02_cell_cycle_proportion.pdf"),
       p_cc_bar, width = 14, height = 6, device = cairo_pdf)

log_message("Cell cycle scoring complete", log_file)


log_message("Analyzing macrophage-associated transcriptional modules...", log_file)


m1_genes <- c("Cd86", "Nos2", "Il1b", "Tnf", "Il6", "Cxcl9", "Cxcl10")
m2_genes <- c("Mrc1", "Arg1", "Cd163", "Retnla", "Chil3", "Il10", "Tgfb1")

m1_present <- m1_genes[m1_genes %in% rownames(seurat_obj)]
m2_present <- m2_genes[m2_genes %in% rownames(seurat_obj)]

message(sprintf("M1 genes: %d/%d found", length(m1_present), length(m1_genes)))
message(sprintf("M2 genes: %d/%d found", length(m2_present), length(m2_genes)))


if (length(m1_present) >= 3) {
  seurat_obj <- AddModuleScore(seurat_obj, features = list(m1_present),
                               name = "M1_score", search = TRUE)
}

if (length(m2_present) >= 3) {
  seurat_obj <- AddModuleScore(seurat_obj, features = list(m2_present),
                               name = "M2_score", search = TRUE)
}


if ("M1_score1" %in% colnames(seurat_obj@meta.data) &&
    "M2_score1" %in% colnames(seurat_obj@meta.data)) {

  seurat_obj$macrophage_module_balance <- seurat_obj$M2_score1 - seurat_obj$M1_score1

  mac_data <- seurat_obj@meta.data %>%
    filter(cell_type == "Macrophage")

  mac_summary <- mac_data %>%
    group_by(condition) %>%
    summarise(
      M1_mean = mean(M1_score1),
      M2_mean = mean(M2_score1),
      M2_minus_M1 = mean(macrophage_module_balance),
      n = n(),
      .groups = "drop"
    )

  write.csv(
    mac_summary,
    file.path(dirs$advanced, "macrophage_module_scores.csv"),
    row.names = FALSE
  )

  p_polar <- ggplot(mac_data, aes(x = sample, y = macrophage_module_balance, fill = sample)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = condition_colors) +
    stat_compare_means(method = "wilcox.test", label = "p.signif",
                       label.x = 1.5,
                       label.y = max(mac_data$macrophage_module_balance, na.rm = TRUE) * 0.9) +
    labs(title = "Macrophage-Associated Module Balance",
         subtitle = "M2-associated minus M1-associated module score",
         x = NULL, y = "Module balance (M2 - M1)") +
    theme_publication(base_size = 14) +
    theme(legend.position = "none")

  ggsave(file.path(dirs$advanced, "Fig03_macrophage_module_balance.pdf"),
         p_polar, width = 8, height = 7, device = cairo_pdf)
  ggsave(file.path(dirs$advanced, "Fig03_macrophage_module_balance.png"),
         p_polar, width = 8, height = 7, dpi = 300)


  seurat_mac <- subset(seurat_obj, cell_type == "Macrophage")

  p_polar_umap <- FeaturePlot(seurat_mac,
                              features = "macrophage_module_balance",
                              reduction = "umap",
                              split.by = "sample",
                              pt.size = 1) &
    scale_color_gradientn(colors = c("#2166AC", "white", "#B2182B"),
                          limits = c(-1, 1) * max(abs(seurat_mac$macrophage_module_balance), na.rm = TRUE)) &
    theme_publication(base_size = 12) &
    labs(x = "UMAP 1", y = "UMAP 2")

  ggsave(file.path(dirs$advanced, "Fig04_macrophage_module_balance_umap.pdf"),
         p_polar_umap, width = 12, height = 5, device = cairo_pdf)

  log_message("Macrophage module analysis complete", log_file)
}


log_message("Analyzing injury-repair axis...", log_file)


injury_genes <- c("Havcr1", "Lcn2", "Sox9", "Krt8", "Krt18", "Cd44", "Ccn1", "Vcam1")


repair_genes <- c("Pcna", "Mki67", "Ccnd1", "Slc34a1", "Lrp2", "Aqp1", "Cubn")

injury_present <- injury_genes[injury_genes %in% rownames(seurat_obj)]
repair_present <- repair_genes[repair_genes %in% rownames(seurat_obj)]

message(sprintf("Injury genes: %d/%d found", length(injury_present), length(injury_genes)))
message(sprintf("Repair genes: %d/%d found", length(repair_present), length(repair_genes)))


if (length(injury_present) >= 3) {
  seurat_obj <- AddModuleScore(seurat_obj, features = list(injury_present),
                               name = "Injury_score", search = TRUE)
}

if (length(repair_present) >= 3) {
  seurat_obj <- AddModuleScore(seurat_obj, features = list(repair_present),
                               name = "Repair_score", search = TRUE)
}


if ("Injury_score1" %in% colnames(seurat_obj@meta.data)) {
  p_injury <- FeaturePlot(seurat_obj,
                          features = "Injury_score1",
                          reduction = "umap",
                          split.by = "sample",
                          order = TRUE,
                          pt.size = 0.3) &
    scale_color_gradientn(colors = gradient_injury) &
    theme_publication(base_size = 12) &
    labs(x = "UMAP 1", y = "UMAP 2", color = "Injury\nScore")

  ggsave(file.path(dirs$advanced, "Fig05_injury_score_umap.pdf"),
         p_injury, width = 12, height = 5, device = cairo_pdf)
  ggsave(file.path(dirs$advanced, "Fig05_injury_score_umap.png"),
         p_injury, width = 12, height = 5, dpi = 300)
}


if ("Repair_score1" %in% colnames(seurat_obj@meta.data)) {
  p_repair <- FeaturePlot(seurat_obj,
                          features = "Repair_score1",
                          reduction = "umap",
                          split.by = "sample",
                          order = TRUE,
                          pt.size = 0.3) &
    scale_color_gradientn(colors = gradient_repair) &
    theme_publication(base_size = 12) &
    labs(x = "UMAP 1", y = "UMAP 2", color = "Repair\nScore")

  ggsave(file.path(dirs$advanced, "Fig06_repair_score_umap.pdf"),
         p_repair, width = 12, height = 5, device = cairo_pdf)
  ggsave(file.path(dirs$advanced, "Fig06_repair_score_umap.png"),
         p_repair, width = 12, height = 5, dpi = 300)
}

log_message("Injury-repair analysis complete", log_file)


log_message("Comparing PTC subtypes...", log_file)

ptc_cells <- subset(seurat_obj, cell_type %in% c("PTC-1", "PTC-2", "PTC-3"))

if (ncol(ptc_cells) > 100) {

  ptc_markers <- c("Slc34a1", "Havcr1", "Lcn2", "Sox9", "Ccn1", "Ghr")
  ptc_markers_present <- ptc_markers[ptc_markers %in% rownames(ptc_cells)]

  p_ptc_violin <- VlnPlot(ptc_cells,
                          features = ptc_markers_present,
                          group.by = "cell_type",
                          cols = celltype_colors[c("PTC-1", "PTC-2", "PTC-3")],
                          pt.size = 0,
                          ncol = 3) &
    theme_publication(base_size = 12)

  ggsave(file.path(dirs$advanced, "Fig07_PTC_markers_violin.pdf"),
         p_ptc_violin, width = 14, height = 8, device = cairo_pdf)
  ggsave(file.path(dirs$advanced, "Fig07_PTC_markers_violin.png"),
         p_ptc_violin, width = 14, height = 8, dpi = 300)


  if ("Injury_score1" %in% colnames(ptc_cells@meta.data)) {
    p_ptc_injury <- ggplot(ptc_cells@meta.data,
                           aes(x = cell_type, y = Injury_score1, fill = cell_type)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
      scale_fill_manual(values = celltype_colors[c("PTC-1", "PTC-2", "PTC-3")]) +
      facet_wrap(~sample) +
      stat_compare_means(comparisons = list(c("PTC-1", "PTC-3"), c("PTC-2", "PTC-3")),
                         method = "wilcox.test", label = "p.signif") +
      labs(title = "Injury Score by PTC Subtype",
           x = NULL, y = "Injury Score") +
      theme_publication(base_size = 14) +
      theme(legend.position = "none")

    ggsave(file.path(dirs$advanced, "Fig08_PTC_injury_comparison.pdf"),
           p_ptc_injury, width = 10, height = 6, device = cairo_pdf)
  }

  log_message("PTC comparison complete", log_file)
}


log_message("Saving data...", log_file)

saveRDS(seurat_obj, file.path(dirs$processed, "seurat_with_advanced_scores.rds"))

log_message("Data saved", log_file)


sink(file.path(dirs$advanced, "Advanced_Analysis_Report.txt"))
cat("===============================================================================\n")
cat("Advanced Analysis Report\n")
cat("CS-NO Hydrogel scRNA-seq Analysis\n")
cat("===============================================================================\n\n")

cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Random seed: 20260614\n\n")

cat("1. ANALYSES PERFORMED:\n")
cat("----------------------\n")
cat("   - Cell cycle scoring\n")
cat("   - Macrophage-associated M1/M2 module scores\n")
cat("   - Injury-repair axis analysis\n")
cat("   - PTC subtype comparison\n")

cat("\n2. MODULE SCORES ADDED:\n")
cat("-----------------------\n")
score_cols <- grep("_score1$", colnames(seurat_obj@meta.data), value = TRUE)
for (col in score_cols) {
  cat(sprintf("   - %s\n", col))
}

cat("\n3. OUTPUT FILES:\n")
cat("-----------------\n")
cat("   - Cell cycle UMAP and proportions (PDF/PNG)\n")
cat("   - Macrophage module-balance plots (PDF/PNG)\n")
cat("   - Injury and repair score UMAPs (PDF/PNG)\n")
cat("   - PTC marker comparisons (PDF/PNG)\n")

cat("\n===============================================================================\n")
cat("Analysis Complete\n")
cat("Next step: Run 07_nichenet_velocity.R or figure scripts\n")
cat("===============================================================================\n")
sink()

log_message("Advanced Analysis completed!", log_file)

sessionInfo()
