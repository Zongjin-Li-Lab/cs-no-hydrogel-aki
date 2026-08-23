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


theme_publication <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5),
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
  library(CellChat)
  library(tidyverse)
  library(patchwork)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})

log_file <- file.path(dirs$cellchat, "cellchat_analysis_log.txt")
log_message("Starting CellChat Analysis", log_file)


log_message("Loading data...", log_file)

seurat_obj <- readRDS(file.path(dirs$processed, "seurat_with_trajectory.rds"))
message(sprintf("Loaded: %d cells, %d genes", ncol(seurat_obj), nrow(seurat_obj)))
message(sprintf("Conditions: %s", paste(names(table(seurat_obj$sample)), collapse = ", ")))
message(sprintf("Cell types (12): %s", paste(names(table(seurat_obj$cell_type)), collapse = ", ")))


log_message("Creating CellChat objects...", log_file)

create_cellchat <- function(seurat_obj, condition_name) {

  seurat_subset <- subset(seurat_obj, subset = sample == condition_name)


  valid_cells <- !is.na(seurat_subset$cell_type)
  if (sum(!valid_cells) > 0) {
    message(sprintf("Removing %d cells with NA cell_type", sum(!valid_cells)))
    seurat_subset <- seurat_subset[, valid_cells]
  }

  message(sprintf("Creating CellChat for %s: %d cells", condition_name, ncol(seurat_subset)))


  data_input <- GetAssayData(seurat_subset, layer = "data")


  meta <- data.frame(
    labels = as.character(seurat_subset$cell_type),
    row.names = colnames(seurat_subset)
  )


  cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")


  CellChatDB <- CellChatDB.mouse
  cellchat@DB <- CellChatDB


  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)


  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
  cellchat <- filterCommunication(cellchat, min.cells = 10)


  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  return(cellchat)
}


cellchat_ctrl <- create_cellchat(seurat_obj, "Ctrl")
cellchat_no <- create_cellchat(seurat_obj, "NO")

log_message("CellChat objects created", log_file)


log_message("Merging CellChat objects...", log_file)

object.list <- list(Ctrl = cellchat_ctrl, NO = cellchat_no)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

log_message("CellChat objects merged", log_file)


log_message("Creating interaction comparison plots...", log_file)


p_compare <- compareInteractions(cellchat_merged, show.legend = TRUE, group = c(1, 2))

ggsave(file.path(dirs$cellchat, "Fig01_interaction_comparison.pdf"),
       p_compare, width = 10, height = 5, device = cairo_pdf)
ggsave(file.path(dirs$cellchat, "Fig01_interaction_comparison.png"),
       p_compare, width = 10, height = 5, dpi = 300)


log_message("Creating circle plots...", log_file)


group_colors <- celltype_colors[levels(factor(cellchat_ctrl@idents))]


pdf(file.path(dirs$cellchat, "Fig02_circle_Ctrl.pdf"), width = 10, height = 10)
netVisual_circle(cellchat_ctrl@net$count,
                 vertex.weight = as.numeric(table(cellchat_ctrl@idents)),
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Number of Interactions (Ctrl)",
                 color.use = group_colors)
dev.off()


pdf(file.path(dirs$cellchat, "Fig02_circle_NO.pdf"), width = 10, height = 10)
netVisual_circle(cellchat_no@net$count,
                 vertex.weight = as.numeric(table(cellchat_no@idents)),
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Number of Interactions (NO)",
                 color.use = group_colors)
dev.off()

log_message("Circle plots saved", log_file)


log_message("Creating differential interaction heatmap...", log_file)

p_diff_heatmap <- netVisual_heatmap(cellchat_merged,
                                    comparison = c(1, 2),
                                    measure = "count",
                                    title.name = "Differential Interactions (NO vs Ctrl)")

pdf(file.path(dirs$cellchat, "Fig03_differential_heatmap.pdf"), width = 12, height = 10)
draw(p_diff_heatmap)
dev.off()


log_message("Analyzing key signaling pathways...", log_file)


key_pathways <- c("CCL", "CXCL", "TNF", "IL1", "VEGF", "TGFb", "WNT", "FGF")


available_pathways_ctrl <- unique(cellchat_ctrl@netP$pathways)
available_pathways_no <- unique(cellchat_no@netP$pathways)

message("Available pathways in Ctrl: ", paste(head(available_pathways_ctrl, 10), collapse = ", "))
message("Available pathways in NO: ", paste(head(available_pathways_no, 10), collapse = ", "))


pathways_to_plot <- key_pathways[key_pathways %in% available_pathways_ctrl &
                                   key_pathways %in% available_pathways_no]

message("Pathways to analyze: ", paste(pathways_to_plot, collapse = ", "))


for (pathway in pathways_to_plot) {
  tryCatch({

    pdf(file.path(dirs$cellchat, sprintf("Fig04_chord_%s_Ctrl.pdf", pathway)), width = 10, height = 10)
    netVisual_chord_gene(cellchat_ctrl,
                         signaling = pathway,
                         title.name = sprintf("%s Signaling (Ctrl)", pathway),
                         color.use = group_colors,
                         show.legend = TRUE)
    dev.off()


    pdf(file.path(dirs$cellchat, sprintf("Fig04_chord_%s_NO.pdf", pathway)), width = 10, height = 10)
    netVisual_chord_gene(cellchat_no,
                         signaling = pathway,
                         title.name = sprintf("%s Signaling (NO)", pathway),
                         color.use = group_colors,
                         show.legend = TRUE)
    dev.off()

    message(sprintf("  ✓ %s chord diagrams saved", pathway))

  }, error = function(e) {
    message(sprintf("  ✗ %s: %s", pathway, e$message))
  })
}


log_message("Analyzing signaling roles...", log_file)


cellchat_ctrl <- netAnalysis_computeCentrality(cellchat_ctrl)
cellchat_no <- netAnalysis_computeCentrality(cellchat_no)


p_role_ctrl <- netAnalysis_signalingRole_heatmap(cellchat_ctrl,
                                                  pattern = "outgoing",
                                                  title = "Outgoing Signaling (Ctrl)")

pdf(file.path(dirs$cellchat, "Fig05_signaling_role_Ctrl.pdf"), width = 12, height = 10)
draw(p_role_ctrl)
dev.off()


p_role_no <- netAnalysis_signalingRole_heatmap(cellchat_no,
                                                pattern = "outgoing",
                                                title = "Outgoing Signaling (NO)")

pdf(file.path(dirs$cellchat, "Fig05_signaling_role_NO.pdf"), width = 12, height = 10)
draw(p_role_no)
dev.off()

log_message("Signaling role analysis saved", log_file)


log_message("Creating bubble plots...", log_file)


sources_of_interest <- c("Macrophage", "PTC-1", "PTC-2", "PTC-3")
targets_of_interest <- c("Macrophage", "PTC-1", "PTC-2", "PTC-3", "VEC", "Fibroblast")


sources_of_interest <- sources_of_interest[sources_of_interest %in% levels(factor(cellchat_ctrl@idents))]
targets_of_interest <- targets_of_interest[targets_of_interest %in% levels(factor(cellchat_ctrl@idents))]

if (length(sources_of_interest) > 0 && length(targets_of_interest) > 0) {

  p_bubble_ctrl <- netVisual_bubble(cellchat_ctrl,
                                    sources.use = sources_of_interest,
                                    targets.use = targets_of_interest,
                                    remove.isolate = TRUE,
                                    title.name = "Communication (Ctrl)")

  ggsave(file.path(dirs$cellchat, "Fig06_bubble_Ctrl.pdf"),
         p_bubble_ctrl, width = 14, height = 10, device = cairo_pdf)


  p_bubble_no <- netVisual_bubble(cellchat_no,
                                  sources.use = sources_of_interest,
                                  targets.use = targets_of_interest,
                                  remove.isolate = TRUE,
                                  title.name = "Communication (NO)")

  ggsave(file.path(dirs$cellchat, "Fig06_bubble_NO.pdf"),
         p_bubble_no, width = 14, height = 10, device = cairo_pdf)

  log_message("Bubble plots saved", log_file)
}


log_message("Creating information flow comparison...", log_file)

p_info_flow <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE,
                       do.stat = TRUE)

ggsave(file.path(dirs$cellchat, "Fig07_information_flow.pdf"),
       p_info_flow, width = 12, height = 10, device = cairo_pdf)
ggsave(file.path(dirs$cellchat, "Fig07_information_flow.png"),
       p_info_flow, width = 12, height = 10, dpi = 300)

log_message("Information flow comparison saved", log_file)


log_message("Applying Benjamini-Hochberg correction within each library...", log_file)

export_interactions <- function(cellchat_object, condition_label) {
  interactions <- subsetCommunication(cellchat_object)
  interactions$FDR <- p.adjust(interactions$pval, method = "BH")
  interactions$condition <- condition_label
  interactions
}

interactions_bh <- rbind(
  export_interactions(cellchat_ctrl, "Control"),
  export_interactions(cellchat_no, "NO_treated")
)

write.csv(
  interactions_bh,
  file.path(dirs$cellchat, "cellchat_all_bh.csv"),
  row.names = FALSE
)

ccl_fdr <- subset(interactions_bh, pathway_name == "CCL")
write.csv(
  ccl_fdr,
  file.path(dirs$cellchat, "cellchat_ccl_fdr.csv"),
  row.names = FALSE
)

significant_ccl <- subset(ccl_fdr, FDR < 0.05)
log_message(
  paste(
    "Significant CCL source-target ligand-receptor combinations:",
    paste(names(table(significant_ccl$condition)),
          as.integer(table(significant_ccl$condition)), collapse = "; ")
  ),
  log_file
)


log_message("Saving CellChat objects...", log_file)

saveRDS(cellchat_ctrl, file.path(dirs$processed, "cellchat_Ctrl.rds"))
saveRDS(cellchat_no, file.path(dirs$processed, "cellchat_NO.rds"))
saveRDS(cellchat_merged, file.path(dirs$processed, "cellchat_merged.rds"))

log_message("CellChat objects saved", log_file)


sink(file.path(dirs$cellchat, "CellChat_Analysis_Report.txt"))
cat("===============================================================================\n")
cat("CellChat Analysis Report\n")
cat("CS-NO Hydrogel scRNA-seq Analysis\n")
cat("===============================================================================\n\n")

cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Random seed: 20260614\n\n")

cat("1. CELL COUNTS:\n")
cat("---------------\n")
cat(sprintf("   Ctrl: %d cells\n", ncol(subset(seurat_obj, sample == "Ctrl"))))
cat(sprintf("   NO: %d cells\n", ncol(subset(seurat_obj, sample == "NO"))))

cat("\n2. INTERACTION COUNTS:\n")
cat("----------------------\n")
cat(sprintf("   Ctrl: %d interactions\n", sum(cellchat_ctrl@net$count)))
cat(sprintf("   NO: %d interactions\n", sum(cellchat_no@net$count)))
cat(sprintf("   Significant CCL combinations after BH correction: %d Ctrl; %d NO\n",
            sum(significant_ccl$condition == "Control"),
            sum(significant_ccl$condition == "NO_treated")))

cat("\n3. KEY PATHWAYS ANALYZED:\n")
cat("--------------------------\n")
cat(sprintf("   %s\n", paste(pathways_to_plot, collapse = ", ")))

cat("\n4. OUTPUT FILES:\n")
cat("-----------------\n")
cat("   - Circle plots (PDF)\n")
cat("   - Chord diagrams for key pathways (PDF)\n")
cat("   - Differential heatmap (PDF)\n")
cat("   - Signaling role heatmaps (PDF)\n")
cat("   - Bubble plots (PDF)\n")
cat("   - Information flow comparison (PDF)\n")

cat("\n===============================================================================\n")
cat("Analysis Complete\n")
cat("Next step: Run 06_advanced_analyses.R\n")
cat("===============================================================================\n")
sink()

log_message("CellChat Analysis completed!", log_file)

sessionInfo()
