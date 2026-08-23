#!/usr/bin/env Rscript


rm(list = ls())
gc()


source("scripts/_common.R")
base_dir <- project_root()
setwd(base_dir)


dirs <- make_output_dirs(base_dir)

output_dir <- file.path(dirs$figures, "Figure5_panels")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Output folder: ", output_dir)


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


suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
})


message("\n=== Loading CellChat Objects ===")

cellchat_ctrl <- readRDS(file.path(dirs$processed, "cellchat_Ctrl.rds"))
cellchat_no <- readRDS(file.path(dirs$processed, "cellchat_NO.rds"))
message("CellChat objects loaded!")


get_celltype_colors <- function(cellchat_obj) {
  cell_types <- levels(cellchat_obj@idents)


  colors <- sapply(cell_types, function(ct) {
    if (ct %in% names(celltype_colors)) {
      celltype_colors[ct]
    } else {

      rainbow(1, s = 0.7, v = 0.8)
    }
  })

  names(colors) <- cell_types
  return(colors)
}

colors_ctrl <- get_celltype_colors(cellchat_ctrl)
colors_no <- get_celltype_colors(cellchat_no)

message("Ctrl cell types: ", paste(names(colors_ctrl), collapse = ", "))
message("NO cell types: ", paste(names(colors_no), collapse = ", "))


message("\n=== Available Signaling Pathways ===")

ctrl_pathways <- cellchat_ctrl@netP$pathways
no_pathways <- cellchat_no@netP$pathways
common_pathways <- intersect(ctrl_pathways, no_pathways)

message("Common pathways: ", paste(common_pathways, collapse = ", "))


key_pathways <- c("CCL", "CXCL", "TNF", "IL1", "MIF", "SPP1", "VEGF", "TGFb", "COLLAGEN")
pathways_to_plot <- intersect(key_pathways, common_pathways)

if (length(pathways_to_plot) == 0) {
  message("No key pathways found, using first 5 common pathways")
  pathways_to_plot <- head(common_pathways, 5)
}

message("Pathways to plot: ", paste(pathways_to_plot, collapse = ", "))


message("\n=== Creating Pathway-Specific Chord Diagrams ===")

for (pathway in pathways_to_plot) {
  message(sprintf("\nPlotting %s pathway...", pathway))


  tryCatch({
    pdf(file.path(output_dir, sprintf("Fig5C_chord_%s_Ctrl.pdf", pathway)),
        width = 12, height = 12)
    netVisual_chord_gene(cellchat_ctrl,
                         signaling = pathway,
                         color.use = colors_ctrl,
                         lab.cex = 1.0,
                         small.gap = 1,
                         big.gap = 10,
                         title.name = sprintf("%s Signaling (Ctrl)", pathway))
    dev.off()

    png(file.path(output_dir, sprintf("Fig5C_chord_%s_Ctrl.png", pathway)),
        width = 12, height = 12, units = "in", res = 300)
    netVisual_chord_gene(cellchat_ctrl,
                         signaling = pathway,
                         color.use = colors_ctrl,
                         lab.cex = 1.0,
                         small.gap = 1,
                         big.gap = 10,
                         title.name = sprintf("%s Signaling (Ctrl)", pathway))
    dev.off()

    message(sprintf("  ✓ Ctrl %s saved!", pathway))
  }, error = function(e) {
    message(sprintf("  ✗ Ctrl %s error: %s", pathway, e$message))
    try(dev.off(), silent = TRUE)
  })


  tryCatch({
    pdf(file.path(output_dir, sprintf("Fig5C_chord_%s_NO.pdf", pathway)),
        width = 12, height = 12)
    netVisual_chord_gene(cellchat_no,
                         signaling = pathway,
                         color.use = colors_no,
                         lab.cex = 1.0,
                         small.gap = 1,
                         big.gap = 10,
                         title.name = sprintf("%s Signaling (NO)", pathway))
    dev.off()

    png(file.path(output_dir, sprintf("Fig5C_chord_%s_NO.png", pathway)),
        width = 12, height = 12, units = "in", res = 300)
    netVisual_chord_gene(cellchat_no,
                         signaling = pathway,
                         color.use = colors_no,
                         lab.cex = 1.0,
                         small.gap = 1,
                         big.gap = 10,
                         title.name = sprintf("%s Signaling (NO)", pathway))
    dev.off()

    message(sprintf("  ✓ NO %s saved!", pathway))
  }, error = function(e) {
    message(sprintf("  ✗ NO %s error: %s", pathway, e$message))
    try(dev.off(), silent = TRUE)
  })
}


message("\n=== Creating Combined L-R Chord Diagrams ===")


tryCatch({
  pdf(file.path(output_dir, "Fig5C_chord_all_LR_Ctrl.pdf"), width = 14, height = 14)
  netVisual_chord_gene(cellchat_ctrl,
                       slot.name = "netP",
                       color.use = colors_ctrl,
                       lab.cex = 0.8,
                       small.gap = 0.5,
                       big.gap = 5,
                       legend.pos.x = 15,
                       title.name = "All L-R Interactions (Ctrl)")
  dev.off()

  png(file.path(output_dir, "Fig5C_chord_all_LR_Ctrl.png"),
      width = 14, height = 14, units = "in", res = 300)
  netVisual_chord_gene(cellchat_ctrl,
                       slot.name = "netP",
                       color.use = colors_ctrl,
                       lab.cex = 0.8,
                       small.gap = 0.5,
                       big.gap = 5,
                       legend.pos.x = 15,
                       title.name = "All L-R Interactions (Ctrl)")
  dev.off()

  message("✓ Ctrl all L-R chord saved!")
}, error = function(e) {
  message(sprintf("✗ Ctrl all L-R error: %s", e$message))
  try(dev.off(), silent = TRUE)
})


tryCatch({
  pdf(file.path(output_dir, "Fig5C_chord_all_LR_NO.pdf"), width = 14, height = 14)
  netVisual_chord_gene(cellchat_no,
                       slot.name = "netP",
                       color.use = colors_no,
                       lab.cex = 0.8,
                       small.gap = 0.5,
                       big.gap = 5,
                       legend.pos.x = 15,
                       title.name = "All L-R Interactions (NO)")
  dev.off()

  png(file.path(output_dir, "Fig5C_chord_all_LR_NO.png"),
      width = 14, height = 14, units = "in", res = 300)
  netVisual_chord_gene(cellchat_no,
                       slot.name = "netP",
                       color.use = colors_no,
                       lab.cex = 0.8,
                       small.gap = 0.5,
                       big.gap = 5,
                       legend.pos.x = 15,
                       title.name = "All L-R Interactions (NO)")
  dev.off()

  message("✓ NO all L-R chord saved!")
}, error = function(e) {
  message(sprintf("✗ NO all L-R error: %s", e$message))
  try(dev.off(), silent = TRUE)
})


message("\n=== Creating PTC-Focused Chord Diagrams ===")


ptc_types <- grep("PTC", levels(cellchat_ctrl@idents), value = TRUE)
message("PTC sources: ", paste(ptc_types, collapse = ", "))

if (length(ptc_types) > 0) {

  tryCatch({
    pdf(file.path(output_dir, "Fig5C_chord_from_PTC_Ctrl.pdf"), width = 12, height = 12)
    netVisual_chord_gene(cellchat_ctrl,
                         sources.use = ptc_types,
                         color.use = colors_ctrl,
                         lab.cex = 0.9,
                         small.gap = 1,
                         big.gap = 8,
                         title.name = "Signals FROM PTCs (Ctrl)")
    dev.off()
    message("✓ Ctrl from PTC chord saved!")
  }, error = function(e) {
    message(sprintf("✗ Ctrl from PTC error: %s", e$message))
    try(dev.off(), silent = TRUE)
  })


  tryCatch({
    pdf(file.path(output_dir, "Fig5C_chord_from_PTC_NO.pdf"), width = 12, height = 12)
    netVisual_chord_gene(cellchat_no,
                         sources.use = ptc_types,
                         color.use = colors_no,
                         lab.cex = 0.9,
                         small.gap = 1,
                         big.gap = 8,
                         title.name = "Signals FROM PTCs (NO)")
    dev.off()
    message("✓ NO from PTC chord saved!")
  }, error = function(e) {
    message(sprintf("✗ NO from PTC error: %s", e$message))
    try(dev.off(), silent = TRUE)
  })


  tryCatch({
    pdf(file.path(output_dir, "Fig5C_chord_to_PTC_Ctrl.pdf"), width = 12, height = 12)
    netVisual_chord_gene(cellchat_ctrl,
                         targets.use = ptc_types,
                         color.use = colors_ctrl,
                         lab.cex = 0.9,
                         small.gap = 1,
                         big.gap = 8,
                         title.name = "Signals TO PTCs (Ctrl)")
    dev.off()
    message("✓ Ctrl to PTC chord saved!")
  }, error = function(e) {
    message(sprintf("✗ Ctrl to PTC error: %s", e$message))
    try(dev.off(), silent = TRUE)
  })


  tryCatch({
    pdf(file.path(output_dir, "Fig5C_chord_to_PTC_NO.pdf"), width = 12, height = 12)
    netVisual_chord_gene(cellchat_no,
                         targets.use = ptc_types,
                         color.use = colors_no,
                         lab.cex = 0.9,
                         small.gap = 1,
                         big.gap = 8,
                         title.name = "Signals TO PTCs (NO)")
    dev.off()
    message("✓ NO to PTC chord saved!")
  }, error = function(e) {
    message(sprintf("✗ NO to PTC error: %s", e$message))
    try(dev.off(), silent = TRUE)
  })
}


message("\n========================================")
message("CELLCHAT CHORD DIAGRAMS COMPLETE!")
message("========================================")
message("\nOutput folder: ", output_dir)
message("\nGenerated files:")
message("  - Pathway-specific chords (CCL, CXCL, TNF, etc.)")
message("  - All L-R pairs chord (combined view)")
message("  - PTC-focused chords (from/to PTCs)")
message("\nUse the 'all_LR' versions for comprehensive view!")
message("========================================")

sessionInfo()
