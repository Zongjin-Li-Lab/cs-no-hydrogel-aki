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


theme_publication <- function(base_size = 14, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
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
      legend.background = element_blank(),
      panel.background = element_blank(),
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
  library(Matrix)
  library(scater)
  library(DropletUtils)
  library(future)
})


n_workers <- min(4L, future::availableCores())
if (n_workers > 1L) future::plan("multisession", workers = n_workers) else future::plan("sequential")
options(future.globals.maxSize = 8000 * 1024^2)
options(future.rng.onMisuse = "ignore")


log_file <- file.path(dirs$qc, "qc_analysis_log.txt")
log_message("Starting Quality Control Analysis", log_file)
log_message(sprintf("Base directory: %s", base_dir), log_file)


log_message("Loading raw data...", log_file)


sample_paths <- list(
  Ctrl = file.path(dirs$raw_data, "Ctrl/filtered_feature_bc_matrix"),
  NO = file.path(dirs$raw_data, "NO/filtered_feature_bc_matrix")
)


seurat_list <- list()

for (sample_name in names(sample_paths)) {
  log_message(sprintf("Processing sample: %s", sample_name), log_file)


  if (!dir.exists(sample_paths[[sample_name]])) {
    stop(sprintf("Data directory not found: %s", sample_paths[[sample_name]]))
  }


  counts <- Read10X(data.dir = sample_paths[[sample_name]])


  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )


  seurat_obj$sample <- sample_name
  seurat_obj$condition <- ifelse(sample_name == "Ctrl", "Control", "NO_treated")

  seurat_list[[sample_name]] <- seurat_obj

  log_message(sprintf("  %s: %d cells, %d genes",
                      sample_name, ncol(seurat_obj), nrow(seurat_obj)), log_file)
}


log_message("Calculating QC metrics...", log_file)

calculate_qc_metrics <- function(seurat_obj) {

  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")


  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")


  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb[^(p)]")


  seurat_obj[["log10GenesPerUMI"]] <- log10(seurat_obj$nFeature_RNA) /
    log10(seurat_obj$nCount_RNA)

  return(seurat_obj)
}


seurat_list <- lapply(seurat_list, calculate_qc_metrics)


log_message("Creating initial QC visualizations...", log_file)


seurat_combined <- merge(seurat_list[[1]], y = seurat_list[[2]],
                         add.cell.ids = names(seurat_list))


p_violin_features <- VlnPlot(seurat_combined, features = "nFeature_RNA",
                             group.by = "sample", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Genes per Cell", x = NULL, y = "Number of Genes") +
  theme_publication(base_size = 14) +
  theme(legend.position = "none")

p_violin_counts <- VlnPlot(seurat_combined, features = "nCount_RNA",
                           group.by = "sample", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "UMIs per Cell", x = NULL, y = "Number of UMIs") +
  theme_publication(base_size = 14) +
  theme(legend.position = "none")

p_violin_mito <- VlnPlot(seurat_combined, features = "percent.mt",
                         group.by = "sample", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Mitochondrial Content", x = NULL, y = "% Mitochondrial") +
  theme_publication(base_size = 14) +
  theme(legend.position = "none")

p_violin_ribo <- VlnPlot(seurat_combined, features = "percent.ribo",
                         group.by = "sample", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Ribosomal Content", x = NULL, y = "% Ribosomal") +
  theme_publication(base_size = 14) +
  theme(legend.position = "none")


p_qc_violin <- (p_violin_features | p_violin_counts) / (p_violin_mito | p_violin_ribo) +
  plot_annotation(
    title = "Quality Control Metrics (Pre-filtering)",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  )

ggsave(file.path(dirs$qc, "01_QC_violin_prefilter.pdf"),
       p_qc_violin, width = 12, height = 10, device = cairo_pdf)
ggsave(file.path(dirs$qc, "01_QC_violin_prefilter.png"),
       p_qc_violin, width = 12, height = 10, dpi = 300)


scatter_data <- seurat_combined@meta.data

p_scatter <- ggplot(scatter_data, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_gradientn(colors = c("#2E5A87", "#FEE08B", "#B2182B"),
                        name = "% Mito") +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~sample, ncol = 2) +
  labs(title = "Gene-UMI Relationship Colored by Mitochondrial Content",
       x = "Total UMIs (log10)", y = "Number of Genes (log10)") +
  theme_publication(base_size = 14) +
  theme(legend.position = "right")

ggsave(file.path(dirs$qc, "02_QC_scatter_genes_umi.pdf"),
       p_scatter, width = 14, height = 7, device = cairo_pdf)
ggsave(file.path(dirs$qc, "02_QC_scatter_genes_umi.png"),
       p_scatter, width = 14, height = 7, dpi = 300)


p_complexity <- ggplot(scatter_data, aes(x = log10GenesPerUMI, fill = sample)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = condition_colors) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Cell Complexity Distribution",
       subtitle = "Red line indicates minimum threshold (0.8)",
       x = "log10(Genes) / log10(UMIs)", y = "Density", fill = "Condition") +
  theme_publication(base_size = 14)

ggsave(file.path(dirs$qc, "03_QC_complexity.pdf"),
       p_complexity, width = 9, height = 6, device = cairo_pdf)
ggsave(file.path(dirs$qc, "03_QC_complexity.png"),
       p_complexity, width = 9, height = 6, dpi = 300)


log_message("Applying QC filters...", log_file)


qc_thresholds <- list(
  nFeature_RNA_min = 500,
  nFeature_RNA_max = 6000,
  nCount_RNA_min = 1000,
  nCount_RNA_max = 40000,
  percent_mt_max = 20,
  percent_ribo_max = 50,
  complexity_min = 0.8
)


log_message("QC thresholds applied:", log_file)
for (name in names(qc_thresholds)) {
  log_message(sprintf("  %s: %s", name, qc_thresholds[[name]]), log_file)
}


save(qc_thresholds, file = file.path(dirs$qc, "qc_thresholds.RData"))


filter_cells <- function(seurat_obj, thresholds) {
  n_before <- ncol(seurat_obj)

  cells_keep <- seurat_obj$nFeature_RNA >= thresholds$nFeature_RNA_min &
    seurat_obj$nFeature_RNA <= thresholds$nFeature_RNA_max &
    seurat_obj$nCount_RNA >= thresholds$nCount_RNA_min &
    seurat_obj$nCount_RNA <= thresholds$nCount_RNA_max &
    seurat_obj$percent.mt <= thresholds$percent_mt_max &
    seurat_obj$log10GenesPerUMI >= thresholds$complexity_min

  seurat_obj <- seurat_obj[, cells_keep]
  n_after <- ncol(seurat_obj)

  message(sprintf("  Filtered: %d -> %d cells (%.1f%% retained)",
                  n_before, n_after, 100 * n_after / n_before))

  return(seurat_obj)
}


log_message("Filtering cells:", log_file)
for (name in names(seurat_list)) {
  log_message(sprintf("  %s:", name), log_file)
  seurat_list[[name]] <- filter_cells(seurat_list[[name]], qc_thresholds)
}


log_message("Creating post-filtering QC visualizations...", log_file)


seurat_filtered <- merge(seurat_list[[1]], y = seurat_list[[2]],
                         add.cell.ids = names(seurat_list))


qc_summary <- seurat_filtered@meta.data %>%
  group_by(sample) %>%
  summarise(
    n_cells = n(),
    median_genes = median(nFeature_RNA),
    median_umi = median(nCount_RNA),
    median_mito = median(percent.mt),
    mean_genes = mean(nFeature_RNA),
    mean_umi = mean(nCount_RNA),
    mean_mito = mean(percent.mt),
    .groups = "drop"
  )

write.csv(qc_summary, file.path(dirs$qc, "qc_summary_statistics.csv"), row.names = FALSE)
log_message("QC summary saved", log_file)


p_violin_post_features <- VlnPlot(seurat_filtered, features = "nFeature_RNA",
                                  group.by = "sample", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Genes per Cell", x = NULL, y = "Number of Genes") +
  theme_publication(base_size = 14) +
  theme(legend.position = "none")

p_violin_post_counts <- VlnPlot(seurat_filtered, features = "nCount_RNA",
                                group.by = "sample", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "UMIs per Cell", x = NULL, y = "Number of UMIs") +
  theme_publication(base_size = 14) +
  theme(legend.position = "none")

p_violin_post_mito <- VlnPlot(seurat_filtered, features = "percent.mt",
                              group.by = "sample", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Mitochondrial Content", x = NULL, y = "% Mitochondrial") +
  theme_publication(base_size = 14) +
  theme(legend.position = "none")

p_violin_post_ribo <- VlnPlot(seurat_filtered, features = "percent.ribo",
                              group.by = "sample", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Ribosomal Content", x = NULL, y = "% Ribosomal") +
  theme_publication(base_size = 14) +
  theme(legend.position = "none")


p_qc_violin_post <- (p_violin_post_features | p_violin_post_counts) /
  (p_violin_post_mito | p_violin_post_ribo) +
  plot_annotation(
    title = "Quality Control Metrics (Post-filtering)",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  )

ggsave(file.path(dirs$qc, "04_QC_violin_postfilter.pdf"),
       p_qc_violin_post, width = 12, height = 10, device = cairo_pdf)
ggsave(file.path(dirs$qc, "04_QC_violin_postfilter.png"),
       p_qc_violin_post, width = 12, height = 10, dpi = 300)


cell_counts <- data.frame(
  sample = c("Ctrl", "NO"),
  cells = c(ncol(seurat_list$Ctrl), ncol(seurat_list$NO))
)

p_cell_counts <- ggplot(cell_counts, aes(x = sample, y = cells, fill = sample)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = scales::comma(cells)), vjust = -0.5, size = 6, fontface = "bold") +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Cell Counts After QC Filtering",
       x = NULL, y = "Number of Cells") +
  theme_publication(base_size = 14) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::comma)

ggsave(file.path(dirs$qc, "05_cell_counts_postfilter.pdf"),
       p_cell_counts, width = 7, height = 6, device = cairo_pdf)
ggsave(file.path(dirs$qc, "05_cell_counts_postfilter.png"),
       p_cell_counts, width = 7, height = 6, dpi = 300)


prepost_data <- data.frame(
  sample = rep(c("Ctrl", "NO"), 2),
  stage = rep(c("Pre-filter", "Post-filter"), each = 2),
  cells = c(
    ncol(Read10X(sample_paths$Ctrl)),
    ncol(Read10X(sample_paths$NO)),
    ncol(seurat_list$Ctrl),
    ncol(seurat_list$NO)
  )
)
prepost_data$stage <- factor(prepost_data$stage, levels = c("Pre-filter", "Post-filter"))

p_prepost <- ggplot(prepost_data, aes(x = sample, y = cells, fill = stage)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = scales::comma(cells)), position = position_dodge(0.7),
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Pre-filter" = "#E0E0E0", "Post-filter" = "#4A7BB0")) +
  labs(title = "Cell Counts Before and After QC Filtering",
       x = NULL, y = "Number of Cells", fill = "Stage") +
  theme_publication(base_size = 14) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::comma)

ggsave(file.path(dirs$qc, "06_prepost_filtering_comparison.pdf"),
       p_prepost, width = 9, height = 6, device = cairo_pdf)
ggsave(file.path(dirs$qc, "06_prepost_filtering_comparison.png"),
       p_prepost, width = 9, height = 6, dpi = 300)


log_message("Saving processed data...", log_file)


for (name in names(seurat_list)) {
  filepath <- file.path(dirs$processed, sprintf("seurat_%s_qc_filtered.rds", name))
  saveRDS(seurat_list[[name]], filepath)
  log_message(sprintf("Saved: %s", filepath), log_file)
}


saveRDS(seurat_filtered, file.path(dirs$processed, "seurat_combined_qc_filtered.rds"))
log_message("Saved: seurat_combined_qc_filtered.rds", log_file)


saveRDS(seurat_list, file.path(dirs$processed, "seurat_list_qc_filtered.rds"))
log_message("Saved: seurat_list_qc_filtered.rds", log_file)


log_message("Generating QC report...", log_file)

sink(file.path(dirs$qc, "QC_Analysis_Report.txt"))
cat("===============================================================================\n")
cat("Quality Control Analysis Report\n")
cat("CS-NO Hydrogel scRNA-seq Analysis\n")
cat("===============================================================================\n\n")

cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Random seed: 20260614\n")
cat("Base directory:", base_dir, "\n\n")

cat("1. QC THRESHOLDS APPLIED:\n")
cat("-------------------------\n")
for (name in names(qc_thresholds)) {
  cat(sprintf("   %s: %s\n", name, qc_thresholds[[name]]))
}

cat("\n2. CELL COUNTS:\n")
cat("---------------\n")
cat(sprintf("   Ctrl: %d cells\n", ncol(seurat_list$Ctrl)))
cat(sprintf("   NO: %d cells\n", ncol(seurat_list$NO)))
cat(sprintf("   Total: %d cells\n", sum(sapply(seurat_list, ncol))))

cat("\n3. QC SUMMARY STATISTICS:\n")
cat("-------------------------\n")
print(qc_summary)

cat("\n4. OUTPUT FILES:\n")
cat("----------------\n")
cat("   Data files:\n")
cat("   - seurat_Ctrl_qc_filtered.rds\n")
cat("   - seurat_NO_qc_filtered.rds\n")
cat("   - seurat_combined_qc_filtered.rds\n")
cat("   - seurat_list_qc_filtered.rds\n")
cat("   - qc_summary_statistics.csv\n")
cat("   - qc_thresholds.RData\n")
cat("\n   Figure files:\n")
cat("   - 01_QC_violin_prefilter.pdf/png\n")
cat("   - 02_QC_scatter_genes_umi.pdf/png\n")
cat("   - 03_QC_complexity.pdf/png\n")
cat("   - 04_QC_violin_postfilter.pdf/png\n")
cat("   - 05_cell_counts_postfilter.pdf/png\n")
cat("   - 06_prepost_filtering_comparison.pdf/png\n")

cat("\n===============================================================================\n")
cat("QC Analysis Complete\n")
cat("Next step: Run 02_integration_clustering.R\n")
cat("===============================================================================\n")
sink()

log_message("QC Analysis completed successfully!", log_file)


message("\n")
message("===============================================================================")
message("QUALITY CONTROL ANALYSIS COMPLETE")
message("===============================================================================")
message(sprintf("Total cells after QC: %d", ncol(seurat_filtered)))
message(sprintf("  Ctrl: %d cells", ncol(seurat_list$Ctrl)))
message(sprintf("  NO: %d cells", ncol(seurat_list$NO)))
message("\nOutput saved to:")
message(sprintf("  %s", dirs$qc))
message("\nNext step: Run 02_integration_clustering.R")
message("===============================================================================")

cat("\n=== Session Info ===\n")
sessionInfo()
