#!/usr/bin/env Rscript


source("scripts/_common.R")
base_dir <- project_root()
dirs <- make_output_dirs(base_dir)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

input_candidates <- c(
  file.path(dirs$processed, "seurat_with_trajectory.rds"),
  file.path(dirs$processed, "seurat_integrated_annotated.rds")
)
input_path <- input_candidates[file.exists(input_candidates)][1]
if (is.na(input_path)) {
  stop("No annotated Seurat object found. Run scripts 01-04 or extract the processed-object release asset.")
}

seurat_obj <- readRDS(input_path)

signatures <- list(
  Kirita_FR_PTC = c(
    "Vcam1", "Havcr1", "Lcn2", "Spp1", "Ccl2", "Cxcl1", "Cxcl10",
    "Ccl5", "Sox9", "Cd44", "Krt8", "Krt18", "Tgfb1", "Pdgfb", "Il33"
  ),
  Kirita_Healthy_PTC = c(
    "Slc34a1", "Lrp2", "Cubn", "Slc5a2", "Aqp1", "Slc22a6", "Slc22a8",
    "Slc7a13", "Miox", "Gatm"
  ),
  Severe_AKI = c(
    "Havcr1", "Lcn2", "Spp1", "Clu", "Cd44", "Vim", "Fn1", "Timp1",
    "Lgals3", "S100a6"
  ),
  Inflammatory_PTC = c(
    "Vcam1", "Icam1", "Ccl2", "Cxcl1", "Cxcl2", "Cxcl10", "Ccl5",
    "Il1b", "Tnf", "Csf1"
  ),
  MAPK_Activity = c(
    "Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Egr2", "Dusp1",
    "Dusp6", "Atf3"
  ),
  Maladaptive_Repair = c(
    "Tgfb1", "Tgfb2", "Ctgf", "Col1a1", "Col3a1", "Fn1", "Acta2",
    "Pdgfb", "Pdgfrb", "Timp1"
  ),
  Regeneration = c(
    "Pcna", "Mki67", "Top2a", "Ccnd1", "Ccne1", "Myc", "E2f1",
    "Mcm2", "Mcm6", "Cdk1"
  )
)

set.seed(20260614)
for (signature_name in names(signatures)) {
  genes <- intersect(signatures[[signature_name]], rownames(seurat_obj))
  if (length(genes) < 3L) {
    warning("Skipping ", signature_name, ": fewer than three genes found")
    next
  }
  seurat_obj <- AddModuleScore(
    seurat_obj,
    features = list(genes),
    name = paste0(signature_name, "_"),
    search = TRUE,
    seed = 20260614
  )
}

output_path <- file.path(dirs$processed, "seurat_with_kirita_scores.rds")
saveRDS(seurat_obj, output_path)

ptc <- subset(seurat_obj, subset = cell_type %in% c("PTC-1", "PTC-2", "PTC-3"))
score_columns <- grep("_1$", colnames(ptc@meta.data), value = TRUE)

score_summary <- ptc@meta.data %>%
  group_by(cell_type, condition) %>%
  summarise(across(all_of(score_columns), list(mean = mean, sd = sd)),
            n = n(), .groups = "drop")

write.csv(
  score_summary,
  file.path(dirs$advanced, "ptc_module_score_summary.csv"),
  row.names = FALSE
)

if ("Kirita_FR_PTC_1" %in% colnames(ptc@meta.data)) {
  plot_data <- ptc@meta.data
  p <- ggplot(plot_data, aes(x = cell_type, y = Kirita_FR_PTC_1, fill = condition)) +
    geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = TRUE) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.12,
                 outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = c("Control" = "#6B9CD3", "NO_treated" = "#8B5AA3")) +
    labs(
      x = NULL,
      y = "Kirita failed-repair PTC module score",
      fill = "Condition"
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(dirs$figures, "kirita_failed_repair_ptc_scores.pdf"),
    p,
    width = 7,
    height = 4.8
  )
}

cat("Saved annotated object with module scores: ", output_path, "\n", sep = "")
