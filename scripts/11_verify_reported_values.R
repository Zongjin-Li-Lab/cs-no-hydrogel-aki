#!/usr/bin/env Rscript

source("scripts/_common.R")
root <- project_root()

composition <- read.csv(file.path(root, "results", "key_tables", "celltype_composition.csv"),
                        check.names = FALSE)
bootstrap <- read.csv(file.path(root, "results", "key_tables", "nr_ptc_bootstrap.csv"),
                      check.names = FALSE)
qc <- read.csv(file.path(root, "results", "key_tables", "library_qc_metrics.csv"),
               check.names = FALSE)
ccl <- read.csv(file.path(root, "results", "key_tables", "cellchat_ccl_fdr.csv"),
                check.names = FALSE)

ptc2 <- composition[composition$cell_type == "PTC-2", ]
stopifnot(nrow(ptc2) == 1L)

control_fraction <- ptc2$pct_Ctrl / 100
treated_fraction <- ptc2$pct_NO / 100
relative_reduction <- 100 * (1 - treated_fraction / control_fraction)

stopifnot(abs(relative_reduction - 68.33) < 0.1)
stopifnot(ptc2$n_Ctrl == 2790, ptc2$n_NO == 464)
stopifnot(sum(qc$analysis_cells) == 25375)
stopifnot(bootstrap$CI95_low[1] > 65, bootstrap$CI95_high[1] < 72)

significant_ccl <- subset(ccl, FDR < 0.05)
ccl_counts <- table(significant_ccl$condition)
stopifnot(unname(ccl_counts["Control"]) == 359)
stopifnot(unname(ccl_counts["NO_treated"]) == 165)

cat("Verified reported values:\n")
cat(sprintf("- Retained cells: %s\n", format(sum(qc$analysis_cells), big.mark = ",")))
cat(sprintf("- PTC-2 counts: %d control, %d CS-NO\n", ptc2$n_Ctrl, ptc2$n_NO))
cat(sprintf("- PTC-2 relative-representation reduction: %.2f%%\n", relative_reduction))
cat(sprintf("- Bootstrap 95%% interval: %.2f%% to %.2f%%\n",
            bootstrap$CI95_low[1], bootstrap$CI95_high[1]))
cat(sprintf("- Significant CCL combinations: %d control, %d CS-NO\n",
            ccl_counts["Control"], ccl_counts["NO_treated"]))
