#!/usr/bin/env Rscript


source("scripts/_common.R")
base_dir <- project_root()
dirs <- make_output_dirs(base_dir)

metadata_path <- file.path(base_dir, "data", "processed", "cell_metadata.csv.gz")
assert_file(metadata_path)

metadata <- read.csv(gzfile(metadata_path), row.names = 1, check.names = FALSE)
required_columns <- c("condition", "cell_type")
if (!all(required_columns %in% colnames(metadata))) {
  stop("Metadata must contain: ", paste(required_columns, collapse = ", "))
}

conditions <- c("Control", "NO_treated")
if (!all(conditions %in% metadata$condition)) {
  stop("Expected condition labels: ", paste(conditions, collapse = ", "))
}

counts <- table(metadata$cell_type, metadata$condition)
composition <- data.frame(
  cell_type = rownames(counts),
  n_Ctrl = counts[, "Control"],
  n_NO = counts[, "NO_treated"],
  pct_Ctrl = 100 * counts[, "Control"] / sum(counts[, "Control"]),
  pct_NO = 100 * counts[, "NO_treated"] / sum(counts[, "NO_treated"]),
  row.names = NULL,
  check.names = FALSE
)

fisher_results <- lapply(seq_len(nrow(composition)), function(i) {
  in_type_control <- composition$n_Ctrl[i]
  in_type_treated <- composition$n_NO[i]
  contingency <- matrix(
    c(
      in_type_control,
      sum(composition$n_Ctrl) - in_type_control,
      in_type_treated,
      sum(composition$n_NO) - in_type_treated
    ),
    nrow = 2,
    byrow = TRUE
  )

  test <- fisher.test(contingency)
  data.frame(
    cell_type = composition$cell_type[i],
    log2_fraction_ratio = log2(
      (in_type_treated / sum(composition$n_NO)) /
        (in_type_control / sum(composition$n_Ctrl))
    ),
    odds_ratio = unname(test$estimate),
    p_value = test$p.value
  )
})

fisher_results <- do.call(rbind, fisher_results)
fisher_results$FDR <- p.adjust(fisher_results$p_value, method = "BH")
fisher_results$significance <- cut(
  fisher_results$FDR,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns")
)

write.csv(
  composition,
  file.path(dirs$results, "celltype_composition.csv"),
  row.names = FALSE
)
write.csv(
  fisher_results,
  file.path(dirs$results, "compositional_fisher_bh.csv"),
  row.names = FALSE
)

cell_types <- split(metadata$cell_type, metadata$condition)
observed_control <- mean(cell_types$Control == "PTC-2")
observed_treated <- mean(cell_types$NO_treated == "PTC-2")
observed_ratio <- observed_treated / observed_control
observed_reduction <- 100 * (1 - observed_ratio)

set.seed(42)
n_boot <- 2000L
bootstrap_values <- replicate(n_boot, {
  control_fraction <- mean(
    sample(cell_types$Control, length(cell_types$Control), replace = TRUE) == "PTC-2"
  )
  treated_fraction <- mean(
    sample(cell_types$NO_treated, length(cell_types$NO_treated), replace = TRUE) == "PTC-2"
  )
  ratio <- treated_fraction / control_fraction
  c(reduction = 100 * (1 - ratio), ratio = ratio)
})

bootstrap_summary <- rbind(
  data.frame(
    metric = "PTC-2 reduction (%)",
    observed = observed_reduction,
    boot_median = median(bootstrap_values["reduction", ]),
    CI95_low = quantile(bootstrap_values["reduction", ], 0.025),
    CI95_high = quantile(bootstrap_values["reduction", ], 0.975)
  ),
  data.frame(
    metric = "PTC-2 proportion ratio NO/Control",
    observed = observed_ratio,
    boot_median = median(bootstrap_values["ratio", ]),
    CI95_low = quantile(bootstrap_values["ratio", ], 0.025),
    CI95_high = quantile(bootstrap_values["ratio", ], 0.975)
  )
)

write.csv(
  bootstrap_summary,
  file.path(dirs$results, "nr_ptc_bootstrap.csv"),
  row.names = FALSE
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_data <- data.frame(reduction = bootstrap_values["reduction", ])
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = reduction)) +
    ggplot2::geom_histogram(bins = 45, fill = "#4A7BB0", color = "white") +
    ggplot2::geom_vline(xintercept = observed_reduction, color = "#B2182B", linewidth = 0.8) +
    ggplot2::labs(
      x = "Reduction in recovered PTC-2 fraction (%)",
      y = "Bootstrap resamples",
      title = "Cell-level bootstrap sensitivity analysis"
    ) +
    ggplot2::theme_classic(base_size = 12)

  ggplot2::ggsave(
    file.path(dirs$figures, "nr_ptc_bootstrap.pdf"),
    p,
    width = 6.5,
    height = 4.5
  )
}

cat(sprintf("PTC-2 control fraction: %.4f\n", observed_control))
cat(sprintf("PTC-2 CS-NO fraction: %.4f\n", observed_treated))
cat(sprintf("Relative reduction: %.2f%%\n", observed_reduction))
cat(sprintf(
  "Bootstrap 95%% interval: %.2f%% to %.2f%%\n",
  bootstrap_summary$CI95_low[1],
  bootstrap_summary$CI95_high[1]
))
