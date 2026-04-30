#!/usr/bin/env Rscript

############################################################
# 02_sample_distribution.R
# Plot sample-wise expression distributions from log2(TPM + 1).
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({ library(ggplot2) })

in_rdata <- "data/Step01_expression_data.RData"
outdir <- "results/sample_distribution"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(in_rdata)) stop("Input RData not found. Please run 01_prepare_expression_data.R first.")

load(in_rdata)
exprSet <- as.matrix(express_log2tpm)

plot_df <- data.frame(
  expression = as.vector(exprSet),
  sample = rep(colnames(exprSet), each = nrow(exprSet)),
  stringsAsFactors = FALSE
)
plot_df$sample <- factor(plot_df$sample, levels = colnames(exprSet))
plot_df <- merge(plot_df, sample_info[, c("sample", "stage")], by = "sample", all.x = TRUE)
plot_df$sample <- factor(plot_df$sample, levels = colnames(exprSet))
plot_df$stage <- factor(plot_df$stage, levels = levels(sample_info$stage))

p_box <- ggplot(plot_df, aes(x = sample, y = expression, fill = stage)) +
  geom_boxplot(outlier.size = 0.2) + theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab(NULL) + ylab("log2(TPM + 1)")

ggsave(file.path(outdir, "sample_boxplot.pdf"), p_box, width = 10, height = 5)
ggsave(file.path(outdir, "sample_boxplot.png"), p_box, width = 10, height = 5, dpi = 300)

p_violin <- ggplot(plot_df, aes(x = sample, y = expression, fill = stage)) +
  geom_violin(scale = "width") + theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab(NULL) + ylab("log2(TPM + 1)")

ggsave(file.path(outdir, "sample_violin.pdf"), p_violin, width = 10, height = 5)
ggsave(file.path(outdir, "sample_violin.png"), p_violin, width = 10, height = 5, dpi = 300)

p_density <- ggplot(plot_df, aes(x = expression, color = sample, group = sample)) +
  geom_density(linewidth = 0.4, alpha = 0.5) + theme_bw(base_size = 12) +
  xlab("log2(TPM + 1)") + ylab("Density")

ggsave(file.path(outdir, "sample_density.pdf"), p_density, width = 9, height = 5)
ggsave(file.path(outdir, "sample_density.png"), p_density, width = 9, height = 5, dpi = 300)
cat("Sample-distribution plots saved to:", outdir, "\n")
