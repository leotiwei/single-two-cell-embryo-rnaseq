#!/usr/bin/env Rscript

############################################################
# 03_pca_correlation_clustering.R
# Hierarchical clustering, PCA, and sample-correlation heatmap.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(FactoMineR)
  library(factoextra)
  library(pheatmap)
  library(ggplot2)
})

in_rdata <- "data/Step01_expression_data.RData"
outdir <- "results/pca_correlation_clustering"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(in_rdata)) stop("Input RData not found. Please run 01_prepare_expression_data.R first.")

load(in_rdata)
dat <- as.matrix(express_log2tpm)
stopifnot(identical(colnames(dat), sample_info$sample))

sample_tree_average <- hclust(dist(t(dat)), method = "average")
sample_tree_ward <- hclust(dist(t(dat)), method = "ward.D2")

pdf(file.path(outdir, "sample_tree_average.pdf"), width = 8, height = 5)
plot(sample_tree_average, hang = -1, xlab = "", ylab = "", main = "Sample clustering: average")
dev.off()
png(file.path(outdir, "sample_tree_average.png"), width = 2400, height = 1500, res = 300)
plot(sample_tree_average, hang = -1, xlab = "", ylab = "", main = "Sample clustering: average")
dev.off()

pdf(file.path(outdir, "sample_tree_wardD2.pdf"), width = 8, height = 5)
plot(sample_tree_ward, hang = -1, xlab = "", ylab = "", main = "Sample clustering: Ward.D2")
dev.off()
png(file.path(outdir, "sample_tree_wardD2.png"), width = 2400, height = 1500, res = 300)
plot(sample_tree_ward, hang = -1, xlab = "", ylab = "", main = "Sample clustering: Ward.D2")
dev.off()

pca_res <- PCA(as.data.frame(t(dat)), graph = FALSE)
p_pca <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  col.ind = sample_info$stage,
  addEllipses = TRUE,
  legend.title = "Stage"
) + theme_bw(base_size = 12)

ggsave(file.path(outdir, "pca_stage.pdf"), p_pca, width = 6, height = 5)
ggsave(file.path(outdir, "pca_stage.png"), p_pca, width = 6, height = 5, dpi = 300)

expr_mad <- dat[names(sort(apply(dat, 1, mad), decreasing = TRUE)), , drop = FALSE]
cor_mat <- cor(expr_mad, method = "pearson")

anno <- data.frame(stage = sample_info$stage, origin = sample_info$origin)
rownames(anno) <- sample_info$sample

pdf(file.path(outdir, "sample_correlation_heatmap.pdf"), width = 7, height = 6)
pheatmap(cor_mat, annotation_col = anno, annotation_row = anno, fontsize = 7,
         cluster_rows = FALSE, cluster_cols = FALSE,
         breaks = unique(seq(0.6, 1, length.out = 100)), main = "Sample correlation")
dev.off()

png(file.path(outdir, "sample_correlation_heatmap.png"), width = 2400, height = 2100, res = 300)
pheatmap(cor_mat, annotation_col = anno, annotation_row = anno, fontsize = 7,
         cluster_rows = FALSE, cluster_cols = FALSE,
         breaks = unique(seq(0.6, 1, length.out = 100)), main = "Sample correlation")
dev.off()

save(pca_res, cor_mat, sample_tree_average, sample_tree_ward, file = file.path(outdir, "pca_correlation_clustering_results.RData"))
cat("PCA, correlation and clustering results saved to:", outdir, "\n")
