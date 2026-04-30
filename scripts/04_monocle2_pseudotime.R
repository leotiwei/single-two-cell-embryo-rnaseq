#!/usr/bin/env Rscript

############################################################
# 04_monocle2_pseudotime.R
# Build a Seurat object from the raw-count matrix and perform Monocle2 pseudotime analysis.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(Matrix)
  library(ggplot2)
  library(viridis)
})

in_rdata <- "data/Step01_expression_data.RData"
outdir <- "results/monocle2"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(in_rdata)) stop("Input RData not found. Please run 01_prepare_expression_data.R first.")

load(in_rdata)

seu <- CreateSeuratObject(counts = filter_count, project = "two_cell_blastomere_embryo", min.cells = 0, min.features = 0)
seu$sample <- colnames(seu)
meta_match <- sample_info[match(seu$sample, sample_info$sample), ]
seu$stage <- meta_match$stage
seu$origin <- meta_match$origin
seu$pair_id <- meta_match$pair_id
seu$sister_blastomere <- meta_match$sister_blastomere

DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 4000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 10)
seu <- FindNeighbors(seu, dims = 1:5)
seu <- FindClusters(seu, resolution = 1.5)
seu <- RunUMAP(seu, dims = 1:5, n.neighbors = 8)
save(seu, file = file.path(outdir, "seurat_object.RData"))

p_umap <- DimPlot(seu, group.by = "stage") + ggtitle("UMAP colored by stage") + theme_bw(base_size = 12)
ggsave(file.path(outdir, "seurat_umap_stage.pdf"), p_umap, width = 5.5, height = 4.5)
ggsave(file.path(outdir, "seurat_umap_stage.png"), p_umap, width = 5.5, height = 4.5, dpi = 300)

expr <- as(as.matrix(GetAssayData(seu, slot = "counts")), "sparseMatrix")
pd <- seu@meta.data
pd$cell <- rownames(pd)
pd <- new("AnnotatedDataFrame", data = pd)
fd <- data.frame(gene_short_name = rownames(expr), row.names = rownames(expr), stringsAsFactors = FALSE)
fd <- new("AnnotatedDataFrame", data = fd)

cds <- newCellDataSet(expr, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
ordering_genes <- VariableFeatures(seu)
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, method = "DDRTree")
cds <- orderCells(cds)

root_state <- names(sort(table(pData(cds)$State[pData(cds)$stage == "8cell"]), decreasing = TRUE))[1]
if (!is.na(root_state) && length(root_state) == 1) {
  cds <- orderCells(cds, root_state = as.numeric(root_state))
}

pData(cds)$stage <- factor(pData(cds)$stage, levels = c("8cell", "2M", "morula", "2B", "blastocyst"))

theme_traj <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          legend.title = element_text(face = "bold"),
          legend.position = "right")
}

save_plot <- function(p, filename_no_ext, w = 5.5, h = 3.8, dpi = 300) {
  ggsave(paste0(filename_no_ext, ".pdf"), p, width = w, height = h)
  ggsave(paste0(filename_no_ext, ".png"), p, width = w, height = h, dpi = dpi)
}

p_state <- plot_cell_trajectory(cds, color_by = "State", cell_size = 3, show_branch_points = TRUE) +
  ggtitle("Trajectory colored by State") + theme_traj()

p_pseudotime <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 3, show_branch_points = TRUE) +
  scale_color_viridis_c(option = "D", end = 0.95) + ggtitle("Trajectory colored by Pseudotime") + theme_traj()

p_stage <- plot_cell_trajectory(cds, color_by = "stage", cell_size = 3, show_branch_points = TRUE) +
  ggtitle("Trajectory colored by Stage") + theme_traj()

save_plot(p_state, file.path(outdir, "monocle2_trajectory_state"))
save_plot(p_pseudotime, file.path(outdir, "monocle2_trajectory_pseudotime"))
save_plot(p_stage, file.path(outdir, "monocle2_trajectory_stage"))
save(cds, file = file.path(outdir, "monocle2_cds.RData"))
cat("Monocle2 trajectory results saved to:", outdir, "\n")
