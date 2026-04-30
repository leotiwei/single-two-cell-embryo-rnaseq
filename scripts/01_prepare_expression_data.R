#!/usr/bin/env Rscript

############################################################
# 01_prepare_expression_data.R
#
# Start from processed mouse embryo expression matrices.
# Read raw count and filtered TPM matrices, build sample metadata,
# calculate log2(TPM + 1), and save an RData object for downstream analyses.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

count_file <- "data/mouse_embryo_gene_rawcount.txt.gz"
tpm_file   <- "data/mouse_embryo_gene_filter_tpm.txt.gz"
out_rdata  <- "data/Step01_expression_data.RData"
outdir_qc  <- "results/qc"

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_qc, showWarnings = FALSE, recursive = TRUE)

sample_order <- c(
  paste0("8cell_", 1:4),
  "2M_1_1", "2M_1_2", "2M_2_1", "2M_2_2", "2M_3_1", "2M_3_2",
  paste0("morula_", 1:5),
  "2B_1_1", "2B_1_2", "2B_2_1", "2B_2_2", "2B_3_1", "2B_3_2",
  paste0("blastocyst_", 1:3)
)

read_expression_matrix <- function(file) {
  if (!file.exists(file)) stop("Input file not found: ", file)
  mat <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  if (!names(mat)[1] %in% sample_order) {
    rownames(mat) <- make.unique(as.character(mat[[1]]))
    mat <- mat[, -1, drop = FALSE]
  }
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat
}

filter_count <- read_expression_matrix(count_file)
filter_tpm   <- read_expression_matrix(tpm_file)

missing_count <- setdiff(sample_order, colnames(filter_count))
missing_tpm   <- setdiff(sample_order, colnames(filter_tpm))
if (length(missing_count) > 0) stop("Missing samples in count matrix: ", paste(missing_count, collapse = ", "))
if (length(missing_tpm) > 0) stop("Missing samples in TPM matrix: ", paste(missing_tpm, collapse = ", "))

filter_count <- filter_count[, sample_order, drop = FALSE]
filter_tpm   <- filter_tpm[, sample_order, drop = FALSE]

stopifnot(identical(colnames(filter_count), colnames(filter_tpm)))
stopifnot(identical(rownames(filter_count), rownames(filter_tpm)))

get_stage <- function(x) {
  if (grepl("^8cell_", x)) return("8cell")
  if (grepl("^2M_", x)) return("2M")
  if (grepl("^morula_", x)) return("morula")
  if (grepl("^2B_", x)) return("2B")
  if (grepl("^blastocyst_", x)) return("blastocyst")
  NA_character_
}
get_origin <- function(x) ifelse(grepl("^2M_|^2B_", x), "two-cell split", "intact control")
get_pair_id <- function(x) ifelse(grepl("^2M_|^2B_", x), paste0("pair ", strsplit(x, "_")[[1]][2]), "not applicable")
get_sister <- function(x) ifelse(grepl("^2M_|^2B_", x), paste0("blastomere ", strsplit(x, "_")[[1]][3]), "not applicable")

sample_info <- data.frame(
  sample = sample_order,
  stage = vapply(sample_order, get_stage, character(1)),
  origin = vapply(sample_order, get_origin, character(1)),
  pair_id = vapply(sample_order, get_pair_id, character(1)),
  sister_blastomere = vapply(sample_order, get_sister, character(1)),
  stringsAsFactors = FALSE
)
sample_info$stage <- factor(sample_info$stage, levels = c("8cell", "2M", "morula", "2B", "blastocyst"))

group_list <- as.character(sample_info$stage)
names(group_list) <- sample_info$sample

express_log2tpm <- log2(filter_tpm + 1)

gene_number <- colSums(filter_count > 1)
seq_depth   <- colSums(filter_count)

qc_df <- data.frame(
  sample = factor(sample_order, levels = sample_order),
  stage = sample_info$stage,
  detected_genes = as.numeric(gene_number[sample_order]),
  library_size = as.numeric(seq_depth[sample_order])
)

p_gene <- ggplot(qc_df, aes(x = sample, y = detected_genes, fill = stage)) +
  geom_col() + theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab(NULL) + ylab("Detected genes")

p_depth <- ggplot(qc_df, aes(x = sample, y = library_size, fill = stage)) +
  geom_col() + theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab(NULL) + ylab("Library size")

ggsave(file.path(outdir_qc, "detected_genes.pdf"), p_gene, width = 9, height = 4)
ggsave(file.path(outdir_qc, "detected_genes.png"), p_gene, width = 9, height = 4, dpi = 300)
ggsave(file.path(outdir_qc, "library_size.pdf"), p_depth, width = 9, height = 4)
ggsave(file.path(outdir_qc, "library_size.png"), p_depth, width = 9, height = 4, dpi = 300)

save(filter_count, filter_tpm, express_log2tpm, sample_info, group_list, gene_number, seq_depth, file = out_rdata)
cat("Saved:", out_rdata, "\n")
cat("Genes:", nrow(filter_count), "\n")
cat("Samples:", ncol(filter_count), "\n")
