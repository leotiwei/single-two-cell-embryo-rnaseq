# Workflow overview

This repository starts from processed expression matrices rather than raw FASTQ files.

## Input

Download the processed matrices from GEO and place them in a local `data/` folder:

```text
data/mouse_embryo_gene_rawcount.txt.gz
data/mouse_embryo_gene_filter_tpm.txt.gz
```

## Analysis steps

1. `01_prepare_expression_data.R`
   - reads the raw-count and TPM matrices
   - checks sample names and order
   - builds sample metadata
   - calculates `log2(TPM + 1)`
   - saves `data/Step01_expression_data.RData`

2. `02_sample_distribution.R`
   - plots sample-wise expression distributions

3. `03_pca_correlation_clustering.R`
   - performs sample clustering
   - performs PCA
   - computes Pearson sample correlations
   - generates the correlation heatmap

4. `04_monocle2_pseudotime.R`
   - creates a Seurat object from the count matrix
   - identifies highly variable genes
   - constructs a Monocle2 trajectory
   - plots trajectory by State, Pseudotime and Stage
