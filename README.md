# Two-cell blastomere RNA-seq analysis

This repository contains R scripts used to analyze stage-resolved Smart-seq2 RNA-seq data from mouse preimplantation embryos after late two-cell blastomere separation.

The repository starts from processed expression matrices, not from raw FASTQ files. Raw sequencing data and processed matrices should be downloaded from GEO.

## Required input files

Create a local `data/` directory and place the following files inside it:

```text
data/
  mouse_embryo_gene_rawcount.txt.gz
  mouse_embryo_gene_filter_tpm.txt.gz
```

Both matrices should contain genes as rows and samples as columns. The scripts expect the following sample column names:

```text
8cell_1 8cell_2 8cell_3 8cell_4
2M_1_1 2M_1_2 2M_2_1 2M_2_2 2M_3_1 2M_3_2
morula_1 morula_2 morula_3 morula_4 morula_5
2B_1_1 2B_1_2 2B_2_1 2B_2_2 2B_3_1 2B_3_2
blastocyst_1 blastocyst_2 blastocyst_3
```

The first column may contain gene IDs or gene symbols. If extra gene annotation columns are present, the scripts will keep only the expected sample columns for analysis.

## How to run

Run the scripts in order:

```bash
Rscript scripts/01_prepare_expression_data.R
Rscript scripts/02_sample_distribution.R
Rscript scripts/03_pca_correlation_clustering.R
Rscript scripts/04_monocle2_pseudotime.R
```

## Output

Intermediate R objects are saved in `data/`. Figures and result objects are saved in `results/`.

## Required R packages

```text
ggplot2
FactoMineR
factoextra
pheatmap
Seurat
monocle
Matrix
viridis
```

## Data availability

Raw sequencing data and processed expression matrices are deposited in GEO under accession number `GSEXXXXXX`.

## Code availability

No new algorithm was developed in this study. These scripts reproduce expression-matrix processing, sample-distribution plots, PCA, sample-correlation analysis, hierarchical clustering and Monocle2 pseudotime analysis used in the study.

## License

This code is released under the MIT License.
