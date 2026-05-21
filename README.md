# GATA2 Pan-Cancer Multi-Omics Analysis

R/bash pipeline for "Pan-cancer multi-omics analysis of GATA2 reveals a stromal-immune versus proliferative dichotomy across 33 tumour types" (Wu, 2026).

## Pipeline
Scripts in `scripts/` run sequentially (00–14):
- 00–04: setup, data download, QC, sample cleaning, GATA2 extraction
- 05–14: Figures 1–10 and multi-omics validation

## Data sources
- TCGA Pan-Cancer Atlas: UCSC Xena (https://xena.ucsc.edu/)
- CCLE: DepMap 26Q1 (https://depmap.org/)
- NCI-60: CellMiner v2.10
- CPTAC proteomics: UALCAN
- MSigDB Hallmark v2024.1.Hs

## Requirements

Tested on macOS with **R 4.5.2 (2025-10-31)** and **Bioconductor 3.22**.
Run `Rscript scripts/00_install_packages.R` to install everything.

### R packages (versions used)

CRAN — core:
data.table 1.18.4, tidyverse 2.0.0, ggpubr 0.6.3, ggsci 5.0.0, ggrepel 0.9.8,
patchwork 1.3.2, scales 1.4.0, RColorBrewer 1.1.3, viridis 0.6.5, cowplot 1.2.0

CRAN — stats / survival:
survival 3.8.3, survminer 0.5.2, forestplot 3.2.0, rms 8.1.1, Hmisc 5.2.5, broom 1.0.13

CRAN — heatmap / correlation:
pheatmap 1.0.13, corrplot 0.95, circlize 0.4.18, ggcorrplot 0.1.4.1

CRAN — utilities:
openxlsx 4.2.8.1, readxl 1.5.0, writexl 1.5.4, stringr 1.6.0, glue 1.8.0,
future 1.70.0, future.apply 1.20.2, BiocManager 1.30.27, remotes 2.5.0

Bioconductor 3.22:
ComplexHeatmap 2.26.1, GSVA 2.4.9, GSEABase 1.72.0, clusterProfiler 4.18.4,
enrichplot 1.30.5, org.Hs.eg.db 3.22.0, limma 3.66.0, edgeR 4.8.2,
DESeq2 1.50.2, WGCNA 1.74, maftools 2.26.0, msigdbr 26.1.0, biomaRt 2.66.2

R-Forge:
estimate 1.0.13

External (not an R package):
CIBERSORT.R + LM22 signature matrix — download from https://cibersortx.stanford.edu/
(requires login) and place under `scripts/external/`.

## Contact
xwu76@uiowa.edu
