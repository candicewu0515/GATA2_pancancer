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
R 4.5.2, Bioconductor 3.22. See `scripts/00_install_packages.R`.

## Contact
xwu76@uiowa.edu
