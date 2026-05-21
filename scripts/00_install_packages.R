#!/usr/bin/env Rscript
# =============================================================
# 00_install_packages.R
# Run once to install all R packages used in this project.
#
# Usage:
#   Rscript 00_install_packages.R
#
# On Mac, expect 15-30 min the first time.
# =============================================================

# ---- mirror (faster from US) ---------------------------------
options(
    repos = c(CRAN = "https://cloud.r-project.org"),
    Ncpus = parallel::detectCores() - 1,
    timeout = 600
)

# ---- helper --------------------------------------------------
need_cran <- function(pkgs) {
    new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
    if (length(new)) install.packages(new)
}

need_bioc <- function(pkgs) {
    if (!"BiocManager" %in% installed.packages()[, "Package"])
        install.packages("BiocManager")
    new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
    if (length(new)) BiocManager::install(new, update = FALSE, ask = FALSE)
}

# ---- CRAN: core ----------------------------------------------
cran_core <- c(
    "data.table",    # fast big-file reading
    "tidyverse",     # dplyr, ggplot2, readr, tibble, etc.
    "ggpubr",        # stat_compare_means, theme_pubr
    "ggsci",         # NPG / AAAS / JCO color palettes
    "ggrepel",       # non-overlap labels
    "patchwork",     # combine plots
    "scales",        # axis scaling
    "RColorBrewer",
    "viridis",
    "cowplot"
)
need_cran(cran_core)

# ---- CRAN: stats / survival ----------------------------------
cran_stats <- c(
    "survival",
    "survminer",
    "forestplot",
    "rms",           # nomogram
    "Hmisc",
    "broom"
)
need_cran(cran_stats)

# ---- CRAN: heatmap / correlation -----------------------------
cran_viz <- c(
    "pheatmap",
    "corrplot",
    "circlize",
    "ggcorrplot"
)
need_cran(cran_viz)

# ---- CRAN: misc utilities ------------------------------------
cran_util <- c(
    "openxlsx",      # read/write xlsx
    "readxl",
    "writexl",
    "stringr",
    "glue",
    "future",
    "future.apply",
    "BiocManager",
    "remotes"
)
need_cran(cran_util)

# ---- Bioconductor --------------------------------------------
bioc_pkgs <- c(
    "ComplexHeatmap",     # publication heatmap
    "GSVA",               # GSVA scoring
    "GSEABase",
    "clusterProfiler",    # GSEA / GO / KEGG
    "enrichplot",
    "org.Hs.eg.db",
    "limma",              # DE for high-vs-low groups
    "edgeR",              # voom / TMM
    "DESeq2",             # if doing raw count DE
    "WGCNA",              # weighted co-expression
    "maftools",           # mutation waterfall plot
    "msigdbr",            # MSigDB direct in R
    "estimate",           # ESTIMATE immune score (sometimes archived; see below)
    "biomaRt"
)
need_bioc(bioc_pkgs)

# ---- ESTIMATE (R-Forge, may need manual install) -------------
if (!"estimate" %in% installed.packages()[, "Package"]) {
    tryCatch(
        install.packages("estimate", repos = "http://r-forge.r-project.org",
                         dependencies = TRUE),
        error = function(e) message("Install ESTIMATE manually from R-Forge if needed.")
    )
}

# ---- CIBERSORT (not on CRAN; signature matrix only) ----------
# CIBERSORT is a script file, not a package. Download from:
#   https://cibersortx.stanford.edu/  (requires login)
# Save CIBERSORT.R + LM22.txt into  scripts/external/

# ---- final check ---------------------------------------------
needed <- c(cran_core, cran_stats, cran_viz, cran_util, bioc_pkgs)
missing <- needed[!needed %in% installed.packages()[, "Package"]]

cat("\n========================================\n")
if (length(missing) == 0) {
    cat(" All packages installed.\n")
} else {
    cat(" Missing:", paste(missing, collapse = ", "), "\n")
}
cat("========================================\n")

cat("\nR version:", R.version.string, "\n")
cat("Bioconductor:",
    as.character(BiocManager::version()), "\n")
