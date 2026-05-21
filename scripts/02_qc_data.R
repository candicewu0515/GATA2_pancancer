#!/usr/bin/env Rscript
# =============================================================
# 02_qc_data.R
# Sanity check downloaded files
#   - file size
#   - first 5 rows / columns
#   - sample count
#   - GATA2 row present?
#
# Usage:
#   Rscript 02_qc_data.R
#
# HPC:
#   qrsh -q UI -pe smp 2 -l h_vmem=16G -l h_rt=2:00:00
#   module load R
#   Rscript 02_qc_data.R
# =============================================================

suppressPackageStartupMessages({
    library(data.table)
})

DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
RAW <- file.path(DATA_ROOT, "01_raw")

cat("\n=========================================\n")
cat(" QC report - downloaded raw data\n")
cat(" Data root:", normalizePath(RAW, mustWork = FALSE), "\n")
cat("=========================================\n\n")

# helper -------------------------------------------------------
check_file <- function(path, peek_lines = 3) {
    if (!file.exists(path)) {
        cat(sprintf("  [MISSING] %s\n", path))
        return(invisible(FALSE))
    }
    sz <- file.info(path)$size
    cat(sprintf("  [OK %.1f MB] %s\n", sz / 1e6, basename(path)))

    # binary types: just report size, don't try to peek
    if (grepl("\\.(zip|xls|xlsx|bam|bai|rda|rds|RData)$", path, ignore.case = TRUE)) {
        return(invisible(TRUE))
    }

    # peek first N lines (handles .gz). Wrap in tryCatch so a single bad file
    # doesn't halt the whole QC report.
    tryCatch({
        con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
        on.exit(close(con))
        lines <- readLines(con, n = peek_lines, warn = FALSE)
        for (l in lines) cat("    ", substr(l, 1, 150), "\n")
    }, error = function(e) {
        cat("    (cannot peek as text:", conditionMessage(e), ")\n")
    })
    invisible(TRUE)
}

check_matrix_for_gata2 <- function(path, sep = "\t") {
    if (!file.exists(path)) return(invisible(NULL))
    cat(sprintf("\n--- GATA2 presence in %s ---\n", basename(path)))
    # only read first column to find GATA2 row
    dt <- tryCatch(
        fread(path, select = 1, header = TRUE, sep = sep),
        error = function(e) { cat("  cannot parse:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(dt)) return(invisible(NULL))
    first_col <- dt[[1]]
    hits <- grep("^GATA2$|\\|GATA2$|GATA2\\|", first_col, value = TRUE)
    if (length(hits) > 0) {
        cat("  GATA2 found:", paste(head(hits, 3), collapse = ", "), "\n")
    } else {
        # try ENSG lookup
        ensg_hits <- grep("ENSG00000179348", first_col, value = TRUE)  # GATA2 Ensembl ID
        if (length(ensg_hits) > 0) {
            cat("  GATA2 found (Ensembl):", paste(head(ensg_hits, 3), collapse = ", "), "\n")
        } else {
            cat("  WARNING: GATA2 not found in first column. Maybe transposed?\n")
        }
    }
}

# ---- 1.1 PANCAN ---------------------------------------------
cat("[1.1] TCGA PANCAN\n")
d <- file.path(RAW, "xena_pancan")
files_pancan <- c(
    "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz",
    "TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
    "Survival_SupplementalTable_S1_20171025_xena_sp",
    "mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz",
    "jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz",
    "StemnessScores_RNAexp_20170127.2.tsv.gz"
)
for (f in files_pancan) check_file(file.path(d, f))
check_matrix_for_gata2(file.path(d, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"))

# ---- 1.2 Toil -----------------------------------------------
cat("\n[1.2] TCGA + GTEx Toil\n")
d <- file.path(RAW, "xena_toil")
files_toil <- c(
    "TcgaTargetGtex_gene_expected_count.gz",
    "TcgaTargetGTEX_phenotype.txt.gz",
    "gencode.v23.annotation.gene.probemap"
)
for (f in files_toil) check_file(file.path(d, f))
check_matrix_for_gata2(file.path(d, "TcgaTargetGtex_gene_expected_count.gz"))

# ---- 1.3 CCLE -----------------------------------------------
cat("\n[1.3] CCLE (manual download)\n")
d <- file.path(RAW, "ccle")
# DepMap renamed the expression matrix; check whichever variant exists.
# New (24Q4+): OmicsExpressionTPMLogp1HumanAllGenes.csv  (all human genes)
# Old:         OmicsExpressionProteinCodingGenesTPMLogp1.csv  (coding only)
ccle_expr_candidates <- c(
    "OmicsExpressionTPMLogp1HumanAllGenes.csv",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
)
ccle_expr_found <- ccle_expr_candidates[file.exists(file.path(d, ccle_expr_candidates))]
if (length(ccle_expr_found) > 0) {
    for (f in ccle_expr_found) check_file(file.path(d, f))
} else {
    cat(sprintf("  [MISSING] %s (or older %s)\n",
                ccle_expr_candidates[1], ccle_expr_candidates[2]))
}
check_file(file.path(d, "Model.csv"))

# ---- 1.4 CellMiner ------------------------------------------
cat("\n[1.4] CellMiner\n")
d <- file.path(RAW, "cellminer")
cm_files <- list.files(d, pattern = "\\.(zip|xls|xlsx|tsv|txt)$", full.names = FALSE)
for (f in cm_files) check_file(file.path(d, f))

# ---- 1.5 MSigDB ---------------------------------------------
cat("\n[1.5] MSigDB\n")
d <- file.path(RAW, "msigdb")
msig_files <- list.files(d, pattern = "\\.gmt$", full.names = FALSE)
for (f in msig_files) {
    check_file(file.path(d, f), peek_lines = 1)
}

# ---- 1.6 MC3 mutation MAF (optional but useful) -------------
cat("\n[1.6] MC3 mutation MAF\n")
d <- file.path(RAW, "gdc_mc3")
check_file(file.path(d, "mc3.v0.2.8.PUBLIC.maf.gz"))

# ---- summary ------------------------------------------------
cat("\n=========================================\n")
cat(" QC done.  Next: 03_clean_samples.R\n")
cat("=========================================\n\n")
