#!/usr/bin/env Rscript
# =============================================================
# 13_fig9_estimate_tme.R
# Figure 9: ESTIMATE scores + TME signatures vs GATA2 across pan-cancer.
#   Panel A: ESTIMATE Immune/Stromal/ESTIMATE/Purity score correlations
#            (heatmap: cancer x score, Pearson with GATA2)
#   Panel B: TME-relevant signatures (CD8 effector, antigen processing,
#            MMR, EMT, etc.) -- ssGSEA scores vs GATA2
#
# Inputs:
#   data/01_raw/xena_pancan/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
#   data/02_clean/GATA2_expression.rds
#
# Outputs:
#   data/02_clean/estimate_scores.rds
#   data/02_clean/tme_ssgsea.rds
#   data/03_results/Fig9_estimate_correlations.csv
#   data/03_results/Fig9_tme_correlations.csv
#   figures/Fig9A_estimate_heatmap.pdf/.png
#   figures/Fig9B_tme_signatures_heatmap.pdf/.png
#
# Strategy:
#   - Use the `estimate` R package on the PANCAN matrix per cancer.
#   - For TME signatures, use a 14-signature panel (Mariathasan 2018 style):
#       CD8 T effector, Antigen processing, Immune checkpoint,
#       Mismatch repair, Nucleotide excision repair, DNA damage response,
#       DNA replication, Base excision repair, Pan-fibroblast TGFb,
#       EMT1, EMT2, EMT3, TMEscore_A, TMEscore_B
#     Compute via GSVA ssGSEA.
#
# Usage:
#   Rscript scripts/13_fig9_estimate_tme.R
# =============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(estimate)
    library(GSVA)
    library(ggplot2)
    library(ggpubr)
    library(pheatmap)
    library(RColorBrewer)
})

# ---- paths ---------------------------------------------------
DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
RAW       <- file.path(DATA_ROOT, "01_raw", "xena_pancan")
CLEAN     <- file.path(DATA_ROOT, "02_clean")
RESULTS   <- file.path(DATA_ROOT, "03_results")
FIGDIR    <- "./figures"
for (p in c(CLEAN, RESULTS, FIGDIR)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

EXPR_FILE <- file.path(RAW, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")

# ---- 1. load ------------------------------------------------
cat("[1/6] Loading metadata...\n")
g <- readRDS(file.path(CLEAN, "GATA2_expression.rds"))
g <- g[sample_group == "Tumor" & !is.na(cancer_type)]

cat("[2/6] Reading expression matrix (slow)...\n")
expr <- fread(EXPR_FILE)
setnames(expr, 1, "gene")
sample_cols <- intersect(g$sample, colnames(expr))
m <- as.matrix(expr[, c("gene", sample_cols), with = FALSE], rownames = "gene")
rm(expr); gc()
if (anyDuplicated(rownames(m))) m <- limma::avereps(m)
cat("  matrix:", nrow(m), "genes x", ncol(m), "samples\n")

# ---- 2. ESTIMATE -------------------------------------------
est_file <- file.path(CLEAN, "estimate_scores.rds")
if (file.exists(est_file)) {
    cat("[3/6] Loading cached ESTIMATE scores...\n")
    estscore <- readRDS(est_file)
} else {
    cat("[3/6] Running ESTIMATE...\n")
    tmp_in  <- tempfile(fileext = ".gct")
    tmp_filt <- tempfile(fileext = ".gct")
    tmp_out <- tempfile(fileext = ".gct")

    # Despite the .gct extension, estimate::filterCommonGenes reads the
    # input with read.table(header=TRUE, row.names=1) — it does NOT skip
    # the #1.2 / <nrow>\t<ncol> GCT header lines.  So write a plain TSV
    # (NAME\tDescription\t<sample>...) with no GCT preamble.
    expr_df <- data.frame(NAME = rownames(m), Description = rownames(m), m,
                          check.names = FALSE)
    write.table(
        expr_df, file = tmp_in, sep = "\t",
        quote = FALSE, row.names = FALSE, col.names = TRUE
    )

    filterCommonGenes(input.f = tmp_in, output.f = tmp_filt, id = "GeneSymbol")
    estimateScore(input.ds = tmp_filt, output.ds = tmp_out, platform = "illumina")

    raw <- read.table(tmp_out, skip = 2, header = TRUE, sep = "\t",
                      check.names = FALSE)
    # rows: StromalScore / ImmuneScore / ESTIMATEScore
    rownames(raw) <- raw$NAME
    raw <- raw[, -(1:2)]
    # restore sample IDs (estimate adds "X" prefix to names starting with digits)
    colnames(raw) <- sample_cols[1:ncol(raw)]
    estscore <- as.data.table(t(raw), keep.rownames = "sample")
    saveRDS(estscore, est_file)
}
cat("  ESTIMATE rows:", nrow(estscore), "\n")

# ---- 3. TME signatures (Mariathasan 2018 panel) -------------
TME_SIGS <- list(
    CD_8_T_effector = c("CD8A","GZMA","GZMB","IFNG","CXCL9","CXCL10","PRF1","TBX21"),
    Antigen_processing_machinery = c("HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9","B2M"),
    Immune_Checkpoint = c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT"),
    Mismatch_Repair   = c("MLH1","MSH2","MSH6","PMS2","EPCAM"),
    Nucleotide_excision_repair = c("ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","XPA","XPC"),
    DNA_damage_response = c("ATM","ATR","BRCA1","BRCA2","CHEK1","CHEK2","MDC1","TP53BP1"),
    DNA_replication = c("MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","CDC6","CDC45"),
    Base_excision_repair = c("APEX1","OGG1","XRCC1","POLB","FEN1","LIG3","MUTYH","NEIL1"),
    Pan_F_TBRs = c("ACTA2","ACTG2","ADAM12","COL1A1","COL1A2","COL3A1","COL5A1","COL5A2",
                   "FAP","FN1","MMP2","MMP11","PDGFRA","PDGFRB","TGFB1","TGFB2","TGFB3","VIM"),
    EMT1 = c("VIM","CDH2","FN1","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2"),
    EMT2 = c("CDH1","KRT8","KRT18","KRT19","CLDN1","CLDN4","OCLN"),
    EMT3 = c("MMP2","MMP9","MMP14","CTSK","PLAU","PLAUR","TIMP1"),
    TMEscore_A = c("CD8A","GZMA","GZMB","IFNG","PRF1","CD274","PDCD1","CTLA4","TBX21","CXCL9","CXCL10"),
    TMEscore_B = c("ACTA2","COL1A1","COL3A1","FAP","FN1","MMP2","PDGFRA","TGFB1","VIM")
)

ssg_file <- file.path(CLEAN, "tme_ssgsea.rds")
if (file.exists(ssg_file)) {
    cat("[4/6] Loading cached ssGSEA TME scores...\n")
    tme_scores <- readRDS(ssg_file)
} else {
    cat("[4/6] Computing ssGSEA TME scores...\n")
    # GSVA 2.x uses parameter object; fall back for older GSVA
    tme_scores <- tryCatch({
        par <- ssgseaParam(exprData = m, geneSets = TME_SIGS,
                           minSize = 2, maxSize = 500)
        gsva(par, verbose = FALSE)
    }, error = function(e) {
        # legacy GSVA <= 1.50
        gsva(m, TME_SIGS, method = "ssgsea", verbose = FALSE,
             min.sz = 2, max.sz = 500)
    })
    saveRDS(tme_scores, ssg_file)
}
cat("  TME signatures x samples:", paste(dim(tme_scores), collapse = " x "), "\n")

# ---- 4. correlate ESTIMATE / TME with GATA2 per cancer -----
cat("[5/6] Computing per-cancer correlations...\n")

# merge ESTIMATE
em <- merge(g[, .(sample, cancer_type, GATA2)], estscore, by = "sample")
em[, TumorPurity := cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)]

est_long <- melt(em,
    id.vars = c("sample", "cancer_type", "GATA2"),
    measure.vars = c("ImmuneScore", "StromalScore", "ESTIMATEScore", "TumorPurity"),
    variable.name = "score", value.name = "value")

est_cor <- est_long[, {
    if (sum(!is.na(value)) < 10) {
        .(rho = NA_real_, p = NA_real_, n = sum(!is.na(value)))
    } else {
        ct <- suppressWarnings(cor.test(GATA2, value, method = "pearson"))
        .(rho = as.numeric(ct$estimate),
          p   = as.numeric(ct$p.value),
          n   = sum(!is.na(value)))
    }
}, by = .(cancer_type, score)]
est_cor[, sig := fcase(p < 0.001, "***", p < 0.01, "**", p < 0.05, "*", default = "")]
fwrite(est_cor, file.path(RESULTS, "Fig9_estimate_correlations.csv"))

# merge TME ssGSEA
tme_dt <- as.data.table(t(tme_scores), keep.rownames = "sample")
tm <- merge(g[, .(sample, cancer_type, GATA2)], tme_dt, by = "sample")
tme_long <- melt(tm, id.vars = c("sample", "cancer_type", "GATA2"),
                 variable.name = "signature", value.name = "value")
tme_cor <- tme_long[, {
    if (sum(!is.na(value)) < 10) {
        .(rho = NA_real_, p = NA_real_)
    } else {
        ct <- suppressWarnings(cor.test(GATA2, value, method = "pearson"))
        .(rho = as.numeric(ct$estimate),
          p   = as.numeric(ct$p.value))
    }
}, by = .(cancer_type, signature)]
tme_cor[, sig := fcase(p < 0.001, "***", p < 0.01, "**", p < 0.05, "*", default = "")]
fwrite(tme_cor, file.path(RESULTS, "Fig9_tme_correlations.csv"))

# ---- 5. plot Fig 9A: ESTIMATE heatmap ----------------------
cat("[6/6] Plotting...\n")
make_heatmap <- function(cor_dt, file_base, title, h_factor = 0.5) {
    mat <- dcast(cor_dt, signature ~ cancer_type, value.var = "rho",
                 fill = NA_real_)
    sig <- dcast(cor_dt, signature ~ cancer_type, value.var = "sig",
                 fill = "")
    rn <- mat$signature
    mat <- as.matrix(mat[, -1])
    sig <- as.matrix(sig[, -1])
    rownames(mat) <- rn; rownames(sig) <- rn
    sig[is.na(sig)] <- ""

    pal <- colorRampPalette(c("#3B4992", "white", "#EE0000"))(100)
    pheatmap(
        mat, color = pal,
        breaks = seq(-1, 1, length.out = 101),
        cluster_rows = TRUE, cluster_cols = TRUE,
        display_numbers = sig, number_color = "black",
        fontsize_number = 9, fontsize_row = 9, fontsize_col = 9,
        cellwidth = 16, cellheight = 18,
        border_color = "grey90",
        main = title,
        filename = paste0(file_base, ".pdf"),
        width = 12, height = max(4, nrow(mat) * h_factor + 2)
    )
    pheatmap(
        mat, color = pal,
        breaks = seq(-1, 1, length.out = 101),
        cluster_rows = TRUE, cluster_cols = TRUE,
        display_numbers = sig, number_color = "black",
        fontsize_number = 9, fontsize_row = 9, fontsize_col = 9,
        cellwidth = 16, cellheight = 18,
        border_color = "grey90",
        main = title,
        filename = paste0(file_base, ".png"),
        width = 12, height = max(4, nrow(mat) * h_factor + 2)
    )
}

# Fig 9A
est_cor2 <- copy(est_cor)
setnames(est_cor2, "score", "signature")
make_heatmap(est_cor2,
             file.path(FIGDIR, "Fig9A_estimate_heatmap"),
             "GATA2 vs ESTIMATE / TumorPurity (Pearson r)",
             h_factor = 0.6)

# Fig 9B
make_heatmap(tme_cor,
             file.path(FIGDIR, "Fig9B_tme_signatures_heatmap"),
             "GATA2 vs TME signatures (Pearson r)",
             h_factor = 0.4)

cat("\n[saved]\n")
cat("  ", file.path(FIGDIR, "Fig9A_estimate_heatmap.pdf"), "\n")
cat("  ", file.path(FIGDIR, "Fig9B_tme_signatures_heatmap.pdf"), "\n")
cat("  ", file.path(RESULTS, "Fig9_estimate_correlations.csv"), "\n")
cat("  ", file.path(RESULTS, "Fig9_tme_correlations.csv"), "\n")
