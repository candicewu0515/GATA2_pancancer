#!/usr/bin/env Rscript
# =============================================================
# 08_fig4_cibersort.R
# Figure 4: Immune cell infiltration (CIBERSORT) vs GATA2 expression
# across pan-cancer.
#
# Strategy: Instead of running CIBERSORT locally (slow, needs LM22 +
# the closed-source CIBERSORT.R), we use the precomputed TIMER2.0
# pan-cancer immune infiltration estimates that include CIBERSORT,
# CIBERSORT-ABS, EPIC, MCP-counter, quanTIseq, TIMER, and xCell.
#
# Input:
#   data/02_clean/GATA2_expression.rds      (from 04)
#   data/01_raw/timer2/infiltration_estimation_for_tcga.csv
#
# If the TIMER2 CSV doesn't exist, the script downloads it
# (~100 MB) on first run.
#
# Output:
#   figures/Fig4_cibersort_heatmap.pdf / .png
#   figures/Fig4_cibersort_focus_boxplots.pdf  (high vs low for focus cancers)
#   data/03_results/Fig4_cibersort_correlations.csv
#
# Usage:
#   Rscript scripts/08_fig4_cibersort.R                 # default
#   Rscript scripts/08_fig4_cibersort.R PRAD COAD KIRC  # focus cancers
# =============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
    library(pheatmap)
    library(RColorBrewer)
})

# ---- args ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
FOCUS <- if (length(args) >= 1) args else c("PRAD", "COAD", "KIRC")
cat("Focus cancers:", paste(FOCUS, collapse = ", "), "\n\n")

# ---- paths ---------------------------------------------------
DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
CLEAN     <- file.path(DATA_ROOT, "02_clean")
RESULTS   <- file.path(DATA_ROOT, "03_results")
TIMER_DIR <- file.path(DATA_ROOT, "01_raw", "timer2")
FIGDIR    <- "./figures"
dir.create(RESULTS,   showWarnings = FALSE, recursive = TRUE)
dir.create(TIMER_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR,    showWarnings = FALSE, recursive = TRUE)

# ---- 1. ensure TIMER2 file --------------------------------
TIMER_FILE <- file.path(TIMER_DIR, "infiltration_estimation_for_tcga.csv")
if (!file.exists(TIMER_FILE)) {
    # TIMER2 moved domains: old timer.comp-genomics.org now 301-redirects to
    # compbio.cn/timer3, but the /static/data/ path is gone — the file is
    # served from the timer3 site root.
    cat("[0/5] Downloading TIMER2.0 precomputed estimates (~9 MB gz)...\n")
    url <- "https://compbio.cn/timer3/infiltration_estimation_for_tcga.csv.gz"
    tmp_gz <- paste0(TIMER_FILE, ".gz")
    options(timeout = 1800)
    download.file(url, tmp_gz, mode = "wb")
    cat("  decompressing...\n")
    if (!requireNamespace("R.utils", quietly = TRUE)) {
        install.packages("R.utils", repos = "https://cloud.r-project.org")
    }
    R.utils::gunzip(tmp_gz, destname = TIMER_FILE, remove = TRUE)
}

# ---- 2. load expression + immune ---------------------------
cat("[1/5] Loading data...\n")
expr <- readRDS(file.path(CLEAN, "GATA2_expression.rds"))
expr <- expr[sample_group == "Tumor" & !is.na(cancer_type)]

imm <- fread(TIMER_FILE)
setnames(imm, 1, "sample")
cat("  immune matrix:", nrow(imm), "samples x", ncol(imm) - 1, "features\n")
cat("  expression  :", nrow(expr), "tumor samples\n")

# Match sample barcodes (TIMER2 uses TCGA-XX-XXXX-01 format too,
# but some have trailing letters; trim to 15 char primary tumor ID)
imm[,  sample15 := substr(sample, 1, 15)]
expr[, sample15 := substr(sample, 1, 15)]
merged <- merge(expr[, .(sample15, cancer_type, GATA2)],
                imm, by = "sample15")
cat("  merged    :", nrow(merged), "samples\n")

# Keep only CIBERSORT columns (those ending in "_CIBERSORT", not ABS)
cib_cols <- grep("_CIBERSORT$", colnames(merged), value = TRUE)
cat("  CIBERSORT cell types:", length(cib_cols), "\n")
stopifnot(length(cib_cols) > 0)

# ---- 3. correlation per cancer per cell type ---------------
cat("[2/5] Computing Spearman correlations per cancer x cell type...\n")
cancers <- sort(unique(merged$cancer_type))

cor_dt <- rbindlist(lapply(cancers, function(cc) {
    d <- merged[cancer_type == cc]
    if (nrow(d) < 20) return(NULL)
    rbindlist(lapply(cib_cols, function(cell) {
        v <- d[[cell]]
        if (sum(!is.na(v)) < 10 || sd(v, na.rm = TRUE) == 0) return(NULL)
        ct <- suppressWarnings(cor.test(d$GATA2, v, method = "spearman"))
        data.table(cancer = cc, cell = cell, rho = ct$estimate, p = ct$p.value, n = nrow(d))
    }))
}))

cor_dt[, sig := fcase(
    p < 0.001, "***",
    p < 0.01,  "**",
    p < 0.05,  "*",
    default    = ""
)]
# clean cell labels (drop _CIBERSORT)
cor_dt[, cell_label := sub("_CIBERSORT$", "", cell)]
cor_dt[, cell_label := gsub("\\.", " ", cell_label)]
fwrite(cor_dt, file.path(RESULTS, "Fig4_cibersort_correlations.csv"))

# ---- 4. heatmap (rho + sig stars) ---------------------------
cat("[3/5] Plotting heatmap...\n")
mat_rho <- dcast(cor_dt, cell_label ~ cancer, value.var = "rho")
mat_sig <- dcast(cor_dt, cell_label ~ cancer, value.var = "sig")
rn <- mat_rho$cell_label
mat_rho <- as.matrix(mat_rho[, -1])
mat_sig <- as.matrix(mat_sig[, -1])
rownames(mat_rho) <- rn
rownames(mat_sig) <- rn
mat_sig[is.na(mat_sig)] <- ""

pal <- colorRampPalette(c("#3B4992", "white", "#EE0000"))(100)

pheatmap(
    mat_rho,
    color           = pal,
    breaks          = seq(-0.6, 0.6, length.out = 101),
    cluster_rows    = TRUE,
    cluster_cols    = TRUE,
    display_numbers = mat_sig,
    number_color    = "black",
    fontsize_number = 7,
    fontsize_row    = 9,
    fontsize_col    = 9,
    cellwidth       = 16,
    cellheight      = 14,
    border_color    = "grey90",
    main            = "GATA2 vs CIBERSORT immune infiltration (Spearman rho)",
    filename        = file.path(FIGDIR, "Fig4_cibersort_heatmap.pdf"),
    width           = 11,
    height          = 7
)
pheatmap(
    mat_rho,
    color           = pal,
    breaks          = seq(-0.6, 0.6, length.out = 101),
    cluster_rows    = TRUE,
    cluster_cols    = TRUE,
    display_numbers = mat_sig,
    number_color    = "black",
    fontsize_number = 7,
    fontsize_row    = 9,
    fontsize_col    = 9,
    cellwidth       = 16,
    cellheight      = 14,
    border_color    = "grey90",
    main            = "GATA2 vs CIBERSORT immune infiltration (Spearman rho)",
    filename        = file.path(FIGDIR, "Fig4_cibersort_heatmap.png"),
    width           = 11,
    height          = 7
)

# ---- 5. focus-cancer boxplots: high vs low GATA2 -----------
cat("[4/5] Plotting focus-cancer high-vs-low boxplots...\n")

box_one <- function(cc) {
    d <- merged[cancer_type == cc]
    if (nrow(d) < 20) return(NULL)
    cutoff <- median(d$GATA2, na.rm = TRUE)
    d[, group := factor(ifelse(GATA2 > cutoff, "High", "Low"),
                        levels = c("Low", "High"))]

    long <- melt(d, id.vars = c("sample15", "group"),
                 measure.vars = cib_cols,
                 variable.name = "cell", value.name = "frac")
    long[, cell := sub("_CIBERSORT$", "", cell)]
    long[, cell := gsub("\\.", " ", cell)]

    # only keep cells with non-trivial fraction in this cancer
    keep_cells <- long[, .(med = median(frac, na.rm = TRUE)), by = cell][med > 0.005, cell]
    long <- long[cell %in% keep_cells]

    ggplot(long, aes(x = cell, y = frac, fill = group)) +
        geom_boxplot(outlier.size = 0.3, lwd = 0.3,
                     position = position_dodge(width = 0.8), width = 0.7) +
        stat_compare_means(aes(group = group), label = "p.signif",
                           method = "wilcox.test", size = 2.8,
                           hide.ns = TRUE,
                           label.y.npc = 0.95) +
        scale_fill_manual(values = c("Low" = "#4DBBD5", "High" = "#E64B35"),
                          name = paste0(cc, " — GATA2")) +
        labs(x = NULL, y = "Cell fraction (CIBERSORT)", title = cc) +
        theme_pubr(base_size = 10) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            plot.title  = element_text(face = "bold"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
            plot.background = element_rect(fill = "white", color = NA),
            legend.position = "top"
        )
}

box_list <- lapply(FOCUS, box_one)
box_list <- box_list[!sapply(box_list, is.null)]

if (length(box_list) > 0) {
    combined_box <- wrap_plots(box_list, ncol = 1)
    ggsave(file.path(FIGDIR, "Fig4_cibersort_focus_boxplots.pdf"),
           combined_box, width = 10, height = 4 * length(box_list),
           device = cairo_pdf, limitsize = FALSE)
    ggsave(file.path(FIGDIR, "Fig4_cibersort_focus_boxplots.png"),
           combined_box, width = 10, height = 4 * length(box_list),
           dpi = 300, bg = "white", limitsize = FALSE)
}

# ---- 6. summary --------------------------------------------
cat("\n[5/5] Done.\n")
cat("\n[saved]\n")
cat("  ", file.path(FIGDIR, "Fig4_cibersort_heatmap.pdf"), "\n")
cat("  ", file.path(FIGDIR, "Fig4_cibersort_focus_boxplots.pdf"), "\n")
cat("  ", file.path(RESULTS, "Fig4_cibersort_correlations.csv"), "\n")

cat("\nTop 10 most significant (|rho| ranked, p<0.001):\n")
print(cor_dt[p < 0.001][order(-abs(rho))][1:10, .(cancer, cell_label, rho, p, sig)])
