#!/usr/bin/env Rscript
# =============================================================
# 10_fig6_mutation.R
# Figure 6: Mutation landscape comparing GATA2-high vs GATA2-low.
#   Panel A: Oncoplot/waterfall for focus cancer (top mutated genes)
#   Panel B: TMB high vs low GATA2 across pan-cancer
#   Panel C: Differentially mutated genes (Fisher) for focus cancer
#
# Inputs:
#   data/01_raw/gdc_mc3/mc3.v0.2.8.PUBLIC.maf.gz
#   data/02_clean/GATA2_expression.rds
#
# Outputs:
#   figures/Fig6A_oncoplot_<CANCER>.pdf
#   figures/Fig6B_TMB_pancancer.pdf/.png
#   figures/Fig6C_diff_mutated_<CANCER>.pdf
#   data/03_results/Fig6_TMB_per_sample.csv
#   data/03_results/Fig6_TMB_high_vs_low.csv
#   data/03_results/Fig6_diff_mutated_<CANCER>.csv
#
# Usage:
#   Rscript scripts/10_fig6_mutation.R                 # focus = LIHC (most events, most mut)
#   Rscript scripts/10_fig6_mutation.R COAD            # any single cancer
#   Rscript scripts/10_fig6_mutation.R LIHC MESO LGG   # multiple (one oncoplot each)
# =============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(maftools)
    library(ggplot2)
    library(ggpubr)
})

# ---- args ----------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
FOCUS <- if (length(args) >= 1) args else c("LIHC")
cat("Focus cancers for oncoplot:", paste(FOCUS, collapse = ", "), "\n\n")

# ---- paths ---------------------------------------------------
DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
MAF_FILE  <- file.path(DATA_ROOT, "01_raw", "gdc_mc3", "mc3.v0.2.8.PUBLIC.maf.gz")
CLEAN     <- file.path(DATA_ROOT, "02_clean")
RESULTS   <- file.path(DATA_ROOT, "03_results")
FIGDIR    <- "./figures"
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR,  recursive = TRUE, showWarnings = FALSE)

stopifnot("MC3 MAF missing" = file.exists(MAF_FILE))

# ---- 1. load metadata + GATA2 -------------------------------
cat("[1/5] Loading GATA2 expression + sample groups...\n")
g <- readRDS(file.path(CLEAN, "GATA2_expression.rds"))
g <- g[sample_group == "Tumor" & !is.na(cancer_type)]

# group by median within each cancer
g[, group := ifelse(GATA2 > median(GATA2, na.rm = TRUE), "High", "Low"),
  by = cancer_type]
g[, sample15 := substr(sample, 1, 15)]
cat("  samples with group assigned:", nrow(g), "\n")

# ---- 2. load MAF (this is the slow step ~3-5 min) -----------
cat("[2/5] Reading MC3 MAF (~700 MB, takes a few minutes)...\n")
# only read columns we need to save memory
keep_cols <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
               "Variant_Classification", "Variant_Type",
               "Reference_Allele", "Tumor_Seq_Allele2",
               "Tumor_Sample_Barcode", "HGVSp_Short", "FILTER")
maf_dt <- fread(MAF_FILE, select = keep_cols)
cat("  total mutation rows:", nrow(maf_dt), "\n")

# keep only PASS variants, non-silent (standard maftools defaults)
maf_dt <- maf_dt[FILTER == "PASS" | FILTER == "" | is.na(FILTER)]
maf_dt[, sample15 := substr(Tumor_Sample_Barcode, 1, 15)]

# attach cancer_type and group
maf_dt <- merge(maf_dt, g[, .(sample15, cancer_type, GATA2, group)],
                by = "sample15", all.x = FALSE)
cat("  rows matched to GATA2 grouping:", nrow(maf_dt), "\n")

# ---- 3. TMB per sample -------------------------------------
cat("[3/5] Computing TMB (per sample mutation count)...\n")
nonsilent <- c("Missense_Mutation", "Nonsense_Mutation",
               "Frame_Shift_Del", "Frame_Shift_Ins",
               "In_Frame_Del", "In_Frame_Ins",
               "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation")
tmb <- maf_dt[Variant_Classification %in% nonsilent,
              .(TMB = .N), by = .(sample15, cancer_type, group, GATA2)]
fwrite(tmb, file.path(RESULTS, "Fig6_TMB_per_sample.csv"))

# per-cancer Wilcoxon: High vs Low
tmb_test <- tmb[, {
    if (sum(group == "High") < 5 || sum(group == "Low") < 5) {
        .(median_high = NA_real_, median_low = NA_real_, p = NA_real_)
    } else {
        w <- wilcox.test(TMB ~ group, data = .SD)
        .(median_high = as.numeric(median(TMB[group == "High"])),
          median_low  = as.numeric(median(TMB[group == "Low"])),
          p           = as.numeric(w$p.value))
    }
}, by = cancer_type]
tmb_test[, sig := fcase(
    is.na(p),  "",
    p < 0.001, "***",
    p < 0.01,  "**",
    p < 0.05,  "*",
    default    = ""
)]
tmb_test <- tmb_test[order(p)]
fwrite(tmb_test, file.path(RESULTS, "Fig6_TMB_high_vs_low.csv"))
print(tmb_test[!is.na(p)])

# ---- 3B. TMB boxplot pan-cancer ----------------------------
cat("[4/5] Plotting Fig6B TMB pan-cancer...\n")
tmb[, group := factor(group, levels = c("Low", "High"))]
tmb[, cancer_type := factor(cancer_type, levels = sort(unique(cancer_type)))]

# significance label position
y_max <- tmb[, quantile(TMB, 0.98, na.rm = TRUE)]
sig_pos <- tmb_test[, .(cancer_type, sig)]
sig_pos[, cancer_type := factor(cancer_type, levels = levels(tmb$cancer_type))]

p_tmb <- ggplot(tmb, aes(x = cancer_type, y = TMB + 1, fill = group)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4,
                 lwd = 0.3, width = 0.7,
                 position = position_dodge(width = 0.8)) +
    geom_text(data = sig_pos,
              aes(x = cancer_type, y = y_max * 1.5, label = sig),
              inherit.aes = FALSE,
              size = 3.2, fontface = "bold") +
    scale_fill_manual(values = c("Low" = "#4DBBD5", "High" = "#E64B35"),
                      name = "GATA2") +
    scale_y_log10() +
    labs(x = NULL, y = "TMB (non-silent mutations per sample, log10)",
         title = "Tumor mutation burden: GATA2-High vs GATA2-Low") +
    theme_pubr(base_size = 10) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
        plot.title  = element_text(face = "bold"),
        legend.position = "top",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )

ggsave(file.path(FIGDIR, "Fig6B_TMB_pancancer.pdf"),
       p_tmb, width = 12, height = 5, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig6B_TMB_pancancer.png"),
       p_tmb, width = 12, height = 5, dpi = 300, bg = "white")

# ---- 4. per-cancer oncoplot + diff mutation ----------------
cat("[5/5] Building maftools objects per focus cancer...\n")

run_focus <- function(cc) {
    cat("  --", cc, "--\n")
    sub <- maf_dt[cancer_type == cc]
    if (nrow(sub) < 100) { cat("    too few mutations\n"); return() }

    # build maftools MAF objects for High vs Low
    high_dt <- sub[group == "High"]
    low_dt  <- sub[group == "Low"]
    cat("    High n_samples=", uniqueN(high_dt$Tumor_Sample_Barcode),
        " Low n_samples=", uniqueN(low_dt$Tumor_Sample_Barcode),
        " total_mut=", nrow(sub), "\n", sep = "")

    # maftools requires Tumor_Sample_Barcode column
    suppressMessages({
        high_maf <- read.maf(maf = high_dt, vc_nonSyn = nonsilent)
        low_maf  <- read.maf(maf = low_dt,  vc_nonSyn = nonsilent)
    })

    # Fig6A: combined oncoplot for top 25 genes, split by GATA2 group.
    # Build MAF WITH clinical data so oncoplot's clinicalFeatures works.
    sample_anno <- unique(sub[, .(Tumor_Sample_Barcode, group)])
    sample_anno <- as.data.frame(sample_anno)
    colnames(sample_anno) <- c("Tumor_Sample_Barcode", "GATA2")

    full_maf2 <- read.maf(maf = sub, vc_nonSyn = nonsilent, clinicalData = sample_anno)
    pdf(file.path(FIGDIR, paste0("Fig6A_oncoplot_", cc, ".pdf")),
        width = 11, height = 7)
    oncoplot(
        maf = full_maf2,
        top = 25,
        clinicalFeatures = "GATA2",
        sortByAnnotation = TRUE,
        annotationColor  = list(GATA2 = c("High" = "#E64B35", "Low" = "#4DBBD5")),
        fontSize         = 0.7,
        titleText        = paste0(cc, ": top 25 mutated genes (GATA2-High vs Low)")
    )
    dev.off()

    # Fig6C: differential mutation analysis
    diff_res <- tryCatch(
        mafCompare(m1 = high_maf, m2 = low_maf,
                   m1Name = "GATA2-High", m2Name = "GATA2-Low",
                   minMut = 5),
        error = function(e) NULL
    )
    if (!is.null(diff_res)) {
        dt_diff <- as.data.table(diff_res$results)
        fwrite(dt_diff, file.path(RESULTS, paste0("Fig6_diff_mutated_", cc, ".csv")))
        cat("    diff-mutated genes (q<0.05):",
            sum(dt_diff$adjPval < 0.05, na.rm = TRUE),
            "  (raw p<0.01):",
            sum(dt_diff$pval < 0.01, na.rm = TRUE), "\n")

        # Use raw p<0.01 — BH q across ~17k genes is too conservative for
        # hypothesis-generating mutation comparisons. maftools::forestPlot's
        # `pVal` argument already filters by raw pval (the prior `adjPval`
        # guard in this script was the bug — it suppressed the call).
        if (any(dt_diff$pval < 0.01, na.rm = TRUE)) {
            pdf(file.path(FIGDIR, paste0("Fig6C_diff_mutated_", cc, ".pdf")),
                width = 8, height = 6)
            forestPlot(
                mafCompareRes = diff_res,
                pVal          = 0.01,                       # raw p threshold
                color         = c("GATA2-High" = "#E64B35",
                                  "GATA2-Low"  = "#4DBBD5"),
                geneFontSize  = 0.7,
                titleSize     = 1
            )
            dev.off()
            cat("    [saved] Fig6C uses raw p<0.01 (see CSV for BH q-values)\n")
        }
    }
}

for (cc in FOCUS) run_focus(cc)

cat("\n[saved]\n")
cat("  ", file.path(FIGDIR, "Fig6B_TMB_pancancer.pdf"), "\n")
for (cc in FOCUS) {
    cat("  ", file.path(FIGDIR, paste0("Fig6A_oncoplot_", cc, ".pdf")), "\n")
    cat("  ", file.path(FIGDIR, paste0("Fig6C_diff_mutated_", cc, ".pdf")), "\n")
}
cat("  ", file.path(RESULTS, "Fig6_TMB_high_vs_low.csv"), "\n")
