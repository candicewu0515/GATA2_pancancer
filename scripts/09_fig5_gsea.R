#!/usr/bin/env Rscript
# =============================================================
# 09_fig5_gsea.R
# Figure 5: GSEA of GATA2-high vs GATA2-low tumors across cancers.
#   Panel A: Pan-cancer bubble plot (HALLMARK 50, NES per cancer)
#   Panel B: Classic GSEA enrichment curves for focus cancers
#
# Workflow per cancer:
#   1. Split tumor samples by GATA2 median into High / Low
#   2. limma differential expression -> ranked list by t-statistic
#   3. clusterProfiler::GSEA() against HALLMARK + KEGG
#
# Inputs:
#   data/01_raw/xena_pancan/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
#   data/02_clean/GATA2_expression.rds
#   data/01_raw/msigdb/h.all.v2024.1.Hs.symbols.gmt
#   data/01_raw/msigdb/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt
#
# Outputs:
#   data/02_clean/limma_ranks/<CANCER>.rds         (gene-level t-stat ranks)
#   data/03_results/Fig5_gsea_hallmark.csv
#   data/03_results/Fig5_gsea_kegg.csv
#   figures/Fig5A_hallmark_bubble.pdf/.png
#   figures/Fig5B_focus_gsea_curves.pdf/.png
#
# Usage:
#   Rscript scripts/09_fig5_gsea.R                       # all cancers, focus = MESO LGG LIHC
#   Rscript scripts/09_fig5_gsea.R MESO LGG LIHC COAD    # custom focus
#   Rscript scripts/09_fig5_gsea.R --skip-limma          # reuse cached ranks (if rerunning)
# =============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(limma)
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
})

# ---- args ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
SKIP_LIMMA <- "--skip-limma" %in% args
args <- setdiff(args, "--skip-limma")
FOCUS <- if (length(args) >= 1) args else c("MESO", "LGG", "LIHC")
cat("Focus cancers:", paste(FOCUS, collapse = ", "), "\n")
cat("Skip limma (use cached):", SKIP_LIMMA, "\n\n")

# ---- paths ---------------------------------------------------
DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
RAW       <- file.path(DATA_ROOT, "01_raw", "xena_pancan")
MSIG      <- file.path(DATA_ROOT, "01_raw", "msigdb")
CLEAN     <- file.path(DATA_ROOT, "02_clean")
RANKDIR   <- file.path(CLEAN, "limma_ranks")
RESULTS   <- file.path(DATA_ROOT, "03_results")
FIGDIR    <- "./figures"
dir.create(RANKDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR,  showWarnings = FALSE, recursive = TRUE)

# ---- 1. load metadata + GATA2 -------------------------------
cat("[1/5] Loading sample metadata...\n")
gata2 <- readRDS(file.path(CLEAN, "GATA2_expression.rds"))
gata2 <- gata2[sample_group == "Tumor" & !is.na(cancer_type)]

# ---- 2. limma per cancer (cached) ---------------------------
EXPR_FILE <- file.path(RAW, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")
stopifnot(file.exists(EXPR_FILE))

if (!SKIP_LIMMA) {
    cat("[2/5] Reading full expression matrix (this is the slow step, ~2-3 min)...\n")
    expr <- fread(EXPR_FILE)
    gene_col <- colnames(expr)[1]
    setnames(expr, gene_col, "gene")
    cat("  matrix:", nrow(expr), "genes x", ncol(expr) - 1, "samples\n")

    cancers <- sort(unique(gata2$cancer_type))
    cat("[3/5] Running limma per cancer (", length(cancers), "cancers)...\n")

    for (cc in cancers) {
        out_rank <- file.path(RANKDIR, paste0(cc, ".rds"))
        if (file.exists(out_rank)) {
            cat("  [skip]", cc, "(cached)\n")
            next
        }

        d <- gata2[cancer_type == cc]
        if (nrow(d) < 20) { cat("  [skip]", cc, ": n =", nrow(d), "\n"); next }

        cutoff <- median(d$GATA2, na.rm = TRUE)
        d[, grp := ifelse(GATA2 > cutoff, "High", "Low")]
        # keep only samples in this cancer that are in the expression matrix
        d <- d[sample %in% colnames(expr)]
        if (nrow(d) < 20) { cat("  [skip]", cc, ": no match\n"); next }

        # extract sub-matrix
        m <- as.matrix(expr[, c("gene", d$sample), with = FALSE], rownames = "gene")
        # collapse duplicate gene symbols (median)
        if (anyDuplicated(rownames(m))) {
            m <- limma::avereps(m)
        }

        grp <- factor(d$grp, levels = c("Low", "High"))
        design <- model.matrix(~ grp)
        fit <- lmFit(m, design)
        fit <- eBayes(fit)
        tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

        # ranked vector: t-stat
        ranks <- tt$t
        names(ranks) <- rownames(m)
        ranks <- ranks[!is.na(ranks)]
        ranks <- sort(ranks, decreasing = TRUE)

        saveRDS(ranks, out_rank)
        cat("  [done]", cc, " n=", nrow(d), " (", sum(grp == "High"), "H /", sum(grp == "Low"), "L)",
            " ranks=", length(ranks), "\n", sep = "")
    }

    rm(expr); gc()
} else {
    cat("[2/5] Skipping limma, using cached ranks.\n")
}

# ---- 3. load MSigDB -----------------------------------------
cat("[3/5] Loading MSigDB gene sets...\n")
read_gmt_robust <- function(path) {
    lines <- readLines(path)
    splits <- strsplit(lines, "\t", fixed = TRUE)
    sets <- lapply(splits, function(x) x[-c(1, 2)])
    names(sets) <- sapply(splits, `[`, 1)
    sets
}

hallmark <- read_gmt_robust(file.path(MSIG, "h.all.v2024.1.Hs.symbols.gmt"))
kegg     <- read_gmt_robust(file.path(MSIG, "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"))
cat("  HALLMARK pathways:", length(hallmark), "\n")
cat("  KEGG pathways    :", length(kegg), "\n")

# turn into TERM2GENE for clusterProfiler
to_t2g <- function(lst) {
    rbindlist(Map(function(name, genes) data.table(term = name, gene = genes),
                  names(lst), lst))
}
t2g_hallmark <- to_t2g(hallmark)
t2g_kegg     <- to_t2g(kegg)

# ---- 4. GSEA per cancer per database ------------------------
cat("[4/5] Running GSEA...\n")

run_gsea_one <- function(ranks, t2g, db_name) {
    set.seed(1)
    tryCatch(
        as.data.table(as.data.frame(GSEA(
            geneList    = ranks,
            TERM2GENE   = t2g,
            minGSSize   = 10,
            maxGSSize   = 500,
            pvalueCutoff = 1,        # keep all, filter later
            verbose     = FALSE,
            seed        = TRUE
        ))),
        error = function(e) { cat("    GSEA failed:", conditionMessage(e), "\n"); NULL }
    )
}

rank_files <- list.files(RANKDIR, pattern = "\\.rds$", full.names = TRUE)
cat("  rank files found:", length(rank_files), "\n")

# pan-cancer collection
hm_all <- list()
kg_all <- list()
focus_gsea_obj <- list()

for (rf in rank_files) {
    cc <- sub("\\.rds$", "", basename(rf))
    ranks <- readRDS(rf)
    cat("  -", cc, "(n_genes=", length(ranks), ") ", sep = "")

    # HALLMARK
    rh <- run_gsea_one(ranks, t2g_hallmark, "HALLMARK")
    if (!is.null(rh) && nrow(rh) > 0) {
        rh[, cancer := cc]; rh[, db := "HALLMARK"]
        hm_all[[cc]] <- rh
    }

    # KEGG (don't cache full obj for KEGG, only data)
    rk <- run_gsea_one(ranks, t2g_kegg, "KEGG")
    if (!is.null(rk) && nrow(rk) > 0) {
        rk[, cancer := cc]; rk[, db := "KEGG"]
        kg_all[[cc]] <- rk
    }

    # save full GSEA object for focus cancers (for enrichment curves)
    if (cc %in% FOCUS) {
        set.seed(1)
        focus_gsea_obj[[cc]] <- tryCatch(GSEA(
            geneList    = ranks,
            TERM2GENE   = t2g_hallmark,
            minGSSize   = 10, maxGSSize = 500,
            pvalueCutoff = 0.25,
            verbose = FALSE, seed = TRUE
        ), error = function(e) NULL)
    }

    cat(" H=", if (is.null(rh)) 0 else nrow(rh),
        " K=", if (is.null(rk)) 0 else nrow(rk), "\n", sep = "")
}

hm <- rbindlist(hm_all, fill = TRUE)
kg <- rbindlist(kg_all, fill = TRUE)
fwrite(hm, file.path(RESULTS, "Fig5_gsea_hallmark.csv"))
fwrite(kg, file.path(RESULTS, "Fig5_gsea_kegg.csv"))
cat("  saved HALLMARK rows:", nrow(hm), "\n")
cat("  saved KEGG rows    :", nrow(kg), "\n")

# ---- 5A. pan-cancer bubble plot (HALLMARK) ------------------
cat("[5/5] Plotting...\n")

# select pathways significant in >=3 cancers (p.adjust<0.05) for plotting
sig_paths <- hm[p.adjust < 0.05, .N, by = ID][N >= 3, ID]
pp <- hm[ID %in% sig_paths]
pp[, ID := sub("^HALLMARK_", "", ID)]
pp[, ID := gsub("_", " ", ID)]
# cap NES for color scale visibility
pp[, NES_clip := pmax(pmin(NES, 3), -3)]
pp[, log10p   := -log10(pmax(p.adjust, 1e-10))]

# order
pp[, cancer := factor(cancer, levels = sort(unique(cancer)))]
# order pathways by mean |NES|
order_path <- pp[, .(score = mean(abs(NES))), by = ID][order(-score), ID]
pp[, ID := factor(ID, levels = rev(order_path))]

p_bubble <- ggplot(pp, aes(x = cancer, y = ID)) +
    geom_point(aes(size = log10p, color = NES_clip)) +
    scale_color_gradient2(low = "#3B4992", mid = "white", high = "#EE0000",
                          midpoint = 0, limits = c(-3, 3),
                          name = "NES") +
    scale_size_continuous(name = expression(-log[10] ~ italic(q)),
                          range = c(0.5, 5)) +
    labs(x = NULL, y = NULL,
         title = "Pan-cancer GSEA: HALLMARK pathways (GATA2-High vs Low)") +
    theme_pubr(base_size = 10) +
    theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
        axis.text.y  = element_text(size = 8),
        plot.title   = element_text(face = "bold", size = 11),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "right"
    )

h_bubble <- max(6, min(15, length(unique(pp$ID)) * 0.25 + 3))
ggsave(file.path(FIGDIR, "Fig5A_hallmark_bubble.pdf"),
       p_bubble, width = 12, height = h_bubble, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig5A_hallmark_bubble.png"),
       p_bubble, width = 12, height = h_bubble, dpi = 300, bg = "white")

# ---- 5B. focus GSEA enrichment curves -----------------------
make_focus_curve <- function(cc, obj) {
    if (is.null(obj) || nrow(as.data.frame(obj)) == 0) return(NULL)
    df <- as.data.frame(obj)
    # top 5 pathways by |NES|, q<0.25
    df <- df[df$qvalue < 0.25, ]
    if (nrow(df) == 0) return(NULL)
    df <- df[order(-abs(df$NES)), ]
    top_ids <- head(df$ID, 5)
    p <- gseaplot2(
        obj, geneSetID = top_ids,
        title = paste0(cc, ": top HALLMARK enrichments"),
        color = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F"),
        base_size = 10,
        pvalue_table = FALSE
    )
    p
}

curves <- Map(make_focus_curve, names(focus_gsea_obj), focus_gsea_obj)
curves <- curves[!sapply(curves, is.null)]

# Save each focus cancer's curves on its own (gseaplot2 returns a stacked
# gtable, not a plain ggplot — robust to save individually).
for (cc in names(curves)) {
    ggsave(file.path(FIGDIR, paste0("Fig5B_gsea_", cc, ".pdf")),
           curves[[cc]], width = 7, height = 6, device = cairo_pdf)
    ggsave(file.path(FIGDIR, paste0("Fig5B_gsea_", cc, ".png")),
           curves[[cc]], width = 7, height = 6, dpi = 300, bg = "white")
}

# NOTE: gseaplot2 returns a `gglistlist` (aplot stack), which ggsave handles
# directly but patchwork::wrap_plots / cowplot::plot_grid can't combine.
# We skip the multi-panel combined PDF and rely on the per-cancer files above.

# ---- summary -----------------------------------------------
cat("\n[saved]\n")
cat("  ", file.path(FIGDIR, "Fig5A_hallmark_bubble.pdf"), "\n")
cat("  ", file.path(FIGDIR, "Fig5B_focus_gsea_curves.pdf"), "\n")
cat("  ", file.path(RESULTS, "Fig5_gsea_hallmark.csv"), "\n")
cat("  ", file.path(RESULTS, "Fig5_gsea_kegg.csv"), "\n")

cat("\nTop 10 HALLMARK enrichments (|NES|, p.adjust<0.05):\n")
print(hm[p.adjust < 0.05][order(-abs(NES))][1:10,
    .(cancer, pathway = sub("^HALLMARK_", "", ID), NES = round(NES, 2),
      p.adjust = signif(p.adjust, 2))])
