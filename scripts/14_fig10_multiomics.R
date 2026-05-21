#!/usr/bin/env Rscript
# =============================================================
# 14_fig10_multiomics.R
# Figure 10: Multi-omics validation of GATA2.
#   Panel A: DNA methylation of GATA2 promoter vs mRNA expression
#            (per cancer correlation heatmap + scatter for focus)
#   Panel B: GATA2 mRNA vs methylation tumor-normal difference
#            (volcano/scatter across cancers)
#   Panel C: Protein-level info pointer (link to CPTAC/UALCAN webtools)
#            -- this panel is text + screenshot placeholder.
#            For scripted protein analysis, see methylation only here.
#
# Inputs:
#   data/01_raw/xena_pancan/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz
#   data/02_clean/GATA2_expression.rds
#
# Outputs:
#   data/02_clean/GATA2_methylation.rds
#   data/03_results/Fig10_methylation_correlations.csv
#   figures/Fig10A_methylation_heatmap.pdf/.png
#   figures/Fig10B_methylation_scatter.pdf/.png
#
# Strategy:
#   - The 450K matrix is HUGE (~14 GB uncompressed). Use data.table::fread
#     to scan first column for GATA2 probes, then read only those rows.
#   - GATA2 CpG probes (from Illumina 450K manifest, around the GATA2 TSS):
#     promoter region probes are mostly within chr3:128,200,000-128,205,000.
#     Common GATA2 probes: cg17004290, cg08572970, cg12480923, cg06125883,
#     cg18351549, cg14797808, cg25671438 (TSS1500/TSS200 + 5'UTR).
#
# Usage:
#   Rscript scripts/14_fig10_multiomics.R               # default focus = LIHC
#   Rscript scripts/14_fig10_multiomics.R MESO LGG LIHC
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
args  <- commandArgs(trailingOnly = TRUE)
FOCUS <- if (length(args) >= 1) args else c("LIHC")
cat("Focus cancers:", paste(FOCUS, collapse = ", "), "\n\n")

# ---- paths ---------------------------------------------------
DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
RAW       <- file.path(DATA_ROOT, "01_raw", "xena_pancan")
CLEAN     <- file.path(DATA_ROOT, "02_clean")
RESULTS   <- file.path(DATA_ROOT, "03_results")
FIGDIR    <- "./figures"
for (p in c(CLEAN, RESULTS, FIGDIR)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

METH_FILE <- file.path(RAW, "jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz")
stopifnot("methylation file missing" = file.exists(METH_FILE))

# Known GATA2 CpG probes (Illumina 450K, TSS200/TSS1500/5'UTR)
GATA2_PROBES <- c("cg17004290", "cg08572970", "cg12480923", "cg06125883",
                  "cg18351549", "cg14797808", "cg25671438", "cg00018413",
                  "cg27531412", "cg11484872", "cg05222924", "cg21179742")

# ---- 1. extract GATA2 probes from methylation matrix --------
meth_rds <- file.path(CLEAN, "GATA2_methylation.rds")
if (file.exists(meth_rds)) {
    cat("[1/4] Loading cached GATA2 methylation...\n")
    meth_dt <- readRDS(meth_rds)
} else {
    cat("[1/4] Scanning 450K matrix for GATA2 probes (~3-5 min)...\n")

    # read header for sample columns
    hdr <- fread(METH_FILE, nrows = 0)
    sample_ids <- colnames(hdr)[-1]
    cat("  samples in 450K matrix:", length(sample_ids), "\n")

    # read the probe ID column (first column) only
    probes <- fread(METH_FILE, select = 1, header = TRUE)
    probe_ids <- probes[[1]]
    cat("  total probes:", length(probe_ids), "\n")

    hits <- which(probe_ids %in% GATA2_PROBES)
    cat("  GATA2 probes found:", length(hits), "of", length(GATA2_PROBES),
        "\n  -", paste(probe_ids[hits], collapse = ", "), "\n")

    if (length(hits) == 0) {
        stop("No GATA2 probes found.  Inspect ", METH_FILE, " manually.")
    }

    # read those rows via gz stream (skip approach)
    con <- gzfile(METH_FILE, "r")
    on.exit(close(con))
    invisible(readLines(con, n = 1))  # skip header
    selected <- list()
    max_row <- max(hits)
    cur <- 0L
    while (cur < max_row) {
        line <- readLines(con, n = 1)
        if (length(line) == 0) break
        cur <- cur + 1L
        if (cur %in% hits) {
            parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
            selected[[parts[1]]] <- suppressWarnings(as.numeric(parts[-1]))
        }
    }

    cat("  rows read:", length(selected), "\n")
    meth_dt <- data.table(
        sample = sample_ids,
        as.data.table(do.call(cbind, selected))
    )
    saveRDS(meth_dt, meth_rds)
}
cat("  methylation samples:", nrow(meth_dt), "x", ncol(meth_dt) - 1, "probes\n")

# ---- 2. merge with expression -------------------------------
cat("[2/4] Merging with expression metadata...\n")
g <- readRDS(file.path(CLEAN, "GATA2_expression.rds"))
g <- g[!is.na(cancer_type)]

merged <- merge(g, meth_dt, by = "sample")
probe_cols <- intersect(colnames(meth_dt)[-1], colnames(merged))
cat("  merged samples:", nrow(merged), "\n")

# mean beta across all GATA2 probes per sample (proxy of promoter methylation)
merged[, GATA2_meth := rowMeans(.SD, na.rm = TRUE), .SDcols = probe_cols]

# ---- 3. per-cancer correlation (tumor only) -----------------
cat("[3/4] Computing per-cancer correlations...\n")
tumor <- merged[sample_group == "Tumor"]

cor_dt <- tumor[, {
    out <- list()
    for (pc in c(probe_cols, "GATA2_meth")) {
        v <- get(pc)
        ok <- !is.na(v) & !is.na(GATA2)
        if (sum(ok) < 15) next
        ct <- suppressWarnings(cor.test(GATA2[ok], v[ok], method = "pearson"))
        out[[pc]] <- list(rho = as.numeric(ct$estimate),
                          p   = as.numeric(ct$p.value),
                          n   = sum(ok))
    }
    rbindlist(lapply(names(out), function(nm) {
        data.table(probe = nm, rho = out[[nm]]$rho,
                   p = out[[nm]]$p, n = out[[nm]]$n)
    }))
}, by = cancer_type]

cor_dt[, sig := fcase(p < 0.001, "***", p < 0.01, "**", p < 0.05, "*", default = "")]
fwrite(cor_dt, file.path(RESULTS, "Fig10_methylation_correlations.csv"))

# ---- 4. plots -----------------------------------------------
cat("[4/4] Plotting...\n")

# Fig 10A: heatmap (probes x cancers)
mat <- dcast(cor_dt, probe ~ cancer_type, value.var = "rho", fill = NA_real_)
sig <- dcast(cor_dt, probe ~ cancer_type, value.var = "sig", fill = "")
rn <- mat$probe
mat <- as.matrix(mat[, -1]); rownames(mat) <- rn
sig <- as.matrix(sig[, -1]); rownames(sig) <- rn
sig[is.na(sig)] <- ""

# order: mean probe first, then individual probes
ord <- c("GATA2_meth", setdiff(rownames(mat), "GATA2_meth"))
mat <- mat[ord, ]
sig <- sig[ord, ]

pal <- colorRampPalette(c("#3B4992", "white", "#EE0000"))(100)
pheatmap(
    mat, color = pal,
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = FALSE, cluster_cols = TRUE,
    display_numbers = sig, number_color = "black",
    fontsize_number = 8, fontsize_row = 9, fontsize_col = 9,
    cellwidth = 16, cellheight = 18,
    border_color = "grey90",
    main = "GATA2 mRNA vs promoter methylation (Pearson r, tumor only)",
    filename = file.path(FIGDIR, "Fig10A_methylation_heatmap.pdf"),
    width = 12, height = max(4, nrow(mat) * 0.4 + 2)
)
pheatmap(
    mat, color = pal,
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = FALSE, cluster_cols = TRUE,
    display_numbers = sig, number_color = "black",
    fontsize_number = 8, fontsize_row = 9, fontsize_col = 9,
    cellwidth = 16, cellheight = 18,
    border_color = "grey90",
    main = "GATA2 mRNA vs promoter methylation (Pearson r, tumor only)",
    filename = file.path(FIGDIR, "Fig10A_methylation_heatmap.png"),
    width = 12, height = max(4, nrow(mat) * 0.4 + 2)
)

# Fig 10B: scatter for focus cancers (GATA2 mRNA vs mean beta)
scatter_one <- function(cc) {
    d <- tumor[cancer_type == cc & !is.na(GATA2_meth)]
    if (nrow(d) < 20) return(NULL)
    ct <- cor.test(d$GATA2, d$GATA2_meth, method = "pearson")
    ggplot(d, aes(x = GATA2_meth, y = GATA2)) +
        geom_point(size = 1.6, alpha = 0.6, color = "#E64B35") +
        geom_smooth(method = "lm", formula = y ~ x,
                    color = "grey25", fill = "grey80",
                    linewidth = 0.5, alpha = 0.3, se = TRUE) +
        labs(
            title = cc,
            subtitle = sprintf("r = %.3f, p = %.2g, n = %d",
                               ct$estimate, ct$p.value, nrow(d)),
            x = "GATA2 promoter mean beta",
            y = expression(italic("GATA2") ~ "mRNA  log"[2])
        ) +
        theme_pubr(base_size = 10) +
        theme(plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(size = 9, color = "grey30"),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
              plot.background = element_rect(fill = "white", color = NA))
}

scatters <- lapply(FOCUS, scatter_one)
scatters <- scatters[!sapply(scatters, is.null)]

if (length(scatters) > 0) {
    combined <- wrap_plots(scatters, ncol = length(scatters))
    ggsave(file.path(FIGDIR, "Fig10B_methylation_scatter.pdf"),
           combined, width = 4 * length(scatters), height = 4, device = cairo_pdf)
    ggsave(file.path(FIGDIR, "Fig10B_methylation_scatter.png"),
           combined, width = 4 * length(scatters), height = 4, dpi = 300, bg = "white")
}

# ---- summary -----------------------------------------------
cat("\n[saved]\n")
cat("  ", file.path(FIGDIR, "Fig10A_methylation_heatmap.pdf"), "\n")
cat("  ", file.path(FIGDIR, "Fig10B_methylation_scatter.pdf"), "\n")
cat("  ", file.path(RESULTS, "Fig10_methylation_correlations.csv"), "\n")

cat("\nTop negative correlations (methylation -> low expression, expected):\n")
print(cor_dt[probe == "GATA2_meth" & p < 0.05][order(rho)][1:10])

cat("\n--- Panel C (protein) note ---\n")
cat("  CPTAC protein data is best accessed via UALCAN:\n")
cat("  https://ualcan.path.uab.edu/analysis-prot.html\n")
cat("  Search 'GATA2', screenshot the boxplots for cancers of interest\n")
cat("  (e.g. LIHC, COAD, KIRC, BRCA -- all have CPTAC proteomics).\n")
