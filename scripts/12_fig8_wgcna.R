#!/usr/bin/env Rscript
# =============================================================
# 12_fig8_wgcna.R
# Figure 8: WGCNA in focus cancer.
#   Panel A: soft-threshold power selection (scale-free fit)
#   Panel B: module-trait correlation heatmap (GATA2 + clinical)
#   Panel C: KEGG/GO enrichment of GATA2-most-correlated module
#
# Strategy (memory-conservative for Mac):
#   - Subset to focus cancer tumor samples only
#   - Top 5000 most-variable genes
#   - blockwiseModules with maxBlockSize = 5000 (single block)
#
# Inputs:
#   data/01_raw/xena_pancan/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
#   data/02_clean/GATA2_expression.rds
#
# Outputs:
#   data/02_clean/wgcna_<CANCER>.rds   (full network object)
#   data/03_results/Fig8_modules_<CANCER>.csv  (gene -> module mapping)
#   data/03_results/Fig8_module_trait_<CANCER>.csv
#   data/03_results/Fig8_enrich_<CANCER>.csv
#   figures/Fig8A_sft_<CANCER>.pdf
#   figures/Fig8B_module_trait_<CANCER>.pdf/.png
#   figures/Fig8C_enrichment_<CANCER>.pdf/.png
#
# Usage:
#   Rscript scripts/12_fig8_wgcna.R                  # default = LIHC
#   Rscript scripts/12_fig8_wgcna.R COAD             # any single cancer
#
# Note: WGCNA is the heaviest step. ~10-15 min, ~8-12 GB RAM.
# =============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(WGCNA)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(ggpubr)
    library(pheatmap)
    library(RColorBrewer)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = max(1, parallel::detectCores() - 2))

# ---- args ----------------------------------------------------
args   <- commandArgs(trailingOnly = TRUE)
CANCER <- if (length(args) >= 1) args[1] else "LIHC"
TOP_N  <- 5000
cat("Cancer:", CANCER, " | top-variable genes:", TOP_N, "\n\n")

# ---- paths ---------------------------------------------------
DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
RAW       <- file.path(DATA_ROOT, "01_raw", "xena_pancan")
CLEAN     <- file.path(DATA_ROOT, "02_clean")
RESULTS   <- file.path(DATA_ROOT, "03_results")
FIGDIR    <- "./figures"
for (p in c(CLEAN, RESULTS, FIGDIR)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

EXPR_FILE <- file.path(RAW, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")

# ---- 1. load + subset ---------------------------------------
cat("[1/6] Loading data...\n")
g <- readRDS(file.path(CLEAN, "GATA2_expression.rds"))
g <- g[sample_group == "Tumor" & cancer_type == CANCER]
cat("  ", CANCER, "tumor samples in metadata:", nrow(g), "\n")

cat("[2/6] Reading expression matrix...\n")
expr <- fread(EXPR_FILE)
setnames(expr, 1, "gene")

# subset columns to this cancer's samples
keep_samples <- intersect(g$sample, colnames(expr))
cat("  samples in expression matrix:", length(keep_samples), "\n")
m <- as.matrix(expr[, c("gene", keep_samples), with = FALSE], rownames = "gene")
rm(expr); gc()

# collapse duplicate gene symbols
if (anyDuplicated(rownames(m))) m <- limma::avereps(m)

# top variable genes
v <- apply(m, 1, var, na.rm = TRUE)
top <- names(sort(v, decreasing = TRUE))[1:min(TOP_N, length(v))]
m <- m[top, ]
cat("  matrix for WGCNA:", nrow(m), "genes x", ncol(m), "samples\n")

# transpose: WGCNA wants samples x genes
datExpr <- t(m)

# QC: drop samples with too many NAs / zero variance
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    cat("  after QC:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")
}

# ---- 2. soft-threshold power -------------------------------
cat("[3/6] Picking soft-threshold power...\n")
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0,
                         networkType = "signed")
beta <- sft$powerEstimate
if (is.na(beta)) beta <- 9   # fallback to original paper's choice
cat("  selected beta =", beta, "\n")

# plot SFT
pdf(file.path(FIGDIR, paste0("Fig8A_sft_", CANCER, ".pdf")), width = 9, height = 4.5)
par(mfrow = c(1, 2), mar = c(4, 4.2, 2, 1))
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     type = "b", pch = 19, col = "#E64B35",
     xlab = "Soft-threshold (power)",
     ylab = "Scale-free topology fit (signed R^2)",
     main = paste0(CANCER, ": scale-free fit"))
abline(h = 0.85, lty = 2, col = "grey40")
abline(v = beta,  lty = 2, col = "#4DBBD5")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     type = "b", pch = 19, col = "#4DBBD5",
     xlab = "Soft-threshold (power)",
     ylab = "Mean connectivity",
     main = paste0(CANCER, ": connectivity"))
abline(v = beta, lty = 2, col = "#E64B35")
dev.off()

# ---- 3. network construction --------------------------------
cat("[4/6] Building network (this is the slow step, ~3-8 min)...\n")
# Mask base::cor so WGCNA's internal calls reach WGCNA::cor — without this,
# blockwiseModules crashes with "unused arguments (weights.x, weights.y, cosine)"
# in recent R because base cor() doesn't accept those args.
cor <- WGCNA::cor
on.exit(rm(cor), add = TRUE)
net <- blockwiseModules(
    datExpr,
    power           = beta,
    networkType     = "signed",
    TOMType         = "signed",
    minModuleSize   = 30,
    reassignThreshold = 0,
    mergeCutHeight  = 0.25,
    numericLabels   = FALSE,
    maxBlockSize    = TOP_N + 100,
    saveTOMs        = FALSE,
    verbose         = 0
)
saveRDS(net, file.path(CLEAN, paste0("wgcna_", CANCER, ".rds")))

modules <- data.table(gene = colnames(datExpr), module = net$colors)
modules[, n := .N, by = module]
cat("  modules found:", uniqueN(modules$module), "\n")
print(modules[, .N, by = module][order(-N)])
fwrite(modules, file.path(RESULTS, paste0("Fig8_modules_", CANCER, ".csv")))

# ---- 4. module-trait correlation ----------------------------
cat("[5/6] Computing module-trait correlations...\n")
MEs <- net$MEs
sample_ids <- rownames(datExpr)
gsub_meta <- g[match(sample_ids, sample), ]

traits <- data.frame(
    GATA2 = gsub_meta$GATA2,
    age   = suppressWarnings(as.numeric(gsub_meta$age)),
    OS    = as.integer(gsub_meta$OS),
    OS_time_mo = gsub_meta$OS.time / 30.44
)
if ("gender" %in% colnames(gsub_meta)) {
    traits$gender_M <- as.integer(gsub_meta$gender == "MALE")
}
if ("stage" %in% colnames(gsub_meta)) {
    traits$stage_num <- as.integer(factor(gsub_meta$stage,
                                          levels = c("I", "II", "III", "IV")))
}

# drop trait columns that are all-NA
traits <- traits[, colSums(!is.na(traits)) > nrow(traits) * 0.5, drop = FALSE]
cat("  trait columns:", paste(colnames(traits), collapse = ", "), "\n")

mod_trait_cor <- cor(MEs, traits, use = "pairwise.complete.obs")
mod_trait_p   <- corPvalueStudent(mod_trait_cor, nrow(datExpr))

mt_dt <- data.table(
    module = rownames(mod_trait_cor),
    as.data.table(round(mod_trait_cor, 3))
)
fwrite(mt_dt, file.path(RESULTS, paste0("Fig8_module_trait_", CANCER, ".csv")))

# plot heatmap
disp <- matrix(paste0(round(mod_trait_cor, 2),
                      "\n(", format(mod_trait_p, digits = 1, scientific = TRUE), ")"),
               nrow = nrow(mod_trait_cor))
rownames(disp) <- rownames(mod_trait_cor)
colnames(disp) <- colnames(mod_trait_cor)

pheatmap(
    mod_trait_cor,
    cluster_rows = FALSE, cluster_cols = FALSE,
    display_numbers = disp,
    fontsize_number = 7,
    fontsize_row = 9, fontsize_col = 10,
    color = colorRampPalette(c("#3B4992", "white", "#EE0000"))(100),
    breaks = seq(-1, 1, length.out = 101),
    cellwidth = 50, cellheight = 22,
    main = paste0(CANCER, ": module-trait correlation"),
    filename = file.path(FIGDIR, paste0("Fig8B_module_trait_", CANCER, ".pdf")),
    width = max(6, 1.5 + ncol(mod_trait_cor) * 0.9),
    height = max(5, 1 + nrow(mod_trait_cor) * 0.45)
)
pheatmap(
    mod_trait_cor,
    cluster_rows = FALSE, cluster_cols = FALSE,
    display_numbers = disp,
    fontsize_number = 7,
    fontsize_row = 9, fontsize_col = 10,
    color = colorRampPalette(c("#3B4992", "white", "#EE0000"))(100),
    breaks = seq(-1, 1, length.out = 101),
    cellwidth = 50, cellheight = 22,
    main = paste0(CANCER, ": module-trait correlation"),
    filename = file.path(FIGDIR, paste0("Fig8B_module_trait_", CANCER, ".png")),
    width = max(6, 1.5 + ncol(mod_trait_cor) * 0.9),
    height = max(5, 1 + nrow(mod_trait_cor) * 0.45)
)

# ---- 5. enrichment of top-GATA2-correlated module -----------
cat("[6/6] Enrichment analysis of top GATA2-correlated module...\n")
gata2_cor <- mod_trait_cor[, "GATA2"]
top_mod <- names(sort(abs(gata2_cor), decreasing = TRUE))[1]
top_mod_color <- sub("^ME", "", top_mod)
cat("  top module:", top_mod, "(cor with GATA2 =",
    round(gata2_cor[top_mod], 3), ")\n")

mod_genes <- modules[module == top_mod_color, gene]
cat("  module gene count:", length(mod_genes), "\n")

# map to Entrez
gene_ids <- bitr(mod_genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db, drop = TRUE)

# KEGG enrichment
kegg <- tryCatch(enrichKEGG(gene = gene_ids$ENTREZID, organism = "hsa",
                            pvalueCutoff = 0.05),
                 error = function(e) NULL)
go_bp <- tryCatch(enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db,
                           ont = "BP", pvalueCutoff = 0.05, readable = TRUE),
                  error = function(e) NULL)

enrich_combined <- rbindlist(list(
    if (!is.null(kegg)  && nrow(kegg)  > 0) data.table(db = "KEGG",  as.data.frame(kegg))  else NULL,
    if (!is.null(go_bp) && nrow(go_bp) > 0) data.table(db = "GO_BP", as.data.frame(go_bp)) else NULL
), fill = TRUE)

if (nrow(enrich_combined) > 0) {
    fwrite(enrich_combined, file.path(RESULTS, paste0("Fig8_enrich_", CANCER, ".csv")))
}

# plot top 10 from each
plot_enrich <- function(res, title, color) {
    if (is.null(res) || nrow(res) == 0) return(NULL)
    df <- as.data.frame(res)
    df <- df[order(df$p.adjust), ]
    df <- head(df, 10)
    df$Description <- factor(df$Description, levels = rev(df$Description))
    df$log10p <- -log10(df$p.adjust)
    ggplot(df, aes(x = log10p, y = Description, fill = log10p)) +
        geom_col(width = 0.7) +
        scale_fill_gradient(low = "#A6CEE3", high = color, guide = "none") +
        labs(x = expression(-log[10] ~ "q-value"), y = NULL, title = title) +
        theme_pubr(base_size = 9) +
        theme(plot.title = element_text(face = "bold", size = 10),
              axis.text.y = element_text(size = 8),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
              plot.background = element_rect(fill = "white", color = NA))
}

p1 <- plot_enrich(kegg, "KEGG", "#E64B35")
p2 <- plot_enrich(go_bp, "GO Biological Process", "#3C5488")

plots <- list(p1, p2)
plots <- plots[!sapply(plots, is.null)]
if (length(plots) > 0) {
    combined <- patchwork::wrap_plots(plots, ncol = length(plots))
    ggsave(file.path(FIGDIR, paste0("Fig8C_enrichment_", CANCER, ".pdf")),
           combined, width = 5 * length(plots), height = 5, device = cairo_pdf)
    ggsave(file.path(FIGDIR, paste0("Fig8C_enrichment_", CANCER, ".png")),
           combined, width = 5 * length(plots), height = 5, dpi = 300, bg = "white")
}

cat("\n[saved]\n")
cat("  ", file.path(FIGDIR, paste0("Fig8A_sft_", CANCER, ".pdf")), "\n")
cat("  ", file.path(FIGDIR, paste0("Fig8B_module_trait_", CANCER, ".pdf")), "\n")
cat("  ", file.path(FIGDIR, paste0("Fig8C_enrichment_", CANCER, ".pdf")), "\n")
cat("  ", file.path(RESULTS, paste0("Fig8_modules_", CANCER, ".csv")), "\n")
cat("  ", file.path(RESULTS, paste0("Fig8_module_trait_", CANCER, ".csv")), "\n")
cat("  ", file.path(RESULTS, paste0("Fig8_enrich_", CANCER, ".csv")), "\n")
