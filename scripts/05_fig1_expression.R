#!/usr/bin/env Rscript
# =============================================================
# 05_fig1_expression.R
# Figure 1: GATA2 expression Tumor vs Normal across TCGA cancers.
#
# Input:
#   data/02_clean/GATA2_expression.rds
#
# Output:
#   figures/Fig1_expression.pdf
#   figures/Fig1_expression.png
#   data/03_results/Fig1_stats.csv   <- Wilcoxon p-value per cancer
#
# Usage:
#   Rscript scripts/05_fig1_expression.R                # all 33 cancers
#   Rscript scripts/05_fig1_expression.R COAD           # single cancer (saves a separate file)
# =============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(ggsci)        # Science/Nature/JCO/NPG palettes
})

# ---- args ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
ONE_CANCER <- if (length(args) >= 1) args[1] else NA_character_

# ---- paths ---------------------------------------------------
DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
CLEAN     <- file.path(DATA_ROOT, "02_clean")
RESULTS   <- file.path(DATA_ROOT, "03_results")
FIGDIR    <- "./figures"
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR,  recursive = TRUE, showWarnings = FALSE)

# ---- load ----------------------------------------------------
cat("[1/4] Loading GATA2 expression...\n")
dt <- readRDS(file.path(CLEAN, "GATA2_expression.rds"))
cat("  total rows:", nrow(dt), "\n")

# keep Tumor + Normal only (drop Metastatic/Recurrent for clean comparison)
dt <- dt[sample_group %in% c("Tumor", "Normal") & !is.na(cancer_type)]
dt[, sample_group := factor(sample_group, levels = c("Normal", "Tumor"))]
cat("  after filter (Tumor/Normal only):", nrow(dt), "\n")

# optional: subset to a single cancer
if (!is.na(ONE_CANCER)) {
    if (!ONE_CANCER %in% dt$cancer_type) stop("Cancer type not found: ", ONE_CANCER)
    dt <- dt[cancer_type == ONE_CANCER]
    cat("  filtered to:", ONE_CANCER, " n =", nrow(dt), "\n")
}

# ---- compute Wilcoxon p per cancer ---------------------------
cat("[2/4] Computing per-cancer Wilcoxon p-values...\n")
stats <- dt[, {
    n_tum <- sum(sample_group == "Tumor")
    n_nor <- sum(sample_group == "Normal")
    if (n_nor >= 3 && n_tum >= 3) {
        pv <- tryCatch(
            wilcox.test(GATA2 ~ sample_group, data = .SD)$p.value,
            error = function(e) NA_real_
        )
    } else {
        pv <- NA_real_
    }
    .(n_normal = n_nor, n_tumor = n_tum, p_value = pv)
}, by = cancer_type][order(cancer_type)]

stats[, sig_label := fcase(
    is.na(p_value),    "ns",
    p_value < 0.0001,  "****",
    p_value < 0.001,   "***",
    p_value < 0.01,    "**",
    p_value < 0.05,    "*",
    default            = "ns"
)]
print(stats)
fwrite(stats, file.path(RESULTS, "Fig1_stats.csv"))
cat("  saved:", file.path(RESULTS, "Fig1_stats.csv"), "\n")

# ---- order cancers alphabetically (or by n) ------------------
dt[, cancer_type := factor(cancer_type, levels = sort(unique(cancer_type)))]

# ---- build per-cancer x-axis labels showing sample sizes -----
# format: "BLCA\nn=19/408"  (Normal / Tumor)
x_lab_map <- setNames(
    sprintf("%s\nn=%d/%d", stats$cancer_type, stats$n_normal, stats$n_tumor),
    stats$cancer_type
)

# ---- plot ----------------------------------------------------
cat("[3/4] Plotting...\n")

# significance label positions (top of each cancer's data)
y_max <- max(dt$GATA2, na.rm = TRUE)
label_y <- y_max * 1.08

stats_lab <- stats[, .(cancer_type, sig_label)]
stats_lab[, cancer_type := factor(cancer_type, levels = levels(dt$cancer_type))]

# Nature Publishing Group palette (blue/yellow contrast)
pal_2 <- c("Normal" = "#4DBBD5", "Tumor" = "#E64B35")  # NPG-like

p <- ggplot(dt, aes(x = cancer_type, y = GATA2, fill = sample_group)) +
    geom_boxplot(
        outlier.size  = 0.4,
        outlier.alpha = 0.5,
        lwd           = 0.3,
        width         = 0.7,
        position      = position_dodge(width = 0.8)
    ) +
    geom_text(
        data = stats_lab,
        aes(x = cancer_type, y = label_y, label = sig_label),
        inherit.aes = FALSE,
        size = 3.2, fontface = "bold", color = "black"
    ) +
    scale_fill_manual(values = pal_2, name = "Tissue") +
    scale_x_discrete(labels = x_lab_map) +
    scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.12)),
        breaks = scales::pretty_breaks(n = 6)
    ) +
    labs(
        x = NULL,
        y = expression(italic("GATA2") ~ "expression  " ~ log[2] * "(norm count + 1)"),
        title = NULL
    ) +
    theme_pubr(base_size = 11) +
    theme(
        axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                        size = 8.5, face = "bold",
                                        lineheight = 0.9),
        axis.text.y      = element_text(size = 9),
        axis.title.y     = element_text(size = 11, face = "bold"),
        legend.position  = "top",
        legend.title     = element_text(face = "bold", size = 10),
        legend.text      = element_text(size = 10),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        panel.grid.major.x = element_blank()
    )

# ---- save ----------------------------------------------------
cat("[4/4] Saving...\n")
suffix <- if (is.na(ONE_CANCER)) "" else paste0("_", ONE_CANCER)
pdf_out <- file.path(FIGDIR, paste0("Fig1_expression", suffix, ".pdf"))
png_out <- file.path(FIGDIR, paste0("Fig1_expression", suffix, ".png"))

w <- if (is.na(ONE_CANCER)) 11 else 4
h <- 5

ggsave(pdf_out, p, width = w, height = h, device = cairo_pdf)
ggsave(png_out, p, width = w, height = h, dpi = 300, bg = "white")

cat("\n[saved]\n")
cat("  ", pdf_out, "\n")
cat("  ", png_out, "\n")
cat("\nSig legend: ns=p>0.05, *=p<0.05, **=p<0.01, ***=p<0.001, ****=p<0.0001\n")
