#!/usr/bin/env Rscript
# =============================================================
# 07_fig3_survival.R
# Figure 3: Survival analysis of GATA2 across pan-cancer.
#   Panel A: Univariate Cox forest plot for OS across all 33 cancers
#   Panels B-D: KM curves for selected focus cancers (default PRAD/COAD/KIRC)
#
# Input:
#   data/02_clean/GATA2_expression.rds   (from 04_extract_GATA2.R)
#
# Output:
#   figures/Fig3A_cox_forest.pdf / .png
#   figures/Fig3BCD_km_curves.pdf / .png
#   figures/Fig3_combined.pdf
#   data/03_results/Fig3_cox_results.csv
#   data/03_results/Fig3_km_logrank.csv
#
# Usage:
#   Rscript scripts/07_fig3_survival.R
#   Rscript scripts/07_fig3_survival.R PRAD COAD KIRC LIHC   # custom focus cancers
#   Rscript scripts/07_fig3_survival.R --endpoint=PFI         # use PFI instead of OS
# =============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(survival)
    library(survminer)
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
})

# ---- args ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# parse --endpoint=XXX
endpoint_arg <- grep("^--endpoint=", args, value = TRUE)
ENDPOINT <- if (length(endpoint_arg) > 0)
    sub("^--endpoint=", "", endpoint_arg[1]) else "OS"
args <- setdiff(args, endpoint_arg)

FOCUS <- if (length(args) >= 1) args else c("PRAD", "COAD", "KIRC")
stopifnot(ENDPOINT %in% c("OS", "PFI", "DSS", "DFI"))
cat("Endpoint:", ENDPOINT, "\n")
cat("Focus cancers:", paste(FOCUS, collapse = ", "), "\n\n")

# ---- paths ---------------------------------------------------
DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
CLEAN   <- file.path(DATA_ROOT, "02_clean")
RESULTS <- file.path(DATA_ROOT, "03_results")
FIGDIR  <- "./figures"
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR,  showWarnings = FALSE, recursive = TRUE)

# ---- load ----------------------------------------------------
cat("[1/5] Loading data...\n")
dt <- readRDS(file.path(CLEAN, "GATA2_expression.rds"))
# tumor samples only with valid survival
time_col  <- paste0(ENDPOINT, ".time")
event_col <- ENDPOINT
dt <- dt[sample_group == "Tumor" &
         !is.na(cancer_type) &
         !is.na(get(time_col)) &
         !is.na(get(event_col))]
# convert time from days to months
dt[, time_mo := get(time_col) / 30.44]
dt[, event   := as.integer(get(event_col))]
cat("  total tumor samples with", ENDPOINT, ":", nrow(dt), "\n")

# ---- 2. Cox per cancer ---------------------------------------
cat("[2/5] Running Cox regression per cancer...\n")

cancers <- sort(unique(dt$cancer_type))
cox_tbl <- rbindlist(lapply(cancers, function(cc) {
    d <- dt[cancer_type == cc]
    if (nrow(d) < 20 || sum(d$event) < 5) return(NULL)
    fit <- tryCatch(
        coxph(Surv(time_mo, event) ~ GATA2, data = d),
        error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    s <- summary(fit)
    data.table(
        cancer = cc,
        n      = nrow(d),
        events = sum(d$event),
        HR     = s$coefficients[1, "exp(coef)"],
        HR_lo  = s$conf.int[1, "lower .95"],
        HR_hi  = s$conf.int[1, "upper .95"],
        p      = s$coefficients[1, "Pr(>|z|)"]
    )
}))
cox_tbl[, sig := fcase(
    p < 0.001, "***",
    p < 0.01,  "**",
    p < 0.05,  "*",
    default    = ""
)]
cox_tbl <- cox_tbl[order(p)]
fwrite(cox_tbl, file.path(RESULTS, "Fig3_cox_results.csv"))
cat("  cancers analyzed:", nrow(cox_tbl), "\n")
cat("  significant (p<0.05):", sum(cox_tbl$p < 0.05), "\n")
print(cox_tbl[p < 0.05])

# ---- 3. Fig 3A: forest plot ----------------------------------
cat("\n[3/5] Plotting Cox forest...\n")
fp <- copy(cox_tbl)
setorder(fp, cancer)
fp[, cancer := factor(cancer, levels = rev(cancer))]   # top-to-bottom alphabetical

# protective vs risk color
fp[, role := fifelse(HR > 1, "Risk", "Protective")]
fp[p >= 0.05, role := "n.s."]

# clip HR for plotting
fp[, HR_lo_clip := pmax(HR_lo, 0.1)]
fp[, HR_hi_clip := pmin(HR_hi, 10)]

p_forest <- ggplot(fp, aes(y = cancer, x = HR, color = role)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = HR_lo_clip, xmax = HR_hi_clip),
                   height = 0.25, linewidth = 0.5) +
    geom_point(size = 2.5) +
    geom_text(aes(x = max(HR_hi_clip) * 1.5,
                  label = sprintf("%.2f (%.2f-%.2f) %s",
                                  HR, HR_lo, HR_hi, sig)),
              color = "black", hjust = 0, size = 2.8) +
    scale_x_log10(limits = c(0.1, max(fp$HR_hi_clip) * 6),
                  breaks = c(0.25, 0.5, 1, 2, 4)) +
    scale_color_manual(values = c("Risk"        = "#E64B35",
                                  "Protective"  = "#4DBBD5",
                                  "n.s."        = "grey60")) +
    labs(
        x = paste0("Hazard ratio (", ENDPOINT, ", log scale)"),
        y = NULL, color = NULL,
        title = paste0("Univariate Cox regression: GATA2 vs ", ENDPOINT)
    ) +
    theme_pubr(base_size = 10) +
    theme(
        plot.title       = element_text(face = "bold", size = 12, hjust = 0),
        axis.text.y      = element_text(face = "bold", size = 9),
        legend.position  = "top",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.4),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
    )

ggsave(file.path(FIGDIR, "Fig3A_cox_forest.pdf"),
       p_forest, width = 7, height = 8, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig3A_cox_forest.png"),
       p_forest, width = 7, height = 8, dpi = 300, bg = "white")

# ---- 4. KM curves for focus cancers --------------------------
cat("[4/5] Plotting KM curves for", paste(FOCUS, collapse = ", "), "...\n")

km_results <- data.table()

km_plot_one <- function(cc) {
    d <- dt[cancer_type == cc]
    if (nrow(d) < 20) {
        cat("  skip", cc, ": n =", nrow(d), "\n")
        return(NULL)
    }
    # split by median
    cutoff <- median(d$GATA2, na.rm = TRUE)
    d[, group := factor(ifelse(GATA2 > cutoff, "High", "Low"),
                        levels = c("Low", "High"))]

    fit <- survfit(Surv(time_mo, event) ~ group, data = d)
    sd <- survdiff(Surv(time_mo, event) ~ group, data = d)
    p_lr <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)

    km_results <<- rbind(km_results, data.table(
        cancer    = cc,
        n_high    = sum(d$group == "High"),
        n_low     = sum(d$group == "Low"),
        events_hi = sum(d$event[d$group == "High"]),
        events_lo = sum(d$event[d$group == "Low"]),
        logrank_p = p_lr
    ))

    g <- ggsurvplot(
        fit, data = d,
        palette         = c("#4DBBD5", "#E64B35"),
        risk.table      = TRUE,
        risk.table.height = 0.25,
        risk.table.fontsize = 3,
        pval            = sprintf("Log-rank p = %.3g", p_lr),
        pval.coord      = c(2, 0.15),
        pval.size       = 3.5,
        conf.int        = FALSE,
        legend.title    = paste0(cc, " — GATA2"),
        legend.labs     = c("Low", "High"),
        legend          = c(0.8, 0.85),
        xlab            = "Time (months)",
        ylab            = paste0(ENDPOINT, " probability"),
        ggtheme         = theme_pubr(base_size = 10) +
                          theme(plot.background = element_rect(fill = "white", color = NA),
                                panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.4))
    )
    g
}

km_list <- lapply(FOCUS, km_plot_one)
km_list <- km_list[!sapply(km_list, is.null)]

# save individual KM panels
for (i in seq_along(km_list)) {
    cc <- FOCUS[i]
    out_pdf <- file.path(FIGDIR, paste0("Fig3_KM_", cc, ".pdf"))
    pdf(out_pdf, width = 5, height = 5.5)
    print(km_list[[i]])
    dev.off()
}

# combine into one row
combined_km <- arrange_ggsurvplots(
    km_list, ncol = length(km_list), nrow = 1,
    print = FALSE
)
ggsave(file.path(FIGDIR, "Fig3BCD_km_curves.pdf"),
       combined_km, width = 5 * length(km_list), height = 5.5, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig3BCD_km_curves.png"),
       combined_km, width = 5 * length(km_list), height = 5.5, dpi = 300, bg = "white")

fwrite(km_results, file.path(RESULTS, "Fig3_km_logrank.csv"))

# ---- 5. summary ----------------------------------------------
cat("\n[5/5] Done.\n")
cat("\n[saved]\n")
cat("  ", file.path(FIGDIR,  "Fig3A_cox_forest.pdf"),  "\n")
cat("  ", file.path(FIGDIR,  "Fig3BCD_km_curves.pdf"), "\n")
cat("  ", file.path(RESULTS, "Fig3_cox_results.csv"),  "\n")
cat("  ", file.path(RESULTS, "Fig3_km_logrank.csv"),   "\n")

cat("\nLog-rank summary:\n")
print(km_results)
