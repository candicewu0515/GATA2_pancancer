#!/usr/bin/env Rscript
# =============================================================
# 11_fig7_nomogram.R
# Figure 7: Nomogram + calibration curves for OS prediction
#           using GATA2 + clinical variables.
#
# Inputs:
#   data/02_clean/GATA2_expression.rds
#
# Outputs:
#   figures/Fig7A_nomogram_<CANCER>.pdf
#   figures/Fig7B_calibration_<CANCER>.pdf/.png
#   data/03_results/Fig7_cindex_<CANCER>.csv
#
# Strategy:
#   - Build a multivariate Cox model with GATA2 + age + gender + stage
#   - Use rms::nomogram for the nomogram
#   - Bootstrap calibration at 3-yr and 5-yr OS
#
# Usage:
#   Rscript scripts/11_fig7_nomogram.R               # default = LIHC
#   Rscript scripts/11_fig7_nomogram.R LIHC MESO     # any cancer(s)
#
# Note on cancer choice:
#   - LIHC: best (n=370, stage available, GATA2 protective HR=0.85*)
#   - LGG : no stage, use age + grade only
#   - MESO: small n, only n=87, stage often missing
#   - COAD: stage available, GATA2 HR=1.17*
# =============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(survival)
    library(rms)
    library(ggplot2)
    library(ggpubr)
})

# ---- args ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
FOCUS <- if (length(args) >= 1) args else c("LIHC")
cat("Focus cancers:", paste(FOCUS, collapse = ", "), "\n\n")

# ---- paths ---------------------------------------------------
DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
CLEAN     <- file.path(DATA_ROOT, "02_clean")
RESULTS   <- file.path(DATA_ROOT, "03_results")
FIGDIR    <- "./figures"
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR,  recursive = TRUE, showWarnings = FALSE)

# ---- load --------------------------------------------------
cat("[1/3] Loading data...\n")
dt <- readRDS(file.path(CLEAN, "GATA2_expression.rds"))
dt <- dt[sample_group == "Tumor" & !is.na(cancer_type)]
dt[, time_mo := OS.time / 30.44]
dt[, event   := as.integer(OS)]

run_one <- function(cc) {
    cat("\n=== ", cc, " ===\n", sep = "")
    d <- copy(dt[cancer_type == cc])

    # build feature set, drop rows with any missing
    d <- d[!is.na(GATA2) & !is.na(time_mo) & !is.na(event) & time_mo > 0]
    feats <- "GATA2"
    if ("age" %in% colnames(d)) {
        d <- d[!is.na(age)]
        feats <- c(feats, "age")
    }
    if ("gender" %in% colnames(d)) {
        d[, gender := as.factor(gender)]
        d <- d[!is.na(gender) & gender %in% c("MALE", "FEMALE")]
        d[, gender := droplevels(gender)]
        if (nlevels(d$gender) >= 2) feats <- c(feats, "gender")
    }
    if ("stage" %in% colnames(d) && !all(is.na(d$stage))) {
        d[, stage := factor(stage, levels = c("I", "II", "III", "IV"))]
        d <- d[!is.na(stage)]
        if (nlevels(droplevels(d$stage)) >= 2) {
            d[, stage := droplevels(stage)]
            feats <- c(feats, "stage")
        }
    }

    cat("  n =", nrow(d), "; events =", sum(d$event),
        "; features =", paste(feats, collapse = ", "), "\n")
    if (nrow(d) < 30 || sum(d$event) < 10) {
        cat("  [skip] insufficient events\n")
        return(NULL)
    }

    # rms datadist
    assign("dd", datadist(d), envir = .GlobalEnv)
    options(datadist = "dd")
    fml <- as.formula(paste("Surv(time_mo, event) ~", paste(feats, collapse = " + ")))
    fit <- cph(fml, data = d, x = TRUE, y = TRUE, surv = TRUE,
               time.inc = 36)   # 36 mo = 3 yr

    # nomogram
    surv_obj <- Survival(fit)
    surv3 <- function(x) surv_obj(36, x)
    surv5 <- function(x) surv_obj(60, x)

    nom <- nomogram(
        fit, fun = list(surv3, surv5),
        funlabel = c("3-year OS prob", "5-year OS prob"),
        lp = TRUE
    )

    # save nomogram
    pdf(file.path(FIGDIR, paste0("Fig7A_nomogram_", cc, ".pdf")),
        width = 9, height = 6)
    plot(nom, xfrac = 0.25, cex.axis = 0.7, cex.var = 0.9,
         col.grid = c("grey80", "grey95"))
    title(paste0(cc, ": OS nomogram (GATA2 + clinical)"),
          line = 2.5, cex.main = 1)
    dev.off()

    # calibration at 3-yr and 5-yr
    pdf(file.path(FIGDIR, paste0("Fig7B_calibration_", cc, ".pdf")),
        width = 6, height = 6)
    par(mar = c(4.2, 4.2, 2, 1))

    # need refit with time.inc for each horizon
    fit3 <- cph(fml, data = d, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)
    cal3 <- calibrate(fit3, cmethod = "KM", method = "boot",
                      u = 36, m = max(20, floor(nrow(d) / 4)), B = 500)

    fit5 <- cph(fml, data = d, x = TRUE, y = TRUE, surv = TRUE, time.inc = 60)
    cal5 <- calibrate(fit5, cmethod = "KM", method = "boot",
                      u = 60, m = max(20, floor(nrow(d) / 4)), B = 500)

    plot(cal3, lwd = 2, lty = 1, col = "#4DBBD5",
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Nomogram-predicted OS",
         ylab = "Observed OS",
         subtitles = FALSE,
         errbar.col = "#4DBBD5")
    par(new = TRUE)
    plot(cal5, lwd = 2, lty = 1, col = "#E64B35",
         xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "",
         subtitles = FALSE, axes = FALSE,
         errbar.col = "#E64B35")
    abline(0, 1, lty = 2, col = "grey40")
    legend("bottomright",
           legend = c("3-yr", "5-yr", "Ideal"),
           col    = c("#4DBBD5", "#E64B35", "grey40"),
           lty    = c(1, 1, 2), lwd = c(2, 2, 1),
           bty = "n", cex = 0.9)
    title(paste0(cc, ": calibration curve"), cex.main = 1)
    dev.off()

    # also save calibration as png
    png(file.path(FIGDIR, paste0("Fig7B_calibration_", cc, ".png")),
        width = 6 * 300, height = 6 * 300, res = 300, bg = "white")
    par(mar = c(4.2, 4.2, 2, 1))
    plot(cal3, lwd = 2, lty = 1, col = "#4DBBD5",
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Nomogram-predicted OS",
         ylab = "Observed OS",
         subtitles = FALSE, errbar.col = "#4DBBD5")
    par(new = TRUE)
    plot(cal5, lwd = 2, lty = 1, col = "#E64B35",
         xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "",
         subtitles = FALSE, axes = FALSE,
         errbar.col = "#E64B35")
    abline(0, 1, lty = 2, col = "grey40")
    legend("bottomright",
           legend = c("3-yr", "5-yr", "Ideal"),
           col    = c("#4DBBD5", "#E64B35", "grey40"),
           lty    = c(1, 1, 2), lwd = c(2, 2, 1),
           bty = "n", cex = 0.9)
    title(paste0(cc, ": calibration curve"), cex.main = 1)
    dev.off()

    # C-index
    ci <- rcorrcens(Surv(time_mo, event) ~ predict(fit), data = d)
    cind <- 1 - ci[1, "C"]   # rcorrcens reports Somers' Dxy = 2*(C-0.5), and column "C"
    # safer: get from summary
    csum <- summary(fit)
    # use survConcordance
    sc <- survConcordance(Surv(time_mo, event) ~ predict(fit), data = d)
    cind <- sc$concordance

    cidx <- data.table(
        cancer    = cc,
        n         = nrow(d),
        events    = sum(d$event),
        features  = paste(feats, collapse = ","),
        C_index   = round(cind, 3),
        C_se      = round(sc$std.err, 3)
    )
    fwrite(cidx, file.path(RESULTS, paste0("Fig7_cindex_", cc, ".csv")))
    cat("  C-index =", round(cind, 3), "± ", round(sc$std.err, 3), "\n")
    cidx
}

# ---- run all ------------------------------------------------
all_ci <- list()
for (cc in FOCUS) {
    res <- tryCatch(run_one(cc),
                    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
    if (!is.null(res)) all_ci[[cc]] <- res
}

# combined C-index table
if (length(all_ci) > 0) {
    all_dt <- rbindlist(all_ci)
    fwrite(all_dt, file.path(RESULTS, "Fig7_cindex_all.csv"))
    cat("\nC-index summary:\n")
    print(all_dt)
}

cat("\n[saved]\n")
for (cc in FOCUS) {
    cat("  ", file.path(FIGDIR, paste0("Fig7A_nomogram_", cc, ".pdf")), "\n")
    cat("  ", file.path(FIGDIR, paste0("Fig7B_calibration_", cc, ".pdf")), "\n")
}
