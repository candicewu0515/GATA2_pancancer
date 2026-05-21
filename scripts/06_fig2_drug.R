#!/usr/bin/env Rscript
# =============================================================
# 06_fig2_drug.R
# Figure 2: GATA2 expression vs drug sensitivity (NCI-60 CellMiner)
#
# Inputs (from data/01_raw/cellminer/):
#   - RNA__RNA_seq_composite_expression.xls(x)  <- after unzip
#   - DTP_NCI60_ZSCORE.xls(x)
#
# Output:
#   figures/Fig2_drug.pdf
#   figures/Fig2_drug.png
#   data/03_results/Fig2_drug_correlations.csv
#
# Usage:
#   Rscript scripts/06_fig2_drug.R                 # default: top 6 |cor|
#   Rscript scripts/06_fig2_drug.R 9               # show top N drugs (must be 4/6/9/12)
#
# Notes:
#   - CellMiner files are .xls (Excel 97-2003) inside the zips.
#   - This script reads them via `readxl::read_excel`. If readxl fails,
#     manually open the file in Excel and re-save as .xlsx.
# =============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(readxl)
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
    library(ggsci)
})

# ---- args ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
TOP_N <- if (length(args) >= 1) as.integer(args[1]) else 6
stopifnot(TOP_N %in% c(4, 6, 9, 12))

# ---- paths ---------------------------------------------------
DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
CM        <- file.path(DATA_ROOT, "01_raw", "cellminer")
RESULTS   <- file.path(DATA_ROOT, "03_results")
FIGDIR    <- "./figures"
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR,  recursive = TRUE, showWarnings = FALSE)

# ---- locate files (handles .xls / .xlsx, varies by release) --
find_one <- function(dir, pattern) {
    hits <- list.files(dir, pattern = pattern, full.names = TRUE,
                       ignore.case = TRUE, recursive = TRUE)
    if (length(hits) == 0) {
        stop("File not found in ", dir, " matching: ", pattern,
             "\n  Did you unzip the CellMiner downloads?")
    }
    hits[1]
}

f_expr <- find_one(CM, "RNA.*composite.*expression.*\\.(xls|xlsx)$")
f_drug <- find_one(CM, "DTP.*NCI60.*ZSCORE.*\\.(xls|xlsx)$")
cat("[1/5] Found CellMiner files:\n")
cat("  expression:", basename(f_expr), "\n")
cat("  drug zscore:", basename(f_drug), "\n")

# ---- 1. read expression -------------------------------------
cat("[2/5] Reading expression matrix...\n")
# CellMiner expression file format:
#   Row 1-10 (~): metadata headers (skip via `skip`)
#   First data row: column = "Gene name d", then 60 cell line names
# We try a few skip values until we find the GATA2 row.

read_cm <- function(path) {
    # try a few starting rows
    for (s in c(10, 9, 8, 11)) {
        d <- tryCatch(
            suppressMessages(read_excel(path, skip = s, .name_repair = "minimal")),
            error = function(e) NULL
        )
        if (!is.null(d) && ncol(d) >= 50) {
            # the first column should be gene names; check if "GATA2" appears
            first_col <- d[[1]]
            if (any(grepl("^GATA2$", first_col))) {
                cat("  parsed with skip =", s, "\n")
                return(as.data.table(d))
            }
        }
    }
    stop("Could not auto-detect header row in: ", path,
         "\n  Open the file manually to inspect.")
}
expr_dt <- read_cm(f_expr)
cat("  rows:", nrow(expr_dt), "  cols:", ncol(expr_dt), "\n")

# pull GATA2 row
gata2_row <- expr_dt[get(colnames(expr_dt)[1]) == "GATA2"][1]
stopifnot("GATA2 row not found in expression file" = nrow(gata2_row) == 1)

# isolate numeric columns (the 60 cell lines)
cell_cols <- colnames(expr_dt)[sapply(expr_dt, is.numeric)]
if (length(cell_cols) < 50) {
    # try coercion: some columns may have come in as text
    expr_dt_num <- copy(expr_dt)
    for (cn in colnames(expr_dt_num)[-1]) {
        v <- suppressWarnings(as.numeric(expr_dt_num[[cn]]))
        if (sum(!is.na(v)) > 50) expr_dt_num[, (cn) := v]
    }
    cell_cols <- colnames(expr_dt_num)[sapply(expr_dt_num, is.numeric)]
    gata2_row <- expr_dt_num[get(colnames(expr_dt_num)[1]) == "GATA2"][1]
}
gata2_expr <- as.numeric(gata2_row[, ..cell_cols])
names(gata2_expr) <- cell_cols
cat("  GATA2 expressed in", sum(!is.na(gata2_expr)), "cell lines\n")

# ---- 2. read drug Z-scores -----------------------------------
cat("[3/5] Reading drug Z-scores...\n")
read_cm_drug <- function(path) {
    # Drug Z-score file has 6 metadata columns then 60 cell lines.
    # readxl infers cell-line columns as character because the data contains
    # the literal string "na" mixed with numbers — so we read, then coerce.
    for (s in c(8, 9, 10, 7, 11)) {
        d <- tryCatch(
            suppressMessages(read_excel(path, skip = s, .name_repair = "minimal")),
            error = function(e) NULL
        )
        if (is.null(d) || ncol(d) < 50) next
        cn <- colnames(d)
        # the header row we want has "Drug name" (or similar) as col 2
        if (any(grepl("Drug.*name", cn, ignore.case = TRUE))) {
            cat("  parsed with skip =", s, "\n")
            d <- as.data.table(d)
            # coerce cell-line columns (everything from col 7 onward, by convention)
            # to numeric; literal "na" / "-" become NA
            cell_cols <- cn[7:length(cn)]
            for (cc in cell_cols) {
                set(d, j = cc,
                    value = suppressWarnings(as.numeric(d[[cc]])))
            }
            return(d)
        }
    }
    stop("Could not auto-detect header row in drug file.")
}
drug_dt <- read_cm_drug(f_drug)
cat("  rows (drugs):", nrow(drug_dt), "  cols:", ncol(drug_dt), "\n")

# identify drug name column — prefer the literal "Drug name" header over the
# NSC number (which is also non-numeric and would otherwise grab position 1)
drug_name_col <- grep("^Drug\\s*name", colnames(drug_dt),
                      ignore.case = TRUE, value = TRUE)[1]
if (is.na(drug_name_col)) {
    # fallback: first text column that isn't NSC / PubChem ID
    txt_cols <- colnames(drug_dt)[!sapply(drug_dt, is.numeric)]
    txt_cols <- txt_cols[!grepl("NSC|PubChem|SMILES|FDA|Mechanism", txt_cols, ignore.case = TRUE)]
    drug_name_col <- txt_cols[1]
}
stopifnot("could not find drug-name column" = !is.na(drug_name_col))
cat("  drug name column:", drug_name_col, "\n")

# numeric columns = cell line z-scores
drug_cell_cols <- colnames(drug_dt)[sapply(drug_dt, is.numeric)]
# match cell line column names between expression and drug files
shared <- intersect(names(gata2_expr), drug_cell_cols)
cat("  cell lines shared expr/drug:", length(shared), "\n")
if (length(shared) < 40) {
    cat("  WARN: few shared cell lines. Column names may differ.\n")
    cat("  expr cols (head):", paste(head(names(gata2_expr), 5), collapse = ", "), "\n")
    cat("  drug cols (head):", paste(head(drug_cell_cols, 5), collapse = ", "), "\n")
}

# ---- 3. compute Pearson correlation for every drug -----------
cat("[4/5] Computing correlations...\n")
gata2_vec <- gata2_expr[shared]
cor_tbl <- data.table(
    drug = drug_dt[[drug_name_col]],
    cor  = NA_real_,
    p    = NA_real_,
    n    = NA_integer_
)
for (i in seq_len(nrow(drug_dt))) {
    dv <- as.numeric(drug_dt[i, ..shared])
    ok <- !is.na(dv) & !is.na(gata2_vec)
    if (sum(ok) < 10) next
    ct <- suppressWarnings(cor.test(gata2_vec[ok], dv[ok], method = "pearson"))
    cor_tbl[i, c("cor", "p", "n") := list(ct$estimate, ct$p.value, sum(ok))]
}

# attach FDA status so we can filter the figure to characterized drugs only
# (the full table is still saved for later inspection)
fda_col <- grep("FDA.*status", colnames(drug_dt), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(fda_col)) {
    cor_tbl[, fda_status := drug_dt[[fda_col]]]
}

cor_tbl <- cor_tbl[!is.na(cor) & !is.na(p)]
cor_tbl <- cor_tbl[order(p)]
fwrite(cor_tbl, file.path(RESULTS, "Fig2_drug_correlations.csv"))
cat("  drugs tested:", nrow(cor_tbl), "\n")
cat("  saved:", file.path(RESULTS, "Fig2_drug_correlations.csv"), "\n\n")

# pick top N significant by absolute correlation, p < 0.05.
# Filter to named/characterized drugs: FDA approved or in clinical trials.
# Skipping this would let unnamed library compounds (drug == "-") dominate.
named <- cor_tbl[
    p < 0.05 &
    drug != "-" &
    (is.na(fda_col) | fda_status %in% c("FDA approved", "Clinical trial"))
]
cat("Named drugs (FDA approved + Clinical trial) tested:", nrow(named), "\n")
sig <- named[order(-abs(cor))][seq_len(min(TOP_N, .N))]
cat("Top significant drugs (|cor|, p<0.05):\n")
print(sig)

# ---- 4. plot scatter panels ----------------------------------
cat("[5/5] Plotting...\n")
make_panel <- function(drug_name, cor_val, pval, n_val) {
    drug_idx <- which(drug_dt[[drug_name_col]] == drug_name)[1]
    dv <- as.numeric(drug_dt[drug_idx, ..shared])
    pd <- data.frame(GATA2 = gata2_vec, drug = dv)

    cor_color <- if (cor_val > 0) "#E64B35" else "#4DBBD5"

    ggplot(pd, aes(x = GATA2, y = drug)) +
        geom_point(size = 1.8, alpha = 0.7, color = cor_color) +
        geom_smooth(method = "lm", formula = y ~ x,
                    color = "grey25", fill = "grey75",
                    linewidth = 0.5, alpha = 0.3, se = TRUE) +
        labs(
            title = drug_name,
            subtitle = sprintf("cor = %.3f,  p = %.3g,  n = %d",
                               cor_val, pval, n_val),
            x = expression(italic("GATA2") ~ "expression"),
            y = "Drug activity (Z-score)"
        ) +
        theme_pubr(base_size = 10) +
        theme(
            plot.title    = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 9, color = "grey30"),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background  = element_rect(fill = "white", color = NA)
        )
}

plots <- mapply(make_panel,
                sig$drug, sig$cor, sig$p, sig$n,
                SIMPLIFY = FALSE)

ncol_grid <- switch(as.character(TOP_N), "4" = 2, "6" = 3, "9" = 3, "12" = 4)
combined <- wrap_plots(plots, ncol = ncol_grid)

w <- ncol_grid * 3.2
h <- ceiling(TOP_N / ncol_grid) * 3

pdf_out <- file.path(FIGDIR, "Fig2_drug.pdf")
png_out <- file.path(FIGDIR, "Fig2_drug.png")
ggsave(pdf_out, combined, width = w, height = h, device = cairo_pdf)
ggsave(png_out, combined, width = w, height = h, dpi = 300, bg = "white")

cat("\n[saved]\n")
cat("  ", pdf_out, "\n")
cat("  ", png_out, "\n")
