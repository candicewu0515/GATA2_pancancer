#!/usr/bin/env Rscript
# =============================================================
# 04_extract_GATA2.R
# Extract GATA2 expression from the PANCAN matrix and merge
# with sample metadata.
#
# Input:
#   data/01_raw/xena_pancan/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz
#   data/02_clean/sample_metadata.rds
#
# Output:
#   data/02_clean/GATA2_expression.rds   <- one row per sample
#   data/02_clean/GATA2_expression.csv
#
# Notes on the expression matrix:
#   - rows = genes (gene symbols, e.g. "GATA2")
#   - cols = TCGA barcodes (e.g. "TCGA-XX-XXXX-01")
#   - values = log2(norm_count + 1), batch-corrected (Hoadley 2018)
#
# Usage:
#   Rscript scripts/04_extract_GATA2.R
#   Rscript scripts/04_extract_GATA2.R GENE_OF_INTEREST     # any other gene
# =============================================================

suppressPackageStartupMessages({
    library(data.table)
})

# ---- args ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
GENE <- if (length(args) >= 1) args[1] else "GATA2"
cat("Target gene:", GENE, "\n")

# ---- paths ---------------------------------------------------
DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
RAW   <- file.path(DATA_ROOT, "01_raw", "xena_pancan")
CLEAN <- file.path(DATA_ROOT, "02_clean")

f_expr <- file.path(RAW, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")
f_meta <- file.path(CLEAN, "sample_metadata.rds")

stopifnot(
    "expression matrix missing — run 01_download_data.sh first" = file.exists(f_expr),
    "sample metadata missing — run 03_clean_samples.R first"    = file.exists(f_meta)
)

# ---- 1. find GATA2 row without loading the whole matrix ------
cat("[1/4] Scanning expression matrix header...\n")
hdr <- fread(f_expr, nrows = 0)
n_samples <- ncol(hdr) - 1
cat("  samples in matrix:", n_samples, "\n")

cat("[2/4] Reading first column to locate", GENE, "...\n")
# read only the first column (sample IDs are header, gene IDs are first column)
first_col <- fread(f_expr, select = 1, header = TRUE)
gene_ids <- first_col[[1]]
cat("  total genes:", length(gene_ids), "\n")

# matrix can store either "GATA2" or "GATA2|2624" (Entrez) — match both
hit_idx <- which(gene_ids == GENE | grepl(paste0("^", GENE, "\\|"), gene_ids))
if (length(hit_idx) == 0) {
    stop("Gene '", GENE, "' not found. First 5 IDs:\n  ",
         paste(head(gene_ids, 5), collapse = "\n  "))
}
if (length(hit_idx) > 1) {
    cat("  WARNING: multiple rows matched, using first:\n   ",
        paste(gene_ids[hit_idx], collapse = ", "), "\n")
    hit_idx <- hit_idx[1]
}
cat("  found at row", hit_idx, ":", gene_ids[hit_idx], "\n")

# ---- 2. read just that row -----------------------------------
cat("[3/4] Reading row", hit_idx, "...\n")
# trick: read the whole file but skip to row, only one row
# data.table::fread can subset rows via system call to head/tail, but for safety
# we just read the line directly via `readLines` on the gz connection
con <- gzfile(f_expr, "r")
on.exit(close(con))
# read header to align columns
header_line <- readLines(con, n = 1)
header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
sample_ids <- header[-1]   # drop the "sample" / "gene" first label

# skip rows up to hit_idx - 1 (we already read the header, so gene rows start at line 2)
if (hit_idx > 1) {
    invisible(readLines(con, n = hit_idx - 1))
}
gene_line <- readLines(con, n = 1)
parts <- strsplit(gene_line, "\t", fixed = TRUE)[[1]]
stopifnot(parts[1] == gene_ids[hit_idx])

expr_vals <- as.numeric(parts[-1])
stopifnot(length(expr_vals) == length(sample_ids))

cat("  values: n=", length(expr_vals),
    " | min=", round(min(expr_vals, na.rm = TRUE), 2),
    " | median=", round(median(expr_vals, na.rm = TRUE), 2),
    " | max=", round(max(expr_vals, na.rm = TRUE), 2),
    " | NAs=", sum(is.na(expr_vals)),
    sep = "")
cat("\n")

# ---- 3. merge with metadata ---------------------------------
cat("[4/4] Merging with sample metadata...\n")
expr_dt <- data.table(sample = sample_ids, expression = expr_vals)
setnames(expr_dt, "expression", GENE)

meta <- readRDS(f_meta)
merged <- merge(expr_dt, meta, by = "sample", all.x = TRUE)
cat("  merged rows:", nrow(merged), "\n")
cat("  matched to a cancer type:",
    sum(!is.na(merged$cancer_type)), "/", nrow(merged), "\n")

# quick sanity check per cancer type
cat("\n  expression by cancer (top 10 by n):\n")
chk <- merged[!is.na(cancer_type),
    .(n = .N,
      med = round(median(get(GENE), na.rm = TRUE), 2),
      mean = round(mean(get(GENE), na.rm = TRUE), 2)),
    by = cancer_type][order(-n)][1:10]
print(chk)

# ---- 4. save -------------------------------------------------
out_rds <- file.path(CLEAN, paste0(GENE, "_expression.rds"))
out_csv <- file.path(CLEAN, paste0(GENE, "_expression.csv"))
saveRDS(merged, out_rds)
fwrite(merged, out_csv)
cat("\n[saved]\n  ", out_rds, "\n  ", out_csv, "\n")
