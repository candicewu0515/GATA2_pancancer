#!/usr/bin/env Rscript
# =============================================================
# 03_clean_samples.R
# Build a clean per-sample metadata table for TCGA pan-cancer.
#
# Inputs (from data/01_raw/xena_pancan/):
#   - TCGA_phenotype_denseDataOnlyDownload.tsv.gz
#       cols: sample, sample_type, sample_type_id, _primary_disease
#   - Survival_SupplementalTable_S1_20171025_xena_sp
#       cols: sample, _PATIENT, cancer type abbreviation, OS, OS.time,
#             DSS, DSS.time, DFI, DFI.time, PFI, PFI.time,
#             age_at_initial_pathologic_diagnosis, gender, ajcc_pathologic_tumor_stage, ...
#
# Output:
#   data/02_clean/sample_metadata.rds   <- one row per sample
#   data/02_clean/sample_metadata.csv   <- human-readable copy
#
# Usage:
#   Rscript scripts/03_clean_samples.R
# =============================================================

suppressPackageStartupMessages({
    library(data.table)
    library(stringr)
})

# ---- paths ---------------------------------------------------
DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "./data")
RAW   <- file.path(DATA_ROOT, "01_raw", "xena_pancan")
CLEAN <- file.path(DATA_ROOT, "02_clean")
dir.create(CLEAN, recursive = TRUE, showWarnings = FALSE)

f_pheno <- file.path(RAW, "TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
f_surv  <- file.path(RAW, "Survival_SupplementalTable_S1_20171025_xena_sp")

stopifnot(
    "phenotype file missing"  = file.exists(f_pheno),
    "survival file missing"   = file.exists(f_surv)
)

# ---- 1. phenotype (sample type, primary disease) -------------
cat("[1/4] Reading phenotype...\n")
ph <- fread(f_pheno)
cat("  rows:", nrow(ph), "  cols:", ncol(ph), "\n")
cat("  columns:", paste(colnames(ph), collapse = " | "), "\n")
# expected columns: sample, sample_type, sample_type_id, _primary_disease
setnames(ph,
    old = c("sample_type", "sample_type_id", "_primary_disease"),
    new = c("sample_type", "sample_type_id", "primary_disease"),
    skip_absent = TRUE
)

# ---- 2. survival + clinical ---------------------------------
cat("[2/4] Reading survival...\n")
sv <- fread(f_surv)
cat("  rows:", nrow(sv), "  cols:", ncol(sv), "\n")
# column name in this file is literally "cancer type abbreviation" (with spaces)
abbr_col <- grep("cancer.type.abbreviation|cancer type abbreviation",
                 colnames(sv), value = TRUE, ignore.case = TRUE)[1]
setnames(sv, abbr_col, "cancer_type")

# normalize a few common column names
rename_if <- function(dt, from, to) {
    if (from %in% colnames(dt)) setnames(dt, from, to)
}
rename_if(sv, "_PATIENT", "patient_id")
rename_if(sv, "age_at_initial_pathologic_diagnosis", "age")
rename_if(sv, "ajcc_pathologic_tumor_stage", "stage_raw")

# ---- 3. merge ------------------------------------------------
cat("[3/4] Merging phenotype + survival...\n")
meta <- merge(ph, sv, by = "sample", all = TRUE)
cat("  merged rows:", nrow(meta), "\n")

# fill missing cancer_type from primary_disease where possible
# (some samples in phenotype lack survival rows). Use an explicit mapping
# from primary_disease -> TCGA cancer abbreviation; the old toupper(substr(...,1,4))
# produced bogus codes like STOM/LUNG/UTER that don't exist in TCGA.
pd2abbr <- c(
    "brain lower grade glioma"               = "LGG",
    "breast invasive carcinoma"              = "BRCA",
    "cervical & endocervical cancer"         = "CESC",
    "colon adenocarcinoma"                   = "COAD",
    "glioblastoma multiforme"                = "GBM",
    "lung adenocarcinoma"                    = "LUAD",
    "lung squamous cell carcinoma"           = "LUSC",
    "ovarian serous cystadenocarcinoma"      = "OV",
    "rectum adenocarcinoma"                  = "READ",
    "skin cutaneous melanoma"                = "SKCM",
    "stomach adenocarcinoma"                 = "STAD",
    "testicular germ cell tumor"             = "TGCT",
    "uterine corpus endometrioid carcinoma"  = "UCEC"
)
meta[is.na(cancer_type) & !is.na(primary_disease),
     cancer_type := pd2abbr[primary_disease]]

# warn if any primary_disease strings still went unmapped
unmapped <- meta[is.na(cancer_type) & !is.na(primary_disease), unique(primary_disease)]
if (length(unmapped)) {
    cat("  WARN: unmapped primary_disease values (cancer_type left NA):\n")
    for (s in unmapped) cat("    -", s, "\n")
}

# simplify sample group (tumor / normal)
meta[, sample_group := fcase(
    sample_type %in% c("Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood",
                       "Additional - New Primary", "Primary Blood Derived Cancer - Bone Marrow"),
    "Tumor",
    sample_type %in% c("Solid Tissue Normal", "Buccal Cell Normal",
                       "Blood Derived Normal", "Bone Marrow Normal"),
    "Normal",
    sample_type %in% c("Metastatic", "Additional Metastatic"),
    "Metastatic",
    sample_type %in% c("Recurrent Tumor", "Recurrent Blood Derived Cancer - Bone Marrow",
                       "Recurrent Blood Derived Cancer - Peripheral Blood"),
    "Recurrent",
    default = NA_character_
)]

# clean stage: extract "I", "II", "III", "IV" if present
meta[, stage := fcase(
    grepl("Stage IV",  stage_raw), "IV",
    grepl("Stage III", stage_raw), "III",
    grepl("Stage II",  stage_raw), "II",
    grepl("Stage I",   stage_raw), "I",
    default = NA_character_
)]

# ---- 4. summary + save --------------------------------------
cat("[4/4] Summary by cancer type:\n")
summary_tbl <- meta[!is.na(cancer_type),
    .(N_total   = .N,
      N_tumor   = sum(sample_group == "Tumor",   na.rm = TRUE),
      N_normal  = sum(sample_group == "Normal",  na.rm = TRUE),
      N_metast  = sum(sample_group == "Metastatic", na.rm = TRUE)),
    by = cancer_type
][order(-N_total)]
print(summary_tbl)

cat("\n  total cancer types:", uniqueN(meta$cancer_type, na.rm = TRUE), "\n")
cat("  total tumor samples:",  meta[sample_group == "Tumor",  .N], "\n")
cat("  total normal samples:", meta[sample_group == "Normal", .N], "\n")

out_rds <- file.path(CLEAN, "sample_metadata.rds")
out_csv <- file.path(CLEAN, "sample_metadata.csv")
saveRDS(meta, out_rds)
fwrite(meta, out_csv)
cat("\n[saved]\n  ", out_rds, "\n  ", out_csv, "\n")

# also save the summary
fwrite(summary_tbl, file.path(CLEAN, "sample_count_by_cancer.csv"))
cat("  ", file.path(CLEAN, "sample_count_by_cancer.csv"), "\n")
