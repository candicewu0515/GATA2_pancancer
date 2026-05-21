# GATA2 Pan-Cancer — Mac local setup

## Layout

```
~/Documents/GATA2_pancancer/
├── data/01_raw/          ← downloaded raw data (gitignored)
├── data/02_clean/        ← parsed intermediates
├── data/03_results/      ← numeric output tables
├── scripts/              ← all R / bash code
├── figures/              ← final figures (PDF/PNG)
├── manuscript/           ← paper drafts
└── README.md
```

## One-time setup

### 1. System tools

```bash
# Homebrew (if not installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# unzip is built-in on Mac; install git/wget only if you want them
brew install git
```

### 2. R

Install R 4.4+ from https://cran.r-project.org/bin/macosx/
(or `brew install --cask r`)

RStudio is optional but recommended: https://posit.co/download/rstudio-desktop/

### 3. R packages (one shot, ~15-30 min)

```bash
cd ~/Documents/GATA2_pancancer
Rscript scripts/00_install_packages.R
```

### 4. Create project folders

```bash
cd ~/Documents/GATA2_pancancer
mkdir -p data/{01_raw,02_clean,03_results} figures manuscript scripts
```

## Run order

| Step | Script | What it does | Time |
|------|--------|--------------|------|
| 00 | `00_install_packages.R` | Install all R deps | 15–30 min |
| 01 | `01_download_data.sh` | Download TCGA/GTEx/CCLE/etc | 30 min–2 h |
| 02 | `02_qc_data.R` | Verify file integrity | 5 min |
| 03 | `03_clean_samples.R` | Build sample metadata table | 5 min |
| 04 | `04_extract_GATA2.R` | Extract GATA2 expression vector | 2 min |
| 05 | `05_fig1_expression.R` | Figure 1: tumor vs normal boxplot | 5 min |
| 06 | `06_fig2_drug.R` | Figure 2: CellMiner drug correlations | 5 min |
| 07 | `07_fig3_survival.R` | Figure 3: Cox forest + KM curves | 10 min |
| ... | ... | ... | ... |

Each `0N_*.R` should be runnable on its own — never re-runs the previous step.

## Run modes

Most scripts support single-cancer vs all-cancers:

```bash
Rscript scripts/05_fig1_expression.R              # all 33 cancers
Rscript scripts/05_fig1_expression.R COAD         # single cancer
```

## .gitignore (recommended)

```
data/01_raw/
data/02_clean/*.rds
data/03_results/*.csv
.Rhistory
.RData
.DS_Store
*.log
manuscript/~$*.docx
```

## Disk budget

| Folder | Size |
|---|---|
| data/01_raw | ~4-5 GB |
| data/02_clean | ~500 MB |
| data/03_results | <100 MB |
| figures | <50 MB |

## Memory notes for Mac (16-32 GB RAM)

- Most scripts: fine
- WGCNA on full 5000 genes × thousand samples: needs 16+ GB. If it OOMs, drop to top-3000 variable genes.
- DNA methylation 450K matrix: lazy-load row-by-row with `data.table::fread(select=...)`, don't load whole.

## iCloud warning

If `~/Documents` is iCloud-synced, **exclude this project folder** from iCloud:
- Right-click the project folder in Finder → "Remove Download" then "Don't Sync"
- Or move the project to `~/Code/` or `~/Work/` (outside iCloud)

R writes lots of small intermediate files; iCloud sync will fight with that.
