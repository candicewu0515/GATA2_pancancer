#!/bin/bash
# =============================================================
# 01_download_data.sh   (Mac-compatible, uses curl)
# GATA2 pan-cancer analysis - data download
#
# Usage from project root:
#   bash scripts/01_download_data.sh                # all
#   bash scripts/01_download_data.sh xena_pancan    # PANCAN only
#   bash scripts/01_download_data.sh xena_toil      # TCGA+GTEx Toil
#   bash scripts/01_download_data.sh ccle           # CCLE (manual)
#   bash scripts/01_download_data.sh cellminer
#   bash scripts/01_download_data.sh msigdb
#
# Expected total size: ~4-5 GB.  Time: 30 min - 2 h depending on connection.
# =============================================================

set -euo pipefail

# ---- paths (run from project root) ---------------------------
DATA_ROOT="${DATA_ROOT:-./data}"
RAW="$DATA_ROOT/01_raw"

mkdir -p "$RAW"/{xena_pancan,xena_toil,ccle,cellminer,msigdb,gdc_mc3}

MODULE="${1:-all}"

# ---- prefer curl (Mac default); fall back to wget if available
if command -v curl >/dev/null 2>&1; then
    DLCMD="curl -L --fail --retry 3 --progress-bar -C - -o"   # -C - resumes partial
else
    DLCMD="wget -c -O"
fi

dl () {
    local url="$1"
    local outdir="$2"
    local fname
    fname=$(basename "$url")
    local out="$outdir/$fname"
    if [ -s "$out" ]; then
        echo "[skip] $fname exists ($(du -h "$out" | cut -f1))"
    else
        echo "[get ] $fname"
        $DLCMD "$out" "$url"
        echo "       -> $(du -h "$out" | cut -f1)"
    fi
}

# ============================================================
# 1.1 TCGA Pan-Cancer (PANCAN)
# ============================================================
if [[ "$MODULE" == "all" || "$MODULE" == "xena_pancan" ]]; then
    echo "=== [1.1] TCGA PANCAN (Xena) ==="
    BASE="https://pancanatlas.xenahubs.net/download"
    OUT="$RAW/xena_pancan"

    # NOTE: hit S3 directly with %2B-encoded '+' — the xenahubs.net 302 redirect
    # passes literal '+' through to S3, which interprets it as space and 403s.
    # Inlined so we can keep the on-disk filename with literal '+'.
    _eb_url="https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
    _eb_out="$OUT/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
    if [ -s "$_eb_out" ]; then
        echo "[skip] $(basename "$_eb_out") exists ($(du -h "$_eb_out" | cut -f1))"
    else
        echo "[get ] $(basename "$_eb_out")"
        $DLCMD "$_eb_out" "$_eb_url"
        echo "       -> $(du -h "$_eb_out" | cut -f1)"
    fi
    unset _eb_url _eb_out
    dl "$BASE/TCGA_phenotype_denseDataOnlyDownload.tsv.gz" "$OUT"
    dl "$BASE/Survival_SupplementalTable_S1_20171025_xena_sp" "$OUT"
    dl "$BASE/mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz" "$OUT"
    dl "$BASE/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz" "$OUT"
    dl "$BASE/StemnessScores_RNAexp_20170127.2.tsv.gz" "$OUT"
fi

# ============================================================
# 1.2 TCGA + GTEx Toil
# ============================================================
if [[ "$MODULE" == "all" || "$MODULE" == "xena_toil" ]]; then
    echo "=== [1.2] TCGA + GTEx Toil ==="
    BASE="https://toil-xena-hub.s3.us-east-1.amazonaws.com/download"
    OUT="$RAW/xena_toil"

    dl "$BASE/TcgaTargetGtex_gene_expected_count.gz" "$OUT"
    dl "$BASE/TcgaTargetGTEX_phenotype.txt.gz" "$OUT"
    # probemap sits in a probeMap/ subdir, not under download/ root
    dl "$BASE/probeMap/gencode.v23.annotation.gene.probemap" "$OUT"
fi

# ============================================================
# 1.3 CCLE / DepMap (manual)
# ============================================================
if [[ "$MODULE" == "all" || "$MODULE" == "ccle" ]]; then
    echo "=== [1.3] CCLE / DepMap (manual download required) ==="
    OUT="$RAW/ccle"

    cat <<EOF
[manual step] CCLE/DepMap requires browser download:
  Go to:  https://depmap.org/portal/data_page/?tab=allData
  Click the latest release (e.g. DepMap Public 24Q2 or newer).
  Download these two files:
    - OmicsExpressionProteinCodingGenesTPMLogp1.csv
    - Model.csv
  Save them into: $OUT
EOF
fi

# ============================================================
# 1.4 CellMiner drug sensitivity
# ============================================================
if [[ "$MODULE" == "all" || "$MODULE" == "cellminer" ]]; then
    echo "=== [1.4] CellMiner NCI-60 ==="
    OUT="$RAW/cellminer"
    # CellMiner moved processed datasets to /processeddataset/ and prefixed
    # the RNA-seq composite file with "nci60_"
    BASE="https://discover.nci.nih.gov/cellminer/download/processeddataset"

    dl "$BASE/nci60_RNA__RNA_seq_composite_expression.zip" "$OUT"
    dl "$BASE/DTP_NCI60_ZSCORE.zip" "$OUT"

    cd "$OUT"
    for z in *.zip; do
        [ -f "$z" ] && unzip -n "$z" || true
    done
    cd - >/dev/null
fi

# ============================================================
# 1.5 MSigDB gene sets
# ============================================================
if [[ "$MODULE" == "all" || "$MODULE" == "msigdb" ]]; then
    echo "=== [1.5] MSigDB ==="
    OUT="$RAW/msigdb"
    # MSigDB 2024.1.Hs lays .gmt files flat under the release dir (no symbols/ subdir)
    BASE="https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs"

    dl "$BASE/h.all.v2024.1.Hs.symbols.gmt" "$OUT"
    dl "$BASE/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt" "$OUT"
    dl "$BASE/c2.cp.reactome.v2024.1.Hs.symbols.gmt" "$OUT"
    dl "$BASE/c5.go.bp.v2024.1.Hs.symbols.gmt" "$OUT"
    dl "$BASE/c5.go.cc.v2024.1.Hs.symbols.gmt" "$OUT"
    dl "$BASE/c5.go.mf.v2024.1.Hs.symbols.gmt" "$OUT"
fi

# ============================================================
# 1.6 GDC MC3 full MAF (optional, large)
# ============================================================
if [[ "$MODULE" == "all" || "$MODULE" == "gdc_mc3" ]]; then
    echo "=== [1.6] MC3 full MAF (~1.5 GB) ==="
    OUT="$RAW/gdc_mc3"

    dl "https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc" "$OUT"
    if [ -f "$OUT/1c8cfe5f-e52d-41ba-94da-f15ea1337efc" ]; then
        mv "$OUT/1c8cfe5f-e52d-41ba-94da-f15ea1337efc" "$OUT/mc3.v0.2.8.PUBLIC.maf.gz"
    fi
fi

echo ""
echo "[done] Module: $MODULE"
echo "Data root: $(cd "$DATA_ROOT" && pwd)"
echo ""
echo "Total size so far:"
du -sh "$RAW"/* 2>/dev/null || true
echo ""
echo "Next: Rscript scripts/02_qc_data.R"
