#!/usr/bin/env bash
# =============================================================================
# run_all.sh — Execute the entire bulk RNA-seq reproduction pipeline
# =============================================================================
# Usage:
#   bash run_all.sh                     (uses default conda env: rnaseq_reproduction)
#   ENV_NAME=my_env bash run_all.sh     (override environment name)
#
# Prerequisites:
#   conda env create -f environment.yml   OR  set up the three separate envs
#   (see README.md for details)
# =============================================================================
set -euo pipefail

ENV="${ENV_NAME:-rnaseq_reproduction}"
SCRIPTS_DIR="$(cd "$(dirname "$0")/scripts" && pwd)"

echo "============================================================"
echo " Bulk RNA-Seq Reproduction Pipeline"
echo " Environment: ${ENV}"
echo " Scripts dir: ${SCRIPTS_DIR}"
echo "============================================================"

run_stage() {
    local stage="$1"
    local script="$2"
    echo ""
    echo "──────────────────────────────────────────────────────────"
    echo " Stage ${stage}: ${script}"
    echo "──────────────────────────────────────────────────────────"
    conda run --no-capture-output -n "${ENV}" bash -c "cd '${SCRIPTS_DIR}' && bash '${script}'"
}

run_r_stage() {
    local stage="$1"
    local script="$2"
    echo ""
    echo "──────────────────────────────────────────────────────────"
    echo " Stage ${stage}: ${script}"
    echo "──────────────────────────────────────────────────────────"
    conda run --no-capture-output -n "${ENV}" Rscript "${SCRIPTS_DIR}/${script}"
}

# Stage 0: Resolve sample metadata
run_stage 0 "00_resolve_metadata.sh"

# Stage 1: Download raw reads from NCBI SRA
run_stage 1 "01_download_data.sh"

# Stage 2: Quality control (FastQC + MultiQC)
run_stage 2 "02_run_qc.sh"

# Stage 3: Download Ensembl 110 reference transcriptome & genome
run_stage 3 "03_download_reference.sh"

# Stage 4: Build Salmon index
run_stage 4 "04_build_salmon_index.sh"

# Stage 5: Quantify transcript abundances
run_stage 5 "05_run_salmon_quant.sh"

# Stage 6: Differential expression analysis (DESeq2)
run_r_stage 6 "06_deseq2_analysis.R"

# Stage 7: Generate publication-quality figures
run_r_stage 7 "07_make_figures.R"

echo ""
echo "============================================================"
echo " Pipeline complete."
echo " Results:  results/deseq2/"
echo " Figures:  figures/"
echo "============================================================"
