#!/usr/bin/env bash
# =============================================================================
# 02_run_qc.sh
# -----------------------------------------------------------------------------
# Purpose:
#   Quality-control pass over all 6 raw single-end FASTQ files:
#     1) FastQC  -> per-file HTML/zip reports (adapter content, quality
#                   score distributions, GC content, duplication, etc.)
#     2) MultiQC -> aggregates all FastQC reports into one navigable
#                   summary report for cross-sample comparison.
#
#   The output of this script is the evidence base used to decide whether
#   fastp trimming (Stage 3b, conditional) is actually warranted. Per the
#   project's modern Salmon-based design, light 3'-end quality decay and
#   typical Illumina adapter traces are usually NOT worth trimming, since
#   Salmon's selective alignment already soft-clips poor-quality/adapter
#   bases during mapping. Trimming is only added if FastQC/MultiQC reveal
#   something more serious (e.g. >5% adapter content, severe per-base
#   quality drop-off before position ~50, strong rRNA/poly-A contamination
#   signatures).
#
# Execution context:
#   Run from inside scripts/  ->  bash 02_run_qc.sh
#
# Idempotency:
#   - Skips FastQC for any FASTQ whose <sample>_fastqc.zip already exists
#     in ../results/fastqc/.
#   - Re-runs MultiQC unconditionally (cheap, and picks up new reports).
#
# Requirements:
#   - fastqc, multiqc  -> conda env "rnaseq-tools" (or "metrics" env)
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration (relative paths assume execution from scripts/)
# -----------------------------------------------------------------------------
DATA_DIR="../data"
QC_DIR="../results/fastqc"
MULTIQC_DIR="../results/multiqc"
THREADS=12

mkdir -p "${QC_DIR}" "${MULTIQC_DIR}"

# -----------------------------------------------------------------------------
# Sanity checks
# -----------------------------------------------------------------------------
for tool in fastqc multiqc; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "[ERROR] Required tool '${tool}' not found on PATH." >&2
        echo "[ERROR] Activate the appropriate conda environment first." >&2
        exit 1
    fi
done

shopt -s nullglob
FASTQ_FILES=("${DATA_DIR}"/*.fastq.gz)
shopt -u nullglob

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No .fastq.gz files found in ${DATA_DIR}." >&2
    echo "[ERROR] Run 01_download_data.sh (Stage 2) first." >&2
    exit 1
fi

echo "[INFO] Found ${#FASTQ_FILES[@]} FASTQ files for QC."

# -----------------------------------------------------------------------------
# Step 1: FastQC
# -----------------------------------------------------------------------------
# Build the list of files that still need FastQC (idempotency).
TO_RUN=()
for FQ in "${FASTQ_FILES[@]}"; do
    BASENAME=$(basename "${FQ}" .fastq.gz)
    if [[ -f "${QC_DIR}/${BASENAME}_fastqc.zip" ]]; then
        echo "[SKIP] ${BASENAME}: FastQC report already exists."
    else
        TO_RUN+=("${FQ}")
    fi
done

if [[ ${#TO_RUN[@]} -gt 0 ]]; then
    echo "[INFO] Running FastQC on ${#TO_RUN[@]} file(s) (threads=${THREADS}) ..."
    # -t: parallel threads (1 file per thread up to THREADS)
    # -o: output directory for HTML + zip reports
    fastqc \
        --threads "${THREADS}" \
        --outdir "${QC_DIR}" \
        "${TO_RUN[@]}"
else
    echo "[INFO] All FastQC reports already present - skipping FastQC step."
fi

# -----------------------------------------------------------------------------
# Step 2: MultiQC aggregation
# -----------------------------------------------------------------------------
# MultiQC scans QC_DIR for FastQC zip/html outputs and produces a single
# cross-sample summary report (results/multiqc/multiqc_report.html).
echo "[INFO] Running MultiQC aggregation ..."
multiqc \
    "${QC_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --force

echo "============================================================"
echo "[DONE] Stage 3 QC complete."
echo "[INFO] FastQC reports : ${QC_DIR}/*_fastqc.html"
echo "[INFO] MultiQC report : ${MULTIQC_DIR}/multiqc_report.html"
echo "============================================================"
echo "[ACTION REQUIRED] Open ${MULTIQC_DIR}/multiqc_report.html and check:"
echo "  - 'Sequence Quality Histograms'      (per-base Phred scores)"
echo "  - 'Adapter Content'                  (Illumina Universal Adapter %)"
echo "  - 'Per Sequence GC Content'          (rRNA depletion efficacy)"
echo "  - 'Sequence Duplication Levels'"
echo "If all panels look clean (no red flags, adapter content < ~5%),"
echo "proceed directly to Stage 4 (reference + Salmon index) WITHOUT"
echo "fastp trimming. If issues are found, request the optional"
echo "02b_trim_reads.sh (fastp) script before continuing."
echo "============================================================"
