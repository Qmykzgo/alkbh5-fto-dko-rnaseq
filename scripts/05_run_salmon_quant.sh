#!/usr/bin/env bash
# =============================================================================
# 05_run_salmon_quant.sh
# -----------------------------------------------------------------------------
# Purpose:
#   Quantify transcript abundances for all 6 single-end samples against the
#   decoy-aware Salmon index built in Stage 4, producing one quant.sf per
#   sample under results/salmon/<SRR>/.
#
# Key flags:
#   -l A              : auto-detect library type (fr/unstranded - safe
#                        default for single-end cDNA libraries; Salmon
#                        inspects mapping orientation in the first reads).
#   --validateMappings: enables Salmon's selective-alignment scoring model
#                        (rather than pure pseudo-alignment). This computes
#                        a proper alignment score for each mapping, which
#                        substantially improves quantification accuracy by
#                        rejecting spurious quasi-mappings - the single most
#                        impactful accuracy flag for modern Salmon usage.
#   --gcBias          : models and corrects for fragment-level GC-content
#                        bias introduced during library prep (relevant here:
#                        NEBNext Ultra Directional kit + rRNA depletion are
#                        known to introduce mild GC bias). Produces more
#                        accurate relative abundance estimates, especially
#                        for GC-extreme transcripts.
#   --threads 12      : hard cap per environment constraints.
#
# Execution context:
#   Run from inside scripts/  ->  bash 05_run_salmon_quant.sh
#
# Idempotency:
#   Skips any sample whose results/salmon/<SRR>/quant.sf already exists.
#
# Requirements:
#   - salmon (>= 1.9)  -> conda env "rnaseq-tools"
#   - ../ref/salmon_index/  (Stage 4)
#   - ../data/<SRR>.fastq.gz (Stage 2)
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration (relative paths assume execution from scripts/)
# -----------------------------------------------------------------------------
CONFIG_DIR="../config"
DATA_DIR="../data"
INDEX_DIR="../ref/salmon_index"
RESULTS_DIR="../results/salmon"
DESIGN_MATRIX="${CONFIG_DIR}/design_matrix.csv"
THREADS=12

mkdir -p "${RESULTS_DIR}"

# -----------------------------------------------------------------------------
# Sanity checks
# -----------------------------------------------------------------------------
if ! command -v salmon >/dev/null 2>&1; then
    echo "[ERROR] 'salmon' not found on PATH. Activate the rnaseq-tools env." >&2
    exit 1
fi

if [[ ! -f "${INDEX_DIR}/info.json" ]]; then
    echo "[ERROR] Salmon index not found at ${INDEX_DIR}" >&2
    echo "[ERROR] Run 04_build_salmon_index.sh (Stage 4b) first." >&2
    exit 1
fi

if [[ ! -f "${DESIGN_MATRIX}" ]]; then
    echo "[ERROR] Design matrix not found at ${DESIGN_MATRIX}" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Iterate over every SRR accession in the design matrix (skip header row).
# -----------------------------------------------------------------------------
SAMPLES=$(tail -n +2 "${DESIGN_MATRIX}" | cut -d',' -f1)

for SRR in ${SAMPLES}; do
    FASTQ="${DATA_DIR}/${SRR}.fastq.gz"
    SAMPLE_OUT="${RESULTS_DIR}/${SRR}"
    QUANT_SF="${SAMPLE_OUT}/quant.sf"

    # -------------------------------------------------------------------
    # Idempotency check: skip if quantification already completed.
    # -------------------------------------------------------------------
    if [[ -f "${QUANT_SF}" ]]; then
        echo "[SKIP] ${SRR}: ${QUANT_SF} already exists."
        continue
    fi

    if [[ ! -f "${FASTQ}" ]]; then
        echo "[ERROR] Missing FASTQ for ${SRR}: ${FASTQ}" >&2
        echo "[ERROR] Run 01_download_data.sh (Stage 2) first." >&2
        exit 1
    fi

    echo "============================================================"
    echo "[INFO] Quantifying ${SRR} ..."
    echo "============================================================"

    salmon quant \
        --index "${INDEX_DIR}" \
        --libType A \
        --unmatedReads "${FASTQ}" \
        --threads "${THREADS}" \
        --validateMappings \
        --gcBias \
        --output "${SAMPLE_OUT}"

    echo "[OK] ${SRR} -> ${QUANT_SF}"
done

echo "============================================================"
echo "[DONE] Stage 5 quantification complete."
echo "[INFO] Per-sample logs (mapping rate, library type, etc.) are in:"
for SRR in ${SAMPLES}; do
    echo "  ${RESULTS_DIR}/${SRR}/logs/salmon_quant.log"
done
echo "============================================================"
echo "[ACTION REQUIRED] Check each logs/salmon_quant.log for the"
echo "'Mapping rate' line. For cytoplasmic RNA-Seq against a"
echo "decoy-aware index, expect roughly 80-95% mapping rate."
echo "Mapping rates below ~70% would indicate a contamination or"
echo "library-type mismatch issue worth investigating before"
echo "proceeding to Stage 6 (DESeq2)."
echo "============================================================"
