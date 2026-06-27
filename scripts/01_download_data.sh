#!/usr/bin/env bash
# =============================================================================
# 01_download_data.sh
# -----------------------------------------------------------------------------
# Purpose:
#   Download the 6 raw FASTQ files belonging to BioProject PRJNA813529
#   (GSE198050: HEK293T WT vs HEK293T-dALKBH5-dFTO, single-end, 76bp).
#
#   Pipeline per sample:
#     prefetch (SRA) -> fasterq-dump (decompress to .fastq)
#       -> pigz/gzip (-> .fastq.gz in ../data/)
#       -> delete intermediate .sra and uncompressed .fastq immediately
#
#   This keeps peak disk usage low (~150GB SSD budget) since the
#   uncompressed intermediates never coexist with more than one sample
#   at a time, and the SRA cache under ~/ncbi/public/sra/ is wiped after
#   each sample is processed.
#
# Execution context:
#   Run from inside scripts/  ->  bash 01_download_data.sh
#
# Idempotency:
#   If ../data/<SRR>.fastq.gz already exists, the sample is skipped
#   entirely - safe to re-run after interruption.
#
# Requirements:
#   - sra-tools  (prefetch, fasterq-dump)   -> conda env "rnaseq-tools"
#   - pigz (preferred) or gzip
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration (relative paths assume execution from scripts/)
# -----------------------------------------------------------------------------
CONFIG_DIR="../config"
DATA_DIR="../data"
DESIGN_MATRIX="${CONFIG_DIR}/design_matrix.csv"

THREADS=12                     # Hard cap per environment constraints (16 cores total)
SRA_CACHE_DIR="${HOME}/ncbi/public/sra"   # Default prefetch download location

mkdir -p "${DATA_DIR}"

# -----------------------------------------------------------------------------
# Choose the fastest available gzip implementation.
# pigz parallelizes compression across THREADS; falls back to gzip if absent.
# -----------------------------------------------------------------------------
if command -v pigz >/dev/null 2>&1; then
    COMPRESS_CMD="pigz -p ${THREADS}"
    echo "[INFO] Using pigz (-p ${THREADS}) for compression."
else
    COMPRESS_CMD="gzip"
    echo "[WARN] pigz not found - falling back to single-threaded gzip."
    echo "[WARN] Consider: sudo apt install pigz"
fi

# -----------------------------------------------------------------------------
# Sanity checks
# -----------------------------------------------------------------------------
if [[ ! -f "${DESIGN_MATRIX}" ]]; then
    echo "[ERROR] Design matrix not found at ${DESIGN_MATRIX}" >&2
    echo "[ERROR] Run 00_resolve_metadata.sh (Stage 1) first." >&2
    exit 1
fi

for tool in prefetch fasterq-dump; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "[ERROR] Required tool '${tool}' not found on PATH." >&2
        echo "[ERROR] Activate the sra-tools conda environment first." >&2
        exit 1
    fi
done

# -----------------------------------------------------------------------------
# Iterate over every SRR accession in the design matrix (skip header row).
# -----------------------------------------------------------------------------
SAMPLES=$(tail -n +2 "${DESIGN_MATRIX}" | cut -d',' -f1)

for SRR in ${SAMPLES}; do
    FINAL_GZ="${DATA_DIR}/${SRR}.fastq.gz"

    # -------------------------------------------------------------------
    # Idempotency check: if the compressed FASTQ already exists, skip
    # this sample entirely. Allows safe re-runs after interruption.
    # -------------------------------------------------------------------
    if [[ -f "${FINAL_GZ}" ]]; then
        echo "[SKIP] ${SRR}: ${FINAL_GZ} already exists."
        continue
    fi

    echo "============================================================"
    echo "[INFO] Processing ${SRR}"
    echo "============================================================"

    # ---------------------------------------------------------------
    # Step 1: prefetch - downloads the .sra container to the local
    # SRA cache (~/ncbi/public/sra/<SRR>/<SRR>.sra). Resumable on
    # interruption (prefetch supports partial-download resume).
    # ---------------------------------------------------------------
    echo "[INFO] (1/3) prefetch ${SRR} ..."
    prefetch --max-size u "${SRR}" --output-directory "${SRA_CACHE_DIR}"

    # ---------------------------------------------------------------
    # Step 2: fasterq-dump - extracts FASTQ from the .sra container.
    # -e: extraction threads (12-thread cap)
    # -p: show progress
    # -O: output directory for the uncompressed FASTQ
    # --skip-technical: drop technical/barcode reads (not present
    #   here since the library is single-end cDNA, but safe default)
    # ---------------------------------------------------------------
    echo "[INFO] (2/3) fasterq-dump ${SRR} (threads=${THREADS}) ..."
    fasterq-dump \
        --threads "${THREADS}" \
        --progress \
        --skip-technical \
        --outdir "${DATA_DIR}" \
        "${SRA_CACHE_DIR}/${SRR}/${SRR}.sra"

    # ---------------------------------------------------------------
    # Step 3: compress immediately. Single-end library => one
    # output file named <SRR>.fastq from fasterq-dump.
    # ---------------------------------------------------------------
    echo "[INFO] (3/3) Compressing ${SRR}.fastq -> ${SRR}.fastq.gz ..."
    if [[ -f "${DATA_DIR}/${SRR}.fastq" ]]; then
        ${COMPRESS_CMD} "${DATA_DIR}/${SRR}.fastq"
    else
        echo "[ERROR] Expected ${DATA_DIR}/${SRR}.fastq not found after fasterq-dump." >&2
        echo "[ERROR] (Check whether this run is unexpectedly paired-end.)" >&2
        exit 1
    fi

    # ---------------------------------------------------------------
    # Disk hygiene: remove the .sra container and per-run cache
    # directory immediately. With ~6 samples x ~70-110MB compressed
    # FASTQ each, this keeps the working set tiny on the laptop SSD.
    # ---------------------------------------------------------------
    echo "[INFO] Cleaning SRA cache for ${SRR} ..."
    rm -rf "${SRA_CACHE_DIR}/${SRR}"

    echo "[OK] ${SRR} -> ${FINAL_GZ}"
done

echo "============================================================"
echo "[DONE] Stage 2 download complete. Contents of ${DATA_DIR}:"
ls -lh "${DATA_DIR}"/*.fastq.gz
echo "============================================================"
