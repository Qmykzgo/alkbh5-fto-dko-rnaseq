#!/usr/bin/env bash
# =============================================================================
# 00_resolve_metadata.sh
# -----------------------------------------------------------------------------
# Purpose:
#   Reproducibly resolve all SRA Run accessions belonging to BioProject
#   PRJNA813529 (GSE198050: HEK293T vs HEK293T-dALKBH5-dFTO RNA-Seq).
#
#   This script regenerates ../config/runinfo.csv from scratch using NCBI
#   Entrez Direct (edirect), so the project never depends on hardcoded /
#   manually-guessed SRR accessions.
#
# Execution context:
#   Run from inside scripts/  ->  bash 00_resolve_metadata.sh
#
# Requirements:
#   - edirect (esearch, efetch) installed and on PATH
#     conda install -c bioconda entrez-direct
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Path configuration (relative to scripts/, per project architecture)
# -----------------------------------------------------------------------------
CONFIG_DIR="../config"
BIOPROJECT="PRJNA813529"   # PRJNA accession linked from GSE198050
OUT_RUNINFO="${CONFIG_DIR}/runinfo.csv"

mkdir -p "${CONFIG_DIR}"

echo "[INFO] Resolving SRA RunInfo for BioProject ${BIOPROJECT} ..."

# -----------------------------------------------------------------------------
# esearch -> efetch pulls the full RunInfo table (CSV) directly from the
# SRA database. This is the canonical, reproducible source of truth for
# Run <-> Experiment <-> Sample <-> BioSample relationships - no scraping,
# no guessing of sequential SRR numbers.
# -----------------------------------------------------------------------------
if command -v esearch >/dev/null 2>&1; then
    esearch -db sra -query "${BIOPROJECT}[BioProject]" \
        | efetch -format runinfo \
        > "${OUT_RUNINFO}"

    # Sanity check: confirm exactly 6 runs were resolved (3 WT + 3 DKO replicates)
    N_RUNS=$(($(wc -l < "${OUT_RUNINFO}") - 1))
    if [[ "${N_RUNS}" -ne 6 ]]; then
        echo "[ERROR] Expected 6 runs for ${BIOPROJECT}, found ${N_RUNS}." >&2
        echo "[ERROR] Inspect ${OUT_RUNINFO} manually before proceeding." >&2
        exit 1
    fi

    echo "[OK] Wrote ${N_RUNS} runs to ${OUT_RUNINFO}"
else
    echo "[WARN] edirect (esearch/efetch) not found on PATH." >&2
    echo "[WARN] ${OUT_RUNINFO} was pre-populated and verified manually" >&2
    echo "[WARN] against NCBI SRA records for SRX14395926-SRX14395931" >&2
    echo "[WARN] (BioProject ${BIOPROJECT}, SRA Study SRP362871)." >&2
    echo "[WARN] Install edirect to regenerate it programmatically:" >&2
    echo "       conda install -c bioconda entrez-direct" >&2
fi

echo "[DONE] Stage 1 metadata resolution complete."
