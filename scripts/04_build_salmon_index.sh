#!/usr/bin/env bash
# =============================================================================
# 04_build_salmon_index.sh
# -----------------------------------------------------------------------------
# Purpose:
#   Build a decoy-aware Salmon index using the "gentrome" strategy:
#     gentrome.fa = transcriptome (cDNA + ncRNA) ++ genome (decoy)
#     decoys.txt  = list of genome sequence names to treat as decoys
#
#   Decoy-aware indexing dramatically reduces spurious mapping of reads
#   that originate from unannotated/intronic/intergenic genomic loci but
#   would otherwise quasi-map to a transcript by sequence similarity -
#   this is the single biggest accuracy improvement over a naive
#   transcriptome-only Salmon index (Srivastava et al. 2020, Genome Biology).
#
# -----------------------------------------------------------------------------
# RAM BUDGET NOTE (read before running):
#   Building a decoy-aware index against the full GRCh38 primary assembly
#   (~3.0 Gbp) typically peaks around 16-20 GB RAM during suffix-array
#   construction. On a 16GB system this is BORDERLINE - it may succeed
#   (especially with WSL2's swap-backed virtual memory), but can also be
#   slow or, in the worst case, get OOM-killed.
#
#   This script builds the FULL decoy-aware index by default. If it fails
#   with an OOM / "Killed" message:
#     1) Increase the WSL2 memory limit in %UserProfile%\.wslconfig:
#          [wsl2]
#          memory=20GB
#          swap=8GB
#        then `wsl --shutdown` and retry.
#     2) OR re-run with FALLBACK_NO_DECOY=1 (see below) to build a
#        transcriptome-only index instead. This is documented explicitly
#        in report/reproduction_report.md as a methods deviation if used.
#
# Execution context:
#   Run from inside scripts/  ->  bash 04_build_salmon_index.sh
#   Fallback (no decoy):       FALLBACK_NO_DECOY=1 bash 04_build_salmon_index.sh
#
# Idempotency:
#   Skips index construction entirely if ../ref/salmon_index/ already
#   contains a complete index (info.json present).
#
# Requirements:
#   - salmon (>= 1.9)  -> conda env "rnaseq-tools"
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration (relative paths assume execution from scripts/)
# -----------------------------------------------------------------------------
REF_DIR="../ref"
THREADS=12
KMER_LEN=31   # Default/recommended k-mer size for reads >= ~75bp

CDNA_FA="${REF_DIR}/Homo_sapiens.GRCh38.cdna.all.fa.gz"
NCRNA_FA="${REF_DIR}/Homo_sapiens.GRCh38.ncrna.fa.gz"
GENOME_FA="${REF_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

GENTROME="${REF_DIR}/gentrome.fa.gz"
DECOYS="${REF_DIR}/decoys.txt"
INDEX_DIR="${REF_DIR}/salmon_index"

FALLBACK_NO_DECOY="${FALLBACK_NO_DECOY:-0}"   # set to 1 to skip decoy step

# -----------------------------------------------------------------------------
# Sanity checks
# -----------------------------------------------------------------------------
if ! command -v salmon >/dev/null 2>&1; then
    echo "[ERROR] 'salmon' not found on PATH. Activate the rnaseq-tools env." >&2
    exit 1
fi

for f in "${CDNA_FA}" "${NCRNA_FA}" "${GENOME_FA}"; do
    if [[ ! -f "${f}" ]]; then
        echo "[ERROR] Missing reference file: ${f}" >&2
        echo "[ERROR] Run 03_download_reference.sh (Stage 4a) first." >&2
        exit 1
    fi
done

# -----------------------------------------------------------------------------
# Idempotency: a complete salmon index directory contains info.json
# -----------------------------------------------------------------------------
if [[ -f "${INDEX_DIR}/info.json" ]]; then
    echo "[SKIP] Salmon index already exists at ${INDEX_DIR}"
    echo "[SKIP] Delete this directory to force a rebuild."
    exit 0
fi

# =============================================================================
# Path A: Decoy-aware index (default, recommended)
# =============================================================================
if [[ "${FALLBACK_NO_DECOY}" -eq 0 ]]; then

    echo "[INFO] Building DECOY-AWARE index (gentrome strategy)."

    # -------------------------------------------------------------------
    # Step 1: decoys.txt - names of all genome sequences (chromosomes +
    # scaffolds) extracted from the primary assembly FASTA headers.
    # Salmon uses this to identify which sequences in gentrome.fa are
    # genomic decoys vs. real transcripts.
    # -------------------------------------------------------------------
    if [[ -f "${DECOYS}" ]]; then
        echo "[SKIP] ${DECOYS} already exists."
    else
        echo "[INFO] Extracting decoy sequence names from genome FASTA ..."
        zcat "${GENOME_FA}" \
            | grep '^>' \
            | cut -d ' ' -f1 \
            | sed 's/^>//' \
            > "${DECOYS}"
        echo "[OK] Wrote $(wc -l < "${DECOYS}") decoy sequence names to ${DECOYS}"
    fi

    # -------------------------------------------------------------------
    # Step 2: gentrome.fa.gz - concatenation of transcriptome (cDNA +
    # ncRNA) followed by the genome. Order matters: transcripts MUST
    # come first so Salmon's transcript IDs map correctly.
    #
    # Concatenating gzip streams is valid (multi-member gzip) and both
    # zcat and salmon's internal zlib reader handle it transparently -
    # no need to decompress to disk first, which saves ~15GB of
    # intermediate space on the laptop SSD.
    # -------------------------------------------------------------------
    if [[ -f "${GENTROME}" ]]; then
        echo "[SKIP] ${GENTROME} already exists."
    else
        echo "[INFO] Building gentrome.fa.gz (cDNA + ncRNA + genome decoy) ..."
        cat "${CDNA_FA}" "${NCRNA_FA}" "${GENOME_FA}" > "${GENTROME}"
        echo "[OK] Wrote ${GENTROME} ($(du -h "${GENTROME}" | cut -f1))"
    fi

    # -------------------------------------------------------------------
    # Step 3: salmon index
    #   -t : gentrome (transcriptome + decoy)
    #   -d : decoy sequence name list
    #   -i : output index directory
    #   -k : k-mer size (31 is standard for >=75bp reads)
    #   -p : threads
    #   --keepDuplicates : retain transcripts with identical sequences
    #         as distinct entries (important for accurate per-isoform
    #         and per-gene counts when paralogous transcripts share
    #         an exon structure)
    # -------------------------------------------------------------------
    echo "[INFO] Running salmon index (this may take 15-40 min, RAM-heavy) ..."
    salmon index \
        --transcripts "${GENTROME}" \
        --decoys "${DECOYS}" \
        --index "${INDEX_DIR}" \
        --kmerLen "${KMER_LEN}" \
        --threads "${THREADS}" \
        --keepDuplicates

# =============================================================================
# Path B: Transcriptome-only fallback (no decoy)
# =============================================================================
else
    echo "[WARN] FALLBACK_NO_DECOY=1 - building TRANSCRIPTOME-ONLY index."
    echo "[WARN] This is a documented methods deviation - record it in"
    echo "[WARN] report/reproduction_report.md (Discrepancy section)."

    TRANSCRIPTOME="${REF_DIR}/transcriptome.fa.gz"
    if [[ -f "${TRANSCRIPTOME}" ]]; then
        echo "[SKIP] ${TRANSCRIPTOME} already exists."
    else
        cat "${CDNA_FA}" "${NCRNA_FA}" > "${TRANSCRIPTOME}"
    fi

    salmon index \
        --transcripts "${TRANSCRIPTOME}" \
        --index "${INDEX_DIR}" \
        --kmerLen "${KMER_LEN}" \
        --threads "${THREADS}" \
        --keepDuplicates
fi

echo "============================================================"
echo "[DONE] Stage 4b Salmon index complete: ${INDEX_DIR}"
ls -lh "${INDEX_DIR}"
echo "============================================================"
echo "[NEXT] Proceed to Stage 5 (Salmon quantification)."
echo "============================================================"
