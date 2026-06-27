#!/usr/bin/env bash
# =============================================================================
# 03_download_reference.sh
# -----------------------------------------------------------------------------
# Purpose:
#   Download all reference sequences required to build a decoy-aware
#   Salmon index, plus a tx2gene mapping table for tximport.
#
#   Files retrieved (Ensembl release 110, GRCh38.p14):
#     1) Homo_sapiens.GRCh38.cdna.all.fa.gz   - all protein-coding +
#                                               processed transcripts
#     2) Homo_sapiens.GRCh38.ncrna.fa.gz      - lncRNA / ncRNA transcripts
#                                               (needed so lncRNA genes are
#                                               quantifiable - several m6A
#                                               targets are non-coding)
#     3) Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#                                              - full genome, used purely as
#                                               a Salmon "decoy" sequence set
#                                               (Stage 4b)
#
# Why Ensembl release 110:
#   - Released April 2023, built on the GRCh38.p14 patch of the reference
#     assembly. It is the most recent Ensembl release with mature,
#     stable annotation at the time this pipeline was authored, ensuring
#     gene/transcript IDs are consistent with widely-used downstream
#     annotation packages (e.g. org.Hs.eg.db, AnnotationHub snapshots
#     from 2023+). Pinning an explicit release number (rather than
#     "current_fasta") is essential for reproducibility: Ensembl's
#     "current" symlinks move forward over time and would silently
#     change the reference out from under this pipeline.
#
# Execution context:
#   Run from inside scripts/  ->  bash 03_download_reference.sh
#
# Idempotency:
#   Each download/derived file is skipped if it already exists in ../ref/.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration (relative paths assume execution from scripts/)
# -----------------------------------------------------------------------------
REF_DIR="../ref"
ENSEMBL_RELEASE=110
BASE_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens"

CDNA_FA="Homo_sapiens.GRCh38.cdna.all.fa.gz"
NCRNA_FA="Homo_sapiens.GRCh38.ncrna.fa.gz"
GENOME_FA="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

mkdir -p "${REF_DIR}"

# -----------------------------------------------------------------------------
# Helper: download a file only if it doesn't already exist (idempotent),
# with retries for flaky connections.
# -----------------------------------------------------------------------------
download_if_missing() {
    local url="$1"
    local out_path="$2"

    if [[ -f "${out_path}" ]]; then
        echo "[SKIP] $(basename "${out_path}") already present."
        return 0
    fi

    echo "[INFO] Downloading $(basename "${out_path}") ..."
    wget --continue --tries=5 --timeout=60 \
        -O "${out_path}.partial" \
        "${url}"

    mv "${out_path}.partial" "${out_path}"
    echo "[OK] Saved to ${out_path}"
}

# -----------------------------------------------------------------------------
# Step 1: cDNA transcriptome (protein-coding + processed transcripts)
# -----------------------------------------------------------------------------
download_if_missing \
    "${BASE_URL}/cdna/${CDNA_FA}" \
    "${REF_DIR}/${CDNA_FA}"

# -----------------------------------------------------------------------------
# Step 2: ncRNA transcripts (lncRNA, miRNA precursors, etc.)
# -----------------------------------------------------------------------------
download_if_missing \
    "${BASE_URL}/ncrna/${NCRNA_FA}" \
    "${REF_DIR}/${NCRNA_FA}"

# -----------------------------------------------------------------------------
# Step 3: Primary assembly genome (decoy sequences for Stage 4b)
# -----------------------------------------------------------------------------
download_if_missing \
    "${BASE_URL}/dna/${GENOME_FA}" \
    "${REF_DIR}/${GENOME_FA}"

# -----------------------------------------------------------------------------
# Step 4: Build tx2gene.tsv from FASTA headers (no GTF parsing required).
#
# Ensembl cDNA/ncRNA FASTA headers look like:
#   >ENST00000631435.1 cdna chromosome:GRCh38:7:142847306:142847317:1 \
#     gene:ENSG00000282431.1 gene_biotype:TR_C_gene \
#     transcript_biotype:TR_C_gene gene_symbol:TRBD1 description:...
#
# We extract: transcript_id, gene_id, gene_symbol (fallback to gene_id if
# gene_symbol is absent, which happens for some ncRNA entries).
# -----------------------------------------------------------------------------
TX2GENE="${REF_DIR}/tx2gene.tsv"

if [[ -f "${TX2GENE}" ]]; then
    echo "[SKIP] ${TX2GENE} already exists."
else
    echo "[INFO] Building tx2gene.tsv from cDNA + ncRNA headers ..."

    {
        for FA in "${REF_DIR}/${CDNA_FA}" "${REF_DIR}/${NCRNA_FA}"; do
            zcat "${FA}" \
                | grep '^>' \
                | awk '
                    {
                        # Field 1: ">ENST00000631435.1" -> strip leading ">"
                        tx_id = substr($1, 2)

                        gene_id = "NA"
                        gene_symbol = "NA"

                        for (i = 2; i <= NF; i++) {
                            if ($i ~ /^gene:/) {
                                gene_id = substr($i, 6)
                            }
                            if ($i ~ /^gene_symbol:/) {
                                gene_symbol = substr($i, 13)
                            }
                        }

                        # Fall back to gene_id if no symbol annotated
                        if (gene_symbol == "NA") {
                            gene_symbol = gene_id
                        }

                        print tx_id "\t" gene_id "\t" gene_symbol
                    }
                '
        done
    } > "${TX2GENE}"

    N_TX=$(wc -l < "${TX2GENE}")
    echo "[OK] Wrote ${N_TX} transcript->gene mappings to ${TX2GENE}"
fi

echo "============================================================"
echo "[DONE] Stage 4a reference download complete."
ls -lh "${REF_DIR}"
echo "============================================================"
echo "[NEXT] Run 04_build_salmon_index.sh to construct the"
echo "       decoy-aware Salmon index (gentrome strategy)."
echo "============================================================"
