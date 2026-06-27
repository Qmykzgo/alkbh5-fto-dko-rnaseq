#!/usr/bin/env Rscript
# =============================================================================
# 06_deseq2_analysis.R
# -----------------------------------------------------------------------------
# Purpose:
#   1) Import Salmon quant.sf files (gene-level, via tx2gene) with tximport
#   2) Build a DESeqDataSet with design ~ condition (DKO vs WT)
#   3) Run DESeq2 normalization + Wald testing
#   4) Apply apeglm log-fold-change shrinkage (falls back to "normal" if
#      apeglm is unavailable)
#   5) Export normalized counts + full/significant DEG tables
#   6) Produce exploratory PCA, MA, and volcano plots
#   7) Persist dds / res / vsd as .rds for Stage 7 (publication figures)
#
# Execution context:
#   Run from inside scripts/  ->  Rscript 06_deseq2_analysis.R
#
# Idempotency:
#   This script is deterministic given identical inputs and overwrites its
#   own outputs on each run (DESeq2 analyses are cheap relative to upstream
#   alignment/quant steps, so re-computation rather than skip-logic is the
#   appropriate choice here).
#
# Requirements (R / Bioconductor):
#   - tximport, DESeq2, ggplot2, readr, dplyr, tibble
#   - apeglm (optional, recommended for LFC shrinkage)
# =============================================================================

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tibble)
})

set.seed(42)

# -----------------------------------------------------------------------------
# Configuration (relative paths assume execution from scripts/)
# -----------------------------------------------------------------------------
config_dir   <- "../config"
salmon_dir   <- "../results/salmon"
ref_dir      <- "../ref"
out_dir      <- "../results/deseq2"
fig_dir      <- "../figures"

design_path  <- file.path(config_dir, "design_matrix.csv")
tx2gene_path <- file.path(ref_dir, "tx2gene.tsv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Step 1: Load design matrix and tx2gene mapping
# -----------------------------------------------------------------------------
if (!file.exists(design_path)) {
  stop(sprintf("[ERROR] Design matrix not found: %s (run Stage 1)", design_path))
}
if (!file.exists(tx2gene_path)) {
  stop(sprintf("[ERROR] tx2gene mapping not found: %s (run Stage 4a)", tx2gene_path))
}

coldata <- read_csv(design_path, show_col_types = FALSE)

# Reference level MUST be "WT" so that all log2FoldChange values are
# interpreted as DKO vs WT (matching the biological framing of the study:
# "what changes upon loss of m6A demethylation?")
coldata$condition <- factor(coldata$condition, levels = c("WT", "DKO"))
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample

cat("[INFO] Design matrix:\n")
print(coldata)

# tx2gene.tsv has 3 columns: transcript_id, gene_id, gene_symbol.
# tximport only needs the first two for aggregation; we keep the full table
# separately to annotate DEG results with gene symbols later.
tx2gene_full <- read_tsv(
  tx2gene_path,
  col_names = c("transcript_id", "gene_id", "gene_symbol"),
  show_col_types = FALSE
)

# Strip version suffixes to match tximport(ignoreTxVersion = TRUE) and clean downstream annotations
tx2gene_full$transcript_id <- sub("\\..*$", "", tx2gene_full$transcript_id)
tx2gene_full$gene_id       <- sub("\\..*$", "", tx2gene_full$gene_id)

tx2gene <- tx2gene_full[, c("transcript_id", "gene_id")]

# gene_id -> gene_symbol lookup (one row per unique gene)
gene_annotation <- tx2gene_full %>%
  distinct(gene_id, gene_symbol)


# -----------------------------------------------------------------------------
# Step 2: Locate quant.sf files and run tximport
# -----------------------------------------------------------------------------
quant_files <- file.path(salmon_dir, coldata$sample, "quant.sf")
names(quant_files) <- coldata$sample

missing <- quant_files[!file.exists(quant_files)]
if (length(missing) > 0) {
  stop(sprintf(
    "[ERROR] Missing quant.sf for: %s\n[ERROR] Run Stage 5 (05_run_salmon_quant.sh) first.",
    paste(names(missing), collapse = ", ")
  ))
}

cat("[INFO] Importing", length(quant_files), "Salmon quant files via tximport ...\n")

# ignoreTxVersion = TRUE: Ensembl transcript IDs in quant.sf retain version
# suffixes (e.g. ENST00000631435.1). tx2gene was built from the same FASTA
# headers so versions should already match, but ignoreTxVersion guards
# against any edge-case mismatch (e.g. patch-level annotation drift).
txi <- tximport(
  files = quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

cat("[OK] tximport complete. Gene-level matrix dimensions:",
    paste(dim(txi$counts), collapse = " x "), "\n")

# -----------------------------------------------------------------------------
# Step 3: Build DESeqDataSet and pre-filter low-count genes
# -----------------------------------------------------------------------------
dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = coldata,
  design = ~condition
)

# Pre-filtering: remove genes with very low total counts across all samples.
# This is NOT the same as the internal DESeq2 independent-filtering step
# applied to p-values later - it simply reduces noise/multiple-testing
# burden from genes with essentially no signal (Love et al. 2014, the
# DESeq2 vignette, recommend a minimal rowSum filter before DESeq()).
keep <- rowSums(counts(dds)) >= 10
cat(sprintf(
  "[INFO] Pre-filtering: keeping %d / %d genes (rowSum >= 10 counts).\n",
  sum(keep), length(keep)
))
dds <- dds[keep, ]

# -----------------------------------------------------------------------------
# Step 4: Run DESeq2 (estimates size factors, dispersions, fits NB GLM,
# performs Wald test for the condition coefficient)
# -----------------------------------------------------------------------------
cat("[INFO] Running DESeq() ...\n")
dds <- DESeq(dds)

# -----------------------------------------------------------------------------
# Step 5: Extract results (DKO vs WT) and apply LFC shrinkage
# -----------------------------------------------------------------------------
contrast <- c("condition", "DKO", "WT")
res <- results(dds, contrast = contrast, alpha = 0.05)

cat("[INFO] Raw results summary (condition: DKO vs WT):\n")
summary(res)

# apeglm shrinkage produces more reliable LFC estimates for genes with low
# counts/high dispersion, without distorting the p-values used for
# significance calling (Zhu, Ibrahim & Love 2019). Falls back gracefully
# to "normal" shrinkage if the apeglm package is not installed.
coef_name <- "condition_DKO_vs_WT"
apeglm_available <- requireNamespace("apeglm", quietly = TRUE)

dds_coef_names <- tryCatch(colnames(coef(dds)), error = function(e) character(0))
if (apeglm_available && coef_name %in% dds_coef_names) {
  cat("[INFO] Applying apeglm LFC shrinkage ...\n")
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res)
} else if (apeglm_available) {
  # coef name fallback - look it up dynamically from resultsNames()
  rn <- resultsNames(dds)
  coef_match <- rn[grepl("condition", rn)]
  if (length(coef_match) == 1) {
    cat("[INFO] Applying apeglm LFC shrinkage (coef =", coef_match, ") ...\n")
    res_shrunk <- lfcShrink(dds, coef = coef_match, type = "apeglm", res = res)
  } else {
    cat("[WARN] Could not resolve apeglm coefficient name - using 'normal' shrinkage.\n")
    res_shrunk <- lfcShrink(dds, contrast = contrast, type = "normal", res = res)
  }
} else {
  cat("[WARN] apeglm not installed - using 'normal' shrinkage.\n")
  cat("[WARN] Install with: BiocManager::install('apeglm')\n")
  res_shrunk <- lfcShrink(dds, contrast = contrast, type = "normal", res = res)
}

# -----------------------------------------------------------------------------
# Step 6: Export normalized counts + DEG tables (annotated with gene symbols)
# -----------------------------------------------------------------------------
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_annotation, by = "gene_id") %>%
  relocate(gene_symbol, .after = gene_id)

write_csv(norm_counts_df, file.path(out_dir, "normalized_counts.csv"))
cat("[OK] Wrote normalized_counts.csv (", nrow(norm_counts_df), "genes )\n")

res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_annotation, by = "gene_id") %>%
  relocate(gene_symbol, .after = gene_id) %>%
  arrange(padj)

write_csv(res_df, file.path(out_dir, "deg_results_full.csv"))
cat("[OK] Wrote deg_results_full.csv (", nrow(res_df), "genes )\n")

# Significance thresholds: padj < 0.05 AND |log2FC| > 1 (2-fold change) -
# a conventional, defensible cutoff for a portfolio reproduction study.
sig_df <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

write_csv(sig_df, file.path(out_dir, "deg_significant.csv"))
cat(sprintf(
  "[OK] Wrote deg_significant.csv (%d genes; padj<0.05 & |log2FC|>1)\n",
  nrow(sig_df)
))

# -----------------------------------------------------------------------------
# Step 7: Exploratory figures (PCA, MA, volcano)
#
# These are quick diagnostic plots written here for immediate feedback.
# Stage 7 (07_make_figures.R) regenerates publication-quality versions
# (PNG + PDF, consistent theme, top-variable-gene heatmap) from the .rds
# objects saved at the end of this script.
# -----------------------------------------------------------------------------

# --- PCA plot (variance-stabilizing transform) ------------------------------
vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA - VST-transformed counts (WT vs ALKBH5/FTO DKO)") +
  theme_bw()

ggsave(file.path(fig_dir, "pca_plot.png"), p_pca, width = 6, height = 5, dpi = 150)

# --- MA plot ------------------------------------------------------------------
png(file.path(fig_dir, "ma_plot.png"), width = 1200, height = 900, res = 150)
plotMA(res_shrunk, ylim = c(-5, 5), main = "MA plot - DKO vs WT (apeglm-shrunk LFC)")
dev.off()

# --- Volcano plot ---------------------------------------------------------
volcano_df <- res_df %>%
  mutate(
    neg_log10_padj = -log10(padj),
    significant = case_when(
      is.na(padj)            ~ "NA",
      padj < 0.05 & log2FoldChange >  1 ~ "Up in DKO",
      padj < 0.05 & log2FoldChange < -1 ~ "Down in DKO",
      TRUE                    ~ "NS"
    )
  )

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj, color = significant)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c(
    "Up in DKO"   = "#D62728",
    "Down in DKO" = "#1F77B4",
    "NS"          = "grey70",
    "NA"          = "grey90"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(
    title = "Volcano plot - HEK293T dALKBH5dFTO vs WT",
    x = "log2 fold change (DKO vs WT)",
    y = expression(-log[10]~"adjusted p-value"),
    color = NULL
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "volcano_plot.png"), p_volcano, width = 7, height = 6, dpi = 150)

# -----------------------------------------------------------------------------
# Step 8: Persist key R objects for Stage 7
# -----------------------------------------------------------------------------
saveRDS(dds,        file.path(out_dir, "dds.rds"))
saveRDS(res_shrunk, file.path(out_dir, "res_shrunk.rds"))
saveRDS(vsd,        file.path(out_dir, "vsd.rds"))
saveRDS(gene_annotation, file.path(out_dir, "gene_annotation.rds"))

cat("============================================================\n")
cat("[DONE] Stage 6 DESeq2 analysis complete.\n")
cat("[INFO] Tables  -> ", out_dir, "\n")
cat("[INFO] Figures -> ", fig_dir, "\n")
cat(sprintf("[INFO] %d significant DEGs (padj<0.05, |log2FC|>1)\n", nrow(sig_df)))
cat("============================================================\n")
cat("[NEXT] Proceed to Stage 7 (publication-quality figures + heatmap).\n")
cat("============================================================\n")
