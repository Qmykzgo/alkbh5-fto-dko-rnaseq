#!/usr/bin/env Rscript
# =============================================================================
# 07_make_figures.R
# -----------------------------------------------------------------------------
# Purpose:
#   Generate three publication-quality figures from the DESeq2 objects
#   persisted in Stage 6:
#
#     1) PCA plot        (VST-transformed counts, top 500 most variable genes)
#     2) Volcano plot    (apeglm-shrunk LFC vs -log10 adjusted p-value)
#     3) Heatmap         (top 50 most variable genes, sample-level z-scores)
#
#   Each figure is exported as both PNG (150 dpi, for README/web preview)
#   and PDF (vector, for print/manuscript submission).
#
# Execution context:
#   Run from inside scripts/  ->  Rscript 07_make_figures.R
#
# Requirements (R / Bioconductor):
#   - ggplot2, ggrepel, dplyr, tibble, readr (CRAN)
#   - pheatmap, RColorBrewer                 (CRAN)
#   - DESeq2                                 (Bioconductor)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library(readr)
  library(pheatmap)
  library(RColorBrewer)
  library(DESeq2)
  library(MatrixGenerics)   # rowVars() - re-exported via DESeq2/SummarizedExperiment
})

set.seed(42)

# -----------------------------------------------------------------------------
# Configuration (relative paths assume execution from scripts/)
# -----------------------------------------------------------------------------
deseq_dir <- "../results/deseq2"
fig_dir   <- "../figures"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Helper: save a ggplot object to both PNG and PDF with consistent sizing
# -----------------------------------------------------------------------------
save_plot <- function(plot_obj, stem, width, height) {
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  ggsave(png_path, plot_obj, width = width, height = height, dpi = 150)
  ggsave(pdf_path, plot_obj, width = width, height = height)
  cat(sprintf("[OK] Saved: %s (.png + .pdf)\n", stem))
}

# -----------------------------------------------------------------------------
# Load persisted R objects from Stage 6
# -----------------------------------------------------------------------------
rds_files <- c("dds.rds", "res_shrunk.rds", "vsd.rds", "gene_annotation.rds")
for (f in rds_files) {
  full_path <- file.path(deseq_dir, f)
  if (!file.exists(full_path)) {
    stop(sprintf(
      "[ERROR] %s not found.\n[ERROR] Run 06_deseq2_analysis.R (Stage 6) first.",
      full_path
    ))
  }
}

dds            <- readRDS(file.path(deseq_dir, "dds.rds"))
res_shrunk     <- readRDS(file.path(deseq_dir, "res_shrunk.rds"))
vsd            <- readRDS(file.path(deseq_dir, "vsd.rds"))
gene_annotation <- readRDS(file.path(deseq_dir, "gene_annotation.rds"))

cat("[INFO] R objects loaded from", deseq_dir, "\n")

# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 1: PCA PLOT
# ──────────────────────────────────────────────────────────────────────────────
# plotPCA uses the top 500 most variable genes (ntop=500) by default.
# returnData=TRUE gives us the projected coordinates to pass to ggplot
# for full aesthetic control (labels, point size, custom palette).
# ──────────────────────────────────────────────────────────────────────────────
cat("[INFO] Building PCA plot ...\n")

pca_data    <- plotPCA(vsd, intgroup = "condition", ntop = 500, returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)

# Clean sample labels: strip the SRR prefix for readability on the plot.
pca_data$label <- rownames(pca_data)

condition_colors <- c("WT" = "#2166AC", "DKO" = "#D6604D")

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2,
                               color = condition, label = label)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(
    size         = 3.5,
    fontface     = "bold",
    box.padding  = 0.4,
    show.legend  = FALSE,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = condition_colors,
    labels = c(
      "WT"  = "Wild-type (WT)",
      "DKO" = expression(paste(Delta, "ALKBH5", Delta, "FTO (DKO)"))
    )
  ) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  labs(
    title    = "Principal Component Analysis",
    subtitle = "VST-normalized counts, top 500 variable genes",
    color    = "Condition"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 14),
    plot.subtitle   = element_text(size = 11, color = "grey40"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

save_plot(p_pca, "fig1_pca", width = 7, height = 6)

# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 2: VOLCANO PLOT
# ──────────────────────────────────────────────────────────────────────────────
# Uses apeglm-shrunk LFC values from res_shrunk.
# Labels the top 15 most significant genes by adjusted p-value on each
# side of the fold-change axis using ggrepel (non-overlapping labels).
# ──────────────────────────────────────────────────────────────────────────────
cat("[INFO] Building volcano plot ...\n")

volcano_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_annotation, by = "gene_id") %>%
  filter(!is.na(pvalue)) %>%
  mutate(
    neg_log10_padj = -log10(replace(padj, padj == 0, .Machine$double.xmin)),
    direction = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up in DKO",
      padj < 0.05 & log2FoldChange < -1 ~ "Down in DKO",
      TRUE                               ~ "NS"
    )
  )

# Select top genes to label (top 12 up + top 12 down by significance)
top_up <- volcano_df %>%
  filter(direction == "Up in DKO") %>%
  slice_min(padj, n = 12)

top_down <- volcano_df %>%
  filter(direction == "Down in DKO") %>%
  slice_min(padj, n = 12)

label_df <- bind_rows(top_up, top_down) %>%
  mutate(display_label = if_else(
    is.na(gene_symbol) | gene_symbol == gene_id,
    gene_id,
    gene_symbol
  ))

n_up   <- sum(volcano_df$direction == "Up in DKO",   na.rm = TRUE)
n_down <- sum(volcano_df$direction == "Down in DKO",  na.rm = TRUE)

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj,
                                     color = direction)) +
  # Draw non-significant points first (bottom layer), then significant ones
  geom_point(data = filter(volcano_df, direction == "NS"),
             alpha = 0.3, size = 0.8) +
  geom_point(data = filter(volcano_df, direction != "NS"),
             alpha = 0.7, size = 1.2) +
  # Significance threshold lines
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  # Gene labels (non-overlapping)
  geom_text_repel(
    data         = label_df,
    aes(label    = display_label),
    size         = 2.8,
    fontface     = "italic",
    box.padding  = 0.35,
    point.padding = 0.2,
    max.overlaps = 20,
    show.legend  = FALSE,
    segment.color = "grey60",
    segment.size  = 0.3
  ) +
  scale_color_manual(values = c(
    "Up in DKO"   = "#D6604D",
    "Down in DKO" = "#2166AC",
    "NS"          = "grey75"
  )) +
  # Annotation: counts of DEGs per direction
  annotate("text", x =  Inf, y = Inf,
           label = paste0("Up: ", n_up),
           hjust = 1.1, vjust = 1.5, color = "#D6604D", size = 4) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("Down: ", n_down),
           hjust = -0.1, vjust = 1.5, color = "#2166AC", size = 4) +
  labs(
    title    = expression(paste("Differential Expression: ", Delta, "ALKBH5",
                                Delta, "FTO vs Wild-type")),
    subtitle = "padj < 0.05, |log2FC| > 1  |  apeglm-shrunk LFC  |  Ensembl 110",
    x        = expression(log[2]~"fold change (DKO / WT)"),
    y        = expression(-log[10]~"(adjusted p-value)"),
    color    = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "grey40"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

save_plot(p_volcano, "fig2_volcano", width = 8, height = 7)

# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 3: HEATMAP — TOP 50 MOST VARIABLE GENES
# ──────────────────────────────────────────────────────────────────────────────
# Select the top 50 genes by across-sample variance in VST space.
# Each cell is z-scored per row (gene) so that colour encodes relative
# expression change rather than absolute abundance (which would make
# high-expression housekeeping genes dominate the colour scale).
# ──────────────────────────────────────────────────────────────────────────────
cat("[INFO] Building heatmap (top 50 variable genes) ...\n")

vst_mat <- assay(vsd)

# Variance across all samples
row_vars    <- rowVars(vst_mat)
top50_idx   <- order(row_vars, decreasing = TRUE)[1:50]
top50_mat   <- vst_mat[top50_idx, ]

# Annotate row names with gene symbols where available
gene_sym_map <- setNames(gene_annotation$gene_symbol, gene_annotation$gene_id)
rownames(top50_mat) <- ifelse(
  is.na(gene_sym_map[rownames(top50_mat)]) |
    gene_sym_map[rownames(top50_mat)] == rownames(top50_mat),
  rownames(top50_mat),
  gene_sym_map[rownames(top50_mat)]
)

# Z-score per gene (row-wise)
top50_scaled <- t(scale(t(top50_mat)))

# Column (sample) annotation bar
col_annotation <- data.frame(
  Condition = colData(vsd)$condition,
  row.names = colnames(vsd)
)

annotation_colors <- list(
  Condition = c("WT" = "#2166AC", "DKO" = "#D6604D")
)

# Diverging blue-white-red palette (safe for most colorblindness types)
heatmap_colors <- colorRampPalette(
  rev(brewer.pal(11, "RdBu"))
)(100)

# pheatmap writes directly to file (no ggplot wrapper)
for (ext in c("png", "pdf")) {
  out_path <- file.path(fig_dir, paste0("fig3_heatmap.", ext))

  pheatmap(
    mat                = top50_scaled,
    color              = heatmap_colors,
    breaks             = seq(-2.5, 2.5, length.out = 101),
    annotation_col     = col_annotation,
    annotation_colors  = annotation_colors,
    cluster_rows       = TRUE,
    cluster_cols       = TRUE,
    show_rownames      = TRUE,
    show_colnames      = TRUE,
    fontsize_row       = 8,
    fontsize_col       = 10,
    border_color       = NA,
    main               = "Top 50 Most Variable Genes (VST, z-scored)",
    filename           = out_path,
    width              = 9,
    height             = 12
  )
  cat(sprintf("[OK] Saved: fig3_heatmap.%s\n", ext))
}

# ──────────────────────────────────────────────────────────────────────────────
# SUMMARY
# ──────────────────────────────────────────────────────────────────────────────
cat("============================================================\n")
cat("[DONE] Stage 7 visualization complete.\n")
cat("[INFO] Figures written to:", fig_dir, "\n")
cat("  fig1_pca.{png,pdf}\n")
cat("  fig2_volcano.{png,pdf}\n")
cat("  fig3_heatmap.{png,pdf}\n")
cat("============================================================\n")
cat("[NEXT] Proceed to Stage 8 (README.md + reproduction report).\n")
cat("============================================================\n")
