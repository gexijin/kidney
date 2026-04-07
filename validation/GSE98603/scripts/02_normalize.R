## 02_normalize.R — Verify normalization and generate QC visualizations
## GEO data was already quantile-normalized by the submitters
## (Feature Extraction 11.5.1.1 → Normexp bg correction → quantile normalization)
## No additional normalization is applied.

library(limma)
library(tidyverse)
library(pheatmap)

base_dir <- "/Users/ge/Documents/research/gtex/analyses/Soulage_GSE98603"
data_dir <- file.path(base_dir, "data")
plots_dir <- file.path(base_dir, "results", "plots")
tables_dir <- file.path(base_dir, "results", "tables")

qc <- readRDS(file.path(data_dir, "01_qc_data.rds"))
expr_filtered <- qc$expr_filtered
sample_info <- qc$sample_info

# ── 1. Use existing normalization (already quantile-normalized by submitters) ─

expr_norm <- expr_filtered
cat("Expression matrix (pre-normalized by submitters):", nrow(expr_norm), "probes x", ncol(expr_norm), "samples\n")

# ── 2. Post-normalization boxplot ────────────────────────────────────────────

condition_colors <- c(Control = "#4DAF4A", CKD = "#E41A1C")
sample_colors <- condition_colors[sample_info$condition]

png(file.path(plots_dir, "02_normalized_boxplot.png"), width = 1000, height = 600)
par(mar = c(10, 4, 3, 1))
boxplot(expr_norm, las = 2, col = sample_colors, main = "Normalized Expression Distribution",
        ylab = "log2 expression", cex.axis = 0.7, outline = FALSE)
legend("topright", legend = names(condition_colors), fill = condition_colors)
dev.off()

# ── 3. PCA post-normalization ────────────────────────────────────────────────

pca <- prcomp(t(expr_norm), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  condition = sample_info$condition,
  label = sample_info$title
)
var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3.5) +
  geom_text(aes(label = label), vjust = -0.8, size = 2.3) +
  scale_color_manual(values = condition_colors) +
  labs(title = "PCA — GSE98603 (Normalized)",
       x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(plots_dir, "02_PCA_normalized.png"), p_pca, width = 8, height = 6, dpi = 300)

# ── 4. Heatmap of top 2000 most variable probes ─────────────────────────────

rv <- apply(expr_norm, 1, var)
top2000 <- names(sort(rv, decreasing = TRUE))[1:min(2000, length(rv))]
mat <- expr_norm[top2000, ]
mat_centered <- mat - apply(mat, 1, median)
mat_clipped <- pmin(pmax(mat_centered, -3), 3)

ann_col <- data.frame(Condition = sample_info$condition, row.names = sample_info$sample_id)
ann_colors <- list(Condition = condition_colors)

# Order samples by condition
sample_order <- order(sample_info$condition)

png(file.path(plots_dir, "02_heatmap_top2000.png"), width = 800, height = 1000)
pheatmap(mat_clipped[, sample_order],
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         clustering_distance_rows = "correlation",
         clustering_method = "average",
         show_rownames = FALSE,
         color = colorRampPalette(c("green3", "black", "red3"))(100),
         main = "Top 2000 Variable Probes (median-centered)")
dev.off()

# ── 5. Sample dendrogram ────────────────────────────────────────────────────

dist_mat <- as.dist(1 - cor(expr_norm))
hc <- hclust(dist_mat, method = "average")

png(file.path(plots_dir, "02_sample_dendrogram.png"), width = 800, height = 500)
par(mar = c(5, 4, 3, 1))
plot(hc, main = "Sample Clustering (1 - Pearson, average linkage)",
     xlab = "", sub = "", labels = sample_info$title)
dev.off()

# ── 6. Save ─────────────────────────────────────────────────────────────────

saveRDS(list(
  expr_norm = expr_norm,
  sample_info = sample_info,
  eset = qc$eset
), file.path(data_dir, "02_normalized_data.rds"))

cat("\n✓ Normalization complete:", nrow(expr_norm), "probes\n")
