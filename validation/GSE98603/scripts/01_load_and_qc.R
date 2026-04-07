## 01_load_and_qc.R — Download GSE98603 data, QC, and initial exploration
## Agilent 4x44K v2 (GPL13497): CKD vs Control subcutaneous WAT

library(GEOquery)
library(limma)
library(tidyverse)
library(pheatmap)

base_dir <- "/Users/ge/Documents/research/gtex/analyses/Soulage_GSE98603"
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
plots_dir <- file.path(results_dir, "plots")
tables_dir <- file.path(results_dir, "tables")

# ── 1. Download from GEO ─────────────────────────────────────────────────────

gse <- getGEO("GSE98603", destdir = data_dir, GSEMatrix = TRUE, getGPL = TRUE)
eset <- gse[[1]]

cat("ExpressionSet dimensions:", nrow(eset), "probes x", ncol(eset), "samples\n")

# ── 2. Extract phenotype data ────────────────────────────────────────────────

pdata <- pData(eset)
cat("\nPhenotype columns:\n")
print(names(pdata))
cat("\nCharacteristics:\n")
print(pdata[, grep("characteristics|title|source|description", names(pdata), ignore.case = TRUE)])

# Create clean sample metadata
sample_info <- data.frame(
  sample_id = rownames(pdata),
  title = pdata$title,
  condition = ifelse(grepl("CKD", pdata$title, ignore.case = TRUE), "CKD", "Control"),
  stringsAsFactors = FALSE
)
# Extract any additional characteristics
for (col in grep("characteristics_ch1", names(pdata), value = TRUE)) {
  vals <- as.character(pdata[[col]])
  field_name <- sub(":.*", "", vals[1])
  field_vals <- trimws(sub("^[^:]+:", "", vals))
  sample_info[[field_name]] <- field_vals
}

cat("\nSample metadata:\n")
print(sample_info)
write.csv(sample_info, file.path(tables_dir, "01_sample_metadata.csv"), row.names = FALSE)

# ── 3. Extract expression matrix ────────────────────────────────────────────

expr_mat <- exprs(eset)
cat("\nExpression matrix:", nrow(expr_mat), "probes x", ncol(expr_mat), "samples\n")
cat("Value range:", range(expr_mat, na.rm = TRUE), "\n")
cat("Are values log-scale? (max <", 25, "suggests yes):", max(expr_mat, na.rm = TRUE) < 25, "\n")

# Check if already log2 transformed
if (max(expr_mat, na.rm = TRUE) > 25) {
  cat("Applying log2 transformation...\n")
  expr_mat <- log2(expr_mat + 1)
}

# ── 4. Probe tracking ───────────────────────────────────────────────────────

probe_tracking <- data.frame(
  step = "1. Raw data loaded",
  n_probes = nrow(expr_mat),
  n_removed = 0L,
  description = "Total probes from GEO matrix"
)

# ── 5. NA check ──────────────────────────────────────────────────────────────

n_na <- sum(is.na(expr_mat))
cat(sprintf("\nNA values: %d (%.2f%%)\n", n_na, 100 * n_na / length(expr_mat)))

if (n_na > 0) {
  na_by_sample <- colSums(is.na(expr_mat))
  na_by_probe <- rowSums(is.na(expr_mat))
  cat("Samples with NAs:", sum(na_by_sample > 0), "\n")
  cat("Probes with NAs:", sum(na_by_probe > 0), "\n")

  # Remove probes with any NA
  keep <- complete.cases(expr_mat)
  n_before <- nrow(expr_mat)
  expr_mat <- expr_mat[keep, ]
  probe_tracking <- rbind(probe_tracking, data.frame(
    step = "2. Remove NA probes",
    n_probes = nrow(expr_mat),
    n_removed = n_before - nrow(expr_mat),
    description = "Probes with NA values removed"
  ))
}

# ── 6. Filter low-expression probes ─────────────────────────────────────────

# Keep probes expressed above median in at least one group
median_expr <- median(expr_mat)
min_samples <- min(table(sample_info$condition))  # 9
expressed <- rowSums(expr_mat > median_expr) >= min_samples
n_before <- nrow(expr_mat)
expr_filtered <- expr_mat[expressed, ]

probe_tracking <- rbind(probe_tracking, data.frame(
  step = "3. Expression filter",
  n_probes = nrow(expr_filtered),
  n_removed = n_before - nrow(expr_filtered),
  description = sprintf("Probes above median in >= %d samples", min_samples)
))

# ── 7. Probe tracking summary ───────────────────────────────────────────────

probe_tracking$pct_remaining <- round(100 * probe_tracking$n_probes / probe_tracking$n_probes[1], 1)
cat("\n=== Probe Tracking Summary ===\n")
print(probe_tracking)
write.csv(probe_tracking, file.path(tables_dir, "01_probe_tracking.csv"), row.names = FALSE)

# ── 8. QC Plots ─────────────────────────────────────────────────────────────

condition_colors <- c(Control = "#4DAF4A", CKD = "#E41A1C")
sample_colors <- condition_colors[sample_info$condition]

# 8a. Boxplot of expression distributions
png(file.path(plots_dir, "01_expression_boxplot.png"), width = 1000, height = 600)
par(mar = c(10, 4, 3, 1))
boxplot(expr_filtered, las = 2, col = sample_colors, main = "Expression Distribution by Sample",
        ylab = "log2 expression", cex.axis = 0.7, outline = FALSE)
legend("topright", legend = names(condition_colors), fill = condition_colors, cex = 0.9)
dev.off()

# 8b. Density plot
png(file.path(plots_dir, "01_expression_density.png"), width = 800, height = 600)
plot(density(expr_filtered[, 1]), main = "Expression Density by Sample",
     xlab = "log2 expression", ylim = c(0, 0.5), col = sample_colors[1], lwd = 0.5)
for (i in 2:ncol(expr_filtered)) {
  lines(density(expr_filtered[, i]), col = sample_colors[i], lwd = 0.5)
}
legend("topright", legend = names(condition_colors), col = condition_colors, lwd = 2)
dev.off()

# 8c. Sample correlation heatmap
cor_mat <- cor(expr_filtered, method = "pearson")
ann_col <- data.frame(Condition = sample_info$condition, row.names = sample_info$sample_id)
ann_colors <- list(Condition = condition_colors)

png(file.path(plots_dir, "01_sample_correlation_heatmap.png"), width = 800, height = 700)
pheatmap(cor_mat,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         main = "Sample Correlation (Pearson)",
         display_numbers = TRUE, number_format = "%.3f",
         fontsize_number = 7, fontsize = 8)
dev.off()

# 8d. PCA
pca <- prcomp(t(expr_filtered), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  condition = sample_info$condition,
  label = sample_info$title
)
var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = label), vjust = -0.8, size = 2.5) +
  scale_color_manual(values = condition_colors) +
  labs(title = "PCA — GSE98603 (Filtered)",
       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2])) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, "01_PCA_raw.png"), p_pca, width = 8, height = 6, dpi = 300)

# ── 9. Save processed data ──────────────────────────────────────────────────

saveRDS(list(
  expr_filtered = expr_filtered,
  sample_info = sample_info,
  eset = eset,
  probe_tracking = probe_tracking
), file.path(data_dir, "01_qc_data.rds"))

cat("\n✓ QC complete. Saved to", data_dir, "\n")
cat("  Filtered expression:", nrow(expr_filtered), "probes x", ncol(expr_filtered), "samples\n")
