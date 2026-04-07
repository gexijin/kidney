## 03_differential_expression.R — limma DE: CKD vs Control
## Unpaired design, all male, Agilent 4x44K

library(limma)
library(tidyverse)

base_dir <- "/Users/ge/Documents/research/gtex/analyses/Soulage_GSE98603"
data_dir <- file.path(base_dir, "data")
plots_dir <- file.path(base_dir, "results", "plots")
tables_dir <- file.path(base_dir, "results", "tables")

dat <- readRDS(file.path(data_dir, "02_normalized_data.rds"))
expr_norm <- dat$expr_norm
sample_info <- dat$sample_info

# ── 1. Design matrix ────────────────────────────────────────────────────────

condition <- factor(sample_info$condition, levels = c("Control", "CKD"))
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)

cat("Design matrix:\n")
print(design)

# ── 2. Fit limma model ──────────────────────────────────────────────────────

fit <- lmFit(expr_norm, design)

contrast_matrix <- makeContrasts(
  CKD_vs_Control = CKD - Control,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# ── 3. Extract results ──────────────────────────────────────────────────────

results <- topTable(fit2, coef = "CKD_vs_Control", number = Inf, sort.by = "none")
results$ProbeID <- rownames(results)

cat("\n=== DE Summary (CKD vs Control) ===\n")
cat("Total probes tested:", nrow(results), "\n")
for (fdr in c(0.05, 0.10, 0.25)) {
  sig <- sum(results$adj.P.Val < fdr)
  up <- sum(results$adj.P.Val < fdr & results$logFC > 0)
  dn <- sum(results$adj.P.Val < fdr & results$logFC < 0)
  cat(sprintf("FDR < %.2f: %d DEGs (%d up, %d down)\n", fdr, sig, up, dn))
}

# Also report with logFC cutoff
for (lfc in c(0, 0.5, 1.0)) {
  sig <- sum(results$adj.P.Val < 0.05 & abs(results$logFC) > lfc)
  cat(sprintf("FDR<0.05 & |logFC|>%.1f: %d DEGs\n", lfc, sig))
}

write.csv(results, file.path(tables_dir, "03_DE_CKD_vs_Control.csv"), row.names = FALSE)

# ── 4. Diagnostic plots ─────────────────────────────────────────────────────

# P-value histogram
png(file.path(plots_dir, "03_pvalue_histogram.png"), width = 800, height = 600)
hist(results$P.Value, breaks = 50, main = "P-value Distribution (CKD vs Control)",
     xlab = "P-value", col = "steelblue", border = "white")
abline(h = nrow(results) / 50, col = "red", lty = 2)
legend("topright", "Expected under null", col = "red", lty = 2)
dev.off()

# MA plot
png(file.path(plots_dir, "03_MA_plot.png"), width = 800, height = 600)
plotMA(fit2, coef = 1, main = "MA Plot (CKD vs Control)")
dev.off()

# Mean-variance trend
png(file.path(plots_dir, "03_mean_variance_trend.png"), width = 800, height = 600)
plotSA(fit2, main = "Mean-Variance Trend")
dev.off()

# Volcano plot
results$significance <- case_when(
  results$adj.P.Val < 0.05 & results$logFC > 0.5 ~ "Up",
  results$adj.P.Val < 0.05 & results$logFC < -0.5 ~ "Down",
  TRUE ~ "NS"
)
sig_colors <- c(Up = "#E41A1C", Down = "#377EB8", NS = "grey70")

# Label top genes
top_label <- results %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(P.Value) %>%
  head(20)

p_volcano <- ggplot(results, aes(logFC, -log10(adj.P.Val), color = significance)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = sig_colors) +
  ggrepel::geom_text_repel(data = top_label, aes(label = ProbeID),
                            size = 2, max.overlaps = 15, color = "black") +
  labs(title = "Volcano Plot — CKD vs Control",
       x = "log2 Fold Change", y = "-log10(FDR)") +
  theme_minimal(base_size = 12)
ggsave(file.path(plots_dir, "03_volcano_plot.png"), p_volcano, width = 9, height = 7, dpi = 300)

# ── 5. Save fit object ──────────────────────────────────────────────────────

saveRDS(list(fit2 = fit2, results = results, design = design),
        file.path(data_dir, "03_de_fit.rds"))

cat("\n✓ DE analysis complete\n")
