## 06_gtex_comparison.R — Compare GSE98603 (male CKD) with GTEx CKD results
## Key question: Does the male-only Soulage dataset show the same signal as
## GTEx females, GTEx males, or neither?

library(tidyverse)
library(ggrepel)

base_dir <- "/Users/ge/Documents/research/gtex/analyses/Soulage_GSE98603"
kidney_dir <- "/Users/ge/Documents/research/gtex/kidney"
data_dir <- file.path(base_dir, "data")
plots_dir <- file.path(base_dir, "results", "plots")
tables_dir <- file.path(base_dir, "results", "tables")

# ── 1. Load all DE results ──────────────────────────────────────────────────

soulage <- readRDS(file.path(data_dir, "04_de_clean.rds"))
gtex_female_list <- readRDS(file.path(kidney_dir, "results", "de_female.rds"))
gtex_male_list <- readRDS(file.path(kidney_dir, "results", "de_male.rds"))

# Extract topTable from list format
gtex_f <- gtex_female_list$tt
gtex_m <- gtex_male_list$tt

# GTEx uses Ensembl IDs — map to gene symbols via annotation
gene_annot <- readRDS("/Users/ge/Documents/research/gtex/data/annotation/gencode_v47_genes.rds")

gtex_f$gene_symbol <- gene_annot[rownames(gtex_f), "symbol"]
gtex_m$gene_symbol <- gene_annot[rownames(gtex_m), "symbol"]

# Remove unmapped
gtex_f <- gtex_f %>% filter(!is.na(gene_symbol) & gene_symbol != "")
gtex_m <- gtex_m %>% filter(!is.na(gene_symbol) & gene_symbol != "")

# Collapse GTEx to unique symbols (keep most significant)
gtex_f <- gtex_f %>% arrange(P.Value) %>% filter(!duplicated(gene_symbol))
gtex_m <- gtex_m %>% arrange(P.Value) %>% filter(!duplicated(gene_symbol))

cat("GSE98603 (Soulage, males):", nrow(soulage), "genes,", sum(soulage$adj.P.Val < 0.05), "DEGs\n")
cat("GTEx Female:", nrow(gtex_f), "genes,", sum(gtex_f$adj.P.Val < 0.05), "DEGs\n")
cat("GTEx Male:", nrow(gtex_m), "genes,", sum(gtex_m$adj.P.Val < 0.05), "DEGs\n")

# ── 2. Merge datasets ───────────────────────────────────────────────────────

merged_f <- inner_join(
  soulage %>% dplyr::select(gene_symbol, logFC, t, P.Value, adj.P.Val) %>%
    rename(logFC_soulage = logFC, t_soulage = t, pval_soulage = P.Value, fdr_soulage = adj.P.Val),
  gtex_f %>% dplyr::select(gene_symbol, logFC, t, P.Value, adj.P.Val) %>%
    rename(logFC_gtex_f = logFC, t_gtex_f = t, pval_gtex_f = P.Value, fdr_gtex_f = adj.P.Val),
  by = "gene_symbol"
)

merged_m <- inner_join(
  soulage %>% dplyr::select(gene_symbol, logFC, t, P.Value, adj.P.Val) %>%
    rename(logFC_soulage = logFC, t_soulage = t, pval_soulage = P.Value, fdr_soulage = adj.P.Val),
  gtex_m %>% dplyr::select(gene_symbol, logFC, t, P.Value, adj.P.Val) %>%
    rename(logFC_gtex_m = logFC, t_gtex_m = t, pval_gtex_m = P.Value, fdr_gtex_m = adj.P.Val),
  by = "gene_symbol"
)

cat("\nOverlapping genes with GTEx female:", nrow(merged_f), "\n")
cat("Overlapping genes with GTEx male:", nrow(merged_m), "\n")

# ── 3. Genome-wide correlation analysis ─────────────────────────────────────

cor_f <- cor.test(merged_f$t_soulage, merged_f$t_gtex_f, method = "spearman")
cor_m <- cor.test(merged_m$t_soulage, merged_m$t_gtex_m, method = "spearman")

cat("\n=== Genome-wide t-statistic correlations ===\n")
cat(sprintf("Soulage vs GTEx Female: rho = %.3f, p = %.2e\n", cor_f$estimate, cor_f$p.value))
cat(sprintf("Soulage vs GTEx Male:   rho = %.3f, p = %.2e\n", cor_m$estimate, cor_m$p.value))

cor_lfc_f <- cor.test(merged_f$logFC_soulage, merged_f$logFC_gtex_f, method = "spearman")
cor_lfc_m <- cor.test(merged_m$logFC_soulage, merged_m$logFC_gtex_m, method = "spearman")

cat(sprintf("\nlogFC Soulage vs GTEx Female: rho = %.3f, p = %.2e\n", cor_lfc_f$estimate, cor_lfc_f$p.value))
cat(sprintf("logFC Soulage vs GTEx Male:   rho = %.3f, p = %.2e\n", cor_lfc_m$estimate, cor_lfc_m$p.value))

# ── 4. Direction concordance among GTEx female DEGs ─────────────────────────

gtex_f_degs <- gtex_f %>% filter(adj.P.Val < 0.05) %>% pull(gene_symbol)
soulage_at_f_degs <- merged_f %>% filter(gene_symbol %in% gtex_f_degs)

concordance <- NA
if (nrow(soulage_at_f_degs) > 0) {
  concordance <- mean(sign(soulage_at_f_degs$logFC_soulage) == sign(soulage_at_f_degs$logFC_gtex_f))
  cat(sprintf("\nDirection concordance at GTEx female DEGs: %.1f%% (%d genes)\n",
              100 * concordance, nrow(soulage_at_f_degs)))
  binom_test <- binom.test(sum(sign(soulage_at_f_degs$logFC_soulage) == sign(soulage_at_f_degs$logFC_gtex_f)),
                           nrow(soulage_at_f_degs), 0.5)
  cat(sprintf("Binomial test: p = %.2e\n", binom_test$p.value))

  # Also check: among concordant genes, are Soulage p-values enriched for small values?
  soulage_at_f_degs$concordant <- sign(soulage_at_f_degs$logFC_soulage) == sign(soulage_at_f_degs$logFC_gtex_f)
  wt <- wilcox.test(pval_soulage ~ concordant, data = soulage_at_f_degs)
  cat(sprintf("Concordant genes have lower Soulage p-values? Wilcox p = %.2e\n", wt$p.value))
}

# ── 5. Scatter plots ────────────────────────────────────────────────────────

# Soulage vs GTEx Female
p1 <- ggplot(merged_f, aes(t_gtex_f, t_soulage)) +
  geom_point(alpha = 0.15, size = 0.5, color = "grey50") +
  geom_point(data = merged_f %>% filter(fdr_gtex_f < 0.05),
             alpha = 0.3, size = 0.8, color = "#E41A1C") +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("rho = %.3f\np = %.1e\nn = %d", cor_f$estimate, cor_f$p.value, nrow(merged_f)),
           hjust = 1.1, vjust = 1.3, size = 4) +
  labs(title = "GSE98603 (Males, CKD V) vs GTEx Female CKD",
       subtitle = "Red: GTEx female FDR < 0.05",
       x = "GTEx Female t-statistic", y = "Soulage t-statistic") +
  theme_minimal(base_size = 12)
ggsave(file.path(plots_dir, "06_scatter_vs_gtex_female.png"), p1, width = 8, height = 7, dpi = 300)

# Soulage vs GTEx Male
p2 <- ggplot(merged_m, aes(t_gtex_m, t_soulage)) +
  geom_point(alpha = 0.15, size = 0.5, color = "grey50") +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("rho = %.3f\np = %.1e\nn = %d", cor_m$estimate, cor_m$p.value, nrow(merged_m)),
           hjust = 1.1, vjust = 1.3, size = 4) +
  labs(title = "GSE98603 (Males, CKD V) vs GTEx Male CKD",
       x = "GTEx Male t-statistic", y = "Soulage t-statistic") +
  theme_minimal(base_size = 12)
ggsave(file.path(plots_dir, "06_scatter_vs_gtex_male.png"), p2, width = 8, height = 7, dpi = 300)

# Combined logFC panel
merged_all <- inner_join(
  soulage %>% dplyr::select(gene_symbol, logFC) %>% rename(logFC_soulage = logFC),
  gtex_f %>% dplyr::select(gene_symbol, logFC, adj.P.Val) %>%
    rename(logFC_gtex_f = logFC, fdr_gtex_f = adj.P.Val),
  by = "gene_symbol"
) %>%
  inner_join(
    gtex_m %>% dplyr::select(gene_symbol, logFC) %>% rename(logFC_gtex_m = logFC),
    by = "gene_symbol"
  )

long_df <- merged_all %>%
  pivot_longer(cols = c(logFC_gtex_f, logFC_gtex_m),
               names_to = "gtex_sex", values_to = "logFC_gtex") %>%
  mutate(gtex_sex = ifelse(gtex_sex == "logFC_gtex_f", "GTEx Female", "GTEx Male"))

rho_labels <- data.frame(
  gtex_sex = c("GTEx Female", "GTEx Male"),
  label = c(sprintf("rho = %.3f", cor_lfc_f$estimate),
            sprintf("rho = %.3f", cor_lfc_m$estimate))
)

p3 <- ggplot(long_df, aes(logFC_gtex, logFC_soulage)) +
  geom_point(alpha = 0.1, size = 0.5, color = "grey50") +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_text(data = rho_labels, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.3, size = 4, inherit.aes = FALSE) +
  facet_wrap(~gtex_sex) +
  labs(title = "GSE98603 (Soulage, Males, CKD V) vs GTEx CKD",
       x = "GTEx logFC", y = "GSE98603 logFC") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 13))
ggsave(file.path(plots_dir, "06_logFC_scatter_panel.png"), p3, width = 12, height = 6, dpi = 300)

# ── 6. Summary ──────────────────────────────────────────────────────────────

cat("\n\n========================================\n")
cat("=== GSE98603 vs GTEx CKD COMPARISON ===\n")
cat("========================================\n\n")

cat("STUDY DESIGN:\n")
cat("  GSE98603 (Soulage): 9 CKD-V vs 9 controls, males only, age/BMI-matched, living tissue\n")
cat("  GTEx: ~40+ CKD vs controls, sex-stratified, post-mortem subcutaneous adipose\n\n")

cat("GENE-LEVEL:\n")
cat(sprintf("  Soulage DEGs (FDR<0.05): %d\n", sum(soulage$adj.P.Val < 0.05)))
cat(sprintf("  GTEx Female DEGs (FDR<0.05): %d\n", sum(gtex_f$adj.P.Val < 0.05)))
cat(sprintf("  GTEx Male DEGs (FDR<0.05): %d\n", sum(gtex_m$adj.P.Val < 0.05)))

cat(sprintf("\nGENOME-WIDE CORRELATION (Spearman rho of t-statistics):\n"))
cat(sprintf("  Soulage vs GTEx Female: rho = %.3f (p = %.1e)\n", cor_f$estimate, cor_f$p.value))
cat(sprintf("  Soulage vs GTEx Male:   rho = %.3f (p = %.1e)\n", cor_m$estimate, cor_m$p.value))

cat(sprintf("\nLOGFC CORRELATION:\n"))
cat(sprintf("  Soulage vs GTEx Female: rho = %.3f\n", cor_lfc_f$estimate))
cat(sprintf("  Soulage vs GTEx Male:   rho = %.3f\n", cor_lfc_m$estimate))

if (!is.na(concordance)) {
  cat(sprintf("\nDIRECTION CONCORDANCE (at %d GTEx female DEGs found in Soulage):\n", nrow(soulage_at_f_degs)))
  cat(sprintf("  %.1f%%\n", 100 * concordance))
}

cat("\nGSEA PATHWAY THEMES (GSE98603 / CKD males):\n")
cat("  ACTIVATED: immune activation, antigen presentation, MHC class I, T cell activation\n")
cat("  SUPPRESSED: TGF-beta signaling, ECM organization, Wnt signaling,\n")
cat("              growth factor response, cardiac/muscle development\n")

# Save comparison
comparison <- data.frame(
  metric = c("Soulage DEGs (FDR<0.05)", "GTEx Female DEGs", "GTEx Male DEGs",
             "Overlapping genes (F)", "Overlapping genes (M)",
             "t-stat rho vs GTEx Female", "t-stat rho vs GTEx Male",
             "logFC rho vs GTEx Female", "logFC rho vs GTEx Male",
             "Direction concordance (GTEx F DEGs)"),
  value = c(sum(soulage$adj.P.Val < 0.05),
            sum(gtex_f$adj.P.Val < 0.05),
            sum(gtex_m$adj.P.Val < 0.05),
            nrow(merged_f), nrow(merged_m),
            round(cor_f$estimate, 3), round(cor_m$estimate, 3),
            round(cor_lfc_f$estimate, 3), round(cor_lfc_m$estimate, 3),
            ifelse(!is.na(concordance), sprintf("%.1f%%", 100 * concordance), "N/A"))
)
write.csv(comparison, file.path(tables_dir, "06_gtex_comparison_summary.csv"), row.names = FALSE)

cat("\n✓ GTEx comparison complete\n")
