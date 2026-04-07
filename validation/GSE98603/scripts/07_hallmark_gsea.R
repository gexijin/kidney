## 07_hallmark_gsea.R — Hallmark GSEA for GSE98603 and comparison with GTEx
## Uses msigdbr Hallmark gene sets for direct NES comparison

library(tidyverse)
library(msigdbr)
library(fgsea)
library(pheatmap)

base_dir <- "/Users/ge/Documents/research/gtex/analyses/Soulage_GSE98603"
kidney_dir <- "/Users/ge/Documents/research/gtex/kidney"
data_dir <- file.path(base_dir, "data")
plots_dir <- file.path(base_dir, "results", "plots")
tables_dir <- file.path(base_dir, "results", "tables")

# ── 1. Load Soulage DE results ──────────────────────────────────────────────

de_clean <- readRDS(file.path(data_dir, "04_de_clean.rds"))

# ── 2. Build Hallmark gene sets (Entrez IDs) ────────────────────────────────

hallmark_entrez <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  filter(!is.na(ncbi_gene) & ncbi_gene != "") %>%
  mutate(ncbi_gene = as.character(ncbi_gene))

hallmark_list <- split(hallmark_entrez$ncbi_gene, hallmark_entrez$gs_name)
cat("Hallmark gene sets:", length(hallmark_list), "\n")

# ── 3. Ranked gene list (t-statistic, Entrez IDs) ───────────────────────────

ranks <- de_clean %>%
  filter(!is.na(entrez_id)) %>%
  arrange(desc(t)) %>%
  { setNames(.$t, .$entrez_id) }
ranks <- ranks[!duplicated(names(ranks))]
cat("Ranked genes:", length(ranks), "\n")

# ── 4. Run fgsea ────────────────────────────────────────────────────────────

set.seed(42)
gsea_hallmark <- fgsea(pathways = hallmark_list, stats = ranks,
                       minSize = 15, maxSize = 500, nPermSimple = 10000)

gsea_hallmark <- gsea_hallmark %>% arrange(pval)

cat("\n=== Hallmark GSEA Results ===\n")
cat("Significant (padj < 0.05):", sum(gsea_hallmark$padj < 0.05), "\n")
cat("Activated:", sum(gsea_hallmark$padj < 0.05 & gsea_hallmark$NES > 0), "\n")
cat("Suppressed:", sum(gsea_hallmark$padj < 0.05 & gsea_hallmark$NES < 0), "\n")

cat("\nTop activated:\n")
gsea_hallmark %>%
  filter(NES > 0) %>%
  head(10) %>%
  dplyr::select(pathway, NES, pval, padj, size) %>%
  print(digits = 3)

cat("\nTop suppressed:\n")
gsea_hallmark %>%
  filter(NES < 0) %>%
  head(10) %>%
  dplyr::select(pathway, NES, pval, padj, size) %>%
  print(digits = 3)

# Save Soulage Hallmark results
gsea_out <- gsea_hallmark %>%
  dplyr::select(pathway, pval, padj, ES, NES, size) %>%
  as.data.frame()
write.csv(gsea_out, file.path(tables_dir, "07_hallmark_gsea.csv"), row.names = FALSE)

# ── 5. Load GTEx Hallmark GSEA results ──────────────────────────────────────

stage2 <- readRDS(file.path(kidney_dir, "results", "stage2_results.rds"))

gtex_female <- stage2$gsea_female_full %>%
  as.data.frame() %>%
  dplyr::select(pathway, NES, padj) %>%
  rename(NES_gtex_f = NES, padj_gtex_f = padj)

gtex_male <- stage2$gsea_male %>%
  as.data.frame() %>%
  dplyr::select(pathway, NES, padj) %>%
  rename(NES_gtex_m = NES, padj_gtex_m = padj)

soulage_gsea <- gsea_out %>%
  dplyr::select(pathway, NES, padj) %>%
  rename(NES_soulage = NES, padj_soulage = padj)

# ── 6. Merge and compare NES ────────────────────────────────────────────────

comparison <- soulage_gsea %>%
  inner_join(gtex_female, by = "pathway") %>%
  inner_join(gtex_male, by = "pathway")

cat("\nPathways in all three analyses:", nrow(comparison), "\n")

# NES correlations
cor_f <- cor.test(comparison$NES_soulage, comparison$NES_gtex_f, method = "spearman")
cor_m <- cor.test(comparison$NES_soulage, comparison$NES_gtex_m, method = "spearman")
cor_fm <- cor.test(comparison$NES_gtex_f, comparison$NES_gtex_m, method = "spearman")

cat("\n=== NES Correlations (Spearman) ===\n")
cat(sprintf("Soulage vs GTEx Female: rho = %.3f, p = %.2e\n", cor_f$estimate, cor_f$p.value))
cat(sprintf("Soulage vs GTEx Male:   rho = %.3f, p = %.2e\n", cor_m$estimate, cor_m$p.value))
cat(sprintf("GTEx Female vs Male:    rho = %.3f, p = %.2e\n", cor_fm$estimate, cor_fm$p.value))

# Direction concordance
concord_f <- mean(sign(comparison$NES_soulage) == sign(comparison$NES_gtex_f))
concord_m <- mean(sign(comparison$NES_soulage) == sign(comparison$NES_gtex_m))
cat(sprintf("\nPathway direction concordance:\n"))
cat(sprintf("  Soulage vs GTEx Female: %.0f%% (%d/%d)\n",
            100 * concord_f, sum(sign(comparison$NES_soulage) == sign(comparison$NES_gtex_f)), nrow(comparison)))
cat(sprintf("  Soulage vs GTEx Male:   %.0f%% (%d/%d)\n",
            100 * concord_m, sum(sign(comparison$NES_soulage) == sign(comparison$NES_gtex_m)), nrow(comparison)))

# ── 7. Clean pathway names ──────────────────────────────────────────────────

clean_name <- function(x) {
  x %>%
    str_remove("^HALLMARK_") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_replace("Tnfa", "TNFa") %>%
    str_replace("Nfkb", "NF-kB") %>%
    str_replace("Kras", "KRAS") %>%
    str_replace("Tgf Beta", "TGF-beta") %>%
    str_replace("Il6 Jak", "IL6-JAK") %>%
    str_replace("Il2 Stat5", "IL2-STAT5") %>%
    str_replace("Dna", "DNA") %>%
    str_replace("Uv ", "UV ") %>%
    str_replace("Mtor", "mTOR") %>%
    str_replace("Myc ", "MYC ") %>%
    str_replace("E2f", "E2F") %>%
    str_replace("P53", "p53") %>%
    str_replace("G2m", "G2M") %>%
    str_replace("Pi3k", "PI3K") %>%
    str_replace("Akt", "AKT")
}

comparison$pathway_clean <- clean_name(comparison$pathway)

# ── 8. NES heatmap (all 50 Hallmark pathways) ───────────────────────────────

nes_mat <- comparison %>%
  dplyr::select(pathway_clean, NES_gtex_f, NES_gtex_m, NES_soulage) %>%
  column_to_rownames("pathway_clean") %>%
  as.matrix()

colnames(nes_mat) <- c("GTEx Female\n(n≈218)", "GTEx Male\n(n≈363)", "Soulage Male\n(CKD-V, n=18)")

# Significance markers
padj_mat <- comparison %>%
  dplyr::select(pathway_clean, padj_gtex_f, padj_gtex_m, padj_soulage) %>%
  column_to_rownames("pathway_clean") %>%
  as.matrix()

sig_marks <- matrix("", nrow(padj_mat), ncol(padj_mat))
sig_marks[padj_mat < 0.05] <- "*"
sig_marks[padj_mat < 0.001] <- "**"

# Order by GTEx female NES
row_order <- order(nes_mat[, 1])
nes_mat <- nes_mat[row_order, ]
sig_marks <- sig_marks[row_order, ]

png(file.path(plots_dir, "07_hallmark_nes_heatmap.png"), width = 700, height = 1200, res = 150)
pheatmap(nes_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = sig_marks,
         fontsize_number = 12,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(-3, 3, length.out = 101),
         main = "Hallmark GSEA: CKD vs Control\n(* FDR<0.05, ** FDR<0.001)",
         fontsize_row = 7,
         fontsize_col = 9,
         angle_col = 0,
         cellwidth = 60,
         cellheight = 14)
dev.off()

# ── 9. NES scatter plots ────────────────────────────────────────────────────

# Soulage vs GTEx Female
p1 <- ggplot(comparison, aes(NES_gtex_f, NES_soulage)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_point(aes(size = -log10(padj_soulage),
                 color = padj_soulage < 0.05),
             alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey50"),
                     labels = c("NS", "FDR<0.05"), name = "Soulage FDR") +
  scale_size_continuous(range = c(1, 5), name = "-log10(FDR)") +
  geom_text(aes(label = pathway_clean), size = 2, vjust = -0.8,
            data = comparison %>% filter(padj_soulage < 0.05 | padj_gtex_f < 0.001)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue", alpha = 0.4) +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("rho = %.2f\np = %.1e", cor_f$estimate, cor_f$p.value),
           hjust = 1.1, vjust = -0.5, size = 4) +
  labs(title = "Hallmark NES: Soulage (Males, CKD-V) vs GTEx Female CKD",
       x = "GTEx Female NES", y = "Soulage NES") +
  theme_minimal(base_size = 11) +
  coord_equal(xlim = c(-3, 3), ylim = c(-3, 3))
ggsave(file.path(plots_dir, "07_hallmark_scatter_vs_gtex_female.png"), p1, width = 10, height = 9, dpi = 300)

# Soulage vs GTEx Male
p2 <- ggplot(comparison, aes(NES_gtex_m, NES_soulage)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_point(aes(size = -log10(padj_soulage),
                 color = padj_soulage < 0.05),
             alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey50"),
                     labels = c("NS", "FDR<0.05"), name = "Soulage FDR") +
  scale_size_continuous(range = c(1, 5), name = "-log10(FDR)") +
  geom_text(aes(label = pathway_clean), size = 2, vjust = -0.8,
            data = comparison %>% filter(padj_soulage < 0.05)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue", alpha = 0.4) +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("rho = %.2f\np = %.1e", cor_m$estimate, cor_m$p.value),
           hjust = 1.1, vjust = -0.5, size = 4) +
  labs(title = "Hallmark NES: Soulage (Males, CKD-V) vs GTEx Male CKD",
       x = "GTEx Male NES", y = "Soulage NES") +
  theme_minimal(base_size = 11) +
  coord_equal(xlim = c(-3, 3), ylim = c(-3, 3))
ggsave(file.path(plots_dir, "07_hallmark_scatter_vs_gtex_male.png"), p2, width = 10, height = 9, dpi = 300)

# ── 10. Print full comparison table ─────────────────────────────────────────

cat("\n=== Full Hallmark Comparison ===\n")
comparison %>%
  arrange(NES_gtex_f) %>%
  dplyr::select(pathway_clean, NES_gtex_f, padj_gtex_f, NES_gtex_m, padj_gtex_m,
                NES_soulage, padj_soulage) %>%
  mutate(across(starts_with("NES"), ~round(., 2)),
         across(starts_with("padj"), ~ifelse(is.na(.), 1, round(., 4)))) %>%
  as.data.frame() %>%
  print()

# Save comparison table
write.csv(comparison %>%
            arrange(NES_gtex_f) %>%
            dplyr::select(pathway, pathway_clean,
                          NES_gtex_f, padj_gtex_f,
                          NES_gtex_m, padj_gtex_m,
                          NES_soulage, padj_soulage),
          file.path(tables_dir, "07_hallmark_comparison.csv"), row.names = FALSE)

# ── 11. Summary ─────────────────────────────────────────────────────────────

cat("\n\n====================================\n")
cat("=== HALLMARK GSEA COMPARISON ===\n")
cat("====================================\n\n")

cat("NES CORRELATIONS:\n")
cat(sprintf("  Soulage vs GTEx Female: rho = %.3f (p = %.2e)\n", cor_f$estimate, cor_f$p.value))
cat(sprintf("  Soulage vs GTEx Male:   rho = %.3f (p = %.2e)\n", cor_m$estimate, cor_m$p.value))
cat(sprintf("  GTEx Female vs Male:    rho = %.3f (p = %.2e)\n", cor_fm$estimate, cor_fm$p.value))

cat(sprintf("\nDIRECTION CONCORDANCE:\n"))
cat(sprintf("  Soulage vs GTEx Female: %.0f%%\n", 100 * concord_f))
cat(sprintf("  Soulage vs GTEx Male:   %.0f%%\n", 100 * concord_m))

# Highlight discordant pathways (sig in both, opposite direction)
discordant <- comparison %>%
  filter(padj_soulage < 0.25 & padj_gtex_f < 0.05) %>%
  filter(sign(NES_soulage) != sign(NES_gtex_f))

if (nrow(discordant) > 0) {
  cat("\nDISCORDANT PATHWAYS (sig in Soulage FDR<0.25 & GTEx-F FDR<0.05, opposite direction):\n")
  discordant %>%
    dplyr::select(pathway_clean, NES_gtex_f, NES_soulage) %>%
    print()
}

# Concordant pathways
concordant <- comparison %>%
  filter(padj_soulage < 0.05 & padj_gtex_f < 0.05) %>%
  filter(sign(NES_soulage) == sign(NES_gtex_f))

if (nrow(concordant) > 0) {
  cat("\nCONCORDANT PATHWAYS (sig in both, same direction):\n")
  concordant %>%
    dplyr::select(pathway_clean, NES_gtex_f, NES_soulage) %>%
    print()
}

cat("\n✓ Hallmark GSEA comparison complete\n")
