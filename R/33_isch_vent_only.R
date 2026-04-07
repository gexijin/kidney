#!/usr/bin/env Rscript
# =============================================================================
# 33: Ischemia DE within Ventilator Deaths Only
#
# Hardy-controlled ischemia analysis: continuous ischemic time effect in
# female non-CKD Ventilator deaths only. Compares gene-level and pathway-
# level effects to the CKD signal.
# =============================================================================

library(tidyverse)
library(limma)
library(edgeR)
library(msigdbr)
library(fgsea)
library(ggrepel)
library(patchwork)

source("kidney/R/functions.R")

set.seed(42)
res_dir <- "kidney/results"
fig_dir <- "kidney/figures"

gene_anno <- readRDS("data/annotation/gencode_v47_genes.rds")
gene_map <- setNames(gene_anno$symbol, rownames(gene_anno))

# ── Load data ────────────────────────────────────────────────
dat <- load_sat_data(res_dir = res_dir)
tb <- dat$tb
counts <- dat$counts
tb$MHT2D[is.na(tb$MHT2D)] <- 0L

all_vars <- c("AGE", "RACE", "BMI", "DTHHRDY", "ischemic_hrs", "SMRIN",
              "SMCENTER", "SMEXNCRT", "MHABNWBC", "MHLVRDIS", "MHRNLFLR", "MHT2D")
tb <- tb[complete.cases(tb[, all_vars]), ]

# ── Ventilator-only non-CKD controls ─────────────────────────
vent <- tb %>%
  filter(MHRNLFLR == 0, DTHHRDY == "Ventilator")

cat(sprintf("Ventilator non-CKD controls: n=%d\n", nrow(vent)))
cat(sprintf("  Female: %d | Male: %d\n",
            sum(vent$SEX == "Female"),
            sum(vent$SEX == "Male")))
cat(sprintf("  Ischemic time: median=%.1fh, range=%.1f-%.1fh\n",
            median(vent$ischemic_hrs),
            min(vent$ischemic_hrs),
            max(vent$ischemic_hrs)))

# =============================================================================
# DE: Continuous ischemic time (Ventilator only)
# =============================================================================

cat("\n=== DE: Continuous ischemic time (Ventilator only) ===\n")

form <- ~ ischemic_hrs + SEX + AGE + RACE + BMI + SMRIN +
          SMCENTER + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D
design <- build_design(form, vent)
cat(sprintf("Design: %d samples x %d cols\n",
            nrow(design), ncol(design)))

res_isch <- run_voom_de(counts[, vent$SAMPID],
                        design, "ischemic_hrs")
tt_isch <- res_isch$tt

n_deg <- sum(tt_isch$adj.P.Val < 0.05)
cat(sprintf("Ischemia DEGs (FDR<0.05): %d\n", n_deg))
cat(sprintf("  Up: %d | Down: %d\n",
            sum(tt_isch$adj.P.Val < 0.05 & tt_isch$logFC > 0),
            sum(tt_isch$adj.P.Val < 0.05 & tt_isch$logFC < 0)))
cat(sprintf("  Min FDR: %.2e | Median p: %.4f\n",
            min(tt_isch$adj.P.Val), median(tt_isch$P.Value)))

# =============================================================================
# GSEA: Hallmark
# =============================================================================

cat("\n=== GSEA: Hallmark ===\n")

hallmark_list <- load_hallmark_sets()
gsea_isch <- run_gsea(tt_isch, hallmark_list, label = "Ischemia (Vent only)", seed = 42)
cat(sprintf("Significant pathways (FDR<0.05): %d\n", sum(gsea_isch$padj < 0.05)))

# =============================================================================
# Compare to CKD
# =============================================================================

cat("\n=== Compare to CKD ===\n")

de_fem <- readRDS(file.path(res_dir, "de_female.rds"))
tt_ckd <- de_fem$tt
gsea_ckd <- run_gsea(tt_ckd, hallmark_list, label = "CKD (Female)", seed = 42)

# ── Gene-level comparison ────────────────────────────────────
merged_gene <- inner_join(
  tt_ckd %>% select(gene_id, logFC_ckd = logFC, t_ckd = t, FDR_ckd = adj.P.Val),
  tt_isch %>% select(gene_id, logFC_isch = logFC, t_isch = t, FDR_isch = adj.P.Val),
  by = "gene_id"
)
merged_gene$symbol <- gene_map[merged_gene$gene_id]

r_lfc <- cor(merged_gene$logFC_ckd, merged_gene$logFC_isch)
r_t <- cor(merged_gene$t_ckd, merged_gene$t_isch)
cat(sprintf("Gene logFC r (CKD vs isch): %.3f\n", r_lfc))
cat(sprintf("Gene t-stat r (CKD vs isch): %.3f\n", r_t))

# Direction concordance among CKD DEGs
ckd_sig <- merged_gene %>% filter(FDR_ckd < 0.05)
n_conc <- sum(sign(ckd_sig$logFC_ckd) == sign(ckd_sig$logFC_isch))
cat(sprintf(
  "CKD DEGs direction concordance: %d/%d (%.1f%%)\n",
  n_conc, nrow(ckd_sig), 100 * n_conc / nrow(ckd_sig)
))

merged_gene <- merged_gene %>%
  mutate(
    sig_ckd = FDR_ckd < 0.05, sig_isch = FDR_isch < 0.05,
    category = case_when(
      sig_ckd & sig_isch  ~ "Both",
      sig_ckd & !sig_isch ~ "CKD only",
      !sig_ckd & sig_isch ~ "Ischemia only",
      TRUE                ~ "Neither"
    ),
    category = factor(category,
      levels = c("Neither", "CKD only",
                 "Ischemia only", "Both"))
  )

# ── Pathway-level comparison ─────────────────────────────────
gsea_merged <- inner_join(
  gsea_ckd %>% as.data.frame() %>% select(pathway, NES_CKD = NES, padj_CKD = padj),
  gsea_isch %>% as.data.frame() %>% select(pathway, NES_isch = NES, padj_isch = padj),
  by = "pathway"
)
r_nes <- cor(gsea_merged$NES_CKD, gsea_merged$NES_isch)
cat(sprintf("\nPathway NES r (CKD vs isch): %.3f\n", r_nes))

# =============================================================================
# Combined figure: gene scatter (A) + pathway NES scatter (B)
# =============================================================================

cat("\nGenerating combined figure...\n")

# ── Panel A: Gene-level density scatter ──────────────────────
marker_genes <- c("LEP", "ADIPOQ", "PPARG", "PLIN1",
                  "IL6", "CCL2", "FASN", "LPL")
merged_gene$label <- ifelse(
  merged_gene$symbol %in% marker_genes,
  merged_gene$symbol, NA
)

p1 <- ggplot(merged_gene, aes(x = logFC_isch, y = logFC_ckd)) +
  geom_hex(bins = 80, alpha = 0.9) +
  scale_fill_viridis_c(option = "inferno", trans = "log10",
                       guide = "none") +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_text_repel(
    aes(label = label), size = 3, na.rm = TRUE,
    color = "black", fontface = "italic",
    segment.size = 0.3, segment.color = "grey40",
    box.padding = 0.4, max.overlaps = 20, seed = 42
  ) +
  annotate("text", x = -0.14, y = 2.8,
           label = sprintf("italic(r) == %.2f", r_lfc),
           parse = TRUE, size = 3.5, hjust = 0) +
  labs(x = expression(log[2]*FC/hr ~ "(Ischemia)"),
       y = expression(log[2]*FC ~ "(CKD)")) +
  xlim(-0.15, 0.15) + ylim(-3, 3) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())

# ── Panel B: Pathway NES scatter ─────────────────────────────
gsea_plot <- gsea_merged %>%
  filter(
    (padj_CKD < 0.01 | padj_isch < 0.01),
    abs(NES_CKD) > 1.5 | abs(NES_isch) > 1.5
  ) %>%
  mutate(
    pathway = clean_pathway(pathway),
    sig_cat = case_when(
      padj_CKD < 0.01 & padj_isch < 0.01 ~ "Both",
      padj_CKD < 0.01 ~ "CKD only",
      TRUE ~ "Ischemia only"
    ),
    sig_cat = factor(sig_cat,
      levels = c("Both", "CKD only", "Ischemia only"))
  )

p2 <- ggplot(gsea_plot,
             aes(x = NES_isch, y = NES_CKD, color = sig_cat)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 0.4) +
  geom_point(size = 2, alpha = 0.85) +
  geom_text_repel(
    aes(label = pathway), size = 2.5,
    max.overlaps = 50, seed = 123,
    segment.size = 0.15, segment.color = "grey60",
    box.padding = 0.4, point.padding = 0.2,
    min.segment.length = 0, force = 3, force_pull = 0.5,
    max.iter = 10000
  ) +
  scale_color_manual(
    values = c("Both" = "#D6604D", "CKD only" = "#4393C3",
               "Ischemia only" = "#E69F00"),
    name = NULL
  ) +
  annotate("text", x = min(gsea_plot$NES_isch) * 0.95,
           y = max(gsea_plot$NES_CKD) * 0.95,
           label = sprintf("italic(r) == %.2f", r_nes),
           parse = TRUE, size = 3.5, hjust = 0) + {
    lim <- range(c(gsea_plot$NES_isch, gsea_plot$NES_CKD))
    lim <- c(floor(lim[1]), ceiling(lim[2]))
    coord_fixed(xlim = lim, ylim = lim)
  } +
  labs(x = "Enrichment score (Ischemia)",
       y = "Enrichment score (CKD)") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.82, 0.12),
    legend.background = element_rect(
      fill = alpha("white", 0.9),
      linewidth = 0.3, color = "grey80"),
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 8)
  )

combined <- p1 + p2 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave(file.path(fig_dir, "fig3_ckd_vs_isch_vent.png"),
       combined, width = 10, height = 5, dpi = 300)
cat("Saved:", file.path(fig_dir, "fig3_ckd_vs_isch_vent.png"), "\n")

# =============================================================================
# Save
# =============================================================================

saveRDS(list(
  tt_isch = tt_isch, gsea_isch = gsea_isch,
  gsea_ckd = gsea_ckd, merged_gene = merged_gene,
  gsea_merged = gsea_merged,
  r_lfc = r_lfc, r_t = r_t, r_nes = r_nes,
  n_deg = n_deg
), file.path(res_dir, "isch_vent_only.rds"))
cat("Saved:", file.path(res_dir, "isch_vent_only.rds"), "\n")

# ── Summary ──────────────────────────────────────────────────
n_ckd_deg <- sum(merged_gene$sig_ckd)
ckd_conc <- merged_gene %>%
  filter(sig_ckd) %>%
  summarise(
    n = n(),
    same_dir = sum(sign(logFC_ckd) == sign(logFC_isch)),
    pct = 100 * same_dir / n
  )

cat(sprintf("\n=== Summary ===\n"))
cat(sprintf(
  "Gene-level: logFC r=%.3f, t-stat r=%.3f\n",
  r_lfc, r_t
))
cat(sprintf(
  "CKD DEGs with same isch direction: %d/%d (%.1f%%)\n",
  ckd_conc$same_dir, ckd_conc$n, ckd_conc$pct
))
cat(sprintf("Pathway NES r: %.3f\n", r_nes))
cat("\nDone.\n")
