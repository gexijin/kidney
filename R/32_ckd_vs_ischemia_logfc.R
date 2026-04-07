#!/usr/bin/env Rscript
# =============================================================================
# 32: CKD vs Ischemia Log2FC Comparison (Paired Design)
#
# Disentangles CKD transcriptomic signal from ischemic time confound using
# two paired limma/voom contrasts in female subcutaneous adipose:
#   A: Female CKD vs Female Control (both long ischemia) — CKD effect
#   B: Female Lo-Isch Old vs Female Control — ischemia effect
#
# Both contrasts use paired designs (~ pair + group) with 28 samples per arm.
# =============================================================================

library(tidyverse)
library(edgeR)
library(limma)
library(MatchIt)
library(ggrepel)
library(msigdbr)
library(fgsea)
library(pheatmap)
library(RColorBrewer)

source("kidney/R/functions.R")

set.seed(42)

res_dir <- "kidney/results"
fig_dir <- "kidney/figures"

# ── Gene annotation ──────────────────────────────────────────
gene_anno <- readRDS("data/annotation/gencode_v47_genes.rds")
gene_map  <- setNames(gene_anno$symbol, rownames(gene_anno))

# =============================================================================
# Step 1: Data loading & 3-group construction
# =============================================================================

cat("=== Step 1: Data loading & group construction ===\n")

dat <- load_sat_data(res_dir = res_dir)
tb <- dat$tb
counts <- dat$counts
tb$MHT2D[is.na(tb$MHT2D)] <- 0L  # rescue GTEX-1CB4H

all_vars <- c("AGE", "RACE", "BMI", "DTHHRDY", "ischemic_hrs", "SMRIN",
              "SMCENTER", "SMEXNCRT", "MHABNWBC", "MHLVRDIS", "MHRNLFLR", "MHT2D")
tb <- tb[complete.cases(tb[, all_vars]), ]
counts <- counts[, tb$SAMPID]

# ── G1: Female CKD (all available) ──────────────────────────
f_ckd_df  <- tb %>% filter(SEX == "Female", MHRNLFLR == 1)
f_ckd_ids <- f_ckd_df$SAMPID
cat(sprintf("G1 Female CKD: n=%d\n", length(f_ckd_ids)))

# ── G2: Female Control (PSM, ratio=1) ───────────────────────
fem <- tb %>% filter(SEX == "Female", !is.na(MHRNLFLR))
fem <- add_ischemic_bin(fem)

set.seed(42)
psm_ckd_ctrl <- matchit(
  MHRNLFLR ~ AGE + RACE + DTHHRDY + ischemic_hrs + SMRIN,
  data = fem, method = "nearest", ratio = 1,
  exact = ~isch_bin
)
psm_data <- match.data(psm_ckd_ctrl)
f_ctrl_ids <- psm_data %>% filter(MHRNLFLR == 0) %>% pull(SAMPID)
cat(sprintf("G2 Female Control (PSM): n=%d\n", length(f_ctrl_ids)))

# ── G3: Female Lo-Isch Old (Mahalanobis matched to CKD) ─────
used_ids <- c(f_ckd_ids, f_ctrl_ids)
lo_fem_old <- tb %>%
  filter(SEX == "Female", MHRNLFLR == 0, ischemic_hrs <= 3, AGE >= 50,
         !SAMPID %in% used_ids) %>%
  mutate(RACE = factor(RACE))

combined_lo <- bind_rows(
  f_ckd_df %>% mutate(treat = 1, RACE = factor(RACE)),
  lo_fem_old %>% mutate(treat = 0)
)

set.seed(42)
match_lo <- matchit(
  treat ~ BMI + SMRIN,
  data = combined_lo, method = "nearest",
  distance = "mahalanobis",
  exact = ~RACE, ratio = 1, replace = FALSE
)
lo_data <- match.data(match_lo)
f_lo_old_ids <- lo_data %>% filter(treat == 0) %>% pull(SAMPID)
cat(sprintf("G3 Female Lo-Isch Old: n=%d\n", length(f_lo_old_ids)))

stopifnot(length(f_ckd_ids) == 28,
          length(f_ctrl_ids) == 28,
          length(f_lo_old_ids) == 28)

# =============================================================================
# Step 2: Contrast A — CKD effect (paired)
# =============================================================================

cat("\n=== Step 2: Contrast A — CKD vs Control (paired) ===\n")

# Extract pair assignments from PSM subclass
pairs_a <- psm_data %>%
  select(SAMPID, MHRNLFLR, subclass) %>%
  arrange(subclass, desc(MHRNLFLR))

# Verify 1:1 pairing
pair_counts_a <- pairs_a %>% count(subclass)
stopifnot(all(pair_counts_a$n == 2))
cat(sprintf("Pairs: %d (all 1:1)\n", nrow(pair_counts_a)))

# Build sample metadata for contrast A
samp_a <- pairs_a %>%
  mutate(
    group = factor(ifelse(MHRNLFLR == 1, "CKD", "Control"), levels = c("Control", "CKD")),
    pair  = factor(subclass)
  )

# Design matrix with pair blocking
design_a <- model.matrix(~ pair + group, data = samp_a)

# Limma/voom
counts_a <- counts[, samp_a$SAMPID]
dge_a <- DGEList(counts = counts_a)
dge_a <- calcNormFactors(dge_a, method = "TMM")
v_a   <- voom(dge_a, design_a, plot = FALSE)
fit_a <- eBayes(lmFit(v_a, design_a))
tt_a  <- topTable(fit_a, coef = "groupCKD", number = Inf, sort.by = "none")
tt_a$gene_id <- rownames(tt_a)

n_deg_a <- sum(tt_a$adj.P.Val < 0.05)
cat(sprintf("CKD effect: %d DEGs (FDR<0.05)\n", n_deg_a))
cat(sprintf("  Up: %d | Down: %d\n",
            sum(tt_a$adj.P.Val < 0.05 & tt_a$logFC > 0),
            sum(tt_a$adj.P.Val < 0.05 & tt_a$logFC < 0)))

# =============================================================================
# Step 3: Contrast B — Ischemia effect (paired)
# =============================================================================

cat("\n=== Step 3: Contrast B — Control vs Lo-Isch (paired, ischemia effect) ===\n")

# New matching: Lo-Isch Old vs Control on AGE + BMI + SMRIN
pool_b <- bind_rows(
  tb %>% filter(SAMPID %in% f_lo_old_ids) %>% mutate(treat = 1),
  tb %>% filter(SAMPID %in% f_ctrl_ids)   %>% mutate(treat = 0)
)

set.seed(42)
match_b <- matchit(
  treat ~ AGE + BMI + SMRIN,
  data = pool_b, method = "nearest",
  distance = "mahalanobis",
  ratio = 1, replace = FALSE
)
pairs_b_data <- match.data(match_b)

pair_counts_b <- pairs_b_data %>% count(subclass)
stopifnot(all(pair_counts_b$n == 2))
cat(sprintf("Pairs: %d (all 1:1)\n", nrow(pair_counts_b)))

# Build sample metadata for contrast B
samp_b <- pairs_b_data %>%
  select(SAMPID, treat, subclass) %>%
  mutate(
    group = factor(ifelse(treat == 1, "LoIsch", "Control"), levels = c("LoIsch", "Control")),
    pair  = factor(subclass)
  )

# Design matrix with pair blocking
design_b <- model.matrix(~ pair + group, data = samp_b)

# Limma/voom
counts_b <- counts[, samp_b$SAMPID]
dge_b <- DGEList(counts = counts_b)
dge_b <- calcNormFactors(dge_b, method = "TMM")
v_b   <- voom(dge_b, design_b, plot = FALSE)
fit_b <- eBayes(lmFit(v_b, design_b))
tt_b  <- topTable(fit_b, coef = "groupControl", number = Inf, sort.by = "none")
tt_b$gene_id <- rownames(tt_b)

n_deg_b <- sum(tt_b$adj.P.Val < 0.05)
cat(sprintf("Ischemia effect: %d DEGs (FDR<0.05)\n", n_deg_b))
cat(sprintf("  Up: %d | Down: %d\n",
            sum(tt_b$adj.P.Val < 0.05 & tt_b$logFC > 0),
            sum(tt_b$adj.P.Val < 0.05 & tt_b$logFC < 0)))

# =============================================================================
# Step 3b: Contrast C — Moderate ischemia effect (Med-Isch vs Lo-Isch, paired)
# =============================================================================

cat("\n=== Step 3b: Contrast C — Med-Isch vs Lo-Isch (paired) ===\n")

# Pool: female, non-CKD, 3-10h ischemia, age>=45, not already used
used_ids_all <- c(f_ckd_ids, f_ctrl_ids, f_lo_old_ids)
med_isch_pool <- tb %>%
  filter(SEX == "Female", MHRNLFLR == 0,
         ischemic_hrs > 3, ischemic_hrs <= 10,
         AGE >= 45,
         !SAMPID %in% used_ids_all)
cat(sprintf("Medium-ischemia pool: n=%d (median isch=%.1fh)\n",
            nrow(med_isch_pool), median(med_isch_pool$ischemic_hrs)))

# Match Med-Isch to Lo-Isch Old on AGE + BMI + SMRIN
pool_c <- bind_rows(
  tb %>% filter(SAMPID %in% f_lo_old_ids) %>% mutate(treat = 1),
  med_isch_pool %>% mutate(treat = 0)
)

set.seed(42)
match_c <- matchit(
  treat ~ AGE + BMI + SMRIN,
  data = pool_c, method = "nearest",
  distance = "mahalanobis",
  ratio = 1, replace = FALSE
)
pairs_c_data <- match.data(match_c)
f_med_ids <- pairs_c_data %>% filter(treat == 0) %>% pull(SAMPID)
cat(sprintf("G4 Female Med-Isch: n=%d\n", length(f_med_ids)))
cat(sprintf("  Ischemic time: median=%.1fh, range=[%.1f, %.1f]\n",
            median(tb$ischemic_hrs[tb$SAMPID %in% f_med_ids]),
            min(tb$ischemic_hrs[tb$SAMPID %in% f_med_ids]),
            max(tb$ischemic_hrs[tb$SAMPID %in% f_med_ids])))

pair_counts_c <- pairs_c_data %>% count(subclass)
stopifnot(all(pair_counts_c$n == 2))
cat(sprintf("Pairs: %d (all 1:1)\n", nrow(pair_counts_c)))

# Build sample metadata — Med-Isch vs Lo-Isch (Lo-Isch as reference)
samp_c <- pairs_c_data %>%
  select(SAMPID, treat, subclass) %>%
  mutate(
    group = factor(ifelse(treat == 0, "MedIsch", "LoIsch"), levels = c("LoIsch", "MedIsch")),
    pair  = factor(subclass)
  )

design_c <- model.matrix(~ pair + group, data = samp_c)

counts_c <- counts[, samp_c$SAMPID]
dge_c <- DGEList(counts = counts_c)
dge_c <- calcNormFactors(dge_c, method = "TMM")
v_c   <- voom(dge_c, design_c, plot = FALSE)
fit_c <- eBayes(lmFit(v_c, design_c))
tt_c  <- topTable(fit_c, coef = "groupMedIsch", number = Inf, sort.by = "none")
tt_c$gene_id <- rownames(tt_c)

n_deg_c <- sum(tt_c$adj.P.Val < 0.05)
cat(sprintf("Moderate ischemia effect: %d DEGs (FDR<0.05)\n", n_deg_c))
cat(sprintf("  Up: %d | Down: %d\n",
            sum(tt_c$adj.P.Val < 0.05 & tt_c$logFC > 0),
            sum(tt_c$adj.P.Val < 0.05 & tt_c$logFC < 0)))

# =============================================================================
# Step 4: Merge & scatter plot
# =============================================================================

cat("\n=== Step 4: Log2FC comparison ===\n")

merged <- inner_join(
  tt_a %>% select(gene_id, logFC_ckd = logFC, t_ckd = t,
                  P_ckd = P.Value, FDR_ckd = adj.P.Val),
  tt_b %>% select(gene_id, logFC_isch = logFC, t_isch = t,
                  P_isch = P.Value, FDR_isch = adj.P.Val),
  by = "gene_id"
)
merged$symbol <- gene_map[merged$gene_id]

# Significance categories
merged <- merged %>%
  mutate(
    sig_ckd  = FDR_ckd < 0.05,
    sig_isch = FDR_isch < 0.05,
    category = case_when(
      sig_ckd & sig_isch  ~ "Both",
      sig_ckd & !sig_isch ~ "CKD only",
      !sig_ckd & sig_isch ~ "Ischemia only",
      TRUE                ~ "Neither"
    ),
    category = factor(category, levels = c("Neither", "CKD only", "Ischemia only", "Both"))
  )

cat(sprintf("Genes compared: %d\n", nrow(merged)))
cat("\nSignificance breakdown:\n")
print(table(merged$category))

# Correlations
rho_all <- cor(merged$logFC_ckd, merged$logFC_isch, method = "spearman")
cat(sprintf("\nSpearman rho (all genes): %.3f\n", rho_all))

sig_genes <- merged %>% filter(category != "Neither")
if (nrow(sig_genes) > 0) {
  rho_sig <- cor(sig_genes$logFC_ckd, sig_genes$logFC_isch, method = "spearman")
  cat(sprintf("Spearman rho (any DEG):  %.3f\n", rho_sig))
}

# Quadrant analysis for DEGs
both_sig <- merged %>% filter(category == "Both")
if (nrow(both_sig) > 0) {
  concordant <- sum(sign(both_sig$logFC_ckd) == sign(both_sig$logFC_isch))
  cat(sprintf("\nBoth-significant genes: %d\n", nrow(both_sig)))
  cat(sprintf("  Concordant direction: %d (%.1f%%)\n",
              concordant, 100 * concordant / nrow(both_sig)))
  cat(sprintf("  Discordant direction: %d (%.1f%%)\n",
              nrow(both_sig) - concordant,
              100 * (nrow(both_sig) - concordant) / nrow(both_sig)))
}

# ── Scatter plot ─────────────────────────────────────────────
cat("\nGenerating scatter plot...\n")

marker_genes <- c("LEP", "ADIPOQ", "PPARG", "PLIN1", "FASN", "LPL",
                   "IL6", "CCL2", "CXCL10", "EGR1")
merged$label <- ifelse(merged$symbol %in% marker_genes, merged$symbol, NA)

cat_colors <- c(
  "Neither"       = "grey80",
  "CKD only"      = "#D6604D",
  "Ischemia only" = "#4393C3",
  "Both"          = "#7B3294"
)

p <- ggplot(merged, aes(x = logFC_ckd, y = logFC_isch, color = category)) +
  geom_point(data = merged %>% filter(category == "Neither"),
             size = 0.5, alpha = 0.3) +
  geom_point(data = merged %>% filter(category != "Neither"),
             size = 1.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  geom_text_repel(
    aes(label = label),
    size = 3.5, max.overlaps = 20, na.rm = TRUE,
    color = "black", fontface = "italic",
    segment.color = "grey50", segment.size = 0.3
  ) +
  scale_color_manual(values = cat_colors, name = "Significant in") +
  labs(
    x = expression(log[2]*FC ~ "(CKD vs Control)"),
    y = expression(log[2]*FC ~ "(Control vs Lo-Isch)"),
    title = "CKD vs Ischemia: Gene-Level Effect Decomposition",
    subtitle = sprintf("Female subcutaneous adipose (n=28 per group) | rho = %.2f", rho_all)
  ) +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = alpha("white", 0.9)),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig_ckd_vs_ischemia_logfc.png"), p,
       width = 9, height = 9, dpi = 300)
cat("Saved:", file.path(fig_dir, "fig_ckd_vs_ischemia_logfc.png"), "\n")

# =============================================================================
# Step 4b: CKD vs Moderate Ischemia scatter
# =============================================================================

cat("\n=== Step 4b: CKD vs Moderate Ischemia ===\n")

merged2 <- inner_join(
  tt_a %>% select(gene_id, logFC_ckd = logFC, t_ckd = t,
                  P_ckd = P.Value, FDR_ckd = adj.P.Val),
  tt_c %>% select(gene_id, logFC_med = logFC, t_med = t,
                  P_med = P.Value, FDR_med = adj.P.Val),
  by = "gene_id"
)
merged2$symbol <- gene_map[merged2$gene_id]

merged2 <- merged2 %>%
  mutate(
    sig_ckd = FDR_ckd < 0.05,
    sig_med = FDR_med < 0.05,
    category = case_when(
      sig_ckd & sig_med  ~ "Both",
      sig_ckd & !sig_med ~ "CKD only",
      !sig_ckd & sig_med ~ "Ischemia only",
      TRUE               ~ "Neither"
    ),
    category = factor(category, levels = c("Neither", "CKD only", "Ischemia only", "Both"))
  )

cat(sprintf("Genes compared: %d\n", nrow(merged2)))
cat("\nSignificance breakdown:\n")
print(table(merged2$category))

rho_all2 <- cor(merged2$logFC_ckd, merged2$logFC_med, method = "spearman")
cat(sprintf("\nSpearman rho (all genes): %.3f\n", rho_all2))

sig_genes2 <- merged2 %>% filter(category != "Neither")
if (nrow(sig_genes2) > 0) {
  rho_sig2 <- cor(sig_genes2$logFC_ckd, sig_genes2$logFC_med, method = "spearman")
  cat(sprintf("Spearman rho (any DEG):  %.3f\n", rho_sig2))
}

both_sig2 <- merged2 %>% filter(category == "Both")
if (nrow(both_sig2) > 0) {
  concordant2 <- sum(sign(both_sig2$logFC_ckd) == sign(both_sig2$logFC_med))
  cat(sprintf("\nBoth-significant genes: %d\n", nrow(both_sig2)))
  cat(sprintf("  Concordant direction: %d (%.1f%%)\n",
              concordant2, 100 * concordant2 / nrow(both_sig2)))
  cat(sprintf("  Discordant direction: %d (%.1f%%)\n",
              nrow(both_sig2) - concordant2,
              100 * (nrow(both_sig2) - concordant2) / nrow(both_sig2)))
}

# ── Scatter plot 2 ───────────────────────────────────────────
cat("\nGenerating moderate ischemia scatter plot...\n")

merged2$label <- ifelse(merged2$symbol %in% marker_genes, merged2$symbol, NA)

p2 <- ggplot(merged2, aes(x = logFC_ckd, y = logFC_med, color = category)) +
  geom_point(data = merged2 %>% filter(category == "Neither"),
             size = 0.5, alpha = 0.3) +
  geom_point(data = merged2 %>% filter(category != "Neither"),
             size = 1.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  geom_text_repel(
    aes(label = label),
    size = 3.5, max.overlaps = 20, na.rm = TRUE,
    color = "black", fontface = "italic",
    segment.color = "grey50", segment.size = 0.3
  ) +
  scale_color_manual(values = cat_colors, name = "Significant in") +
  labs(
    x = expression(log[2]*FC ~ "(CKD vs Control)"),
    y = expression(log[2]*FC ~ "(Med-Isch vs Lo-Isch)"),
    title = "CKD vs Moderate Ischemia: Gene-Level Effect Decomposition",
    subtitle = sprintf("Female subcutaneous adipose (n=28 per group) | rho = %.2f", rho_all2)
  ) +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = alpha("white", 0.9)),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig_ckd_vs_med_ischemia_logfc.png"), p2,
       width = 9, height = 9, dpi = 300)
cat("Saved:", file.path(fig_dir, "fig_ckd_vs_med_ischemia_logfc.png"), "\n")

# =============================================================================
# Step 4c: Summary comparison across contrasts
# =============================================================================

cat("\n=== Contrast comparison summary ===\n")
cat(sprintf("%-30s %6s %8s\n", "Contrast", "DEGs", "rho_CKD"))
cat(strrep("-", 50), "\n")
cat(sprintf("%-30s %6d %8s\n", "A: CKD vs Control", n_deg_a, "—"))
cat(sprintf("%-30s %6d %8.3f\n", "B: Control vs Lo-Isch (large)", n_deg_b, rho_all))
cat(sprintf("%-30s %6d %8.3f\n", "C: Med-Isch vs Lo-Isch (mod)", n_deg_c, rho_all2))

# =============================================================================
# Step 5: Contrast D — Control vs Med-Isch (long vs medium ischemia, paired)
# =============================================================================

cat("\n=== Step 5: Contrast D — Control vs Med-Isch (paired) ===\n")

# Match Control (n=28) to Med-Isch (n=24) on AGE + BMI + SMRIN → 24 pairs
pool_d <- bind_rows(
  tb %>% filter(SAMPID %in% f_ctrl_ids)  %>% mutate(treat = 1),
  tb %>% filter(SAMPID %in% f_med_ids)   %>% mutate(treat = 0)
)

set.seed(42)
match_d <- matchit(
  treat ~ AGE + BMI + SMRIN,
  data = pool_d, method = "nearest",
  distance = "mahalanobis",
  ratio = 1, replace = FALSE
)
pairs_d_data <- match.data(match_d)

pair_counts_d <- pairs_d_data %>% count(subclass)
stopifnot(all(pair_counts_d$n == 2))
cat(sprintf("Pairs: %d (all 1:1)\n", nrow(pair_counts_d)))

# Control vs Med-Isch (Med-Isch as reference → positive logFC = higher in long ischemia)
samp_d <- pairs_d_data %>%
  select(SAMPID, treat, subclass) %>%
  mutate(
    group = factor(ifelse(treat == 1, "Control", "MedIsch"), levels = c("MedIsch", "Control")),
    pair  = factor(subclass)
  )

design_d <- model.matrix(~ pair + group, data = samp_d)

counts_d <- counts[, samp_d$SAMPID]
dge_d <- DGEList(counts = counts_d)
dge_d <- calcNormFactors(dge_d, method = "TMM")
v_d   <- voom(dge_d, design_d, plot = FALSE)
fit_d <- eBayes(lmFit(v_d, design_d))
tt_d  <- topTable(fit_d, coef = "groupControl", number = Inf, sort.by = "none")
tt_d$gene_id <- rownames(tt_d)

n_deg_d <- sum(tt_d$adj.P.Val < 0.05)
cat(sprintf("Long vs Medium ischemia: %d DEGs (FDR<0.05)\n", n_deg_d))
cat(sprintf("  Up: %d | Down: %d\n",
            sum(tt_d$adj.P.Val < 0.05 & tt_d$logFC > 0),
            sum(tt_d$adj.P.Val < 0.05 & tt_d$logFC < 0)))

# =============================================================================
# Step 6: Pathway GSEA — 3 contrasts
# =============================================================================

cat("\n=== Step 6: Pathway GSEA (Hallmark) ===\n")

hallmark_list <- load_hallmark_sets()
cat(sprintf("Hallmark gene sets: %d\n", length(hallmark_list)))

gsea_ckd  <- run_gsea(tt_a, hallmark_list, label = "CKD vs Control",     seed = 42)
gsea_med  <- run_gsea(tt_c, hallmark_list, label = "Med vs Short Isch",  seed = 42)
gsea_long <- run_gsea(tt_d, hallmark_list, label = "Long vs Med Isch",   seed = 42)

cat(sprintf("CKD vs Control:     %d pathways FDR<0.05\n", sum(gsea_ckd$padj < 0.05)))
cat(sprintf("Med vs Short Isch:  %d pathways FDR<0.05\n", sum(gsea_med$padj < 0.05)))
cat(sprintf("Long vs Med Isch:   %d pathways FDR<0.05\n", sum(gsea_long$padj < 0.05)))

# =============================================================================
# Step 7: Pathway NES heatmap
# =============================================================================

cat("\n=== Step 7: Pathway NES heatmap ===\n")

# Identify pathways significant in at least one contrast
sig_paths <- unique(c(
  gsea_ckd$pathway[gsea_ckd$padj < 0.05],
  gsea_med$pathway[gsea_med$padj < 0.05],
  gsea_long$pathway[gsea_long$padj < 0.05]
))
cat(sprintf("Pathways significant in >= 1 contrast: %d\n", length(sig_paths)))

# Build NES and padj matrices
col_names <- c("CKD vs\nControl", "Med-Isch vs\nLo-Isch", "Long-Isch vs\nMed-Isch")

nes_matrix <- data.frame(
  row.names = sig_paths,
  V1 = gsea_ckd$NES[match(sig_paths, gsea_ckd$pathway)],
  V2 = gsea_med$NES[match(sig_paths, gsea_med$pathway)],
  V3 = gsea_long$NES[match(sig_paths, gsea_long$pathway)]
)
colnames(nes_matrix) <- col_names

padj_matrix <- data.frame(
  row.names = sig_paths,
  V1 = gsea_ckd$padj[match(sig_paths, gsea_ckd$pathway)],
  V2 = gsea_med$padj[match(sig_paths, gsea_med$pathway)],
  V3 = gsea_long$padj[match(sig_paths, gsea_long$pathway)]
)
colnames(padj_matrix) <- col_names

# Clean pathway names
rownames(nes_matrix) <- clean_pathway(rownames(nes_matrix))
rownames(padj_matrix) <- clean_pathway(rownames(padj_matrix))

# Order by CKD NES
row_order <- order(nes_matrix[, 1], decreasing = TRUE)
nes_matrix <- nes_matrix[row_order, ]
padj_matrix <- padj_matrix[row_order, ]

# Significance stars
heatmap_mat <- as.matrix(nes_matrix)
sig_mat <- as.matrix(padj_matrix)
sig_labels <- matrix("", nrow = nrow(sig_mat), ncol = ncol(sig_mat))
sig_labels[sig_mat < 0.05]  <- "*"
sig_labels[sig_mat < 0.01]  <- "**"
sig_labels[sig_mat < 0.001] <- "***"

# Color scale
max_abs <- max(abs(heatmap_mat), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)
colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# Column annotation
col_annotation <- data.frame(
  Type = c("CKD effect", "Ischemia (moderate)", "Ischemia (large)"),
  row.names = col_names
)
ann_colors <- list(Type = c(
  "CKD effect"          = "#D6604D",
  "Ischemia (moderate)" = "#92C5DE",
  "Ischemia (large)"    = "#2166AC"
))

png(file.path(fig_dir, "fig_ckd_ischemia_pathway_heatmap.png"),
    width = 9, height = max(5.5, 0.4 * nrow(heatmap_mat) + 2),
    units = "in", res = 300)
pheatmap(
  heatmap_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = sig_labels,
  number_color = "black",
  fontsize_number = 12,
  color = colors,
  breaks = breaks,
  annotation_col = col_annotation,
  annotation_colors = ann_colors,
  main = "Pathway NES: CKD vs Ischemia Decomposition",
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  cellwidth = 50,
  cellheight = 20
)
dev.off()
cat("Saved:", file.path(fig_dir, "fig_ckd_ischemia_pathway_heatmap.png"), "\n")

# =============================================================================
# Step 8: Save results
# =============================================================================

cat("\n=== Step 8: Saving results ===\n")

results <- list(
  merged_large  = merged,
  merged_mod    = merged2,
  tt_ckd        = tt_a,
  tt_isch_large = tt_b,
  tt_isch_mod   = tt_c,
  tt_isch_long_med = tt_d,
  pairs_a       = pairs_a,
  pairs_b       = pairs_b_data,
  pairs_c       = pairs_c_data,
  pairs_d       = pairs_d_data,
  gsea_ckd      = gsea_ckd,
  gsea_med      = gsea_med,
  gsea_long     = gsea_long,
  rho_large     = rho_all,
  rho_mod       = rho_all2,
  n_deg_ckd     = n_deg_a,
  n_deg_isch_large = n_deg_b,
  n_deg_isch_mod   = n_deg_c,
  n_deg_isch_long_med = n_deg_d
)

saveRDS(results, file.path(res_dir, "ckd_vs_ischemia_logfc.rds"))
cat("Saved:", file.path(res_dir, "ckd_vs_ischemia_logfc.rds"), "\n")
cat("\nDone.\n")
