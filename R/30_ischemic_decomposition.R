#!/usr/bin/env Rscript
# =============================================================================
# 30: Ischemic Time Confound — Decomposition of CKD vs Ischemic Signals
#
# Loads discovery DE results and decomposes the CKD signal into ischemic
# and CKD-specific components.
# =============================================================================

library(tidyverse)
library(edgeR)
library(limma)
library(msigdbr)
library(fgsea)
library(ggrepel)

source("kidney/R/functions.R")

set.seed(42)
N_CORES <- 10

res_dir <- "kidney/results"
fig_dir <- "kidney/figures"

# Load discovery results
dat <- readRDS(file.path(res_dir, "prep_adipose_subcutaneous.rds"))
tb <- dat$tb
counts <- dat$counts
joint_factors <- dat$joint_factors

res_female <- readRDS(file.path(res_dir, "de_female.rds"))

n_obs_f <- res_female$n_deg

# Hallmark gene sets
hallmark_list <- load_hallmark_sets(filter_empty = FALSE)

# ── 1. Bin female controls by ischemic time ──────────────────────────────────
tb_f_ctrl <- tb %>% filter(SEX == "Female", MHRNLFLR == 0) %>%
  add_ischemic_bin()

cat("Female control ischemic time bins:\n")
print(table(tb_f_ctrl$isch_bin))
cat(sprintf("  Short median: %.1f hrs\n", median(tb_f_ctrl$ischemic_hrs[tb_f_ctrl$isch_bin == "Short"])))
cat(sprintf("  Long  median: %.1f hrs\n", median(tb_f_ctrl$ischemic_hrs[tb_f_ctrl$isch_bin == "Long"])))

# ── 2. DE: Long vs Short (drop Medium) ──────────────────────────────────────
cat("\n--- Long vs Short ischemic time DE (female controls only) ---\n")

tb_ls <- tb_f_ctrl %>% filter(isch_bin %in% c("Short", "Long"))

cat(sprintf("  Short: %d | Long: %d\n",
            sum(tb_ls$isch_bin == "Short"), sum(tb_ls$isch_bin == "Long")))

# Covariates: same as CKD model but replace MHRNLFLR with isch_bin, drop ischemic_hrs
# (ischemic_hrs is the treatment here, not a covariate)
isch_covariates <- setdiff(joint_factors, "MHRNLFLR")
form_isch <- as.formula(paste(
  "~ isch_bin + AGE + RACE + DTHHRDY + SMRIN + SMCENTER +",
  paste(isch_covariates, collapse = " + ")))

isch_fit <- run_limma_analysis(
  tb_ls, counts, form_isch, "isch_bin",
  coef_match = "prefix",
  context = "female control ischemic time model"
)
tt_isch <- isch_fit$tt

cat(sprintf("  Design: %d samples x %d columns\n", nrow(isch_fit$design), ncol(isch_fit$design)))
cat(sprintf("  Ischemic time coefficient: %s\n",
            get_required_coef_prefix(isch_fit$design, "isch_bin", "female control ischemic time model")))

n_deg_isch <- isch_fit$n_deg
n_up_isch <- isch_fit$n_up
n_down_isch <- isch_fit$n_down
pi0_isch <- tryCatch(pi0est(tt_isch$P.Value)$pi0, error = function(e) NA)

cat(sprintf("  Ischemic DEGs (FDR<0.05): %d (%d up, %d down)\n",
            n_deg_isch, n_up_isch, n_down_isch))
cat(sprintf("  Pi0: %.3f\n", pi0_isch))

# ── 3. Gene-level correlation: ischemic vs CKD t-statistics ─────────────────
cat("\n--- Gene-level correlation: ischemic vs CKD ---\n")

tt_ckd <- res_female$tt
de_cmp_isch <- compare_de_results(tt_ckd, tt_isch)

cat(sprintf("  Genes compared: %d\n", nrow(de_cmp_isch$merged)))
cat(sprintf("  t-stat rho (ischemic vs CKD): %.3f\n", de_cmp_isch$rho_t))
cat(sprintf("  logFC rho (ischemic vs CKD):  %.3f\n", de_cmp_isch$rho_lfc))
cat(sprintf("  Ischemic DEGs: %d | CKD DEGs: %d | Overlap: %d (%.1f%% of ischemic, %.1f%% of CKD)\n",
            de_cmp_isch$n_deg_new, de_cmp_isch$n_deg_ref, de_cmp_isch$overlap,
            100 * de_cmp_isch$overlap / max(de_cmp_isch$n_deg_new, 1),
            100 * de_cmp_isch$overlap / max(de_cmp_isch$n_deg_ref, 1)))

# Direction concordance among overlap
if (de_cmp_isch$overlap > 0) {
  cat(sprintf("  Direction concordance in overlap: %.1f%%\n",
              100 * de_cmp_isch$direction_concordance))
}

# ── 4. GSEA: ischemic time (Hallmark) ───────────────────────────────────────
cat("\n--- GSEA: ischemic time effect (Hallmark) ---\n")

# Reuse hallmark_list from severity section
gsea_isch <- run_gsea(tt_isch, hallmark_list, "Ischemic (Long vs Short)")
gsea_ckd_ref <- run_gsea(tt_ckd, hallmark_list, "CKD (Female)")

# ── 5. Pathway NES correlation ──────────────────────────────────────────────
cat("\n--- Pathway NES correlation ---\n")

core_paths <- CORE_PATHWAYS[1:7]
gsea_cmp_isch_res <- compare_gsea_results(
  gsea_ckd_ref, gsea_isch,
  label1 = "CKD", label2 = "isch",
  core_paths = core_paths
)
gsea_cmp_isch <- gsea_cmp_isch_res$merged
rho_nes_isch <- gsea_cmp_isch_res$rho_nes
cat(sprintf("  Pathway NES rho (CKD vs ischemic): %.3f\n", rho_nes_isch))

# ── 6. Figure: NES CKD vs Ischemic dot plot ─────────────────────────────────
cat("\nGenerating ischemic confound figure...\n")

gsea_plot <- gsea_cmp_isch %>%
  filter(padj_CKD < 0.05 | padj_isch < 0.05) %>%
  mutate(
    pathway = clean_pathway(pathway),
    sig_cat = case_when(
      padj_CKD < 0.05 & padj_isch < 0.05 ~ "Both",
      padj_CKD < 0.05 ~ "CKD only",
      TRUE ~ "Ischemic only"
    ),
    sig_cat = factor(sig_cat, levels = c("Both", "CKD only", "Ischemic only"))
  )

p_isch <- ggplot(gsea_plot, aes(x = NES_isch, y = NES_CKD, color = sig_cat)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_point(size = 1.8, alpha = 0.85) +
  ggrepel::geom_text_repel(
    aes(label = pathway), size = 1.8, max.overlaps = 50, seed = 123,
    segment.size = 0.15, segment.color = "grey60",
    box.padding = 0.35, point.padding = 0.15,
    min.segment.length = 0, force = 3, force_pull = 0.5,
    max.iter = 10000) +
  scale_color_manual(
    values = c("Both" = "#D6604D", "CKD only" = "#4393C3", "Ischemic only" = "#E69F00"),
    name = NULL) +
  labs(x = "NES (Ischemic time: Long vs Short)",
       y = "NES (CKD vs Control)") +
  annotate("text", x = -3.8, y = 3.4, label = sprintf("rho = %.2f", rho_nes_isch),
           size = 3, fontface = "italic", hjust = 0) +
  coord_fixed(xlim = c(-4.2, 3.5), ylim = c(-4.2, 3.5)) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = c(0.85, 0.10),
    legend.background = element_rect(fill = alpha("white", 0.9), linewidth = 0.3, color = "grey80"),
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 7),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.15, color = "grey93"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    plot.margin = margin(5, 5, 5, 5))

ggsave(file.path(fig_dir, "30_ischemic_confound.png"), p_isch,
       width = 5, height = 5, dpi = 300)
cat("Saved:", file.path(fig_dir, "30_ischemic_confound.png"), "\n")

# ── 7. Save results ─────────────────────────────────────────────────────────
ischemic_confound <- list(
  tt_isch = tt_isch,
  n_deg_isch = n_deg_isch,
  pi0_isch = pi0_isch,
  rho_t_isch_ckd = de_cmp_isch$rho_t,
  rho_lfc_isch_ckd = de_cmp_isch$rho_lfc,
  rho_nes_isch_ckd = rho_nes_isch,
  gsea_isch = gsea_isch,
  gsea_comparison = gsea_cmp_isch,
  n_short = sum(tb_ls$isch_bin == "Short"),
  n_long = sum(tb_ls$isch_bin == "Long"),
  overlap_degs = de_cmp_isch$overlap_genes
)
saveRDS(ischemic_confound, file.path(res_dir, "ischemic_confound.rds"))
write_csv(gsea_cmp_isch, file.path(res_dir, "ischemic_nes_comparison.csv"))

cat(sprintf("\n=== ISCHEMIC CONFOUND SUMMARY ===
  Short (<=3h): %d | Long (>10h): %d
  Ischemic DEGs: %d (pi0=%.3f)
  Gene logFC rho (isch vs CKD): %.3f
  Pathway NES rho (isch vs CKD): %.3f
  DEG overlap: %d / %d ischemic DEGs (%.1f%%)
",
  sum(tb_ls$isch_bin == "Short"), sum(tb_ls$isch_bin == "Long"),
  n_deg_isch, pi0_isch,
  de_cmp_isch$rho_lfc, rho_nes_isch,
  de_cmp_isch$overlap, n_deg_isch,
  100 * de_cmp_isch$overlap / max(n_deg_isch, 1)
))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1g: CKD DE WITH CATEGORICAL ISCHEMIC TIME ADJUSTMENT           ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 1g: CKD DE with Categorical Ischemic Time Adjustment\n")
cat(strrep("=", 70), "\n")

# ── 1. Apply 3-category ischemic binning to ALL females ──────────────────────
tb_f_all <- tb %>% filter(SEX == "Female") %>%
  add_ischemic_bin()

counts_f_all <- counts[, tb_f_all$SAMPID]

cat("Ischemic time bins (all females):\n")
cat("  Controls:\n")
print(table(tb_f_all$isch_bin[tb_f_all$MHRNLFLR == 0]))
cat("  CKD cases:\n")
print(table(tb_f_all$isch_bin[tb_f_all$MHRNLFLR == 1]))

# ── 2. Model: replace ischemic_hrs with isch_bin ────────────────────────────
cat("\n--- CKD DE with isch_bin (3-category) ---\n")

form_f_ischbin <- as.formula(paste(
  "~ AGE + RACE + BMI + DTHHRDY + isch_bin + SMRIN + SMCENTER + SMEXNCRT +",
  paste(joint_factors, collapse = " + ")))

design_f_ischbin <- build_design(form_f_ischbin, tb_f_all)

cat(sprintf("  Design: %d samples x %d columns\n",
            nrow(design_f_ischbin), ncol(design_f_ischbin)))
cat(sprintf("  Columns: %s\n", paste(colnames(design_f_ischbin), collapse = ", ")))

fc_ib <- get_required_coef_prefix(design_f_ischbin, "MHRNLFLR", "CKD isch_bin model")
cat(sprintf("  CKD coefficient: %s\n", fc_ib))

de_ib <- run_voom_de(counts_f_all, design_f_ischbin, fc_ib)
tt_f_ib <- de_ib$tt

deg_ib <- summarize_degs(tt_f_ib)

cat(sprintf("  CKD DEGs (isch_bin adjusted): %d (%d up, %d down)\n",
            deg_ib$n_deg, deg_ib$n_up, deg_ib$n_down))
cat(sprintf("  Pi0: %.3f\n", deg_ib$pi0))
cat(sprintf("  Compare: linear ischemic = %d DEGs → categorical = %d DEGs (%.1f%% retained)\n",
            n_obs_f, deg_ib$n_deg, 100 * deg_ib$n_deg / n_obs_f))

# ── 3. t-stat correlation: linear vs categorical adjustment ──────────────────
de_cmp_linear_vs_cat <- compare_de_results(res_female$tt, tt_f_ib)
rho_linear_vs_cat <- de_cmp_linear_vs_cat$rho_t

cat(sprintf("  t-stat rho (linear vs categorical ischemic): %.3f\n", rho_linear_vs_cat))

# ── 4. GSEA comparison ──────────────────────────────────────────────────────
cat("\n--- GSEA: categorical ischemic adjustment (Hallmark) ---\n")

gsea_ib <- run_gsea(tt_f_ib, hallmark_list, "CKD (isch_bin adjusted)")

gsea_cmp_ib_res <- compare_gsea_results(
  gsea_ckd_ref, gsea_ib,
  label1 = "linear", label2 = "cat",
  core_paths = CORE_PATHWAYS[1:7]
)
gsea_cmp_ib <- gsea_cmp_ib_res$merged
rho_nes_ib <- gsea_cmp_ib_res$rho_nes
cat(sprintf("  Pathway NES rho (linear vs categorical): %.3f\n", rho_nes_ib))

# ── 5. Also check isch_bin coefficients themselves ───────────────────────────
cat("\n--- Ischemic bin effects in the CKD model ---\n")

isch_coefs <- grep("isch_bin", colnames(design_f_ischbin), value = TRUE)
for (ic in isch_coefs) {
  tt_ic <- topTable(de_ib$fit, coef = ic, number = Inf, sort.by = "none")
  n_ic <- sum(tt_ic$adj.P.Val < 0.05)
  cat(sprintf("  %s: %d DEGs\n", ic, n_ic))
}

# ── 6. Save results ─────────────────────────────────────────────────────────
ischemic_adjustment <- list(
  tt_ischbin = tt_f_ib,
  n_deg_ischbin = deg_ib$n_deg,
  pi0_ischbin = deg_ib$pi0,
  n_deg_linear = n_obs_f,
  rho_linear_vs_cat = rho_linear_vs_cat,
  rho_nes_linear_vs_cat = rho_nes_ib,
  gsea_ischbin = gsea_ib,
  gsea_comparison = gsea_cmp_ib
)
saveRDS(ischemic_adjustment, file.path(res_dir, "ischemic_adjustment.rds"))

cat(sprintf("\n=== ISCHEMIC ADJUSTMENT SUMMARY ===
  Linear ischemic: %d DEGs
  Categorical ischemic (3-bin): %d DEGs (%.1f%% retained)
  t-stat rho (linear vs cat): %.3f
  Pathway NES rho: %.3f
\n",
  n_obs_f, deg_ib$n_deg, 100 * deg_ib$n_deg / n_obs_f,
  rho_linear_vs_cat, rho_nes_ib
))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1h: ISCHEMIC vs CKD GENE-LEVEL DECOMPOSITION                   ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 1h: Ischemic vs CKD Gene-Level Decomposition\n")
cat(strrep("=", 70), "\n")

library(ggrepel)

# ── 1. Merge CKD and ischemic DE results ────────────────────────────────────
tt_ckd_h <- res_female$tt
tt_isch_h <- ischemic_confound$tt_isch

df_h <- inner_join(
  tt_ckd_h %>% dplyr::select(gene_id, logFC_CKD = logFC, t_CKD = t,
                              padj_CKD = adj.P.Val),
  tt_isch_h %>% dplyr::select(gene_id, logFC_isch = logFC, t_isch = t,
                               padj_isch = adj.P.Val),
  by = "gene_id"
)
cat(sprintf("Genes compared: %d\n", nrow(df_h)))

# ── 2. Assign core pathway membership ───────────────────────────────────────
core_path_names <- CORE_PATHWAYS[1:6]
path_labels <- PATHWAY_LABELS[core_path_names]

h_genes <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  filter(gs_name %in% core_path_names)

df_h$pathway <- "Other"
for (p in core_path_names) {
  genes_in <- h_genes$ensembl_gene[h_genes$gs_name == p]
  df_h$pathway[df_h$gene_id %in% genes_in & df_h$pathway == "Other"] <- path_labels[p]
}

# ── 3. Per-pathway logFC correlations ───────────────────────────────────────
cat("\nPer-pathway logFC correlations (ischemic vs CKD):\n")
pw_rhos <- list()
for (pw in c(path_labels, "Other")) {
  sub <- df_h %>% filter(pathway == pw)
  if (nrow(sub) > 5) {
    r <- cor(sub$logFC_CKD, sub$logFC_isch, method = "spearman")
    pw_rhos[[pw]] <- r
    cat(sprintf("  %-15s n=%3d  rho=%.3f\n", pw, nrow(sub), r))
  }
}

# ── 4. Quadrant analysis (CKD DEGs) ────────────────────────────────────────
df_sig <- df_h %>% filter(padj_CKD < 0.05)
q_concordant <- sum(sign(df_sig$logFC_CKD) == sign(df_sig$logFC_isch))
cat(sprintf("\nCKD DEGs direction concordance with ischemic: %d/%d (%.1f%%)\n",
            q_concordant, nrow(df_sig), 100 * q_concordant / nrow(df_sig)))

# ── 5. Gene-level scatter with pathway coloring ────────────────────────────
cat("\nGenerating ischemic vs CKD gene scatter...\n")

gene_anno <- readRDS("data/annotation/gencode_v47_genes.rds")
gene_map <- setNames(gene_anno$symbol, rownames(gene_anno))
df_h$gene_symbol <- gene_map[df_h$gene_id]

pathway_colors <- c(
  "Adipogenesis" = "#2166AC", "Oxphos" = "#4393C3",
  "FA Metabolism" = "#92C5DE", "TNFa/NF-kB" = "#D6604D",
  "EMT" = "#F4A582", "Inflammatory" = "#B2182B",
  "Other" = "grey85"
)
df_h$pathway <- factor(df_h$pathway,
                       levels = c("Other", names(pathway_colors)[1:6]))

rho_lfc_h <- cor(df_h$logFC_CKD, df_h$logFC_isch, method = "spearman")

p_gene <- ggplot(df_h, aes(x = logFC_isch, y = logFC_CKD)) +
  geom_point(data = df_h %>% filter(pathway == "Other"),
             size = 0.3, alpha = 0.2, color = "grey85") +
  geom_point(data = df_h %>% filter(pathway != "Other"),
             aes(color = pathway), size = 1.2, alpha = 0.7) +
  scale_color_manual(values = pathway_colors, name = "Hallmark Pathway") +
  geom_smooth(method = "lm", se = FALSE, color = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_text_repel(
    data = df_h %>% filter(pathway != "Other",
                            abs(logFC_CKD) > 1.5 | abs(logFC_isch) > 1.5),
    aes(label = gene_symbol),
    size = 2.5, max.overlaps = 25, seed = 42,
    segment.size = 0.3, segment.alpha = 0.5,
    min.segment.length = 0.2) +
  labs(x = "logFC (Ischemic: Long vs Short, controls only)",
       y = "logFC (CKD vs Control, females)",
       title = "Gene-Level Effect Size: CKD vs Ischemic Time",
       subtitle = sprintf("logFC Spearman rho = %.2f | %d genes", rho_lfc_h, nrow(df_h))) +
  coord_fixed(xlim = c(-3, 3), ylim = c(-3, 3)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(fig_dir, "30_ischemic_gene_scatter.png"), p_gene,
       width = 9, height = 8, dpi = 300)
cat("Saved:", file.path(fig_dir, "30_ischemic_gene_scatter.png"), "\n")

# ── 5b. Within-model CKD vs ischemic_hrs logFC scatter ──────────────────────
cat("\nGenerating within-model CKD vs ischemic_hrs scatter...\n")

tt_ckd_main <- topTable(res_female$fit, coef = "MHRNLFLR", number = Inf, sort.by = "none")
tt_ckd_main$gene_id <- rownames(tt_ckd_main)

tt_isch_main <- topTable(res_female$fit, coef = "ischemic_hrs", number = Inf, sort.by = "none")
tt_isch_main$gene_id <- rownames(tt_isch_main)

df_main <- inner_join(
  tt_ckd_main %>% dplyr::select(gene_id, logFC_CKD = logFC, padj_CKD = adj.P.Val),
  tt_isch_main %>% dplyr::select(gene_id, logFC_isch = logFC, padj_isch = adj.P.Val),
  by = "gene_id"
)

# Pathway annotation
df_main$pathway <- "Other"
for (p in core_path_names) {
  genes_in <- h_genes$ensembl_gene[h_genes$gs_name == p]
  df_main$pathway[df_main$gene_id %in% genes_in & df_main$pathway == "Other"] <- path_labels[p]
}
df_main$pathway <- factor(df_main$pathway,
                           levels = c("Other", names(pathway_colors)[1:6]))

# Gene symbols
df_main$gene_symbol <- gene_map[df_main$gene_id]

rho_lfc_main <- cor(df_main$logFC_CKD, df_main$logFC_isch, method = "spearman")
cat(sprintf("  Within-model logFC rho (CKD vs ischemic_hrs): %.3f (%d genes)\n",
            rho_lfc_main, nrow(df_main)))

p_main <- ggplot(df_main, aes(x = logFC_isch, y = logFC_CKD)) +
  geom_point(data = df_main %>% filter(pathway == "Other"),
             size = 0.3, alpha = 0.2, color = "grey85") +
  geom_point(data = df_main %>% filter(pathway != "Other"),
             aes(color = pathway), size = 1.2, alpha = 0.7) +
  scale_color_manual(values = pathway_colors, name = "Hallmark Pathway") +
  geom_smooth(method = "lm", se = FALSE, color = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_text_repel(
    data = df_main %>% filter(pathway != "Other",
                               abs(logFC_CKD) > 1.5 | abs(logFC_isch) > 1.5),
    aes(label = gene_symbol),
    size = 2.5, max.overlaps = 25, seed = 42,
    segment.size = 0.3, segment.alpha = 0.5,
    min.segment.length = 0.2) +
  labs(x = "logFC (ischemic_hrs, same model)",
       y = "logFC (CKD vs Control, females)",
       title = "Within-Model: CKD vs Ischemic Time Coefficients",
       subtitle = sprintf("logFC Spearman rho = %.2f | %d genes", rho_lfc_main, nrow(df_main))) +
  coord_cartesian(xlim = range(df_main$logFC_isch, na.rm = TRUE),
                  ylim = c(-3, 3)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(fig_dir, "30_ischemic_within_model_scatter.png"), p_main,
       width = 9, height = 8, dpi = 300)
cat("Saved:", file.path(fig_dir, "30_ischemic_within_model_scatter.png"), "\n")

# ── 5c. Coefficient correlation & PCA ────────────────────────────────────────
cat("\nGenerating coefficient correlation and PCA plots...\n")

# Extract all coefficients (genes x covariates matrix)
coef_mat <- res_female$fit$coefficients
# Drop intercept
coef_mat <- coef_mat[, colnames(coef_mat) != "(Intercept)", drop = FALSE]
coef_mat <- coef_mat[complete.cases(coef_mat), ]
cat(sprintf("  Coefficient matrix: %d genes x %d covariates\n", nrow(coef_mat), ncol(coef_mat)))

# Correlation plot
cor_coef <- cor(coef_mat, method = "spearman")

library(corrplot)
png(file.path(fig_dir, "30_coefficient_correlation.png"), width = 8, height = 8,
    units = "in", res = 300)
corrplot(cor_coef, method = "color", type = "lower", order = "hclust",
         tl.col = "black", tl.cex = 0.8, addCoef.col = "black", number.cex = 0.6,
         title = "Spearman Correlation of Model Coefficients Across Genes",
         mar = c(0, 0, 2, 0))
dev.off()
cat("Saved:", file.path(fig_dir, "30_coefficient_correlation.png"), "\n")

# PCA of coefficient profiles
pca_coef <- prcomp(coef_mat, center = TRUE, scale. = TRUE)
var_explained <- summary(pca_coef)$importance[2, ] * 100

# Project covariates (loadings) into PC space
loadings <- as.data.frame(pca_coef$rotation[, 1:2])
loadings$covariate <- rownames(loadings)

p_pca <- ggplot(loadings, aes(x = PC1, y = PC2, label = covariate)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text_repel(size = 3.5, max.overlaps = 20, seed = 42) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  labs(x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
       y = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
       title = "PCA of Model Coefficients (gene-level loadings)") +
  theme_bw(base_size = 11)

ggsave(file.path(fig_dir, "30_coefficient_pca.png"), p_pca,
       width = 8, height = 7, dpi = 300)
cat("Saved:", file.path(fig_dir, "30_coefficient_pca.png"), "\n")

# ── 6. Residual GSEA: CKD after removing ischemic component ────────────────
cat("\n--- Residual GSEA: CKD after removing ischemic component ---\n")

fit_resid <- lm(t_CKD ~ t_isch, data = df_h)
df_h$t_resid <- residuals(fit_resid)

cat(sprintf("  Regression: t_CKD = %.3f + %.3f * t_isch\n",
            coef(fit_resid)[1], coef(fit_resid)[2]))
cat(sprintf("  R² = %.3f (ischemic explains %.1f%% of CKD t-stat variance)\n",
            summary(fit_resid)$r.squared, 100 * summary(fit_resid)$r.squared))

ranks_ckd_h <- setNames(df_h$t_CKD, df_h$gene_id)
ranks_isch_h <- setNames(df_h$t_isch, df_h$gene_id)
ranks_resid_h <- setNames(df_h$t_resid, df_h$gene_id)

gsea_ckd_h <- run_gsea(data.frame(gene_id = names(ranks_ckd_h),
                                   t = unname(ranks_ckd_h)),
                        hallmark_list, "CKD (original)")
gsea_isch_h <- run_gsea(data.frame(gene_id = names(ranks_isch_h),
                                    t = unname(ranks_isch_h)),
                         hallmark_list, "Ischemic")
gsea_resid_h <- run_gsea(data.frame(gene_id = names(ranks_resid_h),
                                     t = unname(ranks_resid_h)),
                          hallmark_list, "CKD residual")

# Three-way comparison
cmp_h <- gsea_ckd_h %>%
  dplyr::select(pathway, NES_CKD = NES, padj_CKD = padj) %>%
  inner_join(gsea_isch_h %>% dplyr::select(pathway, NES_isch = NES, padj_isch = padj),
             by = "pathway") %>%
  inner_join(gsea_resid_h %>% dplyr::select(pathway, NES_resid = NES, padj_resid = padj),
             by = "pathway") %>%
  mutate(
    pct_retained = ifelse(abs(NES_CKD) > 0.01, 100 * NES_resid / NES_CKD, NA),
    category = case_when(
      padj_CKD < 0.05 & padj_resid < 0.05 & padj_isch >= 0.05 ~ "CKD-specific",
      padj_CKD < 0.05 & padj_resid >= 0.05 & padj_isch < 0.05 ~ "Ischemic-driven",
      padj_CKD < 0.05 & padj_resid < 0.05 & padj_isch < 0.05 ~ "Shared (CKD persists)",
      padj_CKD < 0.05 & padj_resid >= 0.05 ~ "Lost after residualization",
      TRUE ~ "Not significant"
    )
  )

cat("\nPathway categories after residualization:\n")
print(table(cmp_h$category))

# Core pathways table
core_resid <- cmp_h %>%
  filter(padj_CKD < 0.05) %>%
  mutate(pathway = gsub("HALLMARK_", "", pathway)) %>%
  arrange(pct_retained)

cat(sprintf("\n%-40s %8s %8s %8s %6s  %s\n",
            "Pathway", "NES_CKD", "NES_isch", "NES_res", "%ret", "Category"))
cat(strrep("-", 105), "\n")
for (i in seq_len(nrow(core_resid))) {
  cat(sprintf("%-40s %+8.2f %+8.2f %+8.2f %5.0f%%  %s\n",
              core_resid$pathway[i],
              core_resid$NES_CKD[i], core_resid$NES_isch[i],
              core_resid$NES_resid[i], core_resid$pct_retained[i],
              core_resid$category[i]))
}

# NES correlations
rho_ckd_isch_nes <- cor(cmp_h$NES_CKD, cmp_h$NES_isch, method = "spearman")
rho_ckd_resid_nes <- cor(cmp_h$NES_CKD, cmp_h$NES_resid, method = "spearman")
rho_isch_resid_nes <- cor(cmp_h$NES_isch, cmp_h$NES_resid, method = "spearman")
cat(sprintf("\nNES correlations:
  CKD vs Ischemic:  %.3f
  CKD vs Residual:  %.3f
  Isch vs Residual: %.3f\n",
  rho_ckd_isch_nes, rho_ckd_resid_nes, rho_isch_resid_nes))

# ── 7. Figure: Three-signal dot plot ───────────────────────────────────────
cat("\nGenerating residual GSEA figure...\n")

plot_resid <- cmp_h %>%
  filter(padj_CKD < 0.05 | padj_resid < 0.05) %>%
  mutate(pathway = gsub("HALLMARK_", "", pathway),
         pathway = gsub("_", " ", pathway)) %>%
  dplyr::select(pathway, NES_CKD, NES_isch, NES_resid, category) %>%
  pivot_longer(cols = starts_with("NES_"),
               names_to = "model", values_to = "NES",
               names_prefix = "NES_") %>%
  mutate(model = factor(model,
    levels = c("CKD", "isch", "resid"),
    labels = c("CKD (original)", "Ischemic (Long vs Short)", "CKD (residual)")))

path_order_r <- plot_resid %>%
  filter(model == "CKD (original)") %>%
  arrange(NES) %>%
  pull(pathway)
plot_resid$pathway <- factor(plot_resid$pathway, levels = path_order_r)

p_resid <- ggplot(plot_resid, aes(x = NES, y = pathway, color = model)) +
  geom_point(size = 2.5, alpha = 0.8,
             position = position_dodge(width = 0.6)) +
  scale_color_manual(
    values = c("CKD (original)" = "firebrick",
               "Ischemic (Long vs Short)" = "orange",
               "CKD (residual)" = "steelblue"),
    name = "Signal") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Normalized Enrichment Score",
       y = NULL,
       title = "Pathway Decomposition: CKD vs Ischemic Time",
       subtitle = sprintf("Residual = CKD after removing ischemic component (R² = %.2f)",
                          summary(fit_resid)$r.squared)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "30_residual_gsea.png"), p_resid,
       width = 10, height = 9, dpi = 300)
cat("Saved:", file.path(fig_dir, "30_residual_gsea.png"), "\n")

# ── 8. Save all ischemic decomposition results ─────────────────────────────
ischemic_decomp <- list(
  gene_scatter_df = df_h,
  per_pathway_rho = pw_rhos,
  concordance_pct = 100 * q_concordant / nrow(df_sig),
  lm_fit = fit_resid,
  r_squared = summary(fit_resid)$r.squared,
  gsea_ckd = gsea_ckd_h,
  gsea_isch = gsea_isch_h,
  gsea_resid = gsea_resid_h,
  pathway_comparison = cmp_h,
  rho_lfc = rho_lfc_h
)
saveRDS(ischemic_decomp, file.path(res_dir, "ischemic_decomposition.rds"))
write_csv(cmp_h, file.path(res_dir, "ischemic_residual_nes_comparison.csv"))
cat("Saved ischemic decomposition results.\n")

cat("\nDone.\n")
