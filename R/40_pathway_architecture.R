#!/usr/bin/env Rscript
# =============================================================================
# Stage 2: Pathway Architecture — Biological Characterization & Stability
#
# Sections:
#   2a — GSEA on SEX×CKD interaction + ORA (female up/down)
#   2b — GSEA on male CKD and PSM-matched female CKD
#   2c — Consolidated pathway stability table
#   2d — Female CKD PCA
#   2e — Figures
# =============================================================================

library(tidyverse)
library(edgeR)
library(limma)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(pheatmap)
library(RColorBrewer)

source("kidney/R/functions.R")

set.seed(42)

res_dir <- "kidney/results"
fig_dir <- "kidney/figures"


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  LOAD DATA & HALLMARK GENE SETS                                         ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("Stage 2: Pathway Architecture\n")
cat(strrep("=", 70), "\n\n")

# Hallmark gene sets (Ensembl IDs, consistent with stage 1)
hallmark_list <- load_hallmark_sets()
cat("Hallmark gene sets:", length(hallmark_list), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2a: GSEA on SEX×CKD INTERACTION T-STATISTICS                   ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 2a: Interaction GSEA\n")
cat(strrep("=", 70), "\n")

ix <- readRDS(file.path(res_dir, "interaction_results.rds"))
gsea_interaction <- run_gsea(ix$tt_interaction, hallmark_list, seed = 42)

n_sig <- sum(gsea_interaction$padj < 0.05)
cat(sprintf("Interaction GSEA: %d pathways FDR < 0.05\n", n_sig))

# Compare interaction NES to female-only NES
de_female <- readRDS(file.path(res_dir, "de_female.rds"))
gsea_female_full <- run_gsea(de_female$tt, hallmark_list, seed = 42)

compare_ix_res <- compare_gsea_results(
  gsea_female_full, gsea_interaction,
  label1 = "female", label2 = "interaction",
  print_table = FALSE
)
compare_ix <- compare_ix_res$merged
rho_ix <- compare_ix_res$rho_nes
cat(sprintf("Female vs Interaction NES rho = %.3f\n", rho_ix))


# ── ORA: Female up/down DEGs over Hallmark sets ─────────────────────────────
cat("\n--- ORA: Female CKD DEGs (Hallmark, clusterProfiler) ---\n")

hallmark_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, ensembl_gene)

tt_female <- de_female$tt
universe <- tt_female$gene_id

up_genes   <- tt_female %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% pull(gene_id)
down_genes <- tt_female %>% filter(adj.P.Val < 0.05 & logFC < 0) %>% pull(gene_id)

cat(sprintf("  Up genes: %d | Down genes: %d | Universe: %d\n",
            length(up_genes), length(down_genes), length(universe)))

ora_up <- enricher(up_genes, TERM2GENE = hallmark_t2g, universe = universe,
                   pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
ora_down <- enricher(down_genes, TERM2GENE = hallmark_t2g, universe = universe,
                     pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)

ora_up_df <- as.data.frame(ora_up) %>% mutate(direction = "Up")
ora_down_df <- as.data.frame(ora_down) %>% mutate(direction = "Down")

cat(sprintf("  Up: %d enriched (FDR<0.05) | Down: %d enriched (FDR<0.05)\n",
            sum(ora_up_df$p.adjust < 0.05), sum(ora_down_df$p.adjust < 0.05)))

if (sum(ora_up_df$p.adjust < 0.05) > 0) {
  cat("\n  Top Up-regulated pathways:\n")
  top_up <- ora_up_df %>% filter(p.adjust < 0.05) %>% arrange(p.adjust) %>% head(10)
  for (i in seq_len(nrow(top_up))) {
    cat(sprintf("    %-50s  q=%.1e  %s\n",
                gsub("HALLMARK_", "", top_up$ID[i]),
                top_up$p.adjust[i], top_up$GeneRatio[i]))
  }
}

if (sum(ora_down_df$p.adjust < 0.05) > 0) {
  cat("\n  Top Down-regulated pathways:\n")
  top_dn <- ora_down_df %>% filter(p.adjust < 0.05) %>% arrange(p.adjust) %>% head(10)
  for (i in seq_len(nrow(top_dn))) {
    cat(sprintf("    %-50s  q=%.1e  %s\n",
                gsub("HALLMARK_", "", top_dn$ID[i]),
                top_dn$p.adjust[i], top_dn$GeneRatio[i]))
  }
}

ora_combined <- bind_rows(ora_up_df, ora_down_df)
write_csv(ora_combined, file.path(res_dir, "ora_female_hallmark.csv"))
saveRDS(list(up = ora_up, down = ora_down), file.path(res_dir, "ora_female_hallmark.rds"))
cat("Saved:", file.path(res_dir, "ora_female_hallmark.csv"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2b: GSEA on MALE CKD and PSM-MATCHED FEMALE CKD               ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 2b: Male CKD and PSM-matched female GSEA\n")
cat(strrep("=", 70), "\n")

# Male CKD GSEA
de_male <- readRDS(file.path(res_dir, "de_male.rds"))
gsea_male <- run_gsea(de_male$tt, hallmark_list, seed = 42)
n_sig_male <- sum(gsea_male$padj < 0.05)
cat(sprintf("Male CKD GSEA: %d pathways FDR < 0.05\n", n_sig_male))

# PSM-matched female CKD GSEA
psm_female <- readRDS(file.path(res_dir, "de_female_psm.rds"))
gsea_psm <- run_gsea(psm_female$res$tt, hallmark_list, seed = 42)
n_sig_psm <- sum(gsea_psm$padj < 0.05)
cat(sprintf("PSM-matched female GSEA: %d pathways FDR < 0.05\n", n_sig_psm))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2c: CONSOLIDATED PATHWAY STABILITY TABLE                       ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 2c: Pathway Stability Table\n")
cat(strrep("=", 70), "\n")

# Load existing GSEA results
# Note: ischemic_decomposition.rds is the canonical file produced by stage1_01.R
isch_decomp <- readRDS(file.path(res_dir, "ischemic_decomposition.rds"))
isch_adj <- readRDS(file.path(res_dir, "ischemic_adjustment.rds"))
severity <- readRDS(file.path(res_dir, "severity_results.rds"))
# Get categorical ischemic adjustment GSEA (handle old/new field names)
gsea_isch_adj <- if (!is.null(isch_adj$gsea_ischbin)) {
  isch_adj$gsea_ischbin
} else if (!is.null(isch_adj$gsea_linear)) {
  warning("Using gsea_linear from ischemic_adjustment.rds (old field name); ",
          "re-run stage1_01.R to regenerate with gsea_ischbin")
  isch_adj$gsea_linear
} else {
  stop("No ischemic adjustment GSEA found in ischemic_adjustment.rds")
}

# Core pathways for the main table
core_pathways <- CORE_PATHWAYS

# Extract NES and padj for each analysis into a unified format
extract_nes <- function(gsea_df, label) {
  if (is.data.frame(gsea_df)) {
    gsea_df %>%
      as_tibble() %>%
      dplyr::select(pathway, NES, padj) %>%
      mutate(analysis = label)
  } else {
    tibble(pathway = character(), NES = numeric(), padj = numeric(), analysis = character())
  }
}

# Assemble all GSEA results
all_gsea <- bind_rows(
  # Main-text columns
  extract_nes(gsea_female_full, "Full female"),
  extract_nes(gsea_psm, "PSM-matched female"),
  extract_nes(isch_decomp$gsea_resid, "Residual (ischemic removed)"),
  extract_nes(gsea_interaction, "SEX × CKD interaction"),
  # Supplementary columns
  extract_nes(gsea_isch_adj, "Ischemic-time adjusted"),
  extract_nes(severity$gsea_mild, "Mild-only CKD"),
  extract_nes(gsea_male, "Male CKD")
)

# Column order
col_order <- c(
  "Full female",
  "PSM-matched female",
  "Residual (ischemic removed)",
  "SEX × CKD interaction",
  "Ischemic-time adjusted",
  "Mild-only CKD",
  "Male CKD"
)

# Build NES matrix
nes_matrix <- all_gsea %>%
  filter(pathway %in% core_pathways) %>%
  dplyr::select(pathway, NES, analysis) %>%
  pivot_wider(names_from = analysis, values_from = NES) %>%
  mutate(pathway_clean = clean_pathway(pathway)) %>%
  column_to_rownames("pathway_clean") %>%
  dplyr::select(all_of(col_order))

# Build padj matrix
padj_matrix <- all_gsea %>%
  filter(pathway %in% core_pathways) %>%
  dplyr::select(pathway, padj, analysis) %>%
  pivot_wider(names_from = analysis, values_from = padj) %>%
  mutate(pathway_clean = clean_pathway(pathway)) %>%
  column_to_rownames("pathway_clean") %>%
  dplyr::select(all_of(col_order))

cat("Pathway stability table dimensions:", nrow(nes_matrix), "x", ncol(nes_matrix), "\n")

# Direction consistency check (main-text columns only)
main_cols <- col_order[1:4]
direction_check <- nes_matrix[, main_cols]
consistent <- apply(direction_check, 1, function(x) {
  x <- x[!is.na(x)]
  all(x > 0) || all(x < 0)
})
cat(sprintf("Direction consistent across 4 main columns: %d/%d pathways\n",
            sum(consistent), length(consistent)))
cat("Inconsistent:", paste(names(consistent[!consistent]), collapse = ", "), "\n")

# Build the full stability table (wide CSV)
stability_wide <- all_gsea %>%
  filter(pathway %in% core_pathways) %>%
  mutate(pathway_clean = clean_pathway(pathway)) %>%
  pivot_wider(
    id_cols = c(pathway, pathway_clean),
    names_from = analysis,
    values_from = c(NES, padj),
    names_glue = "{analysis}_{.value}"
  )

# Tidy version for supplementary
stability_long <- all_gsea %>%
  filter(pathway %in% core_pathways) %>%
  mutate(
    pathway_clean = clean_pathway(pathway),
    sig = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""),
    direction = ifelse(NES > 0, "Up", "Down"),
    main_text = analysis %in% main_cols
  )

# Save CSV with NES and significance
csv_table <- stability_long %>%
  mutate(NES_sig = sprintf("%.2f%s", NES, sig)) %>%
  dplyr::select(pathway_clean, analysis, NES_sig) %>%
  pivot_wider(names_from = analysis, values_from = NES_sig) %>%
  dplyr::select(pathway_clean, all_of(col_order))

write.csv(csv_table, file.path(res_dir, "pathway_stability_table.csv"),
          row.names = FALSE)
cat("Saved:", file.path(res_dir, "pathway_stability_table.csv"), "\n")

# Print the table
cat("\n--- Pathway Stability Table (NES with significance) ---\n")
print(as.data.frame(csv_table), right = FALSE)

# Also build full-pathway table (all 50 Hallmark, for supplement)
full_table <- all_gsea %>%
  mutate(pathway_clean = clean_pathway(pathway)) %>%
  mutate(NES_sig = sprintf("%.2f%s",
                           NES,
                           case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**",
                                     padj < 0.05 ~ "*", TRUE ~ ""))) %>%
  dplyr::select(pathway_clean, analysis, NES_sig) %>%
  pivot_wider(names_from = analysis, values_from = NES_sig) %>%
  dplyr::select(pathway_clean, any_of(col_order))

write.csv(full_table, file.path(res_dir, "pathway_stability_full.csv"),
          row.names = FALSE)
cat("Saved:", file.path(res_dir, "pathway_stability_full.csv"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2d: FEMALE CKD PCA                                             ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 2d: Female CKD PCA\n")
cat(strrep("=", 70), "\n")

# Reload data for PCA
dat <- readRDS(file.path(res_dir, "prep_adipose_subcutaneous.rds"))
tb <- dat$tb
counts <- dat$counts

# Female subset
tb_f <- tb %>% filter(SEX == "Female")
counts_f <- counts[, tb_f$SAMPID]

# Normalize with voom for PCA
dge_f <- DGEList(counts = counts_f)
dge_f <- calcNormFactors(dge_f, method = "TMM")
logcpm <- cpm(dge_f, log = TRUE, prior.count = 1)

# PCA on top 1000 most variable genes
vars <- apply(logcpm, 1, var)
top_genes <- names(sort(vars, decreasing = TRUE))[1:1000]
pca <- prcomp(t(logcpm[top_genes, ]), scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:5])
pca_df$SAMPID <- rownames(pca_df)
pca_df <- left_join(pca_df, tb_f, by = "SAMPID")

var_explained <- summary(pca)$importance[2, 1:5] * 100

cat(sprintf("PC1: %.1f%%, PC2: %.1f%%, PC3: %.1f%%\n",
            var_explained[1], var_explained[2], var_explained[3]))

# Test CKD association with each PC (CKD - Control direction)
pca_df$CKD <- factor(pca_df$MHRNLFLR, levels = c(0, 1), labels = c("Control", "CKD"))
for (i in 1:5) {
  pc_col <- paste0("PC", i)
  tt <- t.test(pca_df[[pc_col]] ~ pca_df$CKD)
  # t.test reports mean(Control) then mean(CKD); CKD effect = mean(CKD) - mean(Control)
  ckd_effect <- tt$estimate[2] - tt$estimate[1]
  cat(sprintf("  %s ~ CKD: p = %.3g, CKD effect = %.2f\n", pc_col, tt$p.value, ckd_effect))
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2e: FIGURES                                                     ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 2e: Figures\n")
cat(strrep("=", 70), "\n")

# --- Figure 1: Female pathway dot plot ---
cat("Figure 1: Female pathway dot plot\n")

gsea_plot <- gsea_female_full %>%
  as_tibble() %>%
  mutate(pathway_clean = clean_pathway(pathway)) %>%
  arrange(NES) %>%
  mutate(
    sig = padj < 0.05,
    direction = ifelse(NES > 0, "Upregulated", "Downregulated")
  ) %>%
  # Top 10 up + top 10 down by |NES|
  {
    up <- filter(., NES > 0) %>% arrange(desc(NES)) %>% head(10)
    down <- filter(., NES < 0) %>% arrange(NES) %>% head(10)
    bind_rows(down, up)
  } %>%
  mutate(pathway_clean = factor(pathway_clean, levels = pathway_clean))

p_dot <- ggplot(gsea_plot, aes(x = NES, y = pathway_clean)) +
  geom_point(aes(size = size, color = padj < 0.05), shape = 16) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("TRUE" = "firebrick3", "FALSE" = "grey60"),
                     labels = c("TRUE" = "FDR < 0.05", "FALSE" = "NS"),
                     name = "") +
  scale_size_continuous(range = c(2, 6), name = "Gene set size") +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "Hallmark Pathway Enrichment — Female CKD"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

ggsave(file.path(fig_dir, "40_pathway_dotplot.png"), p_dot,
       width = 8, height = 6, dpi = 300)
cat("  Saved: 40_pathway_dotplot.png\n")


# --- Figure 2: Pathway stability heatmap ---
cat("Figure 2: Pathway stability heatmap\n")

# NES matrix for main-text columns
heatmap_mat <- as.matrix(nes_matrix[, main_cols])

# Significance annotation
sig_mat <- as.matrix(padj_matrix[, main_cols])
sig_labels <- matrix("", nrow = nrow(sig_mat), ncol = ncol(sig_mat))
sig_labels[sig_mat < 0.05] <- "*"
sig_labels[sig_mat < 0.01] <- "**"
sig_labels[sig_mat < 0.001] <- "***"

# Color palette
max_abs <- max(abs(heatmap_mat), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)
colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# Column annotation: main vs supplementary
col_annotation <- data.frame(
  Type = c("Primary", "Confounder-controlled", "Ischemia-removed", "Sex-specific"),
  row.names = main_cols
)
ann_colors <- list(Type = c(
  "Primary" = "#2166AC",
  "Confounder-controlled" = "#4DAF4A",
  "Ischemia-removed" = "#E7298A",
  "Sex-specific" = "#FF7F00"
))

png(file.path(fig_dir, "40_stability_heatmap.png"),
    width = 8, height = 5.5, units = "in", res = 300)
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
  main = "Pathway NES Stability Across Analyses",
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 8,
  angle_col = 45,
  cellwidth = 40,
  cellheight = 20
)
dev.off()
cat("  Saved: 40_stability_heatmap.png\n")


# --- Figure 2b: Full 8-column stability heatmap (supplement) ---
cat("Figure 2b: Full stability heatmap (supplement)\n")

heatmap_mat_full <- as.matrix(nes_matrix)
sig_mat_full <- as.matrix(padj_matrix)
sig_labels_full <- matrix("", nrow = nrow(sig_mat_full), ncol = ncol(sig_mat_full))
sig_labels_full[sig_mat_full < 0.05] <- "*"
sig_labels_full[sig_mat_full < 0.01] <- "**"
sig_labels_full[sig_mat_full < 0.001] <- "***"

max_abs_full <- max(abs(heatmap_mat_full), na.rm = TRUE)
breaks_full <- seq(-max_abs_full, max_abs_full, length.out = 101)

col_ann_full <- data.frame(
  Section = c(rep("Main text", 4), rep("Supplement", 3)),
  row.names = col_order
)
ann_colors_full <- list(Section = c("Main text" = "#2166AC", "Supplement" = "#B2182B"))

png(file.path(fig_dir, "40_stability_heatmap_full.png"),
    width = 12, height = 5.5, units = "in", res = 300)
pheatmap(
  heatmap_mat_full,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = sig_labels_full,
  number_color = "black",
  fontsize_number = 10,
  color = colors,
  breaks = breaks_full,
  annotation_col = col_ann_full,
  annotation_colors = ann_colors_full,
  main = "Pathway NES Stability — All 8 Analyses",
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 8,
  angle_col = 45,
  cellwidth = 35,
  cellheight = 20
)
dev.off()
cat("  Saved: 40_stability_heatmap_full.png\n")


# --- Figure 3: Female PCA ---
cat("Figure 3: Female PCA\n")

# PC1 vs PC2, colored by CKD
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = factor(MHRNLFLR), shape = factor(MHRNLFLR)),
             alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("0" = "grey60", "1" = "firebrick3"),
    labels = c("0" = "Control", "1" = "CKD"),
    name = ""
  ) +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17),
    labels = c("0" = "Control", "1" = "CKD"),
    name = ""
  ) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "PCA of Female Subcutaneous Adipose"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "40_female_pca.png"), p_pca,
       width = 6, height = 5, dpi = 300)
cat("  Saved: 40_female_pca.png\n")

# PCA with ischemic time overlay
p_pca_isch <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = ischemic_hrs, shape = factor(MHRNLFLR)),
             alpha = 0.7, size = 2) +
  scale_color_viridis_c(name = "Ischemic\ntime (hrs)") +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17),
    labels = c("0" = "Control", "1" = "CKD"),
    name = ""
  ) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "PCA — Ischemic Time Overlay"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

ggsave(file.path(fig_dir, "40_female_pca_ischemic.png"), p_pca_isch,
       width = 7, height = 5, dpi = 300)
cat("  Saved: 40_female_pca_ischemic.png\n")


# --- Figure 4: Interaction vs Female NES scatter ---
cat("Figure 4: Interaction vs Female NES scatter\n")

compare_plot <- compare_ix %>%
  as_tibble() %>%
  mutate(pathway_clean = clean_pathway(pathway)) %>%
  mutate(
    both_sig = padj_female < 0.05 & padj_interaction < 0.05,
    female_only = padj_female < 0.05 & padj_interaction >= 0.05,
    label = ifelse(both_sig | female_only, pathway_clean, "")
  )

p_scatter <- ggplot(compare_plot, aes(x = NES_female, y = NES_interaction)) +
  geom_point(aes(color = both_sig), size = 2.5, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  ggrepel::geom_text_repel(
    aes(label = label), size = 2.8, max.overlaps = 20, seed = 42
  ) +
  scale_color_manual(
    values = c("TRUE" = "firebrick3", "FALSE" = "grey60"),
    labels = c("TRUE" = "Both FDR < 0.05", "FALSE" = "Not both significant"),
    name = ""
  ) +
  labs(
    x = "NES — Full Female CKD Model",
    y = "NES — SEX × CKD Interaction",
    title = sprintf("Female CKD Signal = Sex-Dimorphic Signal (rho = %.2f)", rho_ix)
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11))

ggsave(file.path(fig_dir, "40_interaction_scatter.png"), p_scatter,
       width = 7, height = 6, dpi = 300)
cat("  Saved: 40_interaction_scatter.png\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SAVE RESULTS                                                           ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("Saving stage 2 results\n")
cat(strrep("=", 70), "\n")

stage2 <- list(
  # GSEA results
  gsea_female_full = gsea_female_full,
  gsea_interaction = gsea_interaction,
  gsea_male = gsea_male,
  gsea_psm = gsea_psm,

  # Comparisons
  interaction_vs_female = compare_ix,
  rho_interaction_female = rho_ix,

  # Stability table components
  nes_matrix = nes_matrix,
  padj_matrix = padj_matrix,
  all_gsea = all_gsea,
  core_pathways = core_pathways,
  direction_consistent = consistent,

  # PCA
  pca = pca,
  pca_df = pca_df,
  var_explained = var_explained
)

saveRDS(stage2, file.path(res_dir, "stage2_results.rds"))
cat("Saved:", file.path(res_dir, "stage2_results.rds"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SUMMARY                                                                ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("STAGE 2 SUMMARY\n")
cat(strrep("=", 70), "\n\n")

cat("GSEA results:\n")
cat(sprintf("  Female full model: %d pathways FDR < 0.05\n",
            sum(gsea_female_full$padj < 0.05)))
cat(sprintf("  PSM-matched female: %d pathways FDR < 0.05\n", n_sig_psm))
cat(sprintf("  SEX×CKD interaction: %d pathways FDR < 0.05\n", n_sig))
cat(sprintf("  Male CKD: %d pathways FDR < 0.05\n", n_sig_male))
cat(sprintf("  Interaction vs Female NES rho = %.3f\n\n", rho_ix))

cat("Direction consistency (main-text columns):\n")
for (pw in names(consistent)) {
  nes_vals <- round(as.numeric(nes_matrix[pw, main_cols]), 2)
  status <- ifelse(consistent[pw], "CONSISTENT", "INCONSISTENT")
  cat(sprintf("  %-30s %s  [%s]\n", pw, status,
              paste(nes_vals, collapse = ", ")))
}

cat(sprintf("\nPCA: CKD separates on PC%d (p = %.3g)\n",
            which.min(sapply(1:5, function(i) t.test(pca_df[[paste0("PC", i)]] ~ pca_df$MHRNLFLR)$p.value)),
            min(sapply(1:5, function(i) t.test(pca_df[[paste0("PC", i)]] ~ pca_df$MHRNLFLR)$p.value))))

cat("\nFiles saved:\n")
cat("  results/stage2_results.rds\n")
cat("  results/pathway_stability_table.csv\n")
cat("  results/pathway_stability_full.csv\n")
cat("  figures/40_pathway_dotplot.png\n")
cat("  figures/40_stability_heatmap.png\n")
cat("  figures/40_stability_heatmap_full.png\n")
cat("  figures/40_female_pca.png\n")
cat("  figures/40_female_pca_ischemic.png\n")
cat("  figures/40_interaction_scatter.png\n")

cat("\n", strrep("=", 70), "\n")
cat("Stage 2 complete.\n")
cat(strrep("=", 70), "\n")
