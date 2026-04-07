#!/usr/bin/env Rscript
# =============================================================================
# Vascular Calcification & Osteogenic Program in CKD Subcutaneous Adipose
#
# Tests whether CKD activates a coordinated pro-calcification / vascular
# remodeling program in female SAT, and contrasts with male SAT and female VAT.
#
# Sections:
#   1 — Load DE results, gene annotation, GO:BP gene sets
#   2 — Full GO:BP GSEA for three contrasts (female SAT, male SAT, female VAT)
#   3 — Curated calcification pathway table (three contrasts)
#   4 — Curated gene-level heatmap (functional categories × three contrasts)
#   5 — Top GO:BP pathways in female SAT (contextualize calcification signal)
#   6 — Dot plot and heatmap figures
# =============================================================================

library(tidyverse)
library(fgsea)
library(msigdbr)
library(pheatmap)
library(RColorBrewer)

source("kidney/R/functions.R")

set.seed(42)

res_dir  <- "kidney/results"
fig_dir  <- "kidney/figures"
ann_file <- "data/annotation/gencode_v47_genes.rds"


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1: LOAD DATA                                                   ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("Section 1: Load data")

# Gene annotation (unversioned Ensembl IDs as rownames)
ann <- readRDS(ann_file)
id2sym <- setNames(ann$symbol, rownames(ann))

# DE results — extract tt (topTable) from each
de_female_sat <- readRDS(file.path(res_dir, "de_female.rds"))$tt
de_male_sat   <- readRDS(file.path(res_dir, "de_male.rds"))$tt
de_vat_female <- readRDS(file.path(res_dir, "de_vat_female.rds"))$tt

cat(sprintf("Female SAT: %d genes, %d DEGs (FDR<0.05)\n",
            nrow(de_female_sat), sum(de_female_sat$adj.P.Val < 0.05)))
cat(sprintf("Male SAT:   %d genes, %d DEGs (FDR<0.05)\n",
            nrow(de_male_sat), sum(de_male_sat$adj.P.Val < 0.05)))
cat(sprintf("Female VAT: %d genes, %d DEGs (FDR<0.05)\n",
            nrow(de_vat_female), sum(de_vat_female$adj.P.Val < 0.05)))

# Add gene symbols
add_symbol <- function(tt) {
  tt$symbol <- id2sym[tt$gene_id]
  tt
}
de_female_sat <- add_symbol(de_female_sat)
de_male_sat   <- add_symbol(de_male_sat)
de_vat_female <- add_symbol(de_vat_female)

# GO:BP gene sets from msigdbr (Ensembl IDs)
cat("\nLoading GO:BP gene sets from msigdbr...\n")
gobp_df <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  filter(ensembl_gene != "")

gobp_list <- split(gobp_df$ensembl_gene, gobp_df$gs_name)
cat(sprintf("GO:BP gene sets: %d\n", length(gobp_list)))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2: FULL GO:BP GSEA — THREE CONTRASTS                           ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("Section 2: GO:BP GSEA — three contrasts")

gsea_f_sat <- run_gsea(de_female_sat, gobp_list, label = "Female SAT", seed = 42)
cat(sprintf("Female SAT GO:BP GSEA: %d terms, %d FDR<0.05\n",
            nrow(gsea_f_sat), sum(gsea_f_sat$padj < 0.05)))

gsea_m_sat <- run_gsea(de_male_sat, gobp_list, label = "Male SAT", seed = 42)
cat(sprintf("Male SAT GO:BP GSEA:   %d terms, %d FDR<0.05\n",
            nrow(gsea_m_sat), sum(gsea_m_sat$padj < 0.05)))

gsea_f_vat <- run_gsea(de_vat_female, gobp_list, label = "Female VAT", seed = 42)
cat(sprintf("Female VAT GO:BP GSEA: %d terms, %d FDR<0.05\n",
            nrow(gsea_f_vat), sum(gsea_f_vat$padj < 0.05)))

# Save full GSEA results
saveRDS(gsea_f_sat, file.path(res_dir, "gsea_gobp_female_sat.rds"))
saveRDS(gsea_m_sat, file.path(res_dir, "gsea_gobp_male_sat.rds"))
saveRDS(gsea_f_vat, file.path(res_dir, "gsea_gobp_vat_female.rds"))
cat("Saved full GO:BP GSEA results.\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 3: CURATED CALCIFICATION PATHWAY TABLE                         ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("Section 3: Curated calcification pathway table")

# Essential calcification-related GO:BP terms
calc_terms <- c(
  "GOBP_OSSIFICATION"                    = "Ossification",
  "GOBP_BONE_MINERALIZATION"             = "Bone mineralization",
  "GOBP_OSTEOBLAST_DIFFERENTIATION"      = "Osteoblast differentiation",
  "GOBP_REGULATION_OF_BMP_SIGNALING_PATHWAY" = "BMP signaling",
  "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY" = "TGF-beta signaling",
  "GOBP_WNT_SIGNALING_PATHWAY"           = "Wnt signaling",
  "GOBP_COLLAGEN_FIBRIL_ORGANIZATION"    = "Collagen fibril organization",
  "GOBP_BLOOD_VESSEL_REMODELING"         = "Blood vessel remodeling"
)

# Check which terms exist in the GSEA results
available <- names(calc_terms)[names(calc_terms) %in% gsea_f_sat$pathway]
missing   <- names(calc_terms)[!names(calc_terms) %in% gsea_f_sat$pathway]
cat(sprintf("Curated terms: %d available, %d missing\n",
            length(available), length(missing)))
if (length(missing) > 0) {
  cat("Missing terms:\n")
  for (m in missing) cat(sprintf("  %s\n", m))
}

# Use available terms
calc_terms <- calc_terms[available]

# Build the contrast table
build_contrast_row <- function(term, gsea1, gsea2, gsea3) {
  get_vals <- function(gsea, term) {
    row <- gsea[gsea$pathway == term, ]
    if (nrow(row) == 0) return(c(NES = NA_real_, padj = NA_real_))
    c(NES = row$NES, padj = row$padj)
  }
  v1 <- get_vals(gsea1, term)
  v2 <- get_vals(gsea2, term)
  v3 <- get_vals(gsea3, term)
  tibble(
    pathway       = term,
    label         = calc_terms[term],
    NES_F_SAT     = v1["NES"],  padj_F_SAT     = v1["padj"],
    NES_M_SAT     = v2["NES"],  padj_M_SAT     = v2["padj"],
    NES_F_VAT     = v3["NES"],  padj_F_VAT     = v3["padj"]
  )
}

calc_table <- bind_rows(lapply(names(calc_terms), build_contrast_row,
                               gsea1 = gsea_f_sat,
                               gsea2 = gsea_m_sat,
                               gsea3 = gsea_f_vat))

# Add significance stars
add_stars <- function(padj) {
  case_when(
    is.na(padj) ~ "",
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE          ~ ""
  )
}

calc_table <- calc_table %>%
  mutate(
    F_SAT = sprintf("%+.2f%s", NES_F_SAT, add_stars(padj_F_SAT)),
    M_SAT = sprintf("%+.2f%s", NES_M_SAT, add_stars(padj_M_SAT)),
    F_VAT = sprintf("%+.2f%s", NES_F_VAT, add_stars(padj_F_VAT))
  )

# Print table
cat("\n--- Calcification GO:BP Pathways Across Contrasts ---\n")
cat(sprintf("%-35s %12s %12s %12s\n", "Pathway", "Female SAT", "Male SAT", "Female VAT"))
cat(strrep("-", 73), "\n")
for (i in seq_len(nrow(calc_table))) {
  cat(sprintf("%-35s %12s %12s %12s\n",
              calc_table$label[i],
              calc_table$F_SAT[i],
              calc_table$M_SAT[i],
              calc_table$F_VAT[i]))
}

# Count significant terms per contrast
cat(sprintf("\nSignificant terms (FDR<0.05): Female SAT=%d, Male SAT=%d, Female VAT=%d\n",
            sum(calc_table$padj_F_SAT < 0.05, na.rm = TRUE),
            sum(calc_table$padj_M_SAT < 0.05, na.rm = TRUE),
            sum(calc_table$padj_F_VAT < 0.05, na.rm = TRUE)))

# NES correlation between contrasts (curated terms only)
complete_fm <- complete.cases(calc_table$NES_F_SAT, calc_table$NES_M_SAT)
complete_fv <- complete.cases(calc_table$NES_F_SAT, calc_table$NES_F_VAT)

if (sum(complete_fm) >= 5) {
  rho_fm <- cor(calc_table$NES_F_SAT[complete_fm],
                calc_table$NES_M_SAT[complete_fm], method = "spearman")
  cat(sprintf("NES rho (Female SAT vs Male SAT): %.3f\n", rho_fm))
}
if (sum(complete_fv) >= 5) {
  rho_fv <- cor(calc_table$NES_F_SAT[complete_fv],
                calc_table$NES_F_VAT[complete_fv], method = "spearman")
  cat(sprintf("NES rho (Female SAT vs Female VAT): %.3f\n", rho_fv))
}

# Save table
write_csv(calc_table %>% dplyr::select(pathway, label, starts_with("NES_"), starts_with("padj_")),
          file.path(res_dir, "calcification_pathway_table.csv"))
cat("Saved:", file.path(res_dir, "calcification_pathway_table.csv"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 4: CURATED GENE-LEVEL HEATMAP                                  ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("Section 4: Curated gene-level heatmap")

# Key genes per functional category (significant or biologically critical)
gene_categories <- list(
  "Calcification inhibitors" = c("ABCC6", "ENPP1", "VKORC1L1"),
  "Osteogenic program"       = c("RUNX2", "ALPL", "TGFB2", "TGFB3"),
  "Vascular smooth muscle"   = c("ACTA2", "CNN1", "DES", "TAGLN", "SRF"),
  "Wnt / BMP signaling"      = c("SFRP2", "GREM1", "GREM2", "BMPR1B"),
  "ECM remodeling"           = c("FN1", "FNDC1", "COL4A3", "THBS2", "TIMP1"),
  "Angiogenic imbalance"     = c("VEGFB", "ANGPT1", "HIF1A")
)

# Build symbol-to-ensembl lookup
sym2id <- setNames(rownames(ann), ann$symbol)

# Extract logFC for each gene × contrast
extract_gene_stats <- function(tt, symbols, sym2id) {
  ids <- sym2id[symbols]
  matched <- tt[match(ids, tt$gene_id), ]
  tibble(
    symbol  = symbols,
    gene_id = ids,
    logFC   = matched$logFC,
    pval    = matched$P.Value,
    fdr     = matched$adj.P.Val
  )
}

gene_stats_f_sat <- extract_gene_stats(de_female_sat, unlist(gene_categories), sym2id)
gene_stats_m_sat <- extract_gene_stats(de_male_sat,   unlist(gene_categories), sym2id)
gene_stats_f_vat <- extract_gene_stats(de_vat_female,  unlist(gene_categories), sym2id)

# Combine into wide format
gene_wide <- gene_stats_f_sat %>%
  dplyr::select(symbol, logFC_F_SAT = logFC, fdr_F_SAT = fdr) %>%
  left_join(gene_stats_m_sat %>% dplyr::select(symbol, logFC_M_SAT = logFC, fdr_M_SAT = fdr),
            by = "symbol") %>%
  left_join(gene_stats_f_vat %>% dplyr::select(symbol, logFC_F_VAT = logFC, fdr_F_VAT = fdr),
            by = "symbol")

# Add category
gene_wide$category <- rep(names(gene_categories), times = lengths(gene_categories))

# Filter to genes present in at least one dataset
gene_wide <- gene_wide %>%
  filter(!is.na(logFC_F_SAT) | !is.na(logFC_M_SAT) | !is.na(logFC_F_VAT))

cat(sprintf("Curated genes with data: %d / %d\n",
            nrow(gene_wide), length(unlist(gene_categories))))

# Print gene-level table
cat("\n--- Gene-Level logFC Across Contrasts ---\n")
cat(sprintf("%-8s %-35s %8s %8s %8s %8s %8s %8s\n",
            "Symbol", "Category", "lFC_F", "FDR_F", "lFC_M", "FDR_M", "lFC_V", "FDR_V"))
cat(strrep("-", 115), "\n")
for (i in seq_len(nrow(gene_wide))) {
  row <- gene_wide[i, ]
  cat(sprintf("%-8s %-35s %+8.3f %8.4f %+8.3f %8.4f %+8.3f %8.4f\n",
              row$symbol, row$category,
              ifelse(is.na(row$logFC_F_SAT), 0, row$logFC_F_SAT),
              ifelse(is.na(row$fdr_F_SAT), 1, row$fdr_F_SAT),
              ifelse(is.na(row$logFC_M_SAT), 0, row$logFC_M_SAT),
              ifelse(is.na(row$fdr_M_SAT), 1, row$fdr_M_SAT),
              ifelse(is.na(row$logFC_F_VAT), 0, row$logFC_F_VAT),
              ifelse(is.na(row$fdr_F_VAT), 1, row$fdr_F_VAT)))
}

# Count significant genes per contrast
cat(sprintf("\nSignificant genes (FDR<0.05): Female SAT=%d, Male SAT=%d, Female VAT=%d\n",
            sum(gene_wide$fdr_F_SAT < 0.05, na.rm = TRUE),
            sum(gene_wide$fdr_M_SAT < 0.05, na.rm = TRUE),
            sum(gene_wide$fdr_F_VAT < 0.05, na.rm = TRUE)))

# Save gene table
write_csv(gene_wide, file.path(res_dir, "calcification_gene_table.csv"))
cat("Saved:", file.path(res_dir, "calcification_gene_table.csv"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 5: TOP GO:BP IN FEMALE SAT — CALCIFICATION IN CONTEXT          ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("Section 5: Top GO:BP pathways in Female SAT")

gsea_sig <- gsea_f_sat %>%
  as_tibble() %>%
  filter(padj < 0.05) %>%
  arrange(padj)

cat(sprintf("Total significant GO:BP terms (FDR<0.05): %d\n", nrow(gsea_sig)))

# Top 20 upregulated
top_up <- gsea_sig %>% filter(NES > 0) %>% head(20)
cat("\n--- Top 20 Upregulated GO:BP (Female SAT CKD) ---\n")
for (i in seq_len(nrow(top_up))) {
  is_calc <- top_up$pathway[i] %in% names(calc_terms)
  marker  <- ifelse(is_calc, " [CALC]", "")
  cat(sprintf("  %2d. %-55s NES=%+.2f  q=%.1e%s\n",
              i, gsub("GOBP_", "", top_up$pathway[i]),
              top_up$NES[i], top_up$padj[i], marker))
}

# Top 20 downregulated
top_dn <- gsea_sig %>% filter(NES < 0) %>% arrange(NES) %>% head(20)
cat("\n--- Top 20 Downregulated GO:BP (Female SAT CKD) ---\n")
for (i in seq_len(nrow(top_dn))) {
  is_calc <- top_dn$pathway[i] %in% names(calc_terms)
  marker  <- ifelse(is_calc, " [CALC]", "")
  cat(sprintf("  %2d. %-55s NES=%-.2f  q=%.1e%s\n",
              i, gsub("GOBP_", "", top_dn$pathway[i]),
              top_dn$NES[i], top_dn$padj[i], marker))
}

# How many of our curated calcification terms are in the top 100?
top100_up <- gsea_sig %>% filter(NES > 0) %>% head(100) %>% pull(pathway)
n_calc_top100 <- sum(names(calc_terms) %in% top100_up)
cat(sprintf("\nCurated calcification terms in top 100 upregulated: %d / %d\n",
            n_calc_top100, length(calc_terms)))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 6: FIGURES                                                      ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("Section 6: Figures")

# --- Figure 1: Calcification pathway dot plot (three contrasts) ---
cat("Figure 1: Calcification pathway dot plot\n")

dot_data <- calc_table %>%
  dplyr::select(label, starts_with("NES_"), starts_with("padj_")) %>%
  pivot_longer(
    cols = -label,
    names_to = c(".value", "contrast"),
    names_pattern = "(NES|padj)_(.*)"
  ) %>%
  mutate(
    contrast = factor(contrast,
                      levels = c("F_SAT", "M_SAT", "F_VAT"),
                      labels = c("Female SAT", "Male SAT", "Female VAT")),
    neg_log_padj = -log10(pmax(padj, 1e-10)),
    sig = padj < 0.05
  )

# Order pathways by female SAT NES (descending so highest is at top)
pathway_order <- calc_table %>% arrange(NES_F_SAT) %>% pull(label)
dot_data$label <- factor(dot_data$label, levels = pathway_order)

nes_lim <- max(abs(dot_data$NES), na.rm = TRUE)

p_dot <- ggplot(dot_data, aes(x = contrast, y = label)) +
  geom_point(aes(size = neg_log_padj, fill = NES),
             shape = 21, color = "grey70", stroke = 0.3) +
  geom_point(data = dot_data %>% filter(sig),
             aes(size = neg_log_padj, fill = NES),
             shape = 21, stroke = 0.9, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = "NES", limits = c(-nes_lim, nes_lim)
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name = expression(-log[10](FDR)),
    breaks = c(1, 2, 3)
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x     = element_text(size = 13, face = "bold"),
    axis.text.y     = element_text(size = 12),
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggsave(file.path(fig_dir, "50_calcification_dotplot.png"), p_dot,
       width = 6.5, height = 4.5, dpi = 300)
cat("  Saved: 50_calcification_dotplot.png\n")


# --- Figure 2: Gene-level heatmap (ggplot2 with facet strip labels) ---
cat("Figure 2: Gene-level heatmap\n")

hm_genes <- gene_wide %>% filter(!is.na(logFC_F_SAT))

# Build long-format data for ggplot
hm_long <- hm_genes %>%
  dplyr::select(symbol, category,
                logFC_F_SAT, logFC_M_SAT, logFC_F_VAT,
                fdr_F_SAT, fdr_M_SAT, fdr_F_VAT) %>%
  pivot_longer(
    cols = starts_with("logFC_"),
    names_to = "contrast", names_prefix = "logFC_",
    values_to = "logFC"
  ) %>%
  left_join(
    hm_genes %>%
      dplyr::select(symbol, fdr_F_SAT, fdr_M_SAT, fdr_F_VAT) %>%
      pivot_longer(cols = starts_with("fdr_"),
                   names_to = "contrast", names_prefix = "fdr_",
                   values_to = "fdr"),
    by = c("symbol", "contrast")
  ) %>%
  dplyr::select(symbol, category, contrast, logFC, fdr) %>%
  mutate(
    logFC = replace_na(logFC, 0),
    stars = case_when(
      is.na(fdr)    ~ "",
      fdr < 0.001   ~ "***",
      fdr < 0.01    ~ "**",
      fdr < 0.05    ~ "*",
      TRUE          ~ ""
    ),
    contrast = factor(contrast,
                      levels = c("F_SAT", "M_SAT", "F_VAT"),
                      labels = c("Female\nSAT", "Male\nSAT", "Female\nVAT")),
    category = factor(category, levels = names(gene_categories))
  )

# Order genes within each category (preserve input order)
gene_order <- unlist(lapply(names(gene_categories), function(cat) {
  hm_genes$symbol[hm_genes$category == cat]
}))
gene_order <- gene_order[gene_order %in% hm_genes$symbol]
hm_long$symbol <- factor(hm_long$symbol, levels = rev(gene_order))

# Symmetric color limit
lfc_lim <- max(abs(hm_long$logFC), na.rm = TRUE)

# Short category labels for strip text
cat_labels <- c(
  "Calcification inhibitors" = "Calcification\ninhibitors",
  "Osteogenic program"       = "Osteogenic\nprogram",
  "Vascular smooth muscle"   = "Vascular\nsmooth muscle",
  "Wnt / BMP signaling"      = "Wnt / BMP\nsignaling",
  "ECM remodeling"           = "ECM\nremodeling",
  "Angiogenic imbalance"     = "Angiogenic\nimbalance"
)

p_hm <- ggplot(hm_long, aes(x = contrast, y = symbol, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = stars), size = 4, vjust = 0.75) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y",
             switch = "y",
             labeller = labeller(category = cat_labels)) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-lfc_lim, lfc_lim),
    name = "logFC"
  ) +
  scale_x_discrete(position = "bottom") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 0, hjust = 1, size = 10,
                                       face = "bold", lineheight = 0.9),
    axis.text.y        = element_text(size = 11, face = "italic"),
    axis.text.x        = element_text(size = 12, face = "bold", lineheight = 0.9),
    axis.ticks         = element_blank(),
    panel.grid         = element_blank(),
    panel.spacing.y    = unit(3, "pt"),
    legend.position    = "right",
    legend.title       = element_text(size = 11),
    legend.text        = element_text(size = 10),
    plot.margin        = margin(5, 5, 5, 5)
  )

ggsave(file.path(fig_dir, "50_calcification_gene_heatmap.png"), p_hm,
       width = 5.5, height = 7, dpi = 300)
cat("  Saved: 50_calcification_gene_heatmap.png\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SUMMARY                                                                ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("SUMMARY")

cat("GO:BP GSEA (FDR<0.05 terms):\n")
cat(sprintf("  Female SAT: %d\n", sum(gsea_f_sat$padj < 0.05)))
cat(sprintf("  Male SAT:   %d\n", sum(gsea_m_sat$padj < 0.05)))
cat(sprintf("  Female VAT: %d\n", sum(gsea_f_vat$padj < 0.05)))

cat("\nCalcification pathways enriched (FDR<0.05):\n")
cat(sprintf("  Female SAT: %d / %d\n",
            sum(calc_table$padj_F_SAT < 0.05, na.rm = TRUE), nrow(calc_table)))
cat(sprintf("  Male SAT:   %d / %d\n",
            sum(calc_table$padj_M_SAT < 0.05, na.rm = TRUE), nrow(calc_table)))
cat(sprintf("  Female VAT: %d / %d\n",
            sum(calc_table$padj_F_VAT < 0.05, na.rm = TRUE), nrow(calc_table)))

cat("\nCalcification genes significant (FDR<0.05):\n")
cat(sprintf("  Female SAT: %d / %d\n",
            sum(gene_wide$fdr_F_SAT < 0.05, na.rm = TRUE), nrow(gene_wide)))
cat(sprintf("  Male SAT:   %d / %d\n",
            sum(gene_wide$fdr_M_SAT < 0.05, na.rm = TRUE), nrow(gene_wide)))
cat(sprintf("  Female VAT: %d / %d\n",
            sum(gene_wide$fdr_F_VAT < 0.05, na.rm = TRUE), nrow(gene_wide)))

cat("\nFiles saved:\n")
cat("  results/gsea_gobp_female_sat.rds\n")
cat("  results/gsea_gobp_male_sat.rds\n")
cat("  results/gsea_gobp_vat_female.rds\n")
cat("  results/calcification_pathway_table.csv\n")
cat("  results/calcification_gene_table.csv\n")
cat("  figures/50_calcification_dotplot.png\n")
cat("  figures/50_calcification_gene_heatmap.png\n")

cat("\n", strrep("=", 70), "\n")
cat("Script 50 complete.\n")
cat(strrep("=", 70), "\n")
