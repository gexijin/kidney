## 05_go_enrichment.R — GO BP, KEGG, and GSEA for CKD vs Control
## Uses clusterProfiler with proper background gene list

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(enrichplot)

base_dir <- "/Users/ge/Documents/research/gtex/analyses/Soulage_GSE98603"
data_dir <- file.path(base_dir, "data")
plots_dir <- file.path(base_dir, "results", "plots")
tables_dir <- file.path(base_dir, "results", "tables")

de_clean <- readRDS(file.path(data_dir, "04_de_clean.rds"))

cat("Gene-level DE results:", nrow(de_clean), "genes\n")
cat("With Entrez ID:", sum(!is.na(de_clean$entrez_id)), "\n")

# ── 0. Background genes ─────────────────────────────────────────────────────

background_genes <- de_clean %>%
  filter(!is.na(entrez_id)) %>%
  pull(entrez_id) %>%
  unique() %>%
  as.character()

cat("Background genes:", length(background_genes), "\n")

# Significant gene lists
sig_up <- de_clean %>%
  filter(adj.P.Val < 0.05 & logFC > 0 & !is.na(entrez_id)) %>%
  pull(entrez_id) %>% unique() %>% as.character()

sig_down <- de_clean %>%
  filter(adj.P.Val < 0.05 & logFC < 0 & !is.na(entrez_id)) %>%
  pull(entrez_id) %>% unique() %>% as.character()

cat("Significant up:", length(sig_up), "\n")
cat("Significant down:", length(sig_down), "\n")

# ── Helper: prepare ORA data for combined plots ─────────────────────────────

prepare_ora_data <- function(enrichResult, direction_label) {
  if (is.null(enrichResult) || nrow(as.data.frame(enrichResult)) == 0) return(NULL)
  df <- as.data.frame(enrichResult)
  df %>% mutate(
    gene_ratio_num = as.numeric(sub("/.*", "", GeneRatio)),
    gene_ratio_denom = as.numeric(sub(".*/", "", GeneRatio)),
    bg_ratio_num = as.numeric(sub("/.*", "", BgRatio)),
    bg_ratio_denom = as.numeric(sub(".*/", "", BgRatio)),
    GeneRatioValue = gene_ratio_num / gene_ratio_denom,
    BgRatioValue = bg_ratio_num / bg_ratio_denom,
    FoldEnrichment = GeneRatioValue / BgRatioValue,
    direction = direction_label
  )
}

plot_combined_ora <- function(df, n_show = 10, title = "") {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  top <- df %>%
    group_by(direction) %>%
    slice_min(p.adjust, n = n_show) %>%
    ungroup()

  top$Description <- factor(top$Description, levels = rev(unique(top$Description)))
  dir_colors <- c(UP = "#E41A1C", DOWN = "#377EB8")

  ggplot(top, aes(FoldEnrichment, Description, size = Count, color = -log10(p.adjust))) +
    geom_point() +
    scale_color_viridis_c(name = "-log10(FDR)") +
    scale_size_continuous(range = c(2, 8)) +
    facet_wrap(~direction, scales = "free_y") +
    labs(title = title, x = "Fold Enrichment") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold", size = 12))
}

# ── 1. GO Biological Process ────────────────────────────────────────────────

cat("\n=== GO BP Enrichment ===\n")

ego_up <- ego_down <- NULL

if (length(sig_up) >= 5) {
  ego_up <- enrichGO(gene = sig_up, universe = background_genes,
                     OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                     ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  ego_up <- clusterProfiler::simplify(ego_up, cutoff = 0.7, by = "p.adjust")
  cat("GO BP up:", nrow(as.data.frame(ego_up)), "terms\n")
  write.csv(as.data.frame(ego_up), file.path(tables_dir, "05_GO_BP_upregulated.csv"), row.names = FALSE)
}

if (length(sig_down) >= 5) {
  ego_down <- enrichGO(gene = sig_down, universe = background_genes,
                       OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                       ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  ego_down <- clusterProfiler::simplify(ego_down, cutoff = 0.7, by = "p.adjust")
  cat("GO BP down:", nrow(as.data.frame(ego_down)), "terms\n")
  write.csv(as.data.frame(ego_down), file.path(tables_dir, "05_GO_BP_downregulated.csv"), row.names = FALSE)
}

# Combined GO plot
go_combined <- bind_rows(
  prepare_ora_data(ego_up, "UP"),
  prepare_ora_data(ego_down, "DOWN")
)

if (!is.null(go_combined) && nrow(go_combined) > 0) {
  p_go <- plot_combined_ora(go_combined, n_show = 15, title = "GO BP: CKD vs Control — Subcutaneous WAT")
  if (!is.null(p_go)) {
    ggsave(file.path(plots_dir, "05_GO_BP_combined.png"), p_go, width = 14, height = 10, dpi = 300)
  }
}

# ── 2. KEGG Pathways ────────────────────────────────────────────────────────

cat("\n=== KEGG Enrichment ===\n")

kegg_up <- kegg_down <- NULL

if (length(sig_up) >= 5) {
  kegg_up <- enrichKEGG(gene = sig_up, universe = background_genes,
                        organism = "hsa", keyType = "ncbi-geneid",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05)
  cat("KEGG up:", nrow(as.data.frame(kegg_up)), "pathways\n")
  write.csv(as.data.frame(kegg_up), file.path(tables_dir, "05_KEGG_upregulated.csv"), row.names = FALSE)
}

if (length(sig_down) >= 5) {
  kegg_down <- enrichKEGG(gene = sig_down, universe = background_genes,
                          organism = "hsa", keyType = "ncbi-geneid",
                          pAdjustMethod = "BH", pvalueCutoff = 0.05)
  cat("KEGG down:", nrow(as.data.frame(kegg_down)), "pathways\n")
  write.csv(as.data.frame(kegg_down), file.path(tables_dir, "05_KEGG_downregulated.csv"), row.names = FALSE)
}

kegg_combined <- bind_rows(
  prepare_ora_data(kegg_up, "UP"),
  prepare_ora_data(kegg_down, "DOWN")
)

if (!is.null(kegg_combined) && nrow(kegg_combined) > 0) {
  p_kegg <- plot_combined_ora(kegg_combined, n_show = 15, title = "KEGG: CKD vs Control — Subcutaneous WAT")
  if (!is.null(p_kegg)) {
    ggsave(file.path(plots_dir, "05_KEGG_combined.png"), p_kegg, width = 14, height = 10, dpi = 300)
  }
}

# ── 3. GSEA ──────────────────────────────────────────────────────────────────

cat("\n=== GSEA ===\n")

# Create ranked gene list (t-statistic)
gene_list <- de_clean %>%
  filter(!is.na(entrez_id)) %>%
  arrange(desc(t)) %>%
  { setNames(.$t, .$entrez_id) }

# Remove duplicates
gene_list <- gene_list[!duplicated(names(gene_list))]
cat("GSEA gene list:", length(gene_list), "genes\n")

# GO BP GSEA
gsea_go <- gseGO(geneList = gene_list,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "BP",
                 minGSSize = 15,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE,
                 seed = TRUE)

gsea_go_df <- as.data.frame(gsea_go)
cat("GSEA GO BP terms:", nrow(gsea_go_df), "\n")
cat("  Activated:", sum(gsea_go_df$NES > 0), "\n")
cat("  Suppressed:", sum(gsea_go_df$NES < 0), "\n")

write.csv(gsea_go_df, file.path(tables_dir, "05_GSEA_GO_BP.csv"), row.names = FALSE)

if (nrow(gsea_go_df) > 0) {
  # Dotplot
  p_dot <- dotplot(gsea_go, showCategory = 20, split = ".sign") +
    facet_grid(~.sign) +
    ggtitle("GSEA GO BP: CKD vs Control")
  ggsave(file.path(plots_dir, "05_GSEA_GO_BP_dotplot.png"), p_dot, width = 14, height = 10, dpi = 300)

  # Ridgeplot
  p_ridge <- ridgeplot(gsea_go, showCategory = 20) +
    ggtitle("GSEA GO BP Ridge Plot")
  ggsave(file.path(plots_dir, "05_GSEA_GO_BP_ridgeplot.png"), p_ridge, width = 12, height = 10, dpi = 300)

  # Top enrichment plots
  top_ids <- gsea_go_df %>%
    arrange(pvalue) %>%
    head(6) %>%
    pull(ID)

  for (i in seq_along(top_ids)) {
    p_run <- gseaplot2(gsea_go, geneSetID = top_ids[i], title = gsea_go_df$Description[gsea_go_df$ID == top_ids[i]])
    ggsave(file.path(plots_dir, sprintf("05_GSEA_running_%d.png", i)), p_run, width = 8, height = 6, dpi = 300)
  }
}

# KEGG GSEA
gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = "hsa",
                     keyType = "ncbi-geneid",
                     minGSSize = 15,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose = FALSE,
                     seed = TRUE)

gsea_kegg_df <- as.data.frame(gsea_kegg)
cat("\nGSEA KEGG pathways:", nrow(gsea_kegg_df), "\n")
write.csv(gsea_kegg_df, file.path(tables_dir, "05_GSEA_KEGG.csv"), row.names = FALSE)

if (nrow(gsea_kegg_df) > 0) {
  p_kegg_dot <- dotplot(gsea_kegg, showCategory = 20, split = ".sign") +
    facet_grid(~.sign) +
    ggtitle("GSEA KEGG: CKD vs Control")
  ggsave(file.path(plots_dir, "05_GSEA_KEGG_dotplot.png"), p_kegg_dot, width = 14, height = 10, dpi = 300)
}

# ── 4. Print top GSEA terms for manuscript context ──────────────────────────

cat("\n=== Top Activated GO BP Terms (CKD > Control) ===\n")
gsea_go_df %>%
  filter(NES > 0) %>%
  arrange(pvalue) %>%
  head(15) %>%
  select(Description, NES, pvalue, p.adjust, setSize) %>%
  print(digits = 3)

cat("\n=== Top Suppressed GO BP Terms (CKD < Control) ===\n")
gsea_go_df %>%
  filter(NES < 0) %>%
  arrange(pvalue) %>%
  head(15) %>%
  select(Description, NES, pvalue, p.adjust, setSize) %>%
  print(digits = 3)

cat("\n✓ Enrichment analysis complete\n")
