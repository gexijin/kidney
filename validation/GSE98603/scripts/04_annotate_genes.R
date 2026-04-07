## 04_annotate_genes.R — Map Agilent probes to gene symbols via GPL annotation
## GPL13497 (Agilent 4x44K v2) — use fData from GEO ExpressionSet

library(tidyverse)
library(AnnotationDbi)

base_dir <- "/Users/ge/Documents/research/gtex/analyses/Soulage_GSE98603"
data_dir <- file.path(base_dir, "data")
plots_dir <- file.path(base_dir, "results", "plots")
tables_dir <- file.path(base_dir, "results", "tables")

dat <- readRDS(file.path(data_dir, "02_normalized_data.rds"))
de_results <- read.csv(file.path(tables_dir, "03_DE_CKD_vs_Control.csv"))

# ── 1. Get feature annotation from GEO ──────────────────────────────────────

fdata <- fData(dat$eset)
cat("Feature data dimensions:", nrow(fdata), "x", ncol(fdata), "\n")
cat("Feature data columns:\n")
print(names(fdata))

# Identify symbol and Entrez columns
symbol_cols <- grep("symbol|gene_symbol|gene.symbol", names(fdata), ignore.case = TRUE, value = TRUE)
entrez_cols <- grep("^GENE$|entrez|geneid", names(fdata), ignore.case = TRUE, value = TRUE)
name_cols <- grep("gene_name|genename|description", names(fdata), ignore.case = TRUE, value = TRUE)

cat("\nSymbol columns:", symbol_cols, "\n")
cat("Entrez columns:", entrez_cols, "\n")
cat("Name columns:", name_cols, "\n")

# Show a sample of annotation data
cat("\nSample of feature data:\n")
print(head(fdata[, c(symbol_cols, entrez_cols, name_cols)[1:min(5, length(c(symbol_cols, entrez_cols, name_cols)))]], 10))

# ── 2. Build annotation lookup ──────────────────────────────────────────────

# Agilent GPL13497 typically has GENE_SYMBOL column
# Try to find the best symbol column
if (length(symbol_cols) > 0) {
  sym_col <- symbol_cols[1]
} else {
  # Check for 'Gene Symbol' or similar
  sym_col <- grep("Gene.*Symbol|GENE.*SYMBOL", names(fdata), value = TRUE)[1]
  if (is.na(sym_col)) {
    # Last resort: print all columns and first few rows
    cat("\nAll columns:\n")
    print(names(fdata))
    cat("\nFirst 3 rows:\n")
    print(head(fdata, 3))
    stop("Cannot find gene symbol column")
  }
}

cat("\nUsing symbol column:", sym_col, "\n")

annotation_df <- data.frame(
  ProbeID = rownames(fdata),
  gene_symbol = as.character(fdata[[sym_col]]),
  stringsAsFactors = FALSE
)

# Add Entrez if available
if (length(entrez_cols) > 0) {
  annotation_df$entrez_id <- as.character(fdata[[entrez_cols[1]]])
}

# Add gene name if available
if (length(name_cols) > 0) {
  annotation_df$gene_name <- as.character(fdata[[name_cols[1]]])
}

# Clean up annotations
annotation_df$gene_symbol[annotation_df$gene_symbol == ""] <- NA
if ("entrez_id" %in% names(annotation_df)) {
  annotation_df$entrez_id[annotation_df$entrez_id == ""] <- NA
}

cat("\nAnnotation coverage:\n")
cat("  Total probes:", nrow(annotation_df), "\n")
cat("  With gene symbol:", sum(!is.na(annotation_df$gene_symbol)), "\n")
if ("entrez_id" %in% names(annotation_df)) {
  cat("  With Entrez ID:", sum(!is.na(annotation_df$entrez_id)), "\n")
}

# Save full probe annotation
write.csv(annotation_df, file.path(tables_dir, "04_probe_annotation.csv"), row.names = FALSE)

# ── 3. Join with DE results ─────────────────────────────────────────────────

de_annotated <- de_results %>%
  left_join(annotation_df, by = "ProbeID")

cat("\nDE results annotation:\n")
cat("  Total probes:", nrow(de_annotated), "\n")
cat("  With gene symbol:", sum(!is.na(de_annotated$gene_symbol)), "\n")
cat("  Unique gene symbols:", n_distinct(de_annotated$gene_symbol, na.rm = TRUE), "\n")

write.csv(de_annotated, file.path(tables_dir, "04_DE_CKD_vs_Control_annotated.csv"), row.names = FALSE)

# ── 4. Collapse to unique genes ─────────────────────────────────────────────

# Remove probes mapping to multiple genes (contain "///")
de_annotated <- de_annotated %>%
  filter(!is.na(gene_symbol)) %>%
  filter(!grepl("///", gene_symbol))

# Keep most significant probe per gene
de_clean <- de_annotated %>%
  arrange(P.Value, -abs(logFC)) %>%
  filter(!duplicated(gene_symbol)) %>%
  arrange(adj.P.Val)

cat("\nAfter collapsing to unique genes:", nrow(de_clean), "\n")

# Re-calculate FDR after collapsing
de_clean$adj.P.Val <- p.adjust(de_clean$P.Value, method = "BH")

cat("\n=== Gene-Level DE Summary (CKD vs Control) ===\n")
for (fdr in c(0.05, 0.10, 0.25)) {
  sig <- sum(de_clean$adj.P.Val < fdr)
  up <- sum(de_clean$adj.P.Val < fdr & de_clean$logFC > 0)
  dn <- sum(de_clean$adj.P.Val < fdr & de_clean$logFC < 0)
  cat(sprintf("FDR < %.2f: %d DEGs (%d up, %d down)\n", fdr, sig, up, dn))
}

write.csv(de_clean, file.path(tables_dir, "04_DE_CKD_vs_Control_clean_symbol.csv"), row.names = FALSE)

# ── 5. Top DEGs table ───────────────────────────────────────────────────────

cat("\nTop 30 DEGs by FDR:\n")
top30 <- de_clean %>%
  head(30) %>%
  dplyr::select(gene_symbol, logFC, P.Value, adj.P.Val, AveExpr)
print(top30, digits = 3)

# Save up/down gene lists
sig_genes <- de_clean %>% filter(adj.P.Val < 0.05)
write.csv(sig_genes %>% filter(logFC > 0), file.path(tables_dir, "04_genes_upregulated.csv"), row.names = FALSE)
write.csv(sig_genes %>% filter(logFC < 0), file.path(tables_dir, "04_genes_downregulated.csv"), row.names = FALSE)

# ── 6. Annotated volcano plot ───────────────────────────────────────────────

de_clean$significance <- case_when(
  de_clean$adj.P.Val < 0.05 & de_clean$logFC > 0.5 ~ "Up",
  de_clean$adj.P.Val < 0.05 & de_clean$logFC < -0.5 ~ "Down",
  TRUE ~ "NS"
)

top_label <- de_clean %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(P.Value) %>%
  head(25)

sig_colors <- c(Up = "#E41A1C", Down = "#377EB8", NS = "grey70")

p_volcano <- ggplot(de_clean, aes(logFC, -log10(adj.P.Val), color = significance)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = sig_colors) +
  ggrepel::geom_text_repel(data = top_label, aes(label = gene_symbol),
                            size = 2.5, max.overlaps = 20, color = "black") +
  labs(title = "GSE98603: CKD vs Control — Subcutaneous WAT (Males)",
       subtitle = sprintf("%d DEGs at FDR<0.05 (%d up, %d down)",
                          sum(de_clean$adj.P.Val < 0.05),
                          sum(de_clean$adj.P.Val < 0.05 & de_clean$logFC > 0),
                          sum(de_clean$adj.P.Val < 0.05 & de_clean$logFC < 0)),
       x = "log2 Fold Change", y = "-log10(FDR)") +
  theme_minimal(base_size = 12)
ggsave(file.path(plots_dir, "04_volcano_annotated.png"), p_volcano, width = 10, height = 8, dpi = 300)

# ── 7. Entrez ID mapping via org.Hs.eg.db (always, for consistency) ──────────

# Always use org.Hs.eg.db rather than GPL annotation for Entrez IDs
# GPL IDs may be stale; org.Hs.eg.db is current
cat("\nMapping gene symbols to Entrez IDs via org.Hs.eg.db...\n")
library(org.Hs.eg.db)

symbol_to_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(de_clean$gene_symbol),
  columns = "ENTREZID",
  keytype = "SYMBOL"
)
# Remove 1-to-many mappings
symbol_to_entrez <- symbol_to_entrez %>%
  filter(!is.na(ENTREZID)) %>%
  filter(!duplicated(SYMBOL))

# Drop any existing entrez_id column from GPL, replace with org.Hs.eg.db
if ("entrez_id" %in% names(de_clean)) de_clean$entrez_id <- NULL

de_clean <- de_clean %>%
  left_join(symbol_to_entrez, by = c("gene_symbol" = "SYMBOL")) %>%
  dplyr::rename(entrez_id = ENTREZID)

cat("Mapped Entrez IDs:", sum(!is.na(de_clean$entrez_id)), "/", nrow(de_clean), "\n")

# Save final annotated file
write.csv(de_clean, file.path(tables_dir, "04_DE_CKD_vs_Control_clean_symbol.csv"), row.names = FALSE)

saveRDS(de_clean, file.path(data_dir, "04_de_clean.rds"))
cat("\n✓ Annotation complete\n")
