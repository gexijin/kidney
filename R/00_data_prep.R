#!/usr/bin/env Rscript
# =============================================================================
# Stage 1, Script 00: Data Prep for Adipose Subcutaneous, Visceral, and Ovary
#
# Loads phenotype + protein-coding gene counts, applies expression filter.
# Saves one data object per tissue for downstream stage1 analyses.
#
# Row IDs are Ensembl gene IDs without version (e.g. ENSG00000141510).
#
# Output: explore/renal_fat/paper/results/stage1_data_{tissue_slug}.rds
#   List with: tissue_df, counts, tissue_name
# =============================================================================

library(tidyverse)

tissues <- c(
  "Adipose - Subcutaneous",
  "Adipose - Visceral (Omentum)",
  "Ovary"
)

out_dir <- "kidney/results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load shared resources ────────────────────────────────────────────────
pheno  <- readRDS("data/phenotype/merged/merged_phenotype.rds")
counts_all <- readRDS("data/reads/gene_reads_v11.rds")

cat("Full matrix:", nrow(counts_all), "genes x", ncol(counts_all), "samples\n")

# ── 2. Process each tissue ──────────────────────────────────────────────────
for (tissue_name in tissues) {
  tissue_slug <- tolower(gsub("[^A-Za-z0-9]+", "_", tissue_name))
  tissue_slug <- gsub("^_|_$", "", tissue_slug)

  cat("\n=== Data Prep:", tissue_name, "===\n")

  # Phenotype
  tissue_df <- pheno %>% filter(SMTSD == tissue_name)
  cat("Samples for", tissue_name, ":", nrow(tissue_df), "\n")

  # Subset counts to tissue samples
  tissue_samps <- intersect(tissue_df$SAMPID, colnames(counts_all))
  cat("Samples in counts:", length(tissue_samps), "\n")
  if (length(tissue_samps) == 0) {
    cat("SKIP: no matched samples\n")
    next
  }

  counts <- counts_all[, tissue_samps]
  tissue_df <- tissue_df %>% filter(SAMPID %in% tissue_samps)

  # Expression filter: ≥10 counts in ≥20% of samples
  min_samples <- ceiling(0.20 * ncol(counts))
  keep_expr   <- rowSums(counts >= 10) >= min_samples
  counts <- counts[keep_expr, ]

  # Save
  out <- list(tissue_df = tissue_df, counts = counts, tissue_name = tissue_name)

  out_file <- file.path(out_dir, paste0("stage1_data_", tissue_slug, ".rds"))
  saveRDS(out, out_file)

  cat("Samples:", nrow(tissue_df), "\n")
  cat("Genes:", nrow(counts), "\n")
  cat("Saved:", out_file, "\n")
}

cat("\nDone.\n")
