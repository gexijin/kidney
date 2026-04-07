#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(splines)
})

replicate_dirs <- local({
  root <- "kidney/manuscript/replicate"
  list(
    root = root,
    figures = file.path(root, "figures"),
    tables = file.path(root, "tables"),
    rds = file.path(root, "rds"),
    logs = file.path(root, "logs")
  )
})

dir.create(replicate_dirs$figures, recursive = TRUE, showWarnings = FALSE)
dir.create(replicate_dirs$tables, recursive = TRUE, showWarnings = FALSE)
dir.create(replicate_dirs$rds, recursive = TRUE, showWarnings = FALSE)
dir.create(replicate_dirs$logs, recursive = TRUE, showWarnings = FALSE)

base_covariates <- c(
  "AGE", "RACE", "BMI", "DTHHRDY", "ischemic_hrs",
  "SMRIN", "SMCENTER", "SMEXNCRT"
)

joint_factors <- c("MHABNWBC", "MHLVRDIS", "MHRNLFLR", "MHT2D")

model_vars <- c(base_covariates, joint_factors)

read_inputs <- function() {
  list(
    pheno = readRDS("data/phenotype/merged/merged_phenotype.rds"),
    counts_all = readRDS("data/reads/gene_reads_v11.rds")
  )
}

prepare_tissue <- function(pheno, counts_all, tissue_name, complete_vars = model_vars) {
  tb <- pheno %>% filter(SMTSD == tissue_name)
  sample_ids <- intersect(tb$SAMPID, colnames(counts_all))
  tb <- tb %>% filter(SAMPID %in% sample_ids)
  counts <- counts_all[, sample_ids, drop = FALSE]

  min_samples <- ceiling(0.20 * ncol(counts))
  keep <- rowSums(counts >= 10) >= min_samples
  counts <- counts[keep, , drop = FALSE]

  if (!is.null(complete_vars)) {
    tb <- tb %>% filter(complete.cases(across(all_of(complete_vars))))
    counts <- counts[, tb$SAMPID, drop = FALSE]
  }

  stopifnot(identical(tb$SAMPID, colnames(counts)))

  list(
    tb = tb,
    counts = counts,
    tissue_name = tissue_name
  )
}

fmt_median_iqr <- function(x) {
  x <- x[is.finite(x)]
  sprintf("%.1f [%.1f, %.1f]", median(x), quantile(x, 0.25), quantile(x, 0.75))
}

binary_smd <- function(x, g) {
  x <- as.numeric(x)
  g <- as.integer(g)
  p1 <- mean(x[g == 1], na.rm = TRUE)
  p0 <- mean(x[g == 0], na.rm = TRUE)
  p <- (p1 + p0) / 2
  if (!is.finite(p) || p %in% c(0, 1)) {
    return(NA_real_)
  }
  (p1 - p0) / sqrt(p * (1 - p))
}

numeric_smd <- function(x, g) {
  g <- as.integer(g)
  m1 <- mean(x[g == 1], na.rm = TRUE)
  m0 <- mean(x[g == 0], na.rm = TRUE)
  s1 <- stats::var(x[g == 1], na.rm = TRUE)
  s0 <- stats::var(x[g == 0], na.rm = TRUE)
  s_pool <- sqrt((s1 + s0) / 2)
  if (!is.finite(s_pool) || s_pool == 0) {
    return(NA_real_)
  }
  (m1 - m0) / s_pool
}

factor_level_smds <- function(x, g, prefix) {
  x <- as.factor(x)
  levels_x <- levels(x)
  map_dfr(levels_x, function(level) {
    tibble(
      variable = sprintf("%s=%s", prefix, level),
      smd = binary_smd(as.integer(x == level), g)
    )
  })
}

sex_ckd_counts <- function(tb) {
  tb %>%
    count(SEX, MHRNLFLR, name = "n") %>%
    complete(SEX, MHRNLFLR, fill = list(n = 0)) %>%
    arrange(SEX, MHRNLFLR)
}

fit_limma <- function(tb, counts, formula, coef_name) {
  design <- model.matrix(formula, data = tb)
  qr_obj <- qr(design)
  if (qr_obj$rank < ncol(design)) {
    design <- design[, qr_obj$pivot[seq_len(qr_obj$rank)], drop = FALSE]
  }
  stopifnot(coef_name %in% colnames(design))

  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  v <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none") %>%
    rownames_to_column("gene_id")

  list(
    design = design,
    fit = fit,
    voom = v,
    tt = tt
  )
}

write_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}
