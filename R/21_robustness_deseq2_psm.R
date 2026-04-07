#!/usr/bin/env Rscript
# =============================================================================
# Stage 1b: DESeq2 Propensity-Matched Validation
#
# Independent method validation using DESeq2 on propensity-matched samples.
# Each sex run separately. Concordance assessed vs limma/voom from 01_stage1.R.
#
# Sections:
#   0 — Load data & prep
#   1 — DESeq2 full-model replication (method change only)
#   2 — Propensity matching (both sexes)
#   3 — Limma/voom on matched samples (sample change only)
#   4 — Repeated randomized matching robustness (10 iterations)
#   5 — Figures
#   6 — Summary
# =============================================================================

library(tidyverse)
library(DESeq2)
library(MatchIt)
library(limma)
library(edgeR)
library(parallel)

source("kidney/R/functions.R")

set.seed(42)
N_CORES <- 10

res_dir <- "kidney/results"
fig_dir <- "kidney/figures"


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 0: LOAD DATA & PREP                                            ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 0: Load Data & Prep\n")
cat(strrep("=", 70), "\n")

dat <- readRDS(file.path(res_dir, "prep_adipose_subcutaneous.rds"))
tb <- dat$tb
counts <- dat$counts
joint_factors <- dat$joint_factors

# Load limma results for concordance
res_female_limma <- readRDS(file.path(res_dir, "de_female.rds"))
res_male_limma   <- readRDS(file.path(res_dir, "de_male.rds"))

# Clean DTHHRDY: unordered factor, safe level names for DESeq2
tb$DTHHRDY <- factor(
  gsub("-", "_", as.character(tb$DTHHRDY)),
  ordered = FALSE)

# Ensure comorbidity factors are integer for matching and modeling
for (f in intersect(c("MHABNWBC", "MHLVRDIS", "MHCVD", "MHT2D", "MHHTN", "MHHRTDIS"),
                    names(tb)))
  tb[[f]] <- as.integer(tb[[f]])

# Center continuous variables for DESeq2 numerical stability
tb$AGE_c <- scale(tb$AGE, center = TRUE, scale = FALSE)[, 1]
tb$ischemic_hrs_c <- scale(tb$ischemic_hrs, center = TRUE, scale = FALSE)[, 1]
tb$SMRIN_c <- scale(tb$SMRIN, center = TRUE, scale = FALSE)[, 1]

cat(sprintf("Loaded: %d samples, %d genes\n", nrow(tb), nrow(counts)))
cat(sprintf("Female: %d (CKD=%d) | Male: %d (CKD=%d)\n",
            sum(tb$SEX == "Female"), sum(tb$SEX == "Female" & tb$MHRNLFLR == 1),
            sum(tb$SEX == "Male"), sum(tb$SEX == "Male" & tb$MHRNLFLR == 1)))
cat("Joint factors:", paste(joint_factors, collapse = ", "), "\n")

# Comorbidity covariates for doubly-robust model (exclude MHRNLFLR = treatment)
comorbidity_covars <- setdiff(joint_factors, "MHRNLFLR")
cat("Comorbidity covariates for DR model:", paste(comorbidity_covars, collapse = ", "), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1: DESeq2 FULL-MODEL REPLICATION                               ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 1: DESeq2 Full-Model Replication (same covariates as limma)\n")
cat(strrep("=", 70), "\n")

# Same formula as limma in script 01, but with centered continuous vars
# Limma: ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + {joint_factors}
# DESeq2: same, with centered AGE/ischemic_hrs/SMRIN

# Clean SMCENTER level names for DESeq2
tb$SMCENTER <- factor(gsub("-", "_", as.character(tb$SMCENTER)), ordered = FALSE)

full_model_terms <- c("AGE_c", "RACE", "BMI", "DTHHRDY", "ischemic_hrs_c", "SMRIN_c",
                       "SMCENTER", "SMEXNCRT", comorbidity_covars, "MHRNLFLR")
full_formula <- as.formula(paste("~", paste(full_model_terms, collapse = " + ")))
cat("Formula:", deparse(full_formula), "\n")

run_deseq2_full <- function(tb_sex, counts, formula, sex_label) {
  cat(sprintf("\n--- %s ---\n", sex_label))
  cat(sprintf("  Samples: %d | CKD: %d | Controls: %d\n",
              nrow(tb_sex), sum(tb_sex$MHRNLFLR == 1), sum(tb_sex$MHRNLFLR == 0)))

  tb_ds <- tb_sex
  tb_ds$MHRNLFLR <- factor(tb_ds$MHRNLFLR, levels = c(0, 1))

  # Drop unused factor levels (e.g. SMCENTER, DTHHRDY levels absent in stratum)
  tb_ds <- droplevels(tb_ds)
  counts_sub <- counts[, tb_ds$SAMPID]

  dds <- DESeqDataSetFromMatrix(counts_sub, tb_ds, formula)
  dds <- DESeq(dds, quiet = TRUE)

  coef_name <- grep("MHRNLFLR", resultsNames(dds), value = TRUE)[1]
  cat(sprintf("  CKD coefficient: %s\n", coef_name))
  cat(sprintf("  Design columns: %d\n", ncol(model.matrix(formula, data = tb_ds))))

  res <- as.data.frame(results(dds, name = coef_name))
  res$gene_id <- rownames(res)

  n_deg <- sum(res$padj < 0.05, na.rm = TRUE)
  n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("  DEGs (FDR<0.05): %d (%d up, %d down)\n", n_deg, n_up, n_down))

  list(results = res, n_deg = n_deg, n_up = n_up, n_down = n_down,
       n_samples = nrow(tb_ds), n_cases = sum(tb_ds$MHRNLFLR == 1))
}

tb_f <- tb %>% filter(SEX == "Female")
tb_m <- tb %>% filter(SEX == "Male")

res_f_full <- run_deseq2_full(tb_f, counts, full_formula, "Female")
res_m_full <- run_deseq2_full(tb_m, counts, full_formula, "Male")

# Concordance: DESeq2 full vs limma full (same model, different method)
cat("\n--- Concordance: DESeq2-full vs Limma-full ---\n")

for (sex_info in list(
  list(limma = res_female_limma$tt, deseq = res_f_full$results, label = "Female"),
  list(limma = res_male_limma$tt, deseq = res_m_full$results, label = "Male")
)) {
  limma_df <- sex_info$limma %>%
    dplyr::select(gene_id, logFC_limma = logFC, t_limma = t, padj_limma = adj.P.Val)
  deseq_df <- sex_info$deseq %>%
    dplyr::select(gene_id, logFC_deseq = log2FoldChange, stat_deseq = stat, padj_deseq = padj) %>%
    filter(!is.na(stat_deseq))
  merged <- inner_join(limma_df, deseq_df, by = "gene_id")

  rho_lfc <- cor(merged$logFC_limma, merged$logFC_deseq, method = "spearman")
  rho_stat <- cor(merged$t_limma, merged$stat_deseq, method = "spearman")
  sig_l <- sum(merged$padj_limma < 0.05)
  sig_d <- sum(merged$padj_deseq < 0.05)
  overlap <- sum(merged$padj_limma < 0.05 & merged$padj_deseq < 0.05)

  if (overlap > 0) {
    shared <- merged %>% filter(padj_limma < 0.05 & padj_deseq < 0.05)
    dir_conc <- mean(sign(shared$logFC_limma) == sign(shared$logFC_deseq))
  } else dir_conc <- NA

  cat(sprintf("  %s: logFC rho=%.3f, stat rho=%.3f | Limma=%d, DESeq2=%d, overlap=%d, dir=%.1f%%\n",
              sex_info$label, rho_lfc, rho_stat, sig_l, sig_d, overlap,
              100 * ifelse(is.na(dir_conc), 0, dir_conc)))
}

saveRDS(list(female = res_f_full, male = res_m_full),
        file.path(res_dir, "de_deseq2_full.rds"))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2: PROPENSITY MATCHING                                         ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 2: Propensity Matching\n")
cat(strrep("=", 70), "\n")

# Match on agonal confounders (the primary threat to validity)
# SMCENTER excluded — too many levels relative to n_cases
# Exact matching on isch_bin forces CKD/control into same ischemic time category
psm_formula <- MHRNLFLR ~ AGE + RACE + DTHHRDY + ischemic_hrs + SMRIN

run_matching <- function(tb_sex, sex_label, seed = 42) {
  set.seed(seed)
  # Shuffle row order so nearest-neighbor tie-breaking varies
  tb_sex <- tb_sex[sample(nrow(tb_sex)), ]

  # 3-category ischemic time bin for exact matching
  tb_sex <- add_ischemic_bin(tb_sex)

  cat(sprintf("\n  Ischemic time bins (CKD): %s\n",
              paste(names(table(tb_sex$isch_bin[tb_sex$MHRNLFLR == 1])),
                    table(tb_sex$isch_bin[tb_sex$MHRNLFLR == 1]),
                    sep = "=", collapse = ", ")))
  cat(sprintf("  Ischemic time bins (Ctrl): %s\n",
              paste(names(table(tb_sex$isch_bin[tb_sex$MHRNLFLR == 0])),
                    table(tb_sex$isch_bin[tb_sex$MHRNLFLR == 0]),
                    sep = "=", collapse = ", ")))

  m <- matchit(psm_formula, data = tb_sex, method = "nearest",
               ratio = 3, caliper = 0.25, std.caliper = TRUE,
               m.order = "random", exact = ~isch_bin)

  n_matched <- sum(m$treat == 1 & !is.na(m$subclass))
  n_unmatched <- sum(m$treat == 1 & is.na(m$subclass))
  tb_matched <- match.data(m)

  cat(sprintf("\n--- %s ---\n", sex_label))
  cat(sprintf("  Matched: %d/%d CKD cases → %d total\n",
              n_matched, sum(m$treat == 1), nrow(tb_matched)))
  if (n_unmatched > 0)
    cat(sprintf("  WARNING: %d cases unmatched (caliper too strict)\n", n_unmatched))

  # Print balance
  bal <- summary(m)
  bal_tab <- as.data.frame(bal$sum.matched)
  for (v in rownames(bal_tab)) {
    smd <- bal_tab[v, "Std. Mean Diff."]
    cat(sprintf("  %-25s SMD = %+.3f %s\n", v, smd,
                ifelse(abs(smd) > 0.1, " *", "")))
  }

  # Ischemic time balance after matching
  cat("  Ischemic time balance (exact matched on isch_bin):\n")
  cat(sprintf("    CKD median: %.1f hrs | Ctrl median: %.1f hrs\n",
              median(tb_matched$ischemic_hrs[tb_matched$MHRNLFLR == 1]),
              median(tb_matched$ischemic_hrs[tb_matched$MHRNLFLR == 0])))
  cat(sprintf("    CKD isch_bin: %s\n",
              paste(names(table(tb_matched$isch_bin[tb_matched$MHRNLFLR == 1])),
                    table(tb_matched$isch_bin[tb_matched$MHRNLFLR == 1]),
                    sep = "=", collapse = ", ")))
  cat(sprintf("    Ctrl isch_bin: %s\n",
              paste(names(table(tb_matched$isch_bin[tb_matched$MHRNLFLR == 0])),
                    table(tb_matched$isch_bin[tb_matched$MHRNLFLR == 0]),
                    sep = "=", collapse = ", ")))

  # Check comorbidity balance (not matched on, but report)
  cat("  Comorbidity balance (not matched on):\n")
  for (cv in comorbidity_covars) {
    p_ckd <- mean(tb_matched[[cv]][tb_matched$MHRNLFLR == 1] == 1, na.rm = TRUE)
    p_ctrl <- mean(tb_matched[[cv]][tb_matched$MHRNLFLR == 0] == 1, na.rm = TRUE)
    cat(sprintf("    %-12s CKD=%.0f%% Ctrl=%.0f%%\n", cv, 100*p_ckd, 100*p_ctrl))
  }

  list(match_obj = m, matched_data = tb_matched, n_unmatched = n_unmatched)
}

tb_f <- tb %>% filter(SEX == "Female")
tb_m <- tb %>% filter(SEX == "Male")

psm_f <- run_matching(tb_f, "Female")
psm_m <- run_matching(tb_m, "Male")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 3: LIMMA/VOOM ON MATCHED SAMPLES                               ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 3: Limma/voom on Propensity-Matched Samples\n")
cat(strrep("=", 70), "\n")

# Same method as script 01, different (balanced) samples
# Doubly robust: agonal covariates + comorbidities (matching handled balance,
# but we adjust anyway for residual imbalance)
dr_limma_terms <- c("AGE", "BMI", "DTHHRDY", "ischemic_hrs", "SMRIN", "SMCENTER",
                     "SMEXNCRT", comorbidity_covars, "MHRNLFLR")
dr_limma_formula <- as.formula(paste("~", paste(dr_limma_terms, collapse = " + ")))
cat("DR formula:", deparse(dr_limma_formula), "\n")

run_limma_psm <- function(tb_matched, counts, dr_formula, limma_ref_tt, label) {
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("  Samples: %d | CKD: %d | Controls: %d\n",
              nrow(tb_matched), sum(tb_matched$MHRNLFLR == 1),
              sum(tb_matched$MHRNLFLR == 0)))

  limma_fit <- tryCatch(
    run_limma_analysis(
      tb_matched, counts, dr_formula, "MHRNLFLR",
      coef_match = "prefix",
      context = sprintf("PSM limma model (%s)", label)
    ),
    error = function(e) NULL
  )
  if (is.null(limma_fit)) { cat("  MHRNLFLR not in model!\n"); return(NULL) }

  tt <- limma_fit$tt
  cat(sprintf("  Design: %d columns\n", ncol(limma_fit$design)))
  cat(sprintf("  DEGs (FDR<0.05): %d (%d up, %d down)\n",
              limma_fit$n_deg, limma_fit$n_up, limma_fit$n_down))

  conc <- compare_de_results(limma_ref_tt, tt)
  merged <- conc$merged %>%
    dplyr::rename(logFC_psm = logFC_new, t_psm = t_new, padj_psm = padj_new)

  cat(sprintf("  vs full-sample limma: logFC rho=%.3f, t-stat rho=%.3f\n",
              conc$rho_lfc, conc$rho_t))
  cat(sprintf("  Full-sample DEGs: %d | PSM DEGs: %d | Overlap: %d | Jaccard: %.3f\n",
              conc$n_deg_ref, conc$n_deg_new, conc$overlap, conc$jaccard))

  if (conc$overlap > 0) {
    shared <- merged %>% filter(gene_id %in% conc$overlap_genes)
    concordant <- sum(sign(shared$logFC_ref) == sign(shared$logFC_psm))
    cat(sprintf("  Direction concordance: %d/%d (%.1f%%)\n",
                concordant, conc$overlap, 100 * concordant / conc$overlap))
  }

  list(tt = tt, n_deg = limma_fit$n_deg, n_up = limma_fit$n_up, n_down = limma_fit$n_down,
       n_samples = nrow(tb_matched), n_cases = sum(tb_matched$MHRNLFLR == 1),
       rho_lfc = conc$rho_lfc, rho_t = conc$rho_t, jaccard = conc$jaccard,
       n_overlap = conc$overlap, merged = merged)
}

res_f_psm <- run_limma_psm(psm_f$matched_data, counts, dr_limma_formula,
                            res_female_limma$tt, "Female PSM + Limma")
res_m_psm <- run_limma_psm(psm_m$matched_data, counts, dr_limma_formula,
                            res_male_limma$tt, "Male PSM + Limma")

saveRDS(list(res = res_f_psm, match_obj = psm_f$match_obj,
             matched_data = psm_f$matched_data),
        file.path(res_dir, "de_female_psm.rds"))
saveRDS(list(res = res_m_psm, match_obj = psm_m$match_obj,
             matched_data = psm_m$matched_data),
        file.path(res_dir, "de_male_psm.rds"))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 4: REPEATED RANDOMIZED MATCHING (LIMMA, 50 ITERATIONS)         ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 4: Repeated Randomized Matching (10 iterations)\n")
cat(strrep("=", 70), "\n")

# Use limma/voom for speed. Each iteration:
#   1. Shuffle row order → matchit (nearest, 1:3, caliper 0.25)
#   2. Fit limma/voom with doubly-robust formula on matched set
#   3. Record n_DEGs, logFC correlation with script-01 full-sample limma

N_ITER <- 10

# Limma doubly-robust formula (mirrors DESeq2 but uses original scale vars)
dr_limma_terms <- c("AGE", "BMI", "DTHHRDY", "ischemic_hrs", "SMRIN", "SMCENTER",
                     "SMEXNCRT", comorbidity_covars, "MHRNLFLR")
dr_limma_formula <- as.formula(paste("~", paste(dr_limma_terms, collapse = " + ")))

run_one_iteration <- function(i, tb_sex, counts, limma_ref_tt, dr_formula, psm_formula) {
  # Shuffle and match
  set.seed(i * 1000 + 7)
  tb_shuf <- tb_sex[sample(nrow(tb_sex)), ]

  # Ensure isch_bin exists for exact matching
  tb_shuf <- add_ischemic_bin(tb_shuf)

  m <- tryCatch(
    matchit(psm_formula, data = tb_shuf, method = "nearest",
            ratio = 3, caliper = 0.25, std.caliper = TRUE,
            m.order = "random", exact = ~isch_bin),
    error = function(e) NULL)
  if (is.null(m)) return(list(n_deg = NA, rho = NA, n_matched = 0))

  n_matched_cases <- sum(m$treat == 1 & !is.na(m$subclass))
  if (n_matched_cases < 10) return(list(n_deg = NA, rho = NA, n_matched = n_matched_cases))

  tb_m <- match.data(m)
  iter_fit <- tryCatch(
    run_limma_analysis(
      tb_m, counts, dr_formula, "MHRNLFLR",
      coef_match = "prefix",
      context = sprintf("PSM iteration %d", i)
    ),
    error = function(e) NULL
  )
  if (is.null(iter_fit)) return(list(n_deg = NA, rho = NA, n_matched = n_matched_cases))

  conc <- compare_de_results(
    limma_ref_tt,
    iter_fit$tt %>% dplyr::select(gene_id, logFC, t, adj.P.Val)
  )

  list(n_deg = iter_fit$n_deg, rho = conc$rho_lfc, n_matched = nrow(tb_m))
}

cat("\n--- Female (10 iterations) ---\n")
iter_f <- lapply(seq_len(N_ITER), function(i) {
  run_one_iteration(i, tb_f, counts, res_female_limma$tt,
                    dr_limma_formula, psm_formula)
})
iter_f_df <- tibble(
  iteration = seq_len(N_ITER),
  n_deg = sapply(iter_f, `[[`, "n_deg"),
  rho = sapply(iter_f, `[[`, "rho"),
  n_matched = sapply(iter_f, `[[`, "n_matched")
)
valid_f <- iter_f_df %>% filter(!is.na(n_deg))
cat(sprintf("  Valid iterations: %d/%d\n", nrow(valid_f), N_ITER))
cat(sprintf("  DEGs: median=%d, IQR=[%d, %d], range=[%d, %d]\n",
            as.integer(median(valid_f$n_deg)),
            as.integer(quantile(valid_f$n_deg, 0.25)),
            as.integer(quantile(valid_f$n_deg, 0.75)),
            min(valid_f$n_deg), max(valid_f$n_deg)))
cat(sprintf("  logFC rho vs full-sample: median=%.3f, IQR=[%.3f, %.3f]\n",
            median(valid_f$rho),
            quantile(valid_f$rho, 0.25), quantile(valid_f$rho, 0.75)))

cat("\n--- Male (10 iterations) ---\n")
iter_m <- lapply(seq_len(N_ITER), function(i) {
  run_one_iteration(i, tb_m, counts, res_male_limma$tt,
                    dr_limma_formula, psm_formula)
})
iter_m_df <- tibble(
  iteration = seq_len(N_ITER),
  n_deg = sapply(iter_m, `[[`, "n_deg"),
  rho = sapply(iter_m, `[[`, "rho"),
  n_matched = sapply(iter_m, `[[`, "n_matched")
)
valid_m <- iter_m_df %>% filter(!is.na(n_deg))
cat(sprintf("  Valid iterations: %d/%d\n", nrow(valid_m), N_ITER))
cat(sprintf("  DEGs: median=%d, IQR=[%d, %d], range=[%d, %d]\n",
            as.integer(median(valid_m$n_deg)),
            as.integer(quantile(valid_m$n_deg, 0.25)),
            as.integer(quantile(valid_m$n_deg, 0.75)),
            min(valid_m$n_deg), max(valid_m$n_deg)))
cat(sprintf("  logFC rho vs full-sample: median=%.3f, IQR=[%.3f, %.3f]\n",
            median(valid_m$rho),
            quantile(valid_m$rho, 0.25), quantile(valid_m$rho, 0.75)))

write_csv(iter_f_df, file.path(res_dir, "psm_repeated_female.csv"))
write_csv(iter_m_df, file.path(res_dir, "psm_repeated_male.csv"))


conc_df <- tibble(
  sex = c("Female", "Male"),
  full_degs = c(res_female_limma$n_deg, res_male_limma$n_deg),
  psm_degs = c(res_f_psm$n_deg, res_m_psm$n_deg),
  overlap = c(res_f_psm$n_overlap, res_m_psm$n_overlap),
  jaccard = c(res_f_psm$jaccard, res_m_psm$jaccard),
  logFC_rho = c(res_f_psm$rho_lfc, res_m_psm$rho_lfc),
  t_stat_rho = c(res_f_psm$rho_t, res_m_psm$rho_t)
)
write_csv(conc_df, file.path(res_dir, "psm_concordance.csv"))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 5: FIGURES                                                     ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 5: Figures\n")
cat(strrep("=", 70), "\n")

# ── Figure: Covariate balance (love plot) ──────────────────────────────────

cat("\nGenerating balance plots...\n")

make_balance_df <- function(match_obj, label) {
  s <- summary(match_obj)
  before <- as.data.frame(s$sum.all) %>%
    mutate(variable = rownames(.), timing = "Before", sex = label)
  after <- as.data.frame(s$sum.matched) %>%
    mutate(variable = rownames(.), timing = "After", sex = label)
  bind_rows(before, after) %>%
    dplyr::select(variable, timing, sex, smd = `Std. Mean Diff.`)
}

bal_df <- bind_rows(
  make_balance_df(psm_f$match_obj, "Female"),
  make_balance_df(psm_m$match_obj, "Male")
) %>%
  mutate(timing = factor(timing, levels = c("Before", "After")))

p_balance <- ggplot(bal_df, aes(x = smd, y = reorder(variable, abs(smd)),
                                 shape = timing, color = timing)) +
  geom_point(size = 3) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, color = "black") +
  scale_color_manual(values = c(Before = "firebrick", After = "steelblue")) +
  facet_wrap(~ sex) +
  labs(x = "Standardized Mean Difference", y = NULL,
       title = "Covariate Balance: Before vs After Propensity Matching",
       color = NULL, shape = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "21_psm_balance.png"), p_balance,
       width = 10, height = 6, dpi = 300)
cat("Saved:", file.path(fig_dir, "21_psm_balance.png"), "\n")


# ── Figure: logFC concordance scatter ──────────────────────────────────────

cat("Generating concordance scatter...\n")

plot_concordance <- function(merged, rho, sex_label) {
  plot_df <- merged %>%
    mutate(sig = case_when(
      padj_ref < 0.05 & padj_psm < 0.05 ~ "Both",
      padj_ref < 0.05 ~ "Full only",
      padj_psm < 0.05 ~ "PSM only",
      TRUE ~ "Neither"
    ))

  ggplot(plot_df, aes(x = logFC_ref, y = logFC_psm, color = sig)) +
    geom_point(data = plot_df %>% filter(sig == "Neither"), size = 0.3, alpha = 0.2) +
    geom_point(data = plot_df %>% filter(sig != "Neither"), size = 0.6, alpha = 0.6) +
    scale_color_manual(values = c(Both = "purple", `Full only` = "red",
                                  `PSM only` = "blue", Neither = "grey80"),
                       name = "Significant in") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = "logFC (Limma, full sample)",
         y = "logFC (Limma, propensity-matched)",
         title = sprintf("%s: Full vs PSM (rho = %.3f)", sex_label, rho)) +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}

p_conc_f <- plot_concordance(res_f_psm$merged, res_f_psm$rho_lfc, "Female")
p_conc_m <- plot_concordance(res_m_psm$merged, res_m_psm$rho_lfc, "Male")

p_conc <- cowplot::plot_grid(p_conc_f, p_conc_m, ncol = 2)
ggsave(file.path(fig_dir, "21_psm_concordance.png"), p_conc,
       width = 14, height = 6, dpi = 300)
cat("Saved:", file.path(fig_dir, "21_psm_concordance.png"), "\n")


# ── Figure: Repeated matching robustness ───────────────────────────────────

cat("Generating repeated matching plots...\n")

rep_df <- bind_rows(
  valid_f %>% mutate(sex = "Female"),
  valid_m %>% mutate(sex = "Male")
)

p_rep_deg <- ggplot(rep_df, aes(x = sex, y = n_deg, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.5) +
  scale_fill_manual(values = c(Female = "firebrick", Male = "steelblue"), guide = "none") +
  labs(x = NULL, y = "DEGs (FDR < 0.05)",
       title = "Repeated Randomized Matching (10 iterations, limma/voom DR)") +
  theme_bw(base_size = 12)

p_rep_rho <- ggplot(rep_df, aes(x = sex, y = rho, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.5) +
  scale_fill_manual(values = c(Female = "firebrick", Male = "steelblue"), guide = "none") +
  labs(x = NULL, y = "Spearman rho vs full-sample logFC",
       title = "Effect Size Concordance Across Match Sets") +
  theme_bw(base_size = 12)

p_rep <- cowplot::plot_grid(p_rep_deg, p_rep_rho, ncol = 2)
ggsave(file.path(fig_dir, "21_psm_robustness.png"), p_rep,
       width = 10, height = 5, dpi = 300)
cat("Saved:", file.path(fig_dir, "21_psm_robustness.png"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 6: SUMMARY                                                     ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("STAGE 1b SUMMARY\n")
cat(strrep("=", 70), "\n")
cat(sprintf("
1. Method change (DESeq2 full model, same covariates as limma):
   Female: %d DEGs (limma: %d) | logFC rho vs limma: see Section 1
   Male:   %d DEGs (limma: %d)

2. Sample change (PSM + limma, doubly-robust):
   Female: %d/%d matched → %d DEGs | logFC rho=%.3f, t rho=%.3f, Jaccard=%.3f
   Male:   %d/%d matched → %d DEGs | logFC rho=%.3f, t rho=%.3f, Jaccard=%.3f

3. Repeated randomized matching (10 iterations, limma DR):
   Female: median=%d DEGs, IQR=[%d, %d], rho=%.3f
   Male:   median=%d DEGs, IQR=[%d, %d], rho=%.3f
",
  res_f_full$n_deg, res_female_limma$n_deg,
  res_m_full$n_deg, res_male_limma$n_deg,
  res_f_psm$n_cases, sum(tb_f$MHRNLFLR == 1), res_f_psm$n_deg,
  res_f_psm$rho_lfc, res_f_psm$rho_t, res_f_psm$jaccard,
  res_m_psm$n_cases, sum(tb_m$MHRNLFLR == 1), res_m_psm$n_deg,
  res_m_psm$rho_lfc, res_m_psm$rho_t, res_m_psm$jaccard,
  as.integer(median(valid_f$n_deg)),
  as.integer(quantile(valid_f$n_deg, 0.25)),
  as.integer(quantile(valid_f$n_deg, 0.75)),
  median(valid_f$rho),
  as.integer(median(valid_m$n_deg)),
  as.integer(quantile(valid_m$n_deg, 0.25)),
  as.integer(quantile(valid_m$n_deg, 0.75)),
  median(valid_m$rho)
))

cat("Done.\n")
