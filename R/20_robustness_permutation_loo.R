#!/usr/bin/env Rscript
# =============================================================================
# 20: Robustness — Permutation, Leave-One-Out, Comorbidity, False-Positive
#
# Loads discovery DE results from 10_discovery_de.R and runs four validation
# tests on the female CKD signal.
# =============================================================================

library(tidyverse)
library(edgeR)
library(limma)
library(parallel)
library(msigdbr)
library(fgsea)

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
res_male <- readRDS(file.path(res_dir, "de_male.rds"))

# Female subset objects
tb_f <- tb %>% filter(SEX == "Female")
counts_f <- counts[, tb_f$SAMPID]

form_f <- as.formula(paste(
  "~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
  paste(joint_factors, collapse = " + ")))

design_f <- build_design(form_f, tb_f)

dge_f <- calcNormFactors(DGEList(counts = counts_f), method = "TMM")

n_obs_f <- res_female$n_deg
cat(sprintf("Observed female DEGs: %d\n", n_obs_f))

# Hallmark gene sets and primary GSEA for pathway correlation
hallmark_list <- load_hallmark_sets(filter_empty = TRUE)
primary_gsea <- run_gsea(res_female$tt, hallmark_list, "primary", seed = 42)

# ── 1c-i. Permutation test (100x, 10 cores) ──────────────────────────────────
cat("\n--- 1c-i. Permutation test (100 permutations, female stratum) ---\n")

N_PERM <- 100

# Pre-compute voom on base design to avoid fork segfaults with mclapply
v_f_perm <- voom(dge_f, design_f, plot = FALSE)

perm_degs <- mclapply(seq_len(N_PERM), function(i) {
  tb_perm <- tb_f
  tb_perm$MHRNLFLR <- sample(tb_perm$MHRNLFLR)

  d_perm <- build_design(form_f, tb_perm)

  fit_p <- eBayes(lmFit(v_f_perm, d_perm))

  fc <- tryCatch(
    get_required_coef_prefix(d_perm, "MHRNLFLR", "female permutation model"),
    error = function(e) NULL
  )
  if (is.null(fc)) return(0)

  tt_p <- topTable(fit_p, coef = fc, number = Inf, sort.by = "none")
  sum(tt_p$adj.P.Val < 0.05)
}, mc.cores = N_CORES)

perm_degs <- unlist(perm_degs)

perm_p <- (sum(perm_degs >= n_obs_f) + 1) / (N_PERM + 1)
cat(sprintf("Null: median=%d, max=%d, mean=%.1f\n",
            median(perm_degs), max(perm_degs), mean(perm_degs)))
cat(sprintf("Observed: %d | Permutation p = %.4f\n", n_obs_f, perm_p))

write_csv(data.frame(perm_degs = perm_degs),
          file.path(res_dir, "permutation_null.csv"))

# Pathway correlation for permutations (sequential — fgsea not fork-safe)
cat("Computing permutation pathway correlations...\n")
set.seed(123)
perm_pathway_rhos <- numeric(N_PERM)
for (i in seq_len(N_PERM)) {
  tb_perm <- tb_f
  tb_perm$MHRNLFLR <- sample(tb_perm$MHRNLFLR)
  d_perm <- build_design(form_f, tb_perm)
  fit_p <- eBayes(lmFit(v_f_perm, d_perm))
  fc <- tryCatch(
    get_required_coef_prefix(d_perm, "MHRNLFLR", "perm pathway_cor"),
    error = function(e) NULL
  )
  if (is.null(fc)) { perm_pathway_rhos[i] <- NA_real_; next }
  tt_p <- topTable(fit_p, coef = fc, number = Inf, sort.by = "none")
  tt_p$gene_id <- rownames(tt_p)
  perm_pathway_rhos[i] <- pathway_cor(tt_p, primary_gsea, hallmark_list)
  if (i %% 20 == 0) cat("  ", i, "/", N_PERM, "\n")
}
cat(sprintf("Permutation pathway rho: median=%.3f, mean=%.3f, range=%.3f--%.3f\n",
            median(perm_pathway_rhos, na.rm = TRUE),
            mean(perm_pathway_rhos, na.rm = TRUE),
            min(perm_pathway_rhos, na.rm = TRUE),
            max(perm_pathway_rhos, na.rm = TRUE)))

# ── 1c-ii. Leave-one-out (26 female CKD cases) ───────────────────────────────
cat("\n--- 1c-ii. Leave-one-out (female CKD cases) ---\n")

female_case_ids <- tb_f %>% filter(MHRNLFLR == 1) %>% pull(SAMPID)
cat("Female CKD cases:", length(female_case_ids), "\n")

loo_results <- data.frame(SAMPID = female_case_ids,
                          n_deg = NA_integer_,
                          pathway_rho = NA_real_,
                          stringsAsFactors = FALSE)

for (li in seq_along(female_case_ids)) {
  drop_id <- female_case_ids[li]
  tb_loo <- tb_f %>% filter(SAMPID != drop_id)
  counts_loo <- counts[, tb_loo$SAMPID]

  loo_fit <- tryCatch(
    run_limma_analysis(
      tb_loo, counts, form_f, "MHRNLFLR",
      coef_match = "prefix",
      context = sprintf("female leave-one-out model (%s)", drop_id)
    ),
    error = function(e) NULL
  )
  if (is.null(loo_fit)) next

  loo_results$n_deg[li] <- loo_fit$n_deg
  loo_results$pathway_rho[li] <- pathway_cor(
    loo_fit$tt, primary_gsea, hallmark_list)

  rm(loo_fit, counts_loo); gc(verbose = FALSE)
  if (li %% 5 == 0) cat("  ", li, "/", length(female_case_ids), "\n")
}

loo_degs <- loo_results$n_deg
loo_results$pct_change <- 100 * (n_obs_f - loo_degs) / n_obs_f

cat(sprintf("\nLOO results: range %d–%d, median %d\n",
            min(loo_degs, na.rm = TRUE),
            max(loo_degs, na.rm = TRUE),
            as.integer(median(loo_degs, na.rm = TRUE))))

influential <- which(abs(loo_results$pct_change) > 20)
if (length(influential) > 0) {
  cat("Influential donors (>20% change):\n")
  for (idx in influential)
    cat(sprintf("  %s: %d DEGs (%.1f%% change)\n",
                female_case_ids[idx], loo_degs[idx], loo_results$pct_change[idx]))
} else {
  cat("No single donor changes DEG count by >20%% — signal is stable.\n")
}

cat(sprintf("LOO pathway rho: median=%.3f, range=%.3f--%.3f\n",
            median(loo_results$pathway_rho, na.rm = TRUE),
            min(loo_results$pathway_rho, na.rm = TRUE),
            max(loo_results$pathway_rho, na.rm = TRUE)))

write_csv(loo_results, file.path(res_dir, "loo_results.csv"))

# ── 1c-iii. Comorbidity adjustment ───────────────────────────────────────────
cat("\n--- 1c-iii. Comorbidity adjustment (female stratum) ---\n")

form_f_str <- paste(
  "~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
  paste(joint_factors, collapse = " + "))

comorbidity_models <- list(base = n_obs_f)

# Note: MHT2D is already in joint_factors (base model), so the comorbidity
# test adds MHHTN + MHHRTDIS on top of the base model

# + MHHTN + MHHRTDIS
extra_cm <- setdiff(c("MHHTN", "MHHRTDIS"), joint_factors)
tb_f$MHHTN    <- as.integer(tb_f$MHHTN)
tb_f$MHHRTDIS <- as.integer(tb_f$MHHRTDIS)

if (length(extra_cm) > 0) {
  form_cm <- as.formula(paste(form_f_str, "+", paste(extra_cm, collapse = " + ")))
} else {
  form_cm <- as.formula(form_f_str)
}
tb_fcm <- tb_f[complete.cases(tb_f[, c("MHHTN", "MHHRTDIS")]), ]
counts_fcm <- counts[, tb_fcm$SAMPID]
d_cm <- build_design(form_cm, tb_fcm)
dge_fcm <- DGEList(counts = counts_fcm); dge_fcm <- calcNormFactors(dge_fcm, method = "TMM")
v_cm <- voom(dge_fcm, d_cm, plot = FALSE)
fit_cm <- eBayes(lmFit(v_cm, d_cm))
fc_cm <- get_required_coef_prefix(d_cm, "MHRNLFLR", "female comorbidity model")
n_cm <- sum(topTable(fit_cm, coef = fc_cm, number = Inf, sort.by = "none")$adj.P.Val < 0.05)
cat(sprintf("  + HTN + HRTDIS (MHT2D already in base): %d DEGs (%.1f%% retained)\n",
            n_cm, 100 * n_cm / n_obs_f))
comorbidity_models[["plus_HTN_HRTDIS"]] <- n_cm

# ── 1c-iv. False-positive calibration (20 random phenotypes) ──────────────────
cat("\n--- 1c-iv. False-positive calibration (20 random phenotypes, female stratum) ---\n")

ckd_prev_f <- mean(tb_f$MHRNLFLR == 1)
v_f_base <- voom(dge_f, design_f, plot = FALSE)

N_RANDOM <- 20
random_degs <- mclapply(seq_len(N_RANDOM), function(i) {
  tb_rand <- tb_f
  n_cases <- round(nrow(tb_rand) * ckd_prev_f)
  tb_rand$MHRNLFLR <- 0
  tb_rand$MHRNLFLR[sample(nrow(tb_rand), n_cases)] <- 1

  d_rand <- build_design(form_f, tb_rand)

  fit_r <- eBayes(lmFit(v_f_base, d_rand))
  fc <- tryCatch(
    get_required_coef_prefix(d_rand, "MHRNLFLR", "female random-phenotype model"),
    error = function(e) NULL
  )
  if (is.null(fc)) return(0)

  tt_r <- topTable(fit_r, coef = fc, number = Inf, sort.by = "none")
  sum(tt_r$adj.P.Val < 0.05)
}, mc.cores = N_CORES)

random_degs <- unlist(random_degs)
cat(sprintf("Random null: mean=%.1f, median=%d, max=%d\n",
            mean(random_degs), median(random_degs), max(random_degs)))
cat(sprintf("Observed/max ratio: %.0fx\n",
            n_obs_f / max(max(random_degs), 1)))

write_csv(data.frame(random_degs = random_degs),
          file.path(res_dir, "random_phenotype_null.csv"))

# Pathway correlation for random phenotypes
cat("Computing random phenotype pathway correlations...\n")
set.seed(456)
random_pathway_rhos <- numeric(N_RANDOM)
for (i in seq_len(N_RANDOM)) {
  tb_rand <- tb_f
  n_cases <- round(nrow(tb_rand) * ckd_prev_f)
  tb_rand$MHRNLFLR <- 0
  tb_rand$MHRNLFLR[sample(nrow(tb_rand), n_cases)] <- 1
  d_rand <- build_design(form_f, tb_rand)
  fit_r <- eBayes(lmFit(v_f_base, d_rand))
  fc <- tryCatch(
    get_required_coef_prefix(d_rand, "MHRNLFLR", "random pathway_cor"),
    error = function(e) NULL
  )
  if (is.null(fc)) { random_pathway_rhos[i] <- NA_real_; next }
  tt_r <- topTable(fit_r, coef = fc, number = Inf, sort.by = "none")
  tt_r$gene_id <- rownames(tt_r)
  random_pathway_rhos[i] <- pathway_cor(tt_r, primary_gsea, hallmark_list)
}
cat(sprintf("Random phenotype pathway rho: median=%.3f, mean=%.3f, range=%.3f--%.3f\n",
            median(random_pathway_rhos, na.rm = TRUE),
            mean(random_pathway_rhos, na.rm = TRUE),
            min(random_pathway_rhos, na.rm = TRUE),
            max(random_pathway_rhos, na.rm = TRUE)))

# ── Save all validation results ───────────────────────────────────────────────
sensitivity_summary <- list(
  n_obs = n_obs_f,
  comorbidity_models = comorbidity_models,
  perm_degs = perm_degs,
  perm_p = perm_p,
  perm_pathway_rhos = perm_pathway_rhos,
  loo_results = loo_results,
  random_degs = random_degs,
  random_pathway_rhos = random_pathway_rhos,
  pi0_female = res_female$pi0,
  pi0_male = res_male$pi0
)
saveRDS(sensitivity_summary, file.path(res_dir, "sensitivity_summary.rds"))
