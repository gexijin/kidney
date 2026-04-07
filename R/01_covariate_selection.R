#!/usr/bin/env Rscript
# =============================================================================
# Stage 1, Script 00b: Exploratory DE Analysis — Covariate Structure
#
# Questions addressed:
#   1. Which factors independently associate with gene expression?
#   2. For continuous factors, is the effect linear or should we bin/transform?
#   3. How similar are the effects of different factors (logFC correlation)?
#   4. What is the optimal base model for disease DE (e.g., CKD)?
#
# Input:  explore/renal_fat/paper/results/stage1_data_adipose_subcutaneous.rds
# Output: figures + results in explore/renal_fat/paper/{figures,results}/
# =============================================================================

library(tidyverse)
library(edgeR)
library(limma)
library(splines)
library(pheatmap)
library(cowplot)

source("kidney/R/functions.R")

set.seed(42)

res_dir <- "kidney/results"
fig_dir <- "kidney/figures"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 0: LOAD & PREPARE                                              ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 0: Load & Prepare\n")
cat(strrep("=", 70), "\n")

mh_factors_all <- c("MHABNWBC", "MHARTHTS", "MHASTHMA", "MHBCTINF", "MHBLDDND",
                     "MHCOCAINE5", "MHCOPD", "MHCVD", "MHDLYSIS", "MHDPRSSN",
                     "MHHRTATT", "MHHRTDIS", "MHHTN", "MHLVRDIS", "MHOPNWND",
                     "MHORGNTP", "MHPLLABS", "MHPNMNIA", "MHRNLFLR", "MHSEPSIS",
                     "MHSZRSU", "MHT1D", "MHT2D", "MHTXCEXP")
dat <- load_sat_data(res_dir = res_dir, mh_flags = mh_factors_all)
tb <- dat$tb
counts <- dat$counts

cat(sprintf("Loaded: %d samples, %d genes\n", nrow(tb), nrow(counts)))

# Additional factor types specific to covariate selection
tb$smoker <- factor(tb$smoker, levels = c(0, 1), labels = c("No", "Yes"))
tb$drinker <- factor(tb$drinker, levels = c("None", "Light", "Regular"), ordered = FALSE)
tb$SMNABTCHD_Y <- factor(tb$SMNABTCHD_Y)
tb$SMGEBTCHD_Y <- factor(tb$SMGEBTCHD_Y)

# Medication flags as integer
med_factors <- c("insulin_use", "statin_use", "opioid_use")
for (f in med_factors) tb[[f]] <- as.integer(tb[[f]])


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1: VARIANCE PARTITION — WHICH FACTORS MATTER?                   ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 1: Variance Partition — Which Factors Matter?\n")
cat(strrep("=", 70), "\n")

# Approach: For each factor, fit limma/voom in a minimal model (factor + intercept)
# and count DEGs at FDR < 0.05. This gives a quick screen of effect size.
# Then fit a joint model to see which factors remain significant after adjustment.

# ── 1a. Complete-case subset for all analyses ──────────────────────────────
# Use the widest set of covariates we might want
core_vars <- c("AGE", "SEX", "RACE", "BMI", "DTHHRDY", "ischemic_hrs",
               "SMRIN", "SMCENTER", "SMEXNCRT", "SMMAPRT", "SMRRNART",
               "smoker")
# Drop samples with missing core vars
tb_cc <- tb[complete.cases(tb[, core_vars]), ]
counts_cc <- counts[, tb_cc$SAMPID]
cat(sprintf("Complete cases (core vars): %d / %d samples\n", nrow(tb_cc), nrow(tb)))

# TMM normalize once
dge <- DGEList(counts = counts_cc)
dge <- calcNormFactors(dge, method = "TMM")

# ── 1b. Marginal screening: one factor at a time ───────────────────────────
cat("\n--- Marginal Factor Screening (one-at-a-time) ---\n")

# Define all factors to test
factor_list <- list(
  # Demographic
  AGE         = list(formula = ~ AGE,         coef = "AGE",     type = "demographic"),
  SEX         = list(formula = ~ SEX,         coef = "SEXFemale", type = "demographic"),
  RACE        = list(formula = ~ RACE,        coef = "RACEBlack", type = "demographic"),
  BMI         = list(formula = ~ BMI,         coef = "BMI",     type = "demographic"),
  smoker      = list(formula = ~ smoker,      coef = "smokerYes", type = "demographic"),
  drinker     = list(formula = ~ drinker,     coef = "auto",     type = "demographic"),
  # Technical
  DTHHRDY     = list(formula = ~ DTHHRDY,     coef = "auto",    type = "technical"),
  ischemic_hrs= list(formula = ~ ischemic_hrs,coef = "ischemic_hrs", type = "technical"),
  SMRIN       = list(formula = ~ SMRIN,        coef = "SMRIN",    type = "technical"),
  SMCENTER    = list(formula = ~ SMCENTER,     coef = "auto",    type = "technical"),
  SMEXNCRT    = list(formula = ~ SMEXNCRT,     coef = "SMEXNCRT", type = "technical"),
  SMMAPRT     = list(formula = ~ SMMAPRT,      coef = "SMMAPRT",  type = "technical"),
  SMRRNART    = list(formula = ~ SMRRNART,     coef = "SMRRNART", type = "technical"),
  SMNABTCHD_Y = list(formula = ~ SMNABTCHD_Y, coef = "auto",    type = "technical"),
  SMGEBTCHD_Y = list(formula = ~ SMGEBTCHD_Y, coef = "auto",    type = "technical"),
  # Disease
  MHHTN       = list(formula = ~ MHHTN,       coef = "MHHTN",    type = "disease"),
  MHT2D       = list(formula = ~ MHT2D,       coef = "MHT2D",    type = "disease"),
  MHCOPD      = list(formula = ~ MHCOPD,      coef = "MHCOPD",   type = "disease"),
  MHHRTDIS    = list(formula = ~ MHHRTDIS,     coef = "MHHRTDIS", type = "disease"),
  MHHRTATT    = list(formula = ~ MHHRTATT,     coef = "MHHRTATT", type = "disease"),
  MHCVD       = list(formula = ~ MHCVD,       coef = "MHCVD",    type = "disease"),
  MHRNLFLR    = list(formula = ~ MHRNLFLR,    coef = "MHRNLFLR", type = "disease"),
  MHABNWBC    = list(formula = ~ MHABNWBC,     coef = "MHABNWBC", type = "disease"),
  MHLVRDIS    = list(formula = ~ MHLVRDIS,     coef = "MHLVRDIS", type = "disease"),
  MHDPRSSN    = list(formula = ~ MHDPRSSN,     coef = "MHDPRSSN", type = "disease"),
  MHASTHMA    = list(formula = ~ MHASTHMA,     coef = "MHASTHMA", type = "disease"),
  MHPNMNIA    = list(formula = ~ MHPNMNIA,    coef = "MHPNMNIA", type = "disease"),
  MHORGNTP    = list(formula = ~ MHORGNTP,     coef = "MHORGNTP", type = "disease"),
  MHT1D       = list(formula = ~ MHT1D,        coef = "MHT1D",    type = "disease"),
  MHDLYSIS    = list(formula = ~ MHDLYSIS,     coef = "MHDLYSIS", type = "disease")
)

marginal_results <- list()

for (fname in names(factor_list)) {
  fi <- factor_list[[fname]]

  # Skip if too few cases for binary
  if (fname %in% names(tb_cc) && is.numeric(tb_cc[[fname]]) &&
      all(tb_cc[[fname]] %in% c(0, 1, NA))) {
    n_cases <- sum(tb_cc[[fname]] == 1, na.rm = TRUE)
    if (n_cases < 20) {
      cat(sprintf("  SKIP %s: only %d cases\n", fname, n_cases))
      next
    }
  }

  # Subset to non-missing for this factor
  tb_sub <- tb_cc[complete.cases(tb_cc[, fname]), ]
  counts_sub <- counts_cc[, tb_sub$SAMPID]

  design <- tryCatch(model.matrix(fi$formula, data = tb_sub), error = function(e) NULL)
  if (is.null(design)) { cat(sprintf("  SKIP %s: design error\n", fname)); next }

  # model.matrix drops NA rows; subset counts to match
  if (nrow(design) < nrow(tb_sub)) {
    keep_idx <- as.integer(rownames(design))
    tb_sub <- tb_sub[keep_idx, ]
    counts_sub <- counts_cc[, tb_sub$SAMPID]
  }

  # Check rank
  qr_obj <- qr(design)
  if (qr_obj$rank < ncol(design))
    design <- design[, qr_obj$pivot[seq_len(qr_obj$rank)], drop = FALSE]

  dge_sub <- DGEList(counts = counts_sub)
  dge_sub <- calcNormFactors(dge_sub, method = "TMM")
  v <- voom(dge_sub, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))

  # Determine which coefficients to test
  if (identical(fi$coef, "auto")) {
    # F-test on all non-intercept columns
    test_coefs <- 2:ncol(design)
  } else if (is.character(fi$coef) && !(fi$coef %in% colnames(design))) {
    cat(sprintf("  SKIP %s: coef '%s' not in design\n", fname, fi$coef))
    next
  } else {
    test_coefs <- fi$coef
  }
  tt <- topTable(fit, coef = test_coefs, number = Inf, sort.by = "none")

  n_deg <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
  n_nom <- sum(tt$P.Value < 0.05, na.rm = TRUE)
  # F column for multi-coef tests, t column for single-coef
  stat_col <- if ("F" %in% names(tt)) tt$F else tt$t
  median_f <- median(abs(stat_col), na.rm = TRUE)

  marginal_results[[fname]] <- tibble(
    factor = fname,
    type = fi$type,
    n_samples = nrow(tb_sub),
    n_deg_fdr05 = n_deg,
    n_nom_p05 = n_nom,
    pct_nom = 100 * n_nom / nrow(tt),
    median_abs_F = median_f
  )

  cat(sprintf("  %-15s type=%-12s n=%3d  DEGs(FDR<0.05)=%5d  nom(p<0.05)=%5d (%.1f%%)\n",
              fname, fi$type, nrow(tb_sub), n_deg, n_nom, 100 * n_nom / nrow(tt)))
}

marginal_df <- bind_rows(marginal_results) %>% arrange(desc(n_deg_fdr05))
write_csv(marginal_df, file.path(res_dir, "exploratory_marginal_screening.csv"))

cat("\n--- Top factors by DEG count (FDR<0.05) ---\n")
print(marginal_df %>% head(15), n = 15)


# ── 1c. Figure: Marginal DEG counts by factor ──────────────────────────────
cat("\nGenerating marginal screening barplot...\n")

p_marginal <- ggplot(marginal_df %>% filter(n_deg_fdr05 > 0),
                      aes(x = reorder(factor, n_deg_fdr05),
                          y = n_deg_fdr05, fill = type)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(demographic = "#E69F00", technical = "#56B4E9",
                                disease = "#CC79A7")) +
  labs(x = NULL, y = "DEGs (FDR < 0.05)", fill = "Factor type",
       title = "Marginal Factor Screening: Subcutaneous Adipose") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "01_marginal_degs.png"), p_marginal,
       width = 8, height = 8, dpi = 300)
cat("Saved:", file.path(fig_dir, "01_marginal_degs.png"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2: CONTINUOUS FACTOR LINEARITY — AGE, ISCHEMIC TIME, BMI, RIN  ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 2: Linearity Assessment for Continuous Factors\n")
cat(strrep("=", 70), "\n")

# Strategy: for each continuous variable, compare:
#   (a) linear term
#   (b) natural spline with 3 df
#   (c) binned (optimal cut-points)
# Compare via AIC-like criterion (sum of gene-level deviance or BIC from limma)

continuous_vars <- list(
  ischemic_hrs = list(
    bins = list(
      terciles = quantile(tb_cc$ischemic_hrs, c(0, 1/3, 2/3, 1)),
      clinical = c(0, 3, 10, 25),      # clinical: short (<3h), medium (3-10h), long (>10h)
      fine     = c(0, 2, 4, 8, 14, 25)  # finer bins
    )
  ),
  AGE = list(
    bins = list(
      decades = c(20, 30, 40, 50, 60, 70),
      median_split = c(20, median(tb_cc$AGE), 70)
    )
  ),
  BMI = list(
    bins = list(
      who = c(0, 18.5, 25, 30, Inf),     # underweight, normal, overweight, obese
      terciles = quantile(tb_cc$BMI, c(0, 1/3, 2/3, 1))
    )
  ),
  SMRIN = list(
    bins = list(
      terciles = quantile(tb_cc$SMRIN, c(0, 1/3, 2/3, 1))
    )
  )
)

# For linearity test: fit linear vs spline, compare number of genes where
# spline is significantly better (ANOVA-like comparison)
cat("\n--- Linearity Tests (linear vs ns(df=3)) ---\n")

linearity_results <- list()

for (vname in names(continuous_vars)) {
  cat(sprintf("\n  Variable: %s\n", vname))

  x <- tb_cc[[vname]]
  if (all(is.na(x))) { cat("    SKIP: all NA\n"); next }

  # Designs: linear vs spline
  design_lin <- model.matrix(~ x)
  ns_mat <- ns(x, df = 3)
  design_spl <- model.matrix(~ ns_mat)

  v_lin <- voom(dge, design_lin, plot = FALSE)
  v_spl <- voom(dge, design_spl, plot = FALSE)

  fit_lin <- eBayes(lmFit(v_lin, design_lin))
  fit_spl <- eBayes(lmFit(v_spl, design_spl))

  # Gene-level comparison: for each gene, compute residual variance under both
  # Use sigma^2 as proxy
  sigma_lin <- fit_lin$sigma^2
  sigma_spl <- fit_spl$sigma^2

  # Fraction of genes where spline improves > 5%
  improvement <- (sigma_lin - sigma_spl) / sigma_lin
  pct_improved <- 100 * mean(improvement > 0.05)
  median_improvement <- 100 * median(improvement)

  # DEGs: linear test (coef 2) and nonlinearity test (extra spline terms, coefs 3:4)
  tt_lin <- topTable(fit_lin, coef = 2, number = Inf, sort.by = "none")
  tt_nonlin <- topTable(fit_spl, coef = 3:4, number = Inf, sort.by = "none")
  n_lin <- sum(tt_lin$adj.P.Val < 0.05)
  n_nonlin <- sum(tt_nonlin$adj.P.Val < 0.05)

  linearity_results[[vname]] <- tibble(
    variable = vname,
    n_deg_linear = n_lin,
    n_deg_nonlinear = n_nonlin,
    pct_genes_spline_5pct_better = round(pct_improved, 1),
    median_var_reduction_pct = round(median_improvement, 2)
  )

  cat(sprintf("    Linear DEGs: %d | Nonlinear DEGs (departure from linearity): %d\n",
              n_lin, n_nonlin))
  cat(sprintf("    Genes where spline >5%% better: %.1f%% | Median var reduction: %.2f%%\n",
              pct_improved, median_improvement))

  # ── Binning comparison ────────────────────────────────────────────────────
  for (bin_name in names(continuous_vars[[vname]]$bins)) {
    breaks <- continuous_vars[[vname]]$bins[[bin_name]]
    bin_var <- cut(x, breaks = breaks, include.lowest = TRUE)

    # Drop NAs from binning (values outside breaks)
    keep_bin <- !is.na(bin_var)
    bin_var <- droplevels(bin_var[keep_bin])

    if (any(table(bin_var) < 10)) {
      cat(sprintf("    Bin '%s': some bins < 10 samples, skip DE\n", bin_name))
      next
    }

    design_bin <- model.matrix(~ bin_var)
    qr_bin <- qr(design_bin)
    if (qr_bin$rank < ncol(design_bin))
      design_bin <- design_bin[, qr_bin$pivot[seq_len(qr_bin$rank)], drop = FALSE]

    # Subset DGE to matching samples
    dge_bin <- dge[, keep_bin]
    v_bin <- voom(dge_bin, design_bin, plot = FALSE)
    fit_bin <- eBayes(lmFit(v_bin, design_bin))
    tt_bin <- topTable(fit_bin, coef = 2:ncol(design_bin), number = Inf, sort.by = "none")
    n_bin <- sum(tt_bin$adj.P.Val < 0.05)
    cat(sprintf("    Bin '%s' (%d levels, n=%d): %d DEGs\n",
                bin_name, nlevels(bin_var), sum(keep_bin), n_bin))
  }
}

linearity_df <- bind_rows(linearity_results)
write_csv(linearity_df, file.path(res_dir, "exploratory_linearity.csv"))


# ── 2b. Ischemic time deep-dive: gene-level nonlinearity ───────────────────
cat("\n--- Ischemic Time: Deep Dive ---\n")

# Bin ischemic time more finely and look at how mean expression changes
isch_breaks <- c(0, 1, 2, 3, 4, 6, 8, 10, 14, 18, 25)
tb_cc$isch_fine <- cut(tb_cc$ischemic_hrs, breaks = isch_breaks, include.lowest = TRUE)
cat("Ischemic time fine bins:\n")
print(table(tb_cc$isch_fine))

# For top 20 ischemic-time-responsive genes, plot expression vs bin
design_isch <- model.matrix(~ ischemic_hrs, data = tb_cc)
v_isch <- voom(dge, design_isch, plot = FALSE)
fit_isch <- eBayes(lmFit(v_isch, design_isch))
tt_isch <- topTable(fit_isch, coef = "ischemic_hrs", number = 20, sort.by = "p")

# Get logCPM values
logcpm <- cpm(dge, log = TRUE, prior.count = 4)

cat("\nTop 20 ischemic-time genes:\n")
print(tt_isch[, c("logFC", "AveExpr", "t", "adj.P.Val")])

# Plot mean expression per bin for top genes
top_isch_genes <- rownames(tt_isch)[1:min(6, nrow(tt_isch))]

plot_data_isch <- list()
for (g in top_isch_genes) {
  plot_data_isch[[g]] <- tibble(
    gene = g,
    isch_bin = tb_cc$isch_fine,
    isch_hrs = tb_cc$ischemic_hrs,
    logcpm = logcpm[g, tb_cc$SAMPID]
  )
}
plot_df_isch <- bind_rows(plot_data_isch)

p_isch_genes <- ggplot(plot_df_isch, aes(x = isch_bin, y = logcpm)) +
  geom_boxplot(outlier.size = 0.5, fill = "lightyellow") +
  facet_wrap(~ gene, scales = "free_y") +
  labs(x = "Ischemic Time Bin (hours)", y = "log2 CPM",
       title = "Top Ischemic-Time Genes: Expression vs Time Bin") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(file.path(fig_dir, "01_ischemic_time_genes.png"), p_isch_genes,
       width = 12, height = 8, dpi = 300)
cat("Saved:", file.path(fig_dir, "01_ischemic_time_genes.png"), "\n")

# ── 2c. Ischemic time: genome-wide nonlinearity test ───────────────────────
# Compare: number of DEGs from linear, log-transformed, and clinical bins
cat("\n--- Ischemic Time: Transformation Comparison ---\n")

tb_cc$log_isch <- log1p(tb_cc$ischemic_hrs)
tb_cc$sqrt_isch <- sqrt(tb_cc$ischemic_hrs)

isch_transforms <- list(
  linear   = ~ ischemic_hrs,
  log1p    = ~ log_isch,
  sqrt     = ~ sqrt_isch,
  spline3  = as.formula(sprintf("~ ns(ischemic_hrs, df = 3)")),
  clinical_3bin = ~ cut(ischemic_hrs, breaks = c(0, 3, 10, 25), include.lowest = TRUE)
)

for (tname in names(isch_transforms)) {
  design_t <- tryCatch(
    model.matrix(isch_transforms[[tname]], data = tb_cc),
    error = function(e) NULL
  )
  if (is.null(design_t)) { cat(sprintf("  %s: design error\n", tname)); next }

  qr_t <- qr(design_t)
  if (qr_t$rank < ncol(design_t))
    design_t <- design_t[, qr_t$pivot[seq_len(qr_t$rank)], drop = FALSE]

  v_t <- voom(dge, design_t, plot = FALSE)
  fit_t <- eBayes(lmFit(v_t, design_t))
  tt_t <- topTable(fit_t, coef = 2:ncol(design_t), number = Inf, sort.by = "none")
  n_deg <- sum(tt_t$adj.P.Val < 0.05)
  cat(sprintf("  %-20s: %d DEGs (FDR<0.05), %d cols\n", tname, n_deg, ncol(design_t)))
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 3: FACTOR EFFECT SIMILARITY — logFC CORRELATION MATRIX          ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 3: Factor Effect Similarity (logFC Correlations)\n")
cat(strrep("=", 70), "\n")

# Fit a comprehensive model with all major factors and extract per-factor
# t-statistics. The correlation of t-statistics across genes reveals which
# factors have similar genome-wide effects.

# Build a comprehensive design with all factors
# Use only factors with ≥20 cases and reasonable prevalence
factors_for_joint <- c(
  "AGE", "SEX", "RACE", "BMI", "DTHHRDY", "ischemic_hrs", "SMRIN",
  "SMCENTER", "SMEXNCRT", "smoker", "drinker",
  "MHHTN", "MHT2D", "MHCOPD", "MHHRTDIS", "MHHRTATT", "MHCVD",
  "MHRNLFLR", "MHABNWBC", "MHLVRDIS", "MHDPRSSN", "MHASTHMA", "MHPNMNIA"
)

# Complete cases
tb_joint <- tb_cc[complete.cases(tb_cc[, factors_for_joint]), ]
counts_joint <- counts_cc[, tb_joint$SAMPID]
cat(sprintf("Joint model samples: %d\n", nrow(tb_joint)))

joint_formula <- as.formula(paste(
  "~", paste(factors_for_joint, collapse = " + ")
))

design_joint <- model.matrix(joint_formula, data = tb_joint)
qr_joint <- qr(design_joint)
if (qr_joint$rank < ncol(design_joint))
  design_joint <- design_joint[, qr_joint$pivot[seq_len(qr_joint$rank)], drop = FALSE]

cat(sprintf("Joint design: %d samples x %d columns\n", nrow(design_joint), ncol(design_joint)))

dge_joint <- DGEList(counts = counts_joint)
dge_joint <- calcNormFactors(dge_joint, method = "TMM")
v_joint <- voom(dge_joint, design_joint, plot = FALSE)
fit_joint <- eBayes(lmFit(v_joint, design_joint))

# Extract t-statistics for each factor of interest
# For multi-level factors (DTHHRDY), take the first non-intercept level
coef_names <- colnames(design_joint)
cat("Design columns:", paste(coef_names, collapse = ", "), "\n")

# Map factor names to coefficient columns
factor_coef_map <- list(
  AGE = "AGE",
  SEX = "SEXFemale",
  RACE = "RACEBlack",
  BMI = "BMI",
  DTHHRDY_Vent = "DTHHRDYVentilator",
  DTHHRDY_Slow = "DTHHRDYSlow",
  ischemic_hrs = "ischemic_hrs",
  SMRIN = "SMRIN",
  SMEXNCRT = "SMEXNCRT",
  smoker = "smokerYes",
  drinker_Light = "drinkerLight",
  drinker_Regular = "drinkerRegular",
  MHHTN = "MHHTN",
  MHT2D = "MHT2D",
  MHCOPD = "MHCOPD",
  MHHRTDIS = "MHHRTDIS",
  MHHRTATT = "MHHRTATT",
  MHCVD = "MHCVD",
  MHRNLFLR = "MHRNLFLR",
  MHABNWBC = "MHABNWBC",
  MHLVRDIS = "MHLVRDIS",
  MHDPRSSN = "MHDPRSSN",
  MHASTHMA = "MHASTHMA",
  MHPNMNIA = "MHPNMNIA"
)

# Only keep coefficients that exist in design
factor_coef_map <- factor_coef_map[sapply(factor_coef_map, function(x) x %in% coef_names)]

# Extract t-statistics matrix
t_mat <- matrix(NA, nrow = nrow(fit_joint$t), ncol = length(factor_coef_map),
                dimnames = list(rownames(fit_joint$t), names(factor_coef_map)))
for (flab in names(factor_coef_map)) {
  coef_col <- factor_coef_map[[flab]]
  t_mat[, flab] <- fit_joint$t[, coef_col]
}

# ── 3a. Correlation matrix of t-statistics ──────────────────────────────────
cor_t <- cor(t_mat, method = "spearman", use = "pairwise.complete.obs")

cat("\n--- Factor t-statistic correlations (Spearman) ---\n")
# Print strong correlations
for (i in 1:(ncol(cor_t)-1)) {
  for (j in (i+1):ncol(cor_t)) {
    if (abs(cor_t[i,j]) > 0.15) {
      cat(sprintf("  %s ~ %s: rho = %+.3f\n",
                  colnames(cor_t)[i], colnames(cor_t)[j], cor_t[i,j]))
    }
  }
}

# ── 3b. Figure: Correlation heatmap of factor effects ───────────────────────
# Group by type for better visualization
type_order <- c("AGE", "SEX", "RACE", "BMI", "smoker", "drinker_Light", "drinker_Regular",
                "DTHHRDY_Vent", "DTHHRDY_Slow", "ischemic_hrs", "SMRIN", "SMEXNCRT",
                "MHHTN", "MHT2D", "MHCOPD", "MHHRTDIS", "MHHRTATT", "MHCVD",
                "MHRNLFLR", "MHABNWBC", "MHLVRDIS", "MHDPRSSN", "MHASTHMA", "MHPNMNIA")
type_order <- intersect(type_order, colnames(cor_t))
cor_t_ordered <- cor_t[type_order, type_order]

# Annotation row for factor type
ann_type <- data.frame(
  Type = case_when(
    type_order %in% c("AGE", "SEX", "RACE", "BMI", "smoker",
                       "drinker_Light", "drinker_Regular") ~ "Demographic",
    type_order %in% c("DTHHRDY_Vent", "DTHHRDY_Slow", "ischemic_hrs", "SMRIN", "SMEXNCRT") ~ "Technical",
    TRUE ~ "Disease"
  ),
  row.names = type_order
)
ann_colors_type <- list(
  Type = c(Demographic = "#E69F00", Technical = "#56B4E9", Disease = "#CC79A7")
)

png(file.path(fig_dir, "01_factor_correlation.png"),
    width = 12, height = 10, units = "in", res = 300)
pheatmap(
  cor_t_ordered,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 7,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  annotation_row = ann_type,
  annotation_colors = ann_colors_type,
  cluster_rows = FALSE, cluster_cols = FALSE,
  main = "Genome-Wide Effect Similarity: Spearman rho of t-statistics\n(joint model, subcutaneous adipose)"
)
dev.off()
cat("Saved:", file.path(fig_dir, "01_factor_correlation.png"), "\n")


# ── 3c. Hierarchical clustering of factor effects ──────────────────────────
png(file.path(fig_dir, "01_factor_clustering.png"),
    width = 12, height = 10, units = "in", res = 300)
pheatmap(
  cor_t_ordered,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 7,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  annotation_row = ann_type,
  annotation_colors = ann_colors_type,
  cluster_rows = TRUE, cluster_cols = TRUE,
  clustering_distance_rows = as.dist(sqrt(2 * (1 - cor_t_ordered))),
  clustering_distance_cols = as.dist(sqrt(2 * (1 - cor_t_ordered))),
  main = "Factor Effect Clustering (Euclidean distance from correlation)\n(joint model, subcutaneous adipose)"
)
dev.off()
cat("Saved:", file.path(fig_dir, "01_factor_clustering.png"), "\n")


# ── 3d. Adjusted DEG counts from joint model ───────────────────────────────
cat("\n--- Adjusted DEG counts (from joint model) ---\n")

adjusted_degs <- list()
for (flab in names(factor_coef_map)) {
  coef_col <- factor_coef_map[[flab]]
  tt <- topTable(fit_joint, coef = coef_col, number = Inf, sort.by = "none")
  n_deg <- sum(tt$adj.P.Val < 0.05)
  n_up <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 0)
  n_down <- sum(tt$adj.P.Val < 0.05 & tt$logFC < 0)

  adjusted_degs[[flab]] <- tibble(
    factor = flab,
    coef = coef_col,
    n_deg_adjusted = n_deg,
    n_up = n_up,
    n_down = n_down
  )
  if (n_deg > 0) {
    cat(sprintf("  %-15s: %5d DEGs (%d up, %d down)\n", flab, n_deg, n_up, n_down))
  }
}

adjusted_df <- bind_rows(adjusted_degs) %>% arrange(desc(n_deg_adjusted))
write_csv(adjusted_df, file.path(res_dir, "exploratory_adjusted_degs.csv"))


# ── 3e. Compare marginal vs adjusted (same sample set) ─────────────────────
# Re-run marginal screening on tb_joint for an apples-to-apples comparison
cat(sprintf("
--- Marginal vs Adjusted DEG Counts (same n=%d samples) ---
",
            nrow(tb_joint)))

marginal_joint_degs <- list()
already_done <- c()
for (flab in names(factor_coef_map)) {
  coef_col <- factor_coef_map[[flab]]
  # Map coefficient back to original variable
  orig_var <- dplyr::case_when(
    flab %in% c("DTHHRDY_Vent", "DTHHRDY_Slow") ~ "DTHHRDY",
    flab %in% c("drinker_Light", "drinker_Regular") ~ "drinker",
    flab == "SEX" ~ "SEX",
    flab == "RACE" ~ "RACE",
    flab == "smoker" ~ "smoker",
    TRUE ~ coef_col
  )

  form_marg <- as.formula(paste("~", orig_var))
  design_marg <- tryCatch(model.matrix(form_marg, data = tb_joint), error = function(e) NULL)
  if (is.null(design_marg)) next

  qr_marg <- qr(design_marg)
  if (qr_marg$rank < ncol(design_marg))
    design_marg <- design_marg[, qr_marg$pivot[seq_len(qr_marg$rank)], drop = FALSE]

  v_marg <- voom(dge_joint, design_marg, plot = FALSE)
  fit_marg <- eBayes(lmFit(v_marg, design_marg))

  # Use named coef if present, else F-test on all non-intercept
  if (coef_col %in% colnames(design_marg)) {
    tt_marg <- topTable(fit_marg, coef = coef_col, number = Inf, sort.by = "none")
  } else {
    tt_marg <- topTable(fit_marg, coef = 2:ncol(design_marg), number = Inf, sort.by = "none")
  }
  marginal_joint_degs[[flab]] <- sum(tt_marg$adj.P.Val < 0.05)
}

comparison <- tibble(
  factor = names(factor_coef_map),
  n_deg_marginal = as.integer(marginal_joint_degs[names(factor_coef_map)])
) %>%
  left_join(adjusted_df %>% dplyr::select(factor, n_deg_adjusted), by = "factor") %>%
  mutate(retention_pct = ifelse(n_deg_marginal > 0,
                                 100 * n_deg_adjusted / n_deg_marginal, NA)) %>%
  arrange(desc(n_deg_marginal))

cat(sprintf("%-18s %10s %10s %10s
", "Factor", "Marginal", "Adjusted", "Retained%"))
cat(strrep("-", 60), "
")
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  cat(sprintf("%-18s %10d %10s %10s
",
              r$factor, r$n_deg_marginal,
              ifelse(is.na(r$n_deg_adjusted), "---", as.character(r$n_deg_adjusted)),
              ifelse(is.na(r$retention_pct), "---", sprintf("%.0f%%", r$retention_pct))))
}



# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 4: ISCHEMIC TIME — OPTIMAL BINNING                             ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 4: Ischemic Time — Optimal Binning\n")
cat(strrep("=", 70), "\n")

# The key question: does the effect plateau or change slope at certain times?
# Strategy: fit a model with ischemic time bins and plot the bin coefficients
# relative to reference bin. This gives a "dose-response" curve.

# Fine binning for dose-response visualization
isch_bin_breaks <- c(0, 1, 2, 3, 4, 6, 8, 10, 14, 18, 25)
tb_joint$isch_fine <- cut(tb_joint$ischemic_hrs,
                           breaks = isch_bin_breaks,
                           include.lowest = TRUE)

# Use a model with isch_fine instead of linear ischemic_hrs
form_isch_bin <- as.formula(paste(
  "~ AGE + SEX + RACE + BMI + DTHHRDY + SMRIN + SMCENTER + SMEXNCRT +",
  "smoker + MHHTN + MHT2D + MHCOPD + MHHRTDIS + MHRNLFLR + isch_fine"
))

design_isch_bin <- model.matrix(form_isch_bin, data = tb_joint)
qr_ib <- qr(design_isch_bin)
if (qr_ib$rank < ncol(design_isch_bin))
  design_isch_bin <- design_isch_bin[, qr_ib$pivot[seq_len(qr_ib$rank)], drop = FALSE]

v_ib <- voom(dge_joint, design_isch_bin, plot = FALSE)
fit_ib <- eBayes(lmFit(v_ib, design_isch_bin))

# Extract ischemic bin coefficients
isch_coefs <- grep("isch_fine", colnames(design_isch_bin), value = TRUE)
cat("Ischemic bin coefficients:\n")

bin_effects <- list()
for (ic in isch_coefs) {
  tt_ic <- topTable(fit_ib, coef = ic, number = Inf, sort.by = "none")

  # Median absolute logFC across all genes (global effect size)
  med_abs_lfc <- median(abs(tt_ic$logFC))
  # Number of DEGs
  n_deg <- sum(tt_ic$adj.P.Val < 0.05)

  # Extract bin label
  bin_label <- gsub("isch_fine", "", ic)

  bin_effects[[ic]] <- tibble(
    bin = bin_label,
    coef = ic,
    median_abs_logFC = med_abs_lfc,
    n_deg = n_deg,
    mean_logFC = mean(tt_ic$logFC),
    mean_t = mean(tt_ic$t)
  )

  cat(sprintf("  %s: median|logFC|=%.4f, DEGs=%d\n", bin_label, med_abs_lfc, n_deg))
}

bin_effect_df <- bind_rows(bin_effects)
# Add midpoint of each bin for plotting
midpoints <- (isch_bin_breaks[-length(isch_bin_breaks)] + isch_bin_breaks[-1]) / 2
# Reference bin is first; others are relative to it
bin_effect_df$midpoint <- midpoints[-1]  # skip reference

# Plot: ischemic time dose-response
p_isch_dose <- ggplot(bin_effect_df, aes(x = midpoint, y = median_abs_logFC)) +
  geom_point(size = 3) +
  geom_line() +
  geom_text(aes(label = n_deg), vjust = -1, size = 3) +
  labs(x = "Ischemic Time (hours, bin midpoint)",
       y = "Median |logFC| (all genes, vs reference [0,1])",
       title = "Ischemic Time Dose-Response in Subcutaneous Adipose",
       subtitle = "Numbers = DEGs per bin (FDR<0.05). Adjusted for demographics + disease.") +
  theme_bw(base_size = 12)

ggsave(file.path(fig_dir, "01_ischemic_dose_response.png"), p_isch_dose,
       width = 9, height = 6, dpi = 300)
cat("Saved:", file.path(fig_dir, "01_ischemic_dose_response.png"), "\n")


# ── 4b. Top ischemic genes: expression trajectories ────────────────────────
# Identify genes with strongest nonlinear ischemic response
# Fit spline model, find genes with biggest spline improvement over linear

design_lin_full <- model.matrix(
  ~ AGE + SEX + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER +
    SMEXNCRT + smoker + MHHTN + MHT2D + MHCOPD + MHHRTDIS + MHRNLFLR,
  data = tb_joint)

design_spl_full <- model.matrix(
  ~ AGE + SEX + RACE + BMI + DTHHRDY + ns(ischemic_hrs, df = 3) + SMRIN +
    SMCENTER + SMEXNCRT + smoker + MHHTN + MHT2D + MHCOPD + MHHRTDIS + MHRNLFLR,
  data = tb_joint)

v_lin_full <- voom(dge_joint, design_lin_full, plot = FALSE)
v_spl_full <- voom(dge_joint, design_spl_full, plot = FALSE)
fit_lin_full <- lmFit(v_lin_full, design_lin_full)
fit_spl_full <- lmFit(v_spl_full, design_spl_full)

# Compare residual variances
resid_var_lin <- rowMeans(residuals(fit_lin_full, v_lin_full)^2)
resid_var_spl <- rowMeans(residuals(fit_spl_full, v_spl_full)^2)
nonlin_improvement <- (resid_var_lin - resid_var_spl) / resid_var_lin

cat(sprintf("\nNonlinearity improvement (ischemic time):\n"))
cat(sprintf("  Median: %.3f%%\n", 100 * median(nonlin_improvement)))
cat(sprintf("  95th percentile: %.3f%%\n", 100 * quantile(nonlin_improvement, 0.95)))
cat(sprintf("  Genes with >5%% improvement: %d (%.1f%%)\n",
            sum(nonlin_improvement > 0.05),
            100 * mean(nonlin_improvement > 0.05)))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 5: DISEASE EFFECT SIMILARITY DEEP-DIVE                         ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 5: Disease Effect Similarity Deep-Dive\n")
cat(strrep("=", 70), "\n")

# Focus on disease factors: do they look like each other? Like aging?
# Like technical confounders?

# ── 5a. Pairwise scatter of CKD vs other factors ───────────────────────────
key_comparisons <- c("AGE", "ischemic_hrs", "DTHHRDY_Vent", "DTHHRDY_Slow",
                      "MHHTN", "MHT2D", "MHCVD", "MHLVRDIS", "MHABNWBC",
                      "SEX", "BMI")
key_comparisons <- intersect(key_comparisons, colnames(t_mat))

scatter_data <- list()
for (comp in key_comparisons) {
  scatter_data[[comp]] <- tibble(
    gene_id = rownames(t_mat),
    t_CKD = t_mat[, "MHRNLFLR"],
    t_other = t_mat[, comp],
    comparison = comp
  )
}
scatter_df <- bind_rows(scatter_data)

p_ckd_scatter <- ggplot(scatter_df, aes(x = t_other, y = t_CKD)) +
  geom_point(size = 0.2, alpha = 0.1, color = "grey40") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
  facet_wrap(~ comparison, scales = "free_x", ncol = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  labs(x = "t-statistic (comparison factor)",
       y = "t-statistic (CKD / MHRNLFLR)",
       title = "CKD Effect Similarity with Other Factors\n(gene-level t-statistics from joint model)") +
  theme_bw(base_size = 10)

ggsave(file.path(fig_dir, "01_ckd_factor_scatter.png"), p_ckd_scatter,
       width = 14, height = 10, dpi = 300)
cat("Saved:", file.path(fig_dir, "01_ckd_factor_scatter.png"), "\n")


# ── 5b. CKD specificity: what fraction of CKD signal overlaps with each factor?
cat("\n--- CKD Signal Specificity ---\n")

# For each factor, compute: among CKD DEGs (from joint model), what fraction
# are also DEGs for that factor?
ckd_tt <- topTable(fit_joint, coef = "MHRNLFLR", number = Inf, sort.by = "none")
ckd_degs <- rownames(ckd_tt)[ckd_tt$adj.P.Val < 0.05]
cat(sprintf("CKD DEGs (joint model, FDR<0.05): %d\n", length(ckd_degs)))

if (length(ckd_degs) > 0) {
  cat(sprintf("\n%-15s %8s %8s %10s %10s\n",
              "Factor", "Its DEGs", "Overlap", "Overlap%", "rho(CKD)"))
  cat(strrep("-", 60), "\n")

  for (flab in names(factor_coef_map)) {
    if (flab == "MHRNLFLR") next
    coef_col <- factor_coef_map[[flab]]
    tt_f <- topTable(fit_joint, coef = coef_col, number = Inf, sort.by = "none")
    f_degs <- rownames(tt_f)[tt_f$adj.P.Val < 0.05]
    overlap <- intersect(ckd_degs, f_degs)
    rho <- cor_t["MHRNLFLR", flab]

    cat(sprintf("  %-15s %8d %8d %9.1f%% %+10.3f\n",
                flab, length(f_degs), length(overlap),
                100 * length(overlap) / length(ckd_degs), rho))
  }
} else {
  cat("No CKD DEGs in joint model — skip overlap analysis\n")
  cat("(CKD may need sex-stratified analysis to show signal)\n")
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 6: CONFOUNDER STRUCTURE — WHAT COVARIES WITH WHAT?              ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 6: Confounder Structure\n")
cat(strrep("=", 70), "\n")

# ── 6a. Continuous-continuous correlations ──────────────────────────────────
cont_vars <- c("AGE", "BMI", "ischemic_hrs", "SMRIN", "SMEXNCRT", "SMMAPRT", "SMRRNART")
cont_cor <- cor(tb_cc[, cont_vars], use = "pairwise.complete.obs", method = "spearman")

cat("Continuous variable correlations (Spearman):\n")
for (i in 1:(ncol(cont_cor)-1)) {
  for (j in (i+1):ncol(cont_cor)) {
    if (abs(cont_cor[i,j]) > 0.1) {
      cat(sprintf("  %s ~ %s: rho = %+.3f\n",
                  colnames(cont_cor)[i], colnames(cont_cor)[j], cont_cor[i,j]))
    }
  }
}

# ── 6b. Disease-demographic associations ───────────────────────────────────
cat("\n--- Disease Factor ~ Demographics ---\n")

disease_vars <- c("MHHTN", "MHT2D", "MHCOPD", "MHHRTDIS", "MHHRTATT",
                  "MHCVD", "MHRNLFLR", "MHABNWBC", "MHLVRDIS")

cat(sprintf("%-12s %6s %6s %6s %6s %6s %8s\n",
            "Disease", "N+", "Age+", "Age-", "pAge", "%Male", "pSex"))
cat(strrep("-", 60), "\n")

for (d in disease_vars) {
  n_pos <- sum(tb_cc[[d]] == 1, na.rm = TRUE)
  age_pos <- mean(tb_cc$AGE[tb_cc[[d]] == 1], na.rm = TRUE)
  age_neg <- mean(tb_cc$AGE[tb_cc[[d]] == 0], na.rm = TRUE)
  p_age <- wilcox.test(tb_cc$AGE ~ tb_cc[[d]])$p.value
  pct_male <- 100 * mean(tb_cc$SEX[tb_cc[[d]] == 1] == "Male", na.rm = TRUE)
  p_sex <- tryCatch(fisher.test(tb_cc$SEX, tb_cc[[d]])$p.value, error = function(e) NA)

  cat(sprintf("  %-12s %6d %6.1f %6.1f %6.1e %7.0f%% %8.1e\n",
              d, n_pos, age_pos, age_neg, p_age, pct_male, p_sex))
}

# ── 6c. Disease-disease comorbidity (phi coefficients) ──────────────────────
cat("\n--- Disease Comorbidity (phi coefficients) ---\n")

phi_mat <- matrix(NA, length(disease_vars), length(disease_vars),
                  dimnames = list(disease_vars, disease_vars))

for (i in seq_along(disease_vars)) {
  for (j in seq_along(disease_vars)) {
    a <- tb_cc[[disease_vars[i]]]
    b <- tb_cc[[disease_vars[j]]]
    cc <- complete.cases(a, b)
    if (sum(cc) < 50) next
    tab <- table(a[cc], b[cc])
    if (nrow(tab) == 2 && ncol(tab) == 2) {
      n <- sum(tab)
      phi_mat[i, j] <- (tab[2,2]*tab[1,1] - tab[1,2]*tab[2,1]) /
        sqrt(prod(rowSums(tab)) * prod(colSums(tab)))
    }
  }
}

# Print strong comorbidities
for (i in 1:(nrow(phi_mat)-1)) {
  for (j in (i+1):ncol(phi_mat)) {
    if (!is.na(phi_mat[i,j]) && abs(phi_mat[i,j]) > 0.15) {
      cat(sprintf("  %s ~ %s: phi = %+.3f\n",
                  disease_vars[i], disease_vars[j], phi_mat[i,j]))
    }
  }
}

# ── 6d. Ischemic time vs Hardy scale ───────────────────────────────────────
cat("\n--- Ischemic Time vs Hardy Scale ---\n")
for (h in levels(tb_cc$DTHHRDY)) {
  isch <- tb_cc$ischemic_hrs[tb_cc$DTHHRDY == h]
  cat(sprintf("  %s: n=%d, median=%.1f, IQR=[%.1f, %.1f]\n",
              h, length(isch), median(isch),
              quantile(isch, 0.25), quantile(isch, 0.75)))
}

p_isch_hardy <- ggplot(tb_cc, aes(x = DTHHRDY, y = ischemic_hrs, fill = DTHHRDY)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = "Hardy Scale", y = "Ischemic Time (hours)",
       title = "Ischemic Time Distribution by Hardy Scale") +
  theme_bw(base_size = 12)

ggsave(file.path(fig_dir, "01_isch_hardy.png"), p_isch_hardy,
       width = 6, height = 5, dpi = 300)


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 7: RECOMMENDED BASE MODEL                                      ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 7: Recommended Base Model\n")
cat(strrep("=", 70), "\n")

# Based on the analysis above, test candidate base models and compare
# their performance on a held-out metric (variance explained).

# Candidate models:
models <- list(
  minimal = "~ AGE + SEX + RACE + DTHHRDY + ischemic_hrs + SMRIN",
  standard = "~ AGE + SEX + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER",
  extended = "~ AGE + SEX + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + smoker",
  full_tech = "~ AGE + SEX + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + SMMAPRT + SMRRNART + smoker",
  plus_comorbid = "~ AGE + SEX + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + smoker + MHHTN + MHT2D + MHCOPD + MHHRTDIS + MHCVD + MHABNWBC + MHLVRDIS"
)

cat("\n--- Base Model Comparison ---\n")
cat("Metric: median gene-level residual variance (lower = more variance explained)\n\n")

model_results <- list()

for (mname in names(models)) {
  form <- as.formula(models[[mname]])
  tb_sub <- tb_joint  # already complete cases for all

  design_m <- tryCatch(model.matrix(form, data = tb_sub), error = function(e) NULL)
  if (is.null(design_m)) { cat(sprintf("  %s: SKIP (design error)\n", mname)); next }

  qr_m <- qr(design_m)
  if (qr_m$rank < ncol(design_m))
    design_m <- design_m[, qr_m$pivot[seq_len(qr_m$rank)], drop = FALSE]

  v_m <- voom(dge_joint, design_m, plot = FALSE)
  fit_m <- lmFit(v_m, design_m)

  resid_var <- rowMeans(residuals(fit_m, v_m)^2)
  med_resid <- median(resid_var)
  p90_resid <- quantile(resid_var, 0.90)

  model_results[[mname]] <- tibble(
    model = mname,
    n_terms = ncol(design_m) - 1,
    median_resid_var = round(med_resid, 4),
    p90_resid_var = round(p90_resid, 4)
  )

  cat(sprintf("  %-20s: %2d terms, median resid = %.4f, p90 = %.4f\n",
              mname, ncol(design_m) - 1, med_resid, p90_resid))
}

model_comp_df <- bind_rows(model_results)
write_csv(model_comp_df, file.path(res_dir, "exploratory_model_comparison.csv"))

# ── 7b. Test whether adding disease covariates changes CKD signal ──────────
cat("\n--- CKD Signal With/Without Disease Covariates ---\n")

form_no_disease <- ~ AGE + SEX + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN +
  SMCENTER + SMEXNCRT + smoker + MHRNLFLR
form_with_disease <- ~ AGE + SEX + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN +
  SMCENTER + SMEXNCRT + smoker + MHHTN + MHT2D + MHCOPD + MHHRTDIS + MHCVD +
  MHABNWBC + MHLVRDIS + MHRNLFLR

for (form_info in list(
  list(form = form_no_disease, label = "No disease covariates"),
  list(form = form_with_disease, label = "With disease covariates")
)) {
  design_test <- model.matrix(form_info$form, data = tb_joint)
  qr_test <- qr(design_test)
  if (qr_test$rank < ncol(design_test))
    design_test <- design_test[, qr_test$pivot[seq_len(qr_test$rank)], drop = FALSE]

  v_test <- voom(dge_joint, design_test, plot = FALSE)
  fit_test <- eBayes(lmFit(v_test, design_test))
  ckd_coef <- grep("MHRNLFLR", colnames(design_test), value = TRUE)[1]
  tt_test <- topTable(fit_test, coef = ckd_coef, number = Inf, sort.by = "none")
  n_ckd <- sum(tt_test$adj.P.Val < 0.05)

  cat(sprintf("  %s: %d CKD DEGs\n", form_info$label, n_ckd))
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 8: PCA OF GENE EXPRESSION — FACTOR ASSOCIATIONS                ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 8: PCA Factor Associations\n")
cat(strrep("=", 70), "\n")

# Run PCA on logCPM, test association of top PCs with each factor
logcpm_joint <- cpm(dge_joint, log = TRUE, prior.count = 1)
pca <- prcomp(t(logcpm_joint), center = TRUE, scale. = FALSE)

# Variance explained
var_explained <- summary(pca)$importance[2, 1:20]
cat("Variance explained (top 20 PCs):\n")
cat(paste(sprintf("PC%d: %.1f%%", 1:20, 100 * var_explained), collapse = ", "), "\n")

# Test each factor against top 10 PCs
cat("\n--- Factor ~ PC associations (p-values, Kruskal-Wallis or Spearman) ---\n")

test_factors <- c("AGE", "SEX", "RACE", "BMI", "DTHHRDY", "ischemic_hrs",
                  "SMRIN", "SMCENTER", "SMEXNCRT", "smoker",
                  "MHHTN", "MHT2D", "MHRNLFLR", "MHABNWBC")

pc_scores <- pca$x[, 1:10]
rownames(pc_scores) <- tb_joint$SAMPID

pc_assoc <- matrix(NA, length(test_factors), 10,
                   dimnames = list(test_factors, paste0("PC", 1:10)))

for (fac in test_factors) {
  x <- tb_joint[[fac]]
  for (pc_i in 1:10) {
    y <- pc_scores[, pc_i]
    if (is.factor(x) || is.character(x) || (is.numeric(x) && all(x %in% c(0, 1, NA)))) {
      # Categorical: Kruskal-Wallis
      p <- tryCatch(kruskal.test(y ~ factor(x))$p.value, error = function(e) NA)
    } else {
      # Continuous: Spearman correlation test
      p <- tryCatch(cor.test(x, y, method = "spearman")$p.value, error = function(e) NA)
    }
    pc_assoc[fac, pc_i] <- p
  }
}

# Display as -log10(p), capped at 20
log_p <- -log10(pc_assoc)
log_p[log_p > 20] <- 20

cat("\nTop PC associations (-log10 p, * = p < 0.05/n_tests):\n")
bonf <- 0.05 / (length(test_factors) * 10)
cat(sprintf("Bonferroni threshold: p < %.1e\n", bonf))

for (fac in test_factors) {
  sig_pcs <- which(pc_assoc[fac, ] < bonf)
  if (length(sig_pcs) > 0) {
    cat(sprintf("  %-15s: %s\n", fac,
                paste(sprintf("PC%d(p=%.1e)", sig_pcs, pc_assoc[fac, sig_pcs]),
                      collapse = ", ")))
  }
}

# ── 8b. Figure: PC association heatmap ──────────────────────────────────────
png(file.path(fig_dir, "01_pca_associations.png"),
    width = 10, height = 8, units = "in", res = 300)
pheatmap(
  log_p,
  display_numbers = TRUE,
  number_format = "%.1f",
  fontsize_number = 7,
  color = colorRampPalette(c("white", "#FEE08B", "#F46D43", "#A50026"))(100),
  breaks = seq(0, 20, length.out = 101),
  cluster_rows = FALSE, cluster_cols = FALSE,
  main = "-log10(p-value): Factor ~ PC Association\n(Kruskal-Wallis for categorical, Spearman for continuous)"
)
dev.off()
cat("Saved:", file.path(fig_dir, "01_pca_associations.png"), "\n")


# ── 8c. Figure: PC1 vs PC2 colored by key factors ──────────────────────────
pc_df <- as_tibble(pc_scores[, 1:5]) %>%
  mutate(
    SAMPID = tb_joint$SAMPID,
    SEX = tb_joint$SEX,
    DTHHRDY = tb_joint$DTHHRDY,
    ischemic_hrs = tb_joint$ischemic_hrs,
    AGE = tb_joint$AGE,
    MHRNLFLR = factor(ifelse(tb_joint$MHRNLFLR == 1, "CKD", "No"),
                       levels = c("No", "CKD")),
    SMCENTER = tb_joint$SMCENTER
  )

p_pca1 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = ischemic_hrs)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_viridis_c() +
  labs(title = "PC1 vs PC2: Ischemic Time",
       color = "Ischemic\nhours") +
  theme_bw(base_size = 10)

p_pca2 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = DTHHRDY)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PC1 vs PC2: Hardy Scale", color = "Hardy") +
  theme_bw(base_size = 10)

p_pca3 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = SEX)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = c(Male = "steelblue", Female = "firebrick")) +
  labs(title = "PC1 vs PC2: Sex", color = "Sex") +
  theme_bw(base_size = 10)

p_pca4 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = MHRNLFLR)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = c(No = "grey70", CKD = "red")) +
  labs(title = "PC1 vs PC2: CKD", color = "CKD") +
  theme_bw(base_size = 10)

p_pca_grid <- plot_grid(p_pca1, p_pca2, p_pca3, p_pca4, ncol = 2)
ggsave(file.path(fig_dir, "01_pca_factors.png"), p_pca_grid,
       width = 12, height = 10, dpi = 300)
cat("Saved:", file.path(fig_dir, "01_pca_factors.png"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 9: SAVE & SUMMARY                                              ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n")

cat("
EXPLORATORY DE ANALYSIS — KEY FINDINGS
=======================================

1. MARGINAL FACTOR SCREENING
   See: 01_marginal_degs.png, exploratory_marginal_screening.csv
   Factors ranked by independent effect on gene expression (DEGs at FDR<0.05).

2. LINEARITY ASSESSMENT
   See: exploratory_linearity.csv, 01_ischemic_time_genes.png
   For each continuous variable, linear vs spline(df=3) comparison.
   Key question: does ischemic time need nonlinear modeling?

3. FACTOR EFFECT SIMILARITY
   See: 01_factor_correlation.png, 01_factor_clustering.png
   Spearman rho of genome-wide t-statistics reveals which factors produce
   similar/opposite expression changes. Critical for understanding confounding.

4. ISCHEMIC TIME DOSE-RESPONSE
   See: 01_ischemic_dose_response.png
   Gene expression change as function of ischemic time bins.
   Helps determine optimal binning strategy.

5. CKD SPECIFICITY
   See: 01_ckd_factor_scatter.png
   How much CKD signal overlaps with aging, comorbidities, and technical factors.

6. CONFOUNDER STRUCTURE
   Disease-demographic associations and comorbidity patterns.
   Critical for designing the DE model.

7. BASE MODEL COMPARISON
   See: exploratory_model_comparison.csv
   Residual variance under candidate models.

8. PCA FACTOR ASSOCIATIONS
   See: 01_pca_associations.png, 01_pca_factors.png
   Which factors drive the top PCs of expression variation.

Output files:
  - results/exploratory_marginal_screening.csv
  - results/exploratory_linearity.csv
  - results/exploratory_adjusted_degs.csv
  - results/exploratory_model_comparison.csv
  - figures/01_*.png
")

# Save key objects
saveRDS(list(
  marginal_df = marginal_df,
  linearity_df = linearity_df,
  adjusted_df = adjusted_df,
  model_comp_df = model_comp_df,
  cor_t = cor_t,
  t_mat = t_mat,
  pc_assoc = pc_assoc,
  phi_mat = phi_mat,
  bin_effect_df = bin_effect_df
), file.path(res_dir, "exploratory_de_results.rds"))

cat("\nDone.\n")
