#!/usr/bin/env Rscript
# =============================================================================
# Stage 1f: HTN–CKD Robustness Analyses
#
# Demonstrates that CKD drives the female SAT transcriptomic signal
# independently of hypertension. Four complementary analyses:
#
#   1. HTN-only DE in females (excluding CKD cases)
#   2. CKD+HTN vs HTN-only (within-HTN stratum)
#   3. CKD without HTN (pathway concordance with full model)
#   4. Variance partitioning — CKD vs HTN in joint model
#
# Output: results/htn_robustness_*.rds, figures/23_htn_robustness.png
# =============================================================================

library(tidyverse)
library(edgeR)
library(limma)
library(fgsea)
library(msigdbr)

source("kidney/R/functions.R")


set.seed(42)

res_dir <- "kidney/results"
fig_dir <- "kidney/figures"


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 0: LOAD DATA                                                   ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 0: Load Data\n")
cat(strrep("=", 70), "\n")

dat <- load_sat_data(
  res_dir = res_dir,
  mh_flags = c("MHABNWBC", "MHLVRDIS", "MHRNLFLR", "MHT2D", "MHHTN"),
  complete_vars = c("AGE", "RACE", "BMI", "DTHHRDY", "ischemic_hrs", "SMRIN",
                    "SMCENTER", "SMEXNCRT", "MHABNWBC", "MHLVRDIS", "MHT2D",
                    "MHRNLFLR", "MHHTN")
)
tb <- dat$tb
counts <- dat$counts

# Load existing female CKD results for concordance
res_female_orig <- readRDS(file.path(res_dir, "de_female.rds"))

# Female subset
tbf <- tb %>% filter(SEX == "Female")
cat(sprintf("Female samples: %d\n", nrow(tbf)))
cat(sprintf("  CKD+HTN: %d | CKD-noHTN: %d | noCKD+HTN: %d | noCKD-noHTN: %d\n",
            sum(tbf$MHRNLFLR == 1 & tbf$MHHTN == 1),
            sum(tbf$MHRNLFLR == 1 & tbf$MHHTN == 0),
            sum(tbf$MHRNLFLR == 0 & tbf$MHHTN == 1),
            sum(tbf$MHRNLFLR == 0 & tbf$MHHTN == 0)))

# Hallmark gene sets
hallmark_list <- load_hallmark_sets()

# Helper: limma/voom pipeline
run_limma <- function(tb_sub, counts_all, formula, coef_name) {
  res <- tryCatch(
    run_limma_analysis(
      tb_sub, counts_all, formula, coef_name,
      coef_match = "prefix",
      context = sprintf("23_robustness_htn_ckd: %s", coef_name)
    ),
    error = function(e) NULL
  )
  if (is.null(res)) {
    cat(sprintf("  WARNING: %s not found in design columns: %s\n",
                coef_name, paste(colnames(build_design(formula, tb_sub)), collapse = ", ")))
    return(NULL)
  }
  res
}

# Helper: fGSEA from topTable
run_hallmark_gsea <- function(tt, pathways = hallmark_list) run_gsea(tt, pathways)

# Helper: concordance with full-model female CKD
calc_concordance <- function(tt_new, tt_ref = res_female_orig$tt) {
  conc <- compare_de_results(tt_ref, tt_new)
  conc$dir_concordance <- conc$direction_concordance
  conc
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  ANALYSIS 1: HTN-ONLY DE IN FEMALES (EXCLUDING CKD)                     ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS 1: HTN-Only DE in Females (excluding CKD cases)\n")
cat(strrep("=", 70), "\n")

# Exclude all CKD cases → test HTN among CKD-free females
tb_a1 <- tbf %>% filter(MHRNLFLR == 0)
cat(sprintf("  Samples: %d (HTN=%d, noHTN=%d)\n",
            nrow(tb_a1), sum(tb_a1$MHHTN == 1), sum(tb_a1$MHHTN == 0)))

form_a1 <- ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +
  MHABNWBC + MHLVRDIS + MHT2D + MHHTN

res_a1 <- run_limma(tb_a1, counts, form_a1, "MHHTN")
cat(sprintf("  DEGs (FDR<0.05): %d (%d up, %d down)\n",
            res_a1$n_deg, res_a1$n_up, res_a1$n_down))

gsea_a1 <- run_hallmark_gsea(res_a1$tt)
conc_a1 <- calc_concordance(res_a1$tt)
cat(sprintf("  t-stat rho vs full CKD: %.3f\n", conc_a1$rho_t))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  ANALYSIS 2: CKD+HTN vs HTN-ONLY (WITHIN-HTN STRATUM)                   ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS 2: CKD+HTN vs HTN-Only (within-HTN stratum)\n")
cat(strrep("=", 70), "\n")

# Restrict to females with HTN → test CKD within HTN
tb_a2 <- tbf %>% filter(MHHTN == 1)
cat(sprintf("  Samples: %d (CKD=%d, noCKD=%d)\n",
            nrow(tb_a2), sum(tb_a2$MHRNLFLR == 1), sum(tb_a2$MHRNLFLR == 0)))

# No MHHTN covariate — constant within stratum
form_a2 <- ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +
  MHABNWBC + MHLVRDIS + MHT2D + MHRNLFLR

res_a2 <- run_limma(tb_a2, counts, form_a2, "MHRNLFLR")
cat(sprintf("  DEGs (FDR<0.05): %d (%d up, %d down)\n",
            res_a2$n_deg, res_a2$n_up, res_a2$n_down))

gsea_a2 <- run_hallmark_gsea(res_a2$tt)
conc_a2 <- calc_concordance(res_a2$tt)
cat(sprintf("  t-stat rho vs full CKD: %.3f\n", conc_a2$rho_t))
cat(sprintf("  DEG overlap with full model: %d / %d (full) / %d (within-HTN)\n",
            conc_a2$overlap, conc_a2$n_deg_ref, conc_a2$n_deg_new))
if (!is.na(conc_a2$dir_concordance))
  cat(sprintf("  Direction concordance: %.1f%%\n", 100 * conc_a2$dir_concordance))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  ANALYSIS 3: CKD WITHOUT HTN (PATHWAY CONCORDANCE)                      ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS 3: CKD Without HTN (pathway concordance)\n")
cat(strrep("=", 70), "\n")

# Restrict to females without HTN → test CKD (n=6 vs n=88)
tb_a3 <- tbf %>% filter(MHHTN == 0)
cat(sprintf("  Samples: %d (CKD=%d, noCKD=%d)\n",
            nrow(tb_a3), sum(tb_a3$MHRNLFLR == 1), sum(tb_a3$MHRNLFLR == 0)))

# Reduced covariates — only 6 cases, can't support 13+ parameters
form_a3 <- ~ AGE + RACE + DTHHRDY + ischemic_hrs + SMRIN + MHRNLFLR

res_a3 <- run_limma(tb_a3, counts, form_a3, "MHRNLFLR")
cat(sprintf("  DEGs (FDR<0.05): %d (expected: few, given n=6 cases)\n", res_a3$n_deg))

gsea_a3 <- run_hallmark_gsea(res_a3$tt)
conc_a3 <- calc_concordance(res_a3$tt)
cat(sprintf("  t-stat rho vs full CKD: %.3f\n", conc_a3$rho_t))

# GSEA NES concordance with full model
gsea_full <- run_hallmark_gsea(res_female_orig$tt)
cmp_a3 <- compare_gsea_results(gsea_full, gsea_a3,
                                label1 = "full", label2 = "noHTN",
                                print_table = FALSE)
rho_nes_a3 <- cmp_a3$rho_nes
cat(sprintf("  GSEA NES Spearman rho (full vs CKD-noHTN): %.3f\n", rho_nes_a3))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  ANALYSIS 4: VARIANCE PARTITIONING (CKD vs HTN IN JOINT MODEL)          ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS 4: Variance Partitioning — CKD vs HTN in Joint Model\n")
cat(strrep("=", 70), "\n")

# Full female cohort, both MHRNLFLR and MHHTN in the same model
cat(sprintf("  Full female cohort: %d samples\n", nrow(tbf)))

form_a4 <- ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +
  MHABNWBC + MHLVRDIS + MHT2D + MHRNLFLR + MHHTN

# Fit once, extract both coefficients
counts_f <- counts[, tbf$SAMPID]
design_a4 <- build_design(form_a4, tbf)
de_a4 <- run_voom_de(counts_f, design_a4,
                      get_required_coef_prefix(design_a4, "MHRNLFLR", "Analysis 4 CKD"))
tt_ckd <- de_a4$tt

# HTN coefficient from same fit
htn_coef <- get_required_coef_prefix(design_a4, "MHHTN", "Analysis 4 HTN")
tt_htn <- topTable(de_a4$fit, coef = htn_coef, number = Inf, sort.by = "none")
tt_htn$gene_id <- rownames(tt_htn)

deg_ckd <- summarize_degs(tt_ckd)
deg_htn <- summarize_degs(tt_htn)
cat(sprintf("  CKD DEGs (FDR<0.05): %d (%d up, %d down)\n",
            deg_ckd$n_deg, deg_ckd$n_up, deg_ckd$n_down))
cat(sprintf("  HTN DEGs (FDR<0.05): %d (%d up, %d down)\n",
            deg_htn$n_deg, deg_htn$n_up, deg_htn$n_down))

# t-stat comparison
merged_a4 <- inner_join(
  tt_ckd %>% dplyr::select(gene_id, t_CKD = t, logFC_CKD = logFC, padj_CKD = adj.P.Val),
  tt_htn %>% dplyr::select(gene_id, t_HTN = t, logFC_HTN = logFC, padj_HTN = adj.P.Val),
  by = "gene_id"
)
rho_t_a4 <- cor(merged_a4$t_CKD, merged_a4$t_HTN, method = "spearman")
cat(sprintf("  CKD vs HTN t-stat Spearman rho: %.3f\n", rho_t_a4))

# Median absolute t-statistics
cat(sprintf("  Median |t|: CKD = %.2f, HTN = %.2f\n",
            median(abs(merged_a4$t_CKD)), median(abs(merged_a4$t_HTN))))

# GSEA for HTN coefficient
gsea_htn <- run_hallmark_gsea(tt_htn)
cat(sprintf("  HTN significant pathways (FDR<0.05): %d\n",
            sum(gsea_htn$padj < 0.05)))

# CKD concordance with full model (without MHHTN)
conc_a4_ckd <- calc_concordance(tt_ckd)
cat(sprintf("  CKD (joint) vs CKD (full model) t-stat rho: %.3f\n", conc_a4_ckd$rho_t))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 5: SUMMARY TABLE                                               ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n")

summary_df <- tibble(
  analysis = c("Full model (reference)",
               "1: HTN-only (no CKD)",
               "2: CKD within-HTN",
               "3: CKD without HTN",
               "4: CKD (joint model)",
               "4: HTN (joint model)"),
  predictor = c("CKD", "HTN", "CKD", "CKD", "CKD", "HTN"),
  n_samples = c(res_female_orig$n_samples, res_a1$n_samples, res_a2$n_samples,
                res_a3$n_samples, nrow(tbf), nrow(tbf)),
  n_cases = c(res_female_orig$n_cases,
              sum(tb_a1$MHHTN == 1), sum(tb_a2$MHRNLFLR == 1),
              sum(tb_a3$MHRNLFLR == 1),
              sum(tbf$MHRNLFLR == 1), sum(tbf$MHHTN == 1)),
  n_deg = c(res_female_orig$n_deg, res_a1$n_deg, res_a2$n_deg, res_a3$n_deg,
            deg_ckd$n_deg, deg_htn$n_deg),
  rho_t_vs_full = c(1.0, conc_a1$rho_t, conc_a2$rho_t, conc_a3$rho_t,
                     conc_a4_ckd$rho_t, NA)
)

cat("\n")
cat(sprintf("%-28s %-5s %4s %4s %6s %6s\n",
            "Analysis", "Pred", "n", "case", "DEGs", "rho_t"))
cat(strrep("-", 65), "\n")
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  rho_str <- ifelse(is.na(r$rho_t_vs_full), "  —", sprintf("%.3f", r$rho_t_vs_full))
  cat(sprintf("%-28s %-5s %4d %4d %6d %6s\n",
              r$analysis, r$predictor, r$n_samples, r$n_cases, r$n_deg, rho_str))
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 6: FIGURE                                                      ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 6: Figure\n")
cat(strrep("=", 70), "\n")

deg_bar_df <- tibble(
  analysis = factor(c("Full CKD", "HTN-only", "CKD within-HTN", "CKD no-HTN",
                       "Joint: CKD", "Joint: HTN"),
                     levels = c("HTN-only", "Joint: HTN", "CKD no-HTN",
                                "CKD within-HTN", "Joint: CKD", "Full CKD")),
  n_deg = c(res_female_orig$n_deg, res_a1$n_deg, res_a2$n_deg, res_a3$n_deg,
            deg_ckd$n_deg, deg_htn$n_deg),
  predictor = c("CKD", "HTN", "CKD", "CKD", "CKD", "HTN")
)

p <- ggplot(deg_bar_df, aes(x = analysis, y = n_deg, fill = predictor)) +
  geom_col(width = 0.7) +
  # Color-coded labels for zero-height bars where fill is invisible
  geom_text(aes(label = n_deg, color = predictor), hjust = -0.15, size = 3.5,
            show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = c(CKD = "#B2182B", HTN = "#4393C3"), name = "Predictor") +
  scale_color_manual(values = c(CKD = "#B2182B", HTN = "#4393C3")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "DEGs (FDR < 0.05)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "23_htn_robustness.png"), p,
       width = 6, height = 4, dpi = 300)
cat("Saved:", file.path(fig_dir, "23_htn_robustness.png"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 7: SAVE RESULTS                                                ║
# ╚════════════════════════════════════════════════════════════════════════════╝

saveRDS(list(
  a1_htn_only   = list(tt = res_a1$tt, n_deg = res_a1$n_deg, gsea = gsea_a1,
                        concordance = conc_a1, n_samples = res_a1$n_samples),
  a2_within_htn = list(tt = res_a2$tt, n_deg = res_a2$n_deg, gsea = gsea_a2,
                        concordance = conc_a2, n_samples = res_a2$n_samples),
  a3_ckd_no_htn = list(tt = res_a3$tt, n_deg = res_a3$n_deg, gsea = gsea_a3,
                        concordance = conc_a3, rho_nes = rho_nes_a3,
                        n_samples = res_a3$n_samples),
  a4_joint      = list(tt_ckd = tt_ckd, tt_htn = tt_htn,
                        n_deg_ckd = deg_ckd$n_deg, n_deg_htn = deg_htn$n_deg,
                        rho_t_ckd_htn = rho_t_a4, gsea_htn = gsea_htn,
                        concordance_ckd = conc_a4_ckd, merged = merged_a4),
  gsea_full = gsea_full,
  summary   = summary_df
), file.path(res_dir, "htn_robustness_results.rds"))
cat("Saved:", file.path(res_dir, "htn_robustness_results.rds"), "\n")

cat("\nDone.\n")
