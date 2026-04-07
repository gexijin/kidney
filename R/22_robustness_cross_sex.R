#!/usr/bin/env Rscript
# =============================================================================
# Stage 1e: Male Ischemic Stress Test — Cross-Sex Matched Negative Control
#
# Construct a male CKD subset that is a demographic twin of the female CKD
# cases — same age, BMI, Hardy, ischemic time, RIN — and show it produces
# ~0 DEGs. Same confounders, same power, different sex → sex-specific biology.
#
# Output: results/cross_sex_matched_de.rds, figures/22_cross_sex_match.png
# =============================================================================

library(tidyverse)
library(edgeR)
library(limma)
library(MatchIt)
library(fgsea)
library(msigdbr)
library(cowplot)

source("kidney/R/functions.R")

set.seed(42)

res_dir <- "kidney/results"
fig_dir <- "kidney/figures"


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 0: LOAD DATA                                                    ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 0: Load Data\n")
cat(strrep("=", 70), "\n")

dat <- load_sat_data(res_dir = res_dir)
tb <- dat$tb
counts <- dat$counts

res_female <- readRDS(file.path(res_dir, "de_female.rds"))
res_male   <- readRDS(file.path(res_dir, "de_male.rds"))

# Extract female and male CKD cases
female_ckd <- tb %>% filter(SEX == "Female", MHRNLFLR == 1)
male_ckd   <- tb %>% filter(SEX == "Male", MHRNLFLR == 1)
male_ctrl  <- tb %>% filter(SEX == "Male", MHRNLFLR == 0)

cat(sprintf("Female CKD: n=%d\n", nrow(female_ckd)))
cat(sprintf("Male CKD:   n=%d\n", nrow(male_ckd)))
cat(sprintf("Male ctrl:  n=%d\n", nrow(male_ctrl)))

cat("\nFemale CKD Hardy distribution:\n")
print(table(female_ckd$DTHHRDY))
cat("Male CKD Hardy distribution:\n")
print(table(male_ckd$DTHHRDY))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1: CROSS-SEX MATCHING — SELECT MALE CKD TWINS                  ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 1: Cross-Sex Matching (Mahalanobis + exact Hardy/Race)\n")
cat(strrep("=", 70), "\n")

# Pooled dataframe: female CKD (treatment=1) + male CKD (treatment=0)
# Restrict to samples complete on ALL DE model vars (not just matching vars)
# to ensure the matched set = the DE set
match_vars <- c("AGE", "BMI", "ischemic_hrs", "SMRIN")
exact_vars <- c("DTHHRDY", "RACE")
all_match_vars <- c(match_vars, exact_vars)
model_vars <- c("AGE", "RACE", "BMI", "DTHHRDY", "ischemic_hrs", "SMRIN", "SMCENTER",
                "SMEXNCRT", "MHABNWBC", "MHLVRDIS", "MHT2D", "MHRNLFLR")
de_complete_vars <- unique(c(all_match_vars, model_vars))

pool <- bind_rows(
  female_ckd %>% mutate(treat = 1),
  male_ckd   %>% mutate(treat = 0)
) %>%
  filter(complete.cases(across(all_of(de_complete_vars))))

cat(sprintf("Matching pool: %d female CKD (treat) + %d male CKD (control)\n",
            sum(pool$treat == 1), sum(pool$treat == 0)))

# Check Hardy feasibility
cat("\nHardy feasibility (female needed / male available):\n")
for (h in levels(pool$DTHHRDY)) {
  nf <- sum(pool$treat == 1 & pool$DTHHRDY == h)
  nm <- sum(pool$treat == 0 & pool$DTHHRDY == h)
  cat(sprintf("  %s: %d needed / %d available %s\n",
              h, nf, nm, ifelse(nf > nm, "*** INSUFFICIENT", "OK")))
}

# MatchIt: 1:1 nearest Mahalanobis, exact on Hardy + Race
m_out <- matchit(
  treat ~ AGE + BMI + ischemic_hrs + SMRIN,
  data = pool,
  method = "nearest",
  distance = "mahalanobis",
  exact = ~ DTHHRDY + RACE,
  ratio = 1,
  replace = FALSE
)

cat("\nMatchIt summary:\n")
print(summary(m_out))

# Extract matched male IDs
matched_data <- match.data(m_out)
matched_males <- matched_data %>% filter(treat == 0)
cat(sprintf("\nMatched males selected: n=%d\n", nrow(matched_males)))

# ── Balance table ────────────────────────────────────────────────────────
balance_vars <- c("AGE", "BMI", "ischemic_hrs", "SMRIN")
matched_females <- matched_data %>% filter(treat == 1)

balance_table <- tibble(
  variable = balance_vars,
  female_mean = sapply(balance_vars, function(v) mean(matched_females[[v]], na.rm = TRUE)),
  female_sd   = sapply(balance_vars, function(v) sd(matched_females[[v]], na.rm = TRUE)),
  male_mean   = sapply(balance_vars, function(v) mean(matched_males[[v]], na.rm = TRUE)),
  male_sd     = sapply(balance_vars, function(v) sd(matched_males[[v]], na.rm = TRUE))
) %>%
  mutate(
    pooled_sd = sqrt((female_sd^2 + male_sd^2) / 2),
    smd = (female_mean - male_mean) / pooled_sd
  )

cat("\n--- Covariate Balance: Female CKD vs Matched Male CKD ---\n")
cat(sprintf("%-15s %8s %8s %8s %8s %8s\n",
            "Variable", "F_mean", "F_sd", "M_mean", "M_sd", "SMD"))
cat(strrep("-", 60), "\n")
for (i in seq_len(nrow(balance_table))) {
  r <- balance_table[i, ]
  cat(sprintf("%-15s %8.2f %8.2f %8.2f %8.2f %8.3f %s\n",
              r$variable, r$female_mean, r$female_sd,
              r$male_mean, r$male_sd, r$smd,
              ifelse(abs(r$smd) > 0.2, "!!", ifelse(abs(r$smd) > 0.1, "!", ""))))
}

# Exact-match variables
cat("\nExact-match variables:\n")
cat("  Hardy: "); print(table(matched_males$DTHHRDY))
cat("  Race:  "); print(table(matched_males$RACE))

# Key check: matched males should have similarly elevated ischemic time
cat(sprintf("\nIschemic time (hrs):\n"))
cat(sprintf("  Female CKD: median=%.1f [%.1f-%.1f]\n",
            median(matched_females$ischemic_hrs),
            quantile(matched_females$ischemic_hrs, 0.25),
            quantile(matched_females$ischemic_hrs, 0.75)))
cat(sprintf("  Matched male CKD: median=%.1f [%.1f-%.1f]\n",
            median(matched_males$ischemic_hrs),
            quantile(matched_males$ischemic_hrs, 0.25),
            quantile(matched_males$ischemic_hrs, 0.75)))
cat(sprintf("  All male CKD: median=%.1f [%.1f-%.1f]\n",
            median(male_ckd$ischemic_hrs, na.rm = TRUE),
            quantile(male_ckd$ischemic_hrs, 0.25, na.rm = TRUE),
            quantile(male_ckd$ischemic_hrs, 0.75, na.rm = TRUE)))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2: DE ON MATCHED MALES VS MALE CONTROLS                        ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 2: DE — Matched Male CKD vs All Male Controls\n")
cat(strrep("=", 70), "\n")

# Combine matched male CKD + all male controls
# matched_males are already DE-complete (filtered in pool construction)
tb_male_test <- bind_rows(matched_males, male_ctrl)
tb_male_test <- tb_male_test[complete.cases(tb_male_test[, model_vars]), ]
counts_test <- counts[, tb_male_test$SAMPID]

cat(sprintf("Samples: %d (%d matched CKD + %d controls)\n",
            nrow(tb_male_test),
            sum(tb_male_test$MHRNLFLR == 1),
            sum(tb_male_test$MHRNLFLR == 0)))

# limma/voom — same model as stage1_01
form <- as.formula(paste("~", paste(model_vars[-which(model_vars == "MHRNLFLR")],
                                     collapse = " + "), "+ MHRNLFLR"))
matched_fit <- run_limma_analysis(
  tb_male_test, counts, form, "MHRNLFLR",
  coef_match = "prefix",
  context = "cross-sex matched male CKD model"
)
tt_matched <- matched_fit$tt
n_deg_matched <- matched_fit$n_deg
n_up_matched  <- matched_fit$n_up
n_dn_matched  <- matched_fit$n_down

cat(sprintf("\nMatched male CKD DEGs (FDR<0.05): %d (%d up, %d down)\n",
            n_deg_matched, n_up_matched, n_dn_matched))
cat(sprintf("For comparison — Female CKD: %d DEGs, Full male CKD: %d DEGs\n",
            res_female$n_deg, res_male$n_deg))

# Nominal p < 0.05 count for reference
n_nom <- sum(tt_matched$P.Value < 0.05)
cat(sprintf("Nominal p<0.05: %d (%.1f%% of tested genes — expect ~5%%)\n",
            n_nom, 100 * n_nom / nrow(tt_matched)))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 3: COMPARISON WITH FEMALE RESULTS                               ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 3: Comparison with Female CKD Results\n")
cat(strrep("=", 70), "\n")

# ── t-statistic correlation ──────────────────────────────────────────────
tt_female <- res_female$tt
de_conc <- compare_de_results(tt_female, tt_matched)
merged <- de_conc$merged %>%
  rename(
    t_female = t_ref, logFC_female = logFC_ref, padj_female = padj_ref,
    t_matched = t_new, logFC_matched = logFC_new, padj_matched = padj_new
  )
rho_t <- de_conc$rho_t
rho_lfc <- de_conc$rho_lfc
cat(sprintf("t-statistic correlation (female vs matched male): rho = %.3f\n", rho_t))
cat(sprintf("logFC correlation: rho = %.3f\n", rho_lfc))

# ── Pathway GSEA ─────────────────────────────────────────────────────────
hallmark_list <- load_hallmark_sets()

# Female GSEA
gsea_female <- run_gsea(tt_female, hallmark_list)

# Matched male GSEA
gsea_matched <- run_gsea(tt_matched, hallmark_list)

# Compare NES
gsea_cmp_res <- compare_gsea_results(
  gsea_female, gsea_matched,
  label1 = "female", label2 = "matched",
  core_paths = NULL, print_table = FALSE
)
gsea_compare <- gsea_cmp_res$merged
rho_nes <- gsea_cmp_res$rho_nes
cat(sprintf("\nPathway NES correlation (female vs matched male): rho = %.3f\n", rho_nes))

cat("\nFemale-significant pathways in matched males:\n")
fem_sig <- gsea_compare %>% filter(padj_female < 0.05) %>% arrange(NES_female)
cat(sprintf("%-50s %7s %7s %7s\n", "Pathway", "NES_F", "NES_M", "padj_M"))
cat(strrep("-", 75), "\n")
for (i in seq_len(nrow(fem_sig))) {
  r <- fem_sig[i, ]
  short <- gsub("^HALLMARK_", "", r$pathway)
  cat(sprintf("%-50s %7.2f %7.2f %7.3f %s\n",
              short, r$NES_female, r$NES_matched, r$padj_matched,
              ifelse(r$padj_matched < 0.05, "*", "")))
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 4: ROBUSTNESS — REPEATED MATCHING                               ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 4: Robustness — Leave-One-Out Re-Matching\n")
cat(strrep("=", 70), "\n")

# The Slow Hardy bin forces all 17/17 males to be selected, making the primary
# match deterministic. To genuinely test robustness, we drop 1 random female
# case per iteration and re-match, producing distinct male matched sets.

n_iters <- 10
iter_results <- vector("list", n_iters)

# Get DE-complete female CKD from the pool
pool_females <- pool %>% filter(treat == 1)
pool_males   <- pool %>% filter(treat == 0)

for (iter in seq_len(n_iters)) {
  set.seed(42 + iter)

  # Drop 1 random female CKD case
  drop_idx <- sample(nrow(pool_females), 1)
  pool_iter <- bind_rows(
    pool_females[-drop_idx, ],
    pool_males
  )

  m_iter <- tryCatch(
    matchit(
      treat ~ AGE + BMI + ischemic_hrs + SMRIN,
      data = pool_iter,
      method = "nearest",
      distance = "mahalanobis",
      exact = ~ DTHHRDY + RACE,
      ratio = 1,
      replace = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(m_iter)) {
    cat(sprintf("  Iteration %d: matching failed (dropped %s)\n",
                iter, pool_females$SAMPID[drop_idx]))
    next
  }

  md_iter <- match.data(m_iter)
  matched_m_iter <- md_iter %>% filter(treat == 0)

  # DE
  tb_iter <- bind_rows(matched_m_iter, male_ctrl)
  tb_iter <- tb_iter[complete.cases(tb_iter[, model_vars]), ]
  counts_iter <- counts[, tb_iter$SAMPID]

  design_iter <- build_design(form, tb_iter)
  de_iter <- run_voom_de(counts_iter, design_iter,
                          get_required_coef_prefix(design_iter, "MHRNLFLR",
                                                    sprintf("iter %d", iter)))
  tt_iter <- de_iter$tt

  n_deg_iter <- sum(tt_iter$adj.P.Val < 0.05)
  n_selected <- nrow(matched_m_iter)
  overlap_with_primary <- length(intersect(matched_m_iter$SAMPID, matched_males$SAMPID))

  iter_results[[iter]] <- tibble(
    iteration = iter,
    n_matched = n_selected,
    n_deg = n_deg_iter,
    overlap_primary = overlap_with_primary
  )

  cat(sprintf("  Iteration %2d: %d matched, %d DEGs, %d/%d overlap with primary\n",
              iter, n_selected, n_deg_iter, overlap_with_primary, n_selected))
}

iter_df <- bind_rows(iter_results)
cat(sprintf("\nAcross %d iterations:\n", nrow(iter_df)))
cat(sprintf("  DEGs: median=%d, range=[%d, %d]\n",
            median(iter_df$n_deg), min(iter_df$n_deg), max(iter_df$n_deg)))
cat(sprintf("  Overlap with primary: median=%d/%d\n",
            median(iter_df$overlap_primary), iter_df$n_matched[1]))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 5: FIGURE                                                       ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 5: Figure\n")
cat(strrep("=", 70), "\n")

# ── Panel A: Covariate balance (lollipop plot) ───────────────────────────
balance_plot_df <- balance_table %>%
  mutate(variable = fct_reorder(variable, abs(smd)))

p_balance <- ggplot(balance_plot_df, aes(x = abs(smd), y = variable)) +
  geom_segment(aes(x = 0, xend = abs(smd), y = variable, yend = variable),
               color = "grey50", linewidth = 0.8) +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 0.2, linetype = "dashed", color = "red") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "|Standardized Mean Difference|",
       y = NULL,
       title = "Covariate balance",
       subtitle = "Female CKD vs matched male CKD") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

# ── Panel B: Volcano side-by-side ────────────────────────────────────────
vol_female <- tt_female %>%
  mutate(sig = adj.P.Val < 0.05,
         group = "Female CKD\n(all cases)")

vol_matched <- tt_matched %>%
  mutate(sig = adj.P.Val < 0.05,
         group = sprintf("Matched male CKD\n(n=%d, ischemic-twin)", nrow(matched_males)))

vol_df <- bind_rows(vol_female, vol_matched) %>%
  mutate(neg_log10p = -log10(P.Value))

# Cap -log10p for plotting
vol_df$neg_log10p <- pmin(vol_df$neg_log10p, 20)

p_volcano <- ggplot(vol_df, aes(x = logFC, y = neg_log10p)) +
  geom_point(aes(color = sig), size = 0.4, alpha = 0.3) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                     labels = c("NS", "FDR<0.05"), name = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.3) +
  facet_wrap(~group, scales = "free_x") +
  labs(x = "log2 fold change (CKD vs control)",
       y = expression(-log[10](p)),
       title = "CKD differential expression") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 12))

# ── Panel C: Pathway NES comparison ──────────────────────────────────────
# Show top female pathways
top_paths <- gsea_compare %>%
  filter(padj_female < 0.05) %>%
  arrange(NES_female) %>%
  mutate(pathway_short = gsub("^HALLMARK_", "", pathway),
         pathway_short = gsub("_", " ", pathway_short),
         pathway_short = str_to_title(pathway_short))

path_plot_df <- top_paths %>%
  pivot_longer(cols = c(NES_female, NES_matched),
               names_to = "group", values_to = "NES") %>%
  mutate(group = ifelse(group == "NES_female", "Female CKD", "Matched male CKD"),
         pathway_short = fct_reorder(pathway_short, NES))

p_pathway <- ggplot(path_plot_df, aes(x = NES, y = pathway_short, color = group)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
  scale_color_manual(values = c("Female CKD" = "firebrick", "Matched male CKD" = "steelblue"),
                     name = NULL) +
  labs(x = "Normalized Enrichment Score",
       y = NULL,
       title = "Hallmark pathway enrichment",
       subtitle = sprintf("NES correlation: rho = %.2f", rho_nes)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12))

# ── Combine ──────────────────────────────────────────────────────────────
p_top <- plot_grid(p_balance, p_volcano, labels = c("A", "B"),
                   rel_widths = c(0.35, 0.65), nrow = 1)
p_combined <- plot_grid(p_top, p_pathway, labels = c("", "C"),
                        ncol = 1, rel_heights = c(0.45, 0.55))

ggsave(file.path(fig_dir, "22_cross_sex_match.png"), p_combined,
       width = 14, height = 12, dpi = 300)
cat("Saved:", file.path(fig_dir, "22_cross_sex_match.png"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SAVE RESULTS                                                            ║
# ╚════════════════════════════════════════════════════════════════════════════╝

saveRDS(list(
  matched_males = matched_males,
  matched_females = matched_females,
  balance_table = balance_table,
  tt_matched = tt_matched,
  n_deg_matched = n_deg_matched,
  gsea_female = gsea_female,
  gsea_matched = gsea_matched,
  gsea_compare = gsea_compare,
  rho_t = rho_t,
  rho_nes = rho_nes,
  iter_results = iter_df,
  matchit_obj = m_out
), file.path(res_dir, "cross_sex_matched_de.rds"))

cat("\nSaved:", file.path(res_dir, "cross_sex_matched_de.rds"), "\n")
cat("Done.\n")
