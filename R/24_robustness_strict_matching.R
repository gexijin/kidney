#!/usr/bin/env Rscript
# =============================================================================
# Stage 1g: Strict Matching Robustness and Hallmark NES Stability
#
# Evaluates whether the female SAT CKD signal survives stricter common-support
# constraints than the primary ischemia-bin PSM analysis. Compares:
#   1) Paper-style nearest-neighbor matching exact on ischemic bins only
#   2) Stricter nearest-neighbor matching exact on Hardy death mode + ischemic
#      bins
# across ratios 1:1, 1:2, 1:3, with Hallmark GSEA/NES correlation versus the
# primary female SAT model.
#
# Outputs:
#   kidney/results/strict_matching_*.{rds,csv}
#   kidney/figures/24_*.png
# =============================================================================

library(tidyverse)
library(MatchIt)
library(edgeR)
library(limma)
library(msigdbr)
library(fgsea)
library(cowplot)

source("kidney/R/functions.R")

set.seed(42)

res_dir <- "kidney/results"
fig_dir <- "kidney/figures"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

SEEDS <- 1:10
RATIOS <- c(1, 2, 3)


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 0: LOAD DATA                                                    ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("SECTION 0: Load Data")

dat <- load_sat_data(
  res_dir = res_dir,
  mh_flags = c("MHABNWBC", "MHLVRDIS", "MHRNLFLR", "MHT2D"),
  complete_vars = c("AGE", "RACE", "BMI", "DTHHRDY", "ischemic_hrs", "SMRIN",
                    "SMCENTER", "SMEXNCRT", "MHABNWBC", "MHLVRDIS",
                    "MHRNLFLR", "MHT2D")
)

tb <- dat$tb
counts <- dat$counts
joint_factors <- dat$joint_factors

res_female <- readRDS(file.path(res_dir, "de_female.rds"))
hallmark_list <- load_hallmark_sets(filter_empty = TRUE)

tb_f <- tb %>%
  filter(SEX == "Female") %>%
  add_ischemic_bin()

tb_f$hardy_isch <- interaction(tb_f$DTHHRDY, tb_f$isch_bin, drop = TRUE)

cat(sprintf("Female SAT samples: %d (%d CKD, %d controls)\n",
            nrow(tb_f), sum(tb_f$MHRNLFLR == 1), sum(tb_f$MHRNLFLR == 0)))
cat("Ischemic bins by CKD:\n")
print(table(tb_f$isch_bin, tb_f$MHRNLFLR))
cat("Hardy x ischemic bins by CKD:\n")
print(table(tb_f$DTHHRDY, tb_f$isch_bin, tb_f$MHRNLFLR))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1: PRIMARY FEMALE REFERENCE                                     ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("SECTION 1: Primary Female Reference")

primary_formula <- as.formula(paste(
  "~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
  paste(joint_factors, collapse = " + ")
))

primary_fit <- run_limma_analysis(
  tb_f, counts, primary_formula, "MHRNLFLR",
  coef_match = "exact",
  context = "24 primary female SAT model"
)
primary_tt <- primary_fit$tt
primary_gsea <- run_gsea(primary_tt, hallmark_list, "Primary female SAT", seed = 42)

cat(sprintf("Primary female SAT DEGs: %d\n", primary_fit$n_deg))
cat(sprintf("Primary Hallmark pathways significant (FDR<0.05): %d\n",
            sum(primary_gsea$padj < 0.05)))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2: MATCHING HELPERS                                             ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("SECTION 2: Matching Helpers")

numeric_smd <- function(x, g) {
  g <- as.integer(g)
  m1 <- mean(x[g == 1], na.rm = TRUE)
  m0 <- mean(x[g == 0], na.rm = TRUE)
  s1 <- stats::var(x[g == 1], na.rm = TRUE)
  s0 <- stats::var(x[g == 0], na.rm = TRUE)
  s_pool <- sqrt((s1 + s0) / 2)
  if (!is.finite(s_pool) || s_pool == 0) return(NA_real_)
  (m1 - m0) / s_pool
}

binary_smd <- function(x, g) {
  x <- as.numeric(x)
  g <- as.integer(g)
  p1 <- mean(x[g == 1], na.rm = TRUE)
  p0 <- mean(x[g == 0], na.rm = TRUE)
  p <- (p1 + p0) / 2
  if (!is.finite(p) || p %in% c(0, 1)) return(NA_real_)
  (p1 - p0) / sqrt(p * (1 - p))
}

compute_balance_summary <- function(tb_matched) {
  g <- as.integer(tb_matched$MHRNLFLR == 1)
  tibble(
    smd_age = numeric_smd(tb_matched$AGE, g),
    smd_bmi = numeric_smd(tb_matched$BMI, g),
    smd_ischemic = numeric_smd(tb_matched$ischemic_hrs, g),
    smd_rin = numeric_smd(tb_matched$SMRIN, g),
    smd_exonic = numeric_smd(tb_matched$SMEXNCRT, g),
    smd_slow = binary_smd(tb_matched$DTHHRDY == "Slow", g),
    smd_vent = binary_smd(tb_matched$DTHHRDY == "Ventilator", g)
  )
}

run_match_scenario <- function(tb_base, counts_all, strategy, ratio, seed) {
  set.seed(seed)
  tb_in <- tb_base[sample(nrow(tb_base)), , drop = FALSE]

  if (strategy == "exact_hardy") {
    match_obj <- matchit(
      MHRNLFLR ~ AGE + BMI + ischemic_hrs + SMRIN + SMEXNCRT,
      data = tb_in,
      method = "nearest",
      ratio = ratio,
      replace = FALSE,
      caliper = 0.25,
      std.caliper = TRUE,
      m.order = "random",
      exact = ~DTHHRDY + isch_bin
    )

    model_formula <- as.formula(paste(
      "~ AGE + BMI + DTHHRDY + isch_bin + ischemic_hrs + SMRIN + SMEXNCRT +",
      paste(setdiff(joint_factors, "MHRNLFLR"), collapse = " + "),
      "+ MHRNLFLR"
    ))
  } else {
    match_obj <- matchit(
      MHRNLFLR ~ AGE + RACE + DTHHRDY + ischemic_hrs + SMRIN,
      data = tb_in,
      method = "nearest",
      ratio = ratio,
      replace = FALSE,
      caliper = 0.25,
      std.caliper = TRUE,
      m.order = "random",
      exact = ~isch_bin
    )

    model_formula <- as.formula(paste(
      "~ AGE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
      paste(setdiff(joint_factors, "MHRNLFLR"), collapse = " + "),
      "+ MHRNLFLR"
    ))
  }

  tb_matched <- match.data(match_obj)
  fit <- run_limma_analysis(
    tb_matched, counts_all, model_formula, "MHRNLFLR",
    coef_match = "exact",
    context = sprintf("24 %s ratio=%d seed=%d", strategy, ratio, seed)
  )
  gsea <- run_gsea(
    fit$tt, hallmark_list,
    label = sprintf("%s ratio=%d seed=%d", strategy, ratio, seed),
    seed = seed
  )
  gsea_cmp <- compare_gsea_results(
    primary_gsea, gsea,
    label1 = "primary", label2 = "matched",
    core_paths = CORE_PATHWAYS,
    print_table = FALSE
  )

  list(
    match_obj = match_obj,
    tb_matched = tb_matched,
    fit = fit,
    gsea = gsea,
    gsea_cmp = gsea_cmp,
    balance = compute_balance_summary(tb_matched)
  )
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 3: RATIO SENSITIVITY                                            ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("SECTION 3: Ratio Sensitivity")

scenario_grid <- expand_grid(
  strategy = c("exact_hardy", "paper_style"),
  ratio = RATIOS,
  seed = SEEDS
)

scenario_results <- vector("list", nrow(scenario_grid))
iteration_rows <- vector("list", nrow(scenario_grid))

for (i in seq_len(nrow(scenario_grid))) {
  strategy_i <- scenario_grid$strategy[i]
  ratio_i <- scenario_grid$ratio[i]
  seed_i <- scenario_grid$seed[i]

  cat(sprintf("\nRunning %s ratio=%d seed=%d\n", strategy_i, ratio_i, seed_i))
  out <- run_match_scenario(tb_f, counts, strategy_i, ratio_i, seed_i)
  scenario_results[[i]] <- out

  iteration_rows[[i]] <- bind_cols(
    tibble(
      strategy = strategy_i,
      ratio = ratio_i,
      seed = seed_i,
      n_cases = sum(out$tb_matched$MHRNLFLR == 1),
      n_controls = sum(out$tb_matched$MHRNLFLR == 0),
      degs = out$fit$n_deg,
      n_hallmark_sig = sum(out$gsea$padj < 0.05),
      n_hallmark_sig_both = sum(out$gsea_cmp$merged$padj_primary < 0.05 &
                                  out$gsea_cmp$merged$padj_matched < 0.05),
      rho_nes = out$gsea_cmp$rho_nes
    ),
    out$balance
  )
}

iteration_tbl <- bind_rows(iteration_rows)
summary_tbl <- iteration_tbl %>%
  group_by(strategy, ratio) %>%
  summarise(
    median_cases = median(n_cases),
    median_controls = median(n_controls),
    median_degs = median(degs),
    min_degs = min(degs),
    max_degs = max(degs),
    median_n_hallmark_sig = median(n_hallmark_sig),
    median_n_hallmark_sig_both = median(n_hallmark_sig_both),
    median_rho_nes = median(rho_nes),
    min_rho_nes = min(rho_nes),
    max_rho_nes = max(rho_nes),
    median_abs_smd_bmi = median(abs(smd_bmi)),
    median_abs_smd_ischemic = median(abs(smd_ischemic)),
    median_abs_smd_slow = median(abs(smd_slow)),
    .groups = "drop"
  )

print(summary_tbl)


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 4: REPRESENTATIVE HALLMARK TABLE                                ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("SECTION 4: Representative Hallmark Table")

get_scenario_index <- function(strategy, ratio, seed) {
  which(scenario_grid$strategy == strategy &
          scenario_grid$ratio == ratio &
          scenario_grid$seed == seed)
}

seed_ref <- 1
rep_exact_11 <- scenario_results[[get_scenario_index("exact_hardy", 1, seed_ref)]]
rep_exact_12 <- scenario_results[[get_scenario_index("exact_hardy", 2, seed_ref)]]
rep_exact_13 <- scenario_results[[get_scenario_index("exact_hardy", 3, seed_ref)]]
rep_paper_13 <- scenario_results[[get_scenario_index("paper_style", 3, seed_ref)]]

seed1_nes <- list(
  primary_gsea %>% transmute(pathway, primary_NES = NES, primary_padj = padj),
  rep_exact_11$gsea %>% transmute(pathway, exact11_NES = NES, exact11_padj = padj),
  rep_exact_12$gsea %>% transmute(pathway, exact12_NES = NES, exact12_padj = padj),
  rep_exact_13$gsea %>% transmute(pathway, exact13_NES = NES, exact13_padj = padj),
  rep_paper_13$gsea %>% transmute(pathway, paper13_NES = NES, paper13_padj = padj)
) %>%
  reduce(full_join, by = "pathway") %>%
  arrange(primary_NES)

core_seed1 <- seed1_nes %>%
  filter(pathway %in% CORE_PATHWAYS)

print(core_seed1)


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 5: SAVE TABLES                                                  ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("SECTION 5: Save Tables")

saveRDS(
  list(
    primary_fit = primary_fit,
    primary_gsea = primary_gsea,
    scenario_grid = scenario_grid,
    iteration_tbl = iteration_tbl,
    summary_tbl = summary_tbl,
    seed1_nes = seed1_nes,
    core_seed1 = core_seed1
  ),
  file.path(res_dir, "strict_matching_sensitivity.rds")
)

write_csv(iteration_tbl, file.path(res_dir, "strict_matching_iterations.csv"))
write_csv(summary_tbl, file.path(res_dir, "strict_matching_summary.csv"))
write_csv(seed1_nes, file.path(res_dir, "strict_matching_seed1_hallmark.csv"))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 6: FIGURES                                                      ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("SECTION 6: Figures")

supp_fig_dir <- "kidney/manuscript/supplementary/figures"
dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)

strategy_labels <- c(exact_hardy = "Exact on Hardy + ischemic bin",
                     paper_style = "Exact on ischemic bin only")
strategy_colors <- c("Exact on Hardy + ischemic bin" = "#E15759",
                     "Exact on ischemic bin only" = "#4C78A8")

deg_plot_tbl <- iteration_tbl %>%
  mutate(
    strategy = recode(strategy, !!!strategy_labels),
    ratio = factor(paste0("1:", ratio), levels = c("1:1", "1:2", "1:3"))
  )

# --- Panel A: DEG counts ---

p_deg <- ggplot(deg_plot_tbl, aes(x = ratio, y = degs, fill = strategy)) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.6,
               outlier.shape = NA, alpha = 0.85) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
             size = 1.4, alpha = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = strategy_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
  labs(x = "Matching ratio", y = "DEGs (FDR < 0.05)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

# --- Panel B: Hallmark NES correlation ---

p_nes <- ggplot(deg_plot_tbl, aes(x = ratio, y = rho_nes, fill = strategy)) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.6,
               outlier.shape = NA, alpha = 0.85) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
             size = 1.4, alpha = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = strategy_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Matching ratio",
       y = expression(Spearman~rho~"(Hallmark NES vs primary model)"),
       fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

# --- Combined supplementary figure ---

p_combined <- cowplot::plot_grid(
  p_deg + theme(legend.position = "none"),
  p_nes + theme(legend.position = "none"),
  labels = c("A", "B"), label_size = 14, ncol = 2
)

legend <- cowplot::get_legend(
  p_deg + theme(legend.position = "bottom",
                legend.text = element_text(size = 10))
)

p_final <- cowplot::plot_grid(p_combined, legend, ncol = 1, rel_heights = c(1, 0.08))

ggsave(file.path(supp_fig_dir, "fig_s13_strict_matching.png"), p_final,
       width = 10, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "24_strict_matching.png"), p_final,
       width = 10, height = 5.5, dpi = 300)

cat("Saved supplementary figure:\n")
cat("  ", file.path(supp_fig_dir, "fig_s13_strict_matching.png"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 7: SUPPLEMENTARY TABLE                                          ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("SECTION 7: Supplementary Table")

table_s4 <- iteration_tbl %>%
  group_by(strategy, ratio) %>%
  summarise(
    n_seeds = n(),
    median_cases = median(n_cases),
    median_controls = median(n_controls),
    median_degs = median(degs),
    deg_range = sprintf("%d--%d", min(degs), max(degs)),
    median_rho_nes = sprintf("%.3f", median(rho_nes)),
    rho_range = sprintf("%.3f--%.3f", min(rho_nes), max(rho_nes)),
    median_hallmark_sig_both = median(n_hallmark_sig_both),
    .groups = "drop"
  ) %>%
  mutate(
    strategy = recode(strategy, !!!strategy_labels),
    ratio = paste0("1:", ratio)
  ) %>%
  rename(
    Strategy = strategy,
    Ratio = ratio,
    Seeds = n_seeds,
    `Median cases` = median_cases,
    `Median controls` = median_controls,
    `Median DEGs` = median_degs,
    `DEG range` = deg_range,
    `Median NES rho` = median_rho_nes,
    `NES rho range` = rho_range,
    `Shared sig pathways` = median_hallmark_sig_both
  )

write_csv(table_s4, file.path(res_dir, "table_s4_strict_matching.csv"))
cat("Saved supplementary table:\n")
cat("  ", file.path(res_dir, "table_s4_strict_matching.csv"), "\n")
print(table_s4)


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 8: SUMMARY                                                      ║
# ╚════════════════════════════════════════════════════════════════════════════╝

section_header("SECTION 8: Summary")

cat("Strict matching summary table:\n")
print(summary_tbl)

cat("\nCore Hallmark NES table (seed 1):\n")
print(core_seed1)
