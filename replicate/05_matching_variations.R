#!/usr/bin/env Rscript

source("kidney/manuscript/replicate/00_common.R")

suppressPackageStartupMessages({
  library(MatchIt)
})

run_match_scenario <- function(tb_base, counts, seed, strategy, ratio) {
  set.seed(seed)
  tb <- tb_base[sample(nrow(tb_base)), , drop = FALSE]

  if (strategy == "exact_hardy") {
    match_obj <- matchit(
      MHRNLFLR ~ AGE + BMI + ischemic_hrs + SMRIN + SMEXNCRT,
      data = tb,
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
      data = tb,
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
  fit <- fit_limma(tb_matched, counts[, tb_matched$SAMPID, drop = FALSE], model_formula, "MHRNLFLR")
  g <- as.integer(tb_matched$MHRNLFLR == 1)

  tibble(
    seed = seed,
    strategy = strategy,
    ratio = ratio,
    n_samples = nrow(tb_matched),
    n_cases = sum(tb_matched$MHRNLFLR == 1),
    n_controls = sum(tb_matched$MHRNLFLR == 0),
    degs = sum(fit$tt$adj.P.Val < 0.05),
    smd_age = numeric_smd(tb_matched$AGE, g),
    smd_bmi = numeric_smd(tb_matched$BMI, g),
    smd_ischemic = numeric_smd(tb_matched$ischemic_hrs, g),
    smd_rin = numeric_smd(tb_matched$SMRIN, g),
    smd_slow = binary_smd(tb_matched$DTHHRDY == "Slow", g),
    smd_vent = binary_smd(tb_matched$DTHHRDY == "Ventilator", g)
  )
}

cat("05_matching_variations.R\n")
cat(strrep("=", 80), "\n")

inputs <- read_inputs()
sat <- prepare_tissue(inputs$pheno, inputs$counts_all, "Adipose - Subcutaneous")

female_tb <- sat$tb %>%
  filter(SEX == "Female") %>%
  mutate(
    isch_bin = cut(
      ischemic_hrs,
      breaks = c(-Inf, 3, 10, Inf),
      labels = c("Short", "Medium", "Long"),
      right = TRUE
    )
  )

seeds <- 1:10
scenarios <- expand_grid(
  strategy = c("exact_hardy", "paper_style"),
  ratio = c(1, 2, 3),
  seed = seeds
)

results <- pmap_dfr(
  scenarios,
  function(strategy, ratio, seed) {
    run_match_scenario(female_tb, sat$counts, seed, strategy, ratio)
  }
)

summary_tbl <- results %>%
  group_by(strategy, ratio) %>%
  summarise(
    median_cases = median(n_cases),
    median_controls = median(n_controls),
    median_degs = median(degs),
    min_degs = min(degs),
    max_degs = max(degs),
    median_abs_smd_age = median(abs(smd_age)),
    median_abs_smd_bmi = median(abs(smd_bmi)),
    median_abs_smd_ischemic = median(abs(smd_ischemic)),
    median_abs_smd_rin = median(abs(smd_rin)),
    median_abs_smd_slow = median(abs(smd_slow)),
    median_abs_smd_vent = median(abs(smd_vent)),
    .groups = "drop"
  )

write_csv(results, file.path(replicate_dirs$tables, "05_matching_variations_iterations.csv"))
write_csv(summary_tbl, file.path(replicate_dirs$tables, "05_matching_variations_summary.csv"))

p_deg <- results %>%
  mutate(
    strategy = recode(strategy, exact_hardy = "Exact Hardy + ischemia-bin", paper_style = "Paper-style"),
    ratio = factor(paste0("1:", ratio), levels = c("1:1", "1:2", "1:3"))
  ) %>%
  ggplot(aes(x = ratio, y = degs, fill = strategy)) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7), size = 1.8, alpha = 0.8) +
  scale_fill_manual(values = c("Exact Hardy + ischemia-bin" = "#E15759", "Paper-style" = "#4C78A8")) +
  labs(
    title = "Female SAT DEG counts across matching ratios",
    x = "Matching ratio",
    y = "DEGs at FDR < 0.05",
    fill = NULL
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(replicate_dirs$figures, "05_matching_variations_deg.png"), p_deg, width = 8.2, height = 5, dpi = 300)

p_balance <- results %>%
  transmute(
    strategy,
    ratio = paste0("1:", ratio),
    `BMI` = abs(smd_bmi),
    `Ischemic time` = abs(smd_ischemic),
    `Slow death` = abs(smd_slow)
  ) %>%
  pivot_longer(cols = c("BMI", "Ischemic time", "Slow death"), names_to = "metric", values_to = "abs_smd") %>%
  mutate(strategy = recode(strategy, exact_hardy = "Exact Hardy + ischemia-bin", paper_style = "Paper-style")) %>%
  ggplot(aes(x = ratio, y = abs_smd, fill = strategy)) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7), size = 1.3, alpha = 0.7) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c("Exact Hardy + ischemia-bin" = "#E15759", "Paper-style" = "#4C78A8")) +
  labs(
    title = "Residual imbalance across matching ratios",
    x = "Matching ratio",
    y = "|SMD|",
    fill = NULL
  ) +
  theme_bw(base_size = 10)

ggsave(file.path(replicate_dirs$figures, "05_matching_variations_balance.png"), p_balance, width = 9, height = 5.2, dpi = 300)

saveRDS(
  list(
    iterations = results,
    summary = summary_tbl
  ),
  file.path(replicate_dirs$rds, "05_matching_variations.rds")
)

write_session_info(file.path(replicate_dirs$logs, "05_matching_variations_sessionInfo.txt"))

cat("\nMatching variation summary:\n")
print(summary_tbl)
