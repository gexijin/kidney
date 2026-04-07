#!/usr/bin/env Rscript

source("kidney/manuscript/replicate/00_common.R")

suppressPackageStartupMessages({
  library(MatchIt)
})

deg_from_match <- function(tb_matched, counts, formula) {
  fit <- fit_limma(tb_matched, counts[, tb_matched$SAMPID, drop = FALSE], formula, "MHRNLFLR")
  list(
    fit = fit,
    degs = sum(fit$tt$adj.P.Val < 0.05),
    balance = matching_balance(tb_matched)
  )
}

matching_balance <- function(tb_matched) {
  g <- as.integer(tb_matched$MHRNLFLR == 1)
  tibble(
    smd_ischemic = numeric_smd(tb_matched$ischemic_hrs, g),
    smd_rin = numeric_smd(tb_matched$SMRIN, g),
    smd_age = numeric_smd(tb_matched$AGE, g),
    smd_slow = binary_smd(tb_matched$DTHHRDY == "Slow", g),
    smd_vent = binary_smd(tb_matched$DTHHRDY == "Ventilator", g)
  )
}

run_match_iteration <- function(tb, counts, seed, type = c("paper_style", "exact_hardy")) {
  type <- match.arg(type)
  set.seed(seed)
  tb <- tb[sample(nrow(tb)), , drop = FALSE]

  if (type == "paper_style") {
    m <- matchit(
      MHRNLFLR ~ AGE + RACE + DTHHRDY + ischemic_hrs + SMRIN,
      data = tb,
      method = "nearest",
      ratio = 3,
      caliper = 0.25,
      std.caliper = TRUE,
      m.order = "random",
      exact = ~isch_bin
    )
    matched <- match.data(m)
    form <- as.formula(paste(
      "~ AGE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
      paste(setdiff(joint_factors, "MHRNLFLR"), collapse = " + "),
      "+ MHRNLFLR"
    ))
  } else {
    m <- matchit(
      MHRNLFLR ~ AGE + BMI + ischemic_hrs + SMRIN + SMEXNCRT,
      data = tb,
      method = "nearest",
      ratio = 1,
      caliper = 0.25,
      std.caliper = TRUE,
      m.order = "random",
      exact = ~DTHHRDY + isch_bin
    )
    matched <- match.data(m)
    form <- as.formula(paste(
      "~ AGE + BMI + DTHHRDY + isch_bin + ischemic_hrs + SMRIN + SMEXNCRT +",
      paste(setdiff(joint_factors, "MHRNLFLR"), collapse = " + "),
      "+ MHRNLFLR"
    ))
  }

  de <- deg_from_match(matched, counts, form)
  bind_cols(
    tibble(
      seed = seed,
      match_type = type,
      n_samples = nrow(matched),
      n_cases = sum(matched$MHRNLFLR == 1),
      n_controls = sum(matched$MHRNLFLR == 0),
      degs = de$degs
    ),
    de$balance
  )
}

cat("04_matching_sensitivity.R\n")
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

paper_res <- map_dfr(seeds, ~ run_match_iteration(female_tb, sat$counts, .x, "paper_style"))
hardy_res <- map_dfr(seeds, ~ run_match_iteration(female_tb, sat$counts, .x, "exact_hardy"))
all_res <- bind_rows(paper_res, hardy_res)

summary_tbl <- all_res %>%
  group_by(match_type) %>%
  summarise(
    iterations = n(),
    median_cases = median(n_cases),
    median_controls = median(n_controls),
    median_degs = median(degs),
    min_degs = min(degs),
    max_degs = max(degs),
    median_abs_smd_ischemic = median(abs(smd_ischemic)),
    median_abs_smd_slow = median(abs(smd_slow)),
    median_abs_smd_vent = median(abs(smd_vent)),
    .groups = "drop"
  )

write_csv(all_res, file.path(replicate_dirs$tables, "04_matching_iterations.csv"))
write_csv(summary_tbl, file.path(replicate_dirs$tables, "04_matching_summary.csv"))

p_deg <- all_res %>%
  mutate(match_type = recode(match_type, paper_style = "Paper-style match", exact_hardy = "Exact Hardy + ischemia match")) %>%
  ggplot(aes(x = match_type, y = degs, fill = match_type)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("Paper-style match" = "#4C78A8", "Exact Hardy + ischemia match" = "#E15759"), guide = "none") +
  labs(
    title = "Female SAT CKD DEG counts across matching strategies",
    x = NULL,
    y = "DEGs at FDR < 0.05"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(replicate_dirs$figures, "04_matching_deg_boxplot.png"), p_deg, width = 7.5, height = 4.8, dpi = 300)

p_balance <- all_res %>%
  transmute(
    match_type,
    `Ischemic time SMD` = abs(smd_ischemic),
    `Slow-death SMD` = abs(smd_slow),
    `Ventilator SMD` = abs(smd_vent)
  ) %>%
  pivot_longer(-match_type, names_to = "metric", values_to = "abs_smd") %>%
  mutate(match_type = recode(match_type, paper_style = "Paper-style match", exact_hardy = "Exact Hardy + ischemia match")) %>%
  ggplot(aes(x = match_type, y = abs_smd, fill = match_type)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.08, size = 1.6, alpha = 0.8) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c("Paper-style match" = "#4C78A8", "Exact Hardy + ischemia match" = "#E15759"), guide = "none") +
  labs(
    title = "Residual imbalance across matching strategies",
    x = NULL,
    y = "|SMD|"
  ) +
  theme_bw(base_size = 10)

ggsave(file.path(replicate_dirs$figures, "04_matching_balance_boxplot.png"), p_balance, width = 9, height = 4.8, dpi = 300)

saveRDS(
  list(
    iterations = all_res,
    summary = summary_tbl
  ),
  file.path(replicate_dirs$rds, "04_matching_sensitivity.rds")
)

write_session_info(file.path(replicate_dirs$logs, "04_matching_sensitivity_sessionInfo.txt"))

cat("\nMatching sensitivity summary:\n")
print(summary_tbl)
