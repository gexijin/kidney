#!/usr/bin/env Rscript

source("kidney/manuscript/replicate/00_common.R")

suppressPackageStartupMessages({
  library(MatchIt)
})

run_variant <- function(name, match_expr, tb, counts, formula) {
  match_obj <- eval(match_expr)
  tb_matched <- match.data(match_obj)
  fit <- fit_limma(tb_matched, counts[, tb_matched$SAMPID, drop = FALSE], formula, "MHRNLFLR")
  g <- as.integer(tb_matched$MHRNLFLR == 1)

  tibble(
    variant = name,
    n_samples = nrow(tb_matched),
    n_cases = sum(tb_matched$MHRNLFLR == 1),
    n_controls = sum(tb_matched$MHRNLFLR == 0),
    degs = sum(fit$tt$adj.P.Val < 0.05),
    smd_age = numeric_smd(tb_matched$AGE, g),
    smd_bmi = numeric_smd(tb_matched$BMI, g),
    smd_ischemic = numeric_smd(tb_matched$ischemic_hrs, g),
    smd_rin = numeric_smd(tb_matched$SMRIN, g),
    smd_exonic = numeric_smd(tb_matched$SMEXNCRT, g),
    smd_slow = binary_smd(tb_matched$DTHHRDY == "Slow", g),
    smd_vent = binary_smd(tb_matched$DTHHRDY == "Ventilator", g)
  )
}

cat("06_exact_matching_variants.R\n")
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
    ),
    bmi_bin = cut(
      BMI,
      breaks = quantile(BMI, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      include.lowest = TRUE
    )
  )

set.seed(1)
female_tb <- female_tb[sample(nrow(female_tb)), , drop = FALSE]

exact_formula <- as.formula(paste(
  "~ AGE + BMI + DTHHRDY + isch_bin + ischemic_hrs + SMRIN + SMEXNCRT +",
  paste(setdiff(joint_factors, "MHRNLFLR"), collapse = " + "),
  "+ MHRNLFLR"
))

paper_formula <- as.formula(paste(
  "~ AGE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
  paste(setdiff(joint_factors, "MHRNLFLR"), collapse = " + "),
  "+ MHRNLFLR"
))

variants <- list(
  list(
    name = "Exact Hardy+isch 1:1",
    expr = quote(matchit(
      MHRNLFLR ~ AGE + BMI + ischemic_hrs + SMRIN + SMEXNCRT,
      data = female_tb,
      method = "nearest",
      ratio = 1,
      replace = FALSE,
      caliper = 0.25,
      std.caliper = TRUE,
      exact = ~DTHHRDY + isch_bin
    )),
    formula = exact_formula
  ),
  list(
    name = "Exact Hardy+isch 1:2",
    expr = quote(matchit(
      MHRNLFLR ~ AGE + BMI + ischemic_hrs + SMRIN + SMEXNCRT,
      data = female_tb,
      method = "nearest",
      ratio = 2,
      replace = FALSE,
      caliper = 0.25,
      std.caliper = TRUE,
      exact = ~DTHHRDY + isch_bin
    )),
    formula = exact_formula
  ),
  list(
    name = "Exact Hardy+isch 1:2 Mahalanobis",
    expr = quote(matchit(
      MHRNLFLR ~ AGE + BMI + ischemic_hrs + SMRIN + SMEXNCRT,
      data = female_tb,
      method = "nearest",
      distance = "mahalanobis",
      ratio = 2,
      replace = FALSE,
      exact = ~DTHHRDY + isch_bin
    )),
    formula = exact_formula
  ),
  list(
    name = "Exact Hardy+isch+bmi_bin 1:2",
    expr = quote(matchit(
      MHRNLFLR ~ AGE + BMI + ischemic_hrs + SMRIN + SMEXNCRT,
      data = female_tb,
      method = "nearest",
      ratio = 2,
      replace = FALSE,
      caliper = 0.25,
      std.caliper = TRUE,
      exact = ~DTHHRDY + isch_bin + bmi_bin
    )),
    formula = exact_formula
  ),
  list(
    name = "Paper-style 1:2",
    expr = quote(matchit(
      MHRNLFLR ~ AGE + RACE + DTHHRDY + ischemic_hrs + SMRIN,
      data = female_tb,
      method = "nearest",
      ratio = 2,
      replace = FALSE,
      caliper = 0.25,
      std.caliper = TRUE,
      exact = ~isch_bin
    )),
    formula = paper_formula
  )
)

summary_tbl <- map_dfr(
  variants,
  ~ run_variant(.x$name, .x$expr, female_tb, sat$counts, .x$formula)
)

write_csv(summary_tbl, file.path(replicate_dirs$tables, "06_exact_matching_variants_summary.csv"))

p_deg <- summary_tbl %>%
  mutate(variant = forcats::fct_reorder(variant, degs)) %>%
  ggplot(aes(x = degs, y = variant)) +
  geom_col(fill = "#E15759") +
  labs(
    title = "Female SAT DEGs under targeted matching variants",
    x = "DEGs at FDR < 0.05",
    y = NULL
  ) +
  theme_bw(base_size = 10)

ggsave(file.path(replicate_dirs$figures, "06_exact_matching_variants_deg.png"), p_deg, width = 8, height = 4.8, dpi = 300)

p_balance <- summary_tbl %>%
  select(variant, smd_bmi, smd_ischemic, smd_slow) %>%
  pivot_longer(cols = c(smd_bmi, smd_ischemic, smd_slow), names_to = "metric", values_to = "smd") %>%
  mutate(metric = recode(metric, smd_bmi = "BMI", smd_ischemic = "Ischemic time", smd_slow = "Slow death")) %>%
  ggplot(aes(x = smd, y = variant)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey70") +
  geom_point(color = "#4C78A8", size = 2) +
  facet_wrap(~ metric, scales = "free_x") +
  labs(
    title = "Residual imbalance in targeted matching variants",
    x = "SMD",
    y = NULL
  ) +
  theme_bw(base_size = 10)

ggsave(file.path(replicate_dirs$figures, "06_exact_matching_variants_balance.png"), p_balance, width = 9, height = 5.2, dpi = 300)

saveRDS(
  list(summary = summary_tbl),
  file.path(replicate_dirs$rds, "06_exact_matching_variants.rds")
)

write_session_info(file.path(replicate_dirs$logs, "06_exact_matching_variants_sessionInfo.txt"))

cat("\nTargeted matching variant summary:\n")
print(summary_tbl)
