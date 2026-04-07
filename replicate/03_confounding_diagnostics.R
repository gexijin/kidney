#!/usr/bin/env Rscript

source("kidney/manuscript/replicate/00_common.R")

suppressPackageStartupMessages({
  library(MatchIt)
})

count_deg <- function(tt, fdr = 0.05) {
  sum(tt$adj.P.Val < fdr, na.rm = TRUE)
}

compare_to_primary <- function(primary_tt, new_tt) {
  merged <- inner_join(
    primary_tt %>% select(gene_id, t_primary = t, logFC_primary = logFC, padj_primary = adj.P.Val),
    new_tt %>% select(gene_id, t_new = t, logFC_new = logFC, padj_new = adj.P.Val),
    by = "gene_id"
  )

  tibble(
    rho_t = cor(merged$t_primary, merged$t_new, method = "spearman"),
    rho_logFC = cor(merged$logFC_primary, merged$logFC_new, method = "spearman"),
    overlap_deg = sum(merged$padj_primary < 0.05 & merged$padj_new < 0.05),
    new_deg = sum(merged$padj_new < 0.05),
    primary_deg = sum(merged$padj_primary < 0.05)
  )
}

fit_and_summarize <- function(tb, counts, formula, label, primary_tt) {
  fit <- fit_limma(tb, counts[, tb$SAMPID, drop = FALSE], formula, "MHRNLFLR")
  concord <- compare_to_primary(primary_tt, fit$tt)
  tibble(
    analysis = label,
    n_samples = nrow(tb),
    n_cases = sum(tb$MHRNLFLR == 1),
    n_controls = sum(tb$MHRNLFLR == 0),
    degs = count_deg(fit$tt),
    rho_t = concord$rho_t,
    rho_logFC = concord$rho_logFC,
    overlap_with_primary = concord$overlap_deg
  )
}

matching_balance <- function(tb_matched) {
  g <- as.integer(tb_matched$MHRNLFLR == 1)
  bind_rows(
    tibble(
      variable = c("AGE", "BMI", "ischemic_hrs", "SMRIN", "SMEXNCRT"),
      smd = c(
        numeric_smd(tb_matched$AGE, g),
        numeric_smd(tb_matched$BMI, g),
        numeric_smd(tb_matched$ischemic_hrs, g),
        numeric_smd(tb_matched$SMRIN, g),
        numeric_smd(tb_matched$SMEXNCRT, g)
      )
    ),
    factor_level_smds(tb_matched$DTHHRDY, g, "DTHHRDY"),
    factor_level_smds(tb_matched$isch_bin, g, "isch_bin")
  ) %>%
    mutate(abs_smd = abs(smd)) %>%
    arrange(desc(abs_smd))
}

cat("03_confounding_diagnostics.R\n")
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

primary_formula <- as.formula(paste(
  "~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
  paste(joint_factors, collapse = " + ")
))

primary_fit <- fit_limma(
  female_tb,
  sat$counts[, female_tb$SAMPID, drop = FALSE],
  primary_formula,
  "MHRNLFLR"
)

primary_tt <- primary_fit$tt

model_summary <- bind_rows(
  fit_and_summarize(
    female_tb,
    sat$counts,
    primary_formula,
    "Primary female SAT model",
    primary_tt
  ),
  fit_and_summarize(
    female_tb,
    sat$counts,
    as.formula(paste(
      "~ AGE + RACE + BMI + DTHHRDY + splines::ns(ischemic_hrs, df = 4) + SMRIN + SMCENTER + SMEXNCRT +",
      paste(joint_factors, collapse = " + ")
    )),
    "Spline ischemia adjustment (df = 4)",
    primary_tt
  ),
  fit_and_summarize(
    female_tb %>% filter(DTHHRDY == "Ventilator"),
    sat$counts,
    as.formula(paste(
      "~ AGE + RACE + BMI + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
      paste(joint_factors, collapse = " + ")
    )),
    "Ventilator-only females",
    primary_tt
  ),
  fit_and_summarize(
    female_tb %>% filter(
      ischemic_hrs >= quantile(female_tb$ischemic_hrs[female_tb$MHRNLFLR == 1], 0.10),
      ischemic_hrs <= quantile(female_tb$ischemic_hrs[female_tb$MHRNLFLR == 0], 0.95)
    ),
    sat$counts,
    primary_formula,
    "Empirical ischemia-overlap subset",
    primary_tt
  )
)

write_csv(model_summary, file.path(replicate_dirs$tables, "03_model_restriction_summary.csv"))

match_formula_a <- MHRNLFLR ~ AGE + BMI + ischemic_hrs + SMRIN + SMEXNCRT
match_formula_b <- MHRNLFLR ~ AGE + BMI + SMRIN + SMEXNCRT

match_a <- matchit(
  match_formula_a,
  data = female_tb,
  method = "nearest",
  ratio = 1,
  exact = ~DTHHRDY + isch_bin,
  caliper = 0.25,
  std.caliper = TRUE
)

match_b <- matchit(
  match_formula_b,
  data = female_tb,
  method = "nearest",
  ratio = 1,
  exact = ~DTHHRDY + isch_bin,
  caliper = 0.25,
  std.caliper = TRUE
)

matched_a <- match.data(match_a)
matched_b <- match.data(match_b)

stopifnot(
  sum(matched_a$MHRNLFLR == 1) >= 15,
  sum(matched_b$MHRNLFLR == 1) >= 15
)

matched_formula <- as.formula(paste(
  "~ AGE + BMI + DTHHRDY + isch_bin + ischemic_hrs + SMRIN + SMEXNCRT +",
  paste(setdiff(joint_factors, "MHRNLFLR"), collapse = " + "),
  "+ MHRNLFLR"
))

matched_fit_a <- fit_limma(matched_a, sat$counts[, matched_a$SAMPID, drop = FALSE], matched_formula, "MHRNLFLR")
matched_fit_b <- fit_limma(matched_b, sat$counts[, matched_b$SAMPID, drop = FALSE], matched_formula, "MHRNLFLR")

matching_summary <- bind_rows(
  compare_to_primary(primary_tt, matched_fit_a$tt) %>%
    mutate(
      analysis = "Exact Hardy + ischemia-bin matching with ischemic_hrs in distance",
      n_samples = nrow(matched_a),
      n_cases = sum(matched_a$MHRNLFLR == 1),
      n_controls = sum(matched_a$MHRNLFLR == 0),
      degs = count_deg(matched_fit_a$tt)
    ),
  compare_to_primary(primary_tt, matched_fit_b$tt) %>%
    mutate(
      analysis = "Exact Hardy + ischemia-bin matching without ischemic_hrs in distance",
      n_samples = nrow(matched_b),
      n_cases = sum(matched_b$MHRNLFLR == 1),
      n_controls = sum(matched_b$MHRNLFLR == 0),
      degs = count_deg(matched_fit_b$tt)
    )
) %>%
  select(analysis, n_samples, n_cases, n_controls, degs, rho_t, rho_logFC, overlap_with_primary = overlap_deg)

write_csv(matching_summary, file.path(replicate_dirs$tables, "03_exact_matching_summary.csv"))

balance_a <- matching_balance(matched_a) %>% mutate(match_spec = "with_ischemic_distance")
balance_b <- matching_balance(matched_b) %>% mutate(match_spec = "without_ischemic_distance")
balance_all <- bind_rows(balance_a, balance_b)
write_csv(balance_all, file.path(replicate_dirs$tables, "03_exact_matching_balance.csv"))

ps_diag <- glm(
  MHRNLFLR ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D,
  data = female_tb,
  family = binomial()
)

female_tb$ps <- predict(ps_diag, type = "response")

ps_summary <- tibble(
  group = c("Female CKD", "Female control"),
  n = c(sum(female_tb$MHRNLFLR == 1), sum(female_tb$MHRNLFLR == 0)),
  min_ps = c(min(female_tb$ps[female_tb$MHRNLFLR == 1]), min(female_tb$ps[female_tb$MHRNLFLR == 0])),
  median_ps = c(median(female_tb$ps[female_tb$MHRNLFLR == 1]), median(female_tb$ps[female_tb$MHRNLFLR == 0])),
  max_ps = c(max(female_tb$ps[female_tb$MHRNLFLR == 1]), max(female_tb$ps[female_tb$MHRNLFLR == 0]))
)
write_csv(ps_summary, file.path(replicate_dirs$tables, "03_propensity_summary.csv"))

p_ps <- ggplot(female_tb, aes(x = ps, fill = factor(MHRNLFLR))) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 25) +
  scale_fill_manual(values = c("0" = "#4C78A8", "1" = "#E15759"), labels = c("Control", "CKD"), name = NULL) +
  labs(
    title = "Female SAT propensity overlap",
    x = "Estimated propensity for CKD",
    y = "Donor count"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(replicate_dirs$figures, "03_propensity_overlap.png"), p_ps, width = 7, height = 4.5, dpi = 300)

p_match_balance <- balance_all %>%
  mutate(variable = forcats::fct_reorder(variable, abs_smd)) %>%
  ggplot(aes(x = smd, y = variable, color = match_spec)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey70") +
  geom_point(size = 2) +
  labs(
    title = "Balance after exact Hardy + ischemia-bin matching",
    x = "Standardized mean difference",
    y = NULL,
    color = NULL
  ) +
  theme_bw(base_size = 10)

ggsave(file.path(replicate_dirs$figures, "03_exact_matching_balance.png"), p_match_balance, width = 8, height = 6, dpi = 300)

p_deg <- bind_rows(
  model_summary %>% select(analysis, degs),
  matching_summary %>% select(analysis, degs)
) %>%
  mutate(analysis = forcats::fct_reorder(analysis, degs)) %>%
  ggplot(aes(x = degs, y = analysis)) +
  geom_col(fill = "#E15759") +
  labs(
    title = "Female SAT CKD DEG counts under confounding restrictions",
    x = "DEGs at FDR < 0.05",
    y = NULL
  ) +
  theme_bw(base_size = 10)

ggsave(file.path(replicate_dirs$figures, "03_deg_attenuation.png"), p_deg, width = 8, height = 5.5, dpi = 300)

saveRDS(
  list(
    model_summary = model_summary,
    matching_summary = matching_summary,
    balance = balance_all,
    propensity_summary = ps_summary,
    primary_fit = primary_fit,
    matched_fit_a = matched_fit_a,
    matched_fit_b = matched_fit_b
  ),
  file.path(replicate_dirs$rds, "03_confounding_diagnostics.rds")
)

write_session_info(file.path(replicate_dirs$logs, "03_confounding_diagnostics_sessionInfo.txt"))

cat("\nModel restriction summary:\n")
print(model_summary)
cat("\nExact matching summary:\n")
print(matching_summary)
cat("\nPropensity summary:\n")
print(ps_summary)
