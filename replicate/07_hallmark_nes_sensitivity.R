#!/usr/bin/env Rscript

source("kidney/manuscript/replicate/00_common.R")

suppressPackageStartupMessages({
  library(MatchIt)
  library(msigdbr)
  library(fgsea)
})

load_hallmark_sets_local <- function() {
  sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
    dplyr::select(gs_name, ensembl_gene) %>%
    filter(ensembl_gene != "") %>%
    distinct()

  split(sets$ensembl_gene, sets$gs_name)
}

run_gsea_local <- function(tt, pathways) {
  ranks <- setNames(tt$t, tt$gene_id)
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500, nPermSimple = 10000)
}

compare_nes <- function(gsea_ref, gsea_new) {
  merged <- inner_join(
    gsea_ref %>% select(pathway, NES_ref = NES, padj_ref = padj),
    gsea_new %>% select(pathway, NES_new = NES, padj_new = padj),
    by = "pathway"
    )

  keep_all <- is.finite(merged$NES_ref) & is.finite(merged$NES_new)
  keep_sig <- keep_all & (merged$padj_ref < 0.05 | merged$padj_new < 0.05)

  tibble(
    rho_nes = if (sum(keep_all) >= 2) cor(merged$NES_ref[keep_all], merged$NES_new[keep_all], method = "spearman") else NA_real_,
    rho_sig = if (sum(keep_sig) >= 2) cor(merged$NES_ref[keep_sig], merged$NES_new[keep_sig], method = "spearman") else NA_real_,
    n_sig_either = sum(merged$padj_ref < 0.05 | merged$padj_new < 0.05, na.rm = TRUE),
    n_sig_both = sum(merged$padj_ref < 0.05 & merged$padj_new < 0.05, na.rm = TRUE)
  )
}

run_match_and_gsea <- function(tb_base, counts, pathways, seed, strategy, ratio) {
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
  gsea <- run_gsea_local(fit$tt, pathways)

  list(
    matched = tb_matched,
    fit = fit,
    gsea = gsea
  )
}

cat("07_hallmark_nes_sensitivity.R\n")
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

hallmark_sets <- load_hallmark_sets_local()

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
primary_gsea <- run_gsea_local(primary_fit$tt, hallmark_sets)

scenarios <- expand_grid(
  strategy = c("exact_hardy", "paper_style"),
  ratio = c(1, 2, 3),
  seed = 1:10
)

iter_results <- pmap_dfr(
  scenarios,
  function(strategy, ratio, seed) {
    out <- run_match_and_gsea(female_tb, sat$counts, hallmark_sets, seed, strategy, ratio)
    cmp <- compare_nes(primary_gsea, out$gsea)

    tibble(
      strategy = strategy,
      ratio = ratio,
      seed = seed,
      n_cases = sum(out$matched$MHRNLFLR == 1),
      n_controls = sum(out$matched$MHRNLFLR == 0),
      degs = sum(out$fit$tt$adj.P.Val < 0.05),
      rho_nes = cmp$rho_nes,
      rho_nes_sig = cmp$rho_sig,
      n_sig_either = cmp$n_sig_either,
      n_sig_both = cmp$n_sig_both
    )
  }
)

summary_tbl <- iter_results %>%
  group_by(strategy, ratio) %>%
  summarise(
    median_cases = median(n_cases),
    median_controls = median(n_controls),
    median_degs = median(degs),
    median_rho_nes = median(rho_nes, na.rm = TRUE),
    min_rho_nes = if (any(is.finite(rho_nes))) min(rho_nes, na.rm = TRUE) else NA_real_,
    max_rho_nes = if (any(is.finite(rho_nes))) max(rho_nes, na.rm = TRUE) else NA_real_,
    median_rho_nes_sig = median(rho_nes_sig, na.rm = TRUE),
    median_sig_either = median(n_sig_either),
    median_sig_both = median(n_sig_both),
    .groups = "drop"
  )

write_csv(iter_results, file.path(replicate_dirs$tables, "07_hallmark_nes_iterations.csv"))
write_csv(summary_tbl, file.path(replicate_dirs$tables, "07_hallmark_nes_summary.csv"))

seed1_examples <- list(
  exact_hardy_1_1 = run_match_and_gsea(female_tb, sat$counts, hallmark_sets, 1, "exact_hardy", 1),
  exact_hardy_1_2 = run_match_and_gsea(female_tb, sat$counts, hallmark_sets, 1, "exact_hardy", 2),
  exact_hardy_1_3 = run_match_and_gsea(female_tb, sat$counts, hallmark_sets, 1, "exact_hardy", 3),
  paper_style_1_3 = run_match_and_gsea(female_tb, sat$counts, hallmark_sets, 1, "paper_style", 3)
)

seed1_nes <- list(
  primary_gsea %>% transmute(pathway, primary_NES = NES, primary_padj = padj),
  seed1_examples$exact_hardy_1_1$gsea %>% transmute(pathway, exact11_NES = NES, exact11_padj = padj),
  seed1_examples$exact_hardy_1_2$gsea %>% transmute(pathway, exact12_NES = NES, exact12_padj = padj),
  seed1_examples$exact_hardy_1_3$gsea %>% transmute(pathway, exact13_NES = NES, exact13_padj = padj),
  seed1_examples$paper_style_1_3$gsea %>% transmute(pathway, paper13_NES = NES, paper13_padj = padj)
) %>%
  reduce(full_join, by = "pathway") %>%
  arrange(primary_NES)

write_csv(seed1_nes, file.path(replicate_dirs$tables, "07_seed1_hallmark_nes.csv"))

p_corr <- iter_results %>%
  filter(is.finite(rho_nes)) %>%
  mutate(
    strategy = recode(strategy, exact_hardy = "Exact Hardy + ischemia-bin", paper_style = "Paper-style"),
    ratio = factor(paste0("1:", ratio), levels = c("1:1", "1:2", "1:3"))
  ) %>%
  ggplot(aes(x = ratio, y = rho_nes, fill = strategy)) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7), size = 1.6, alpha = 0.8) +
  scale_fill_manual(values = c("Exact Hardy + ischemia-bin" = "#E15759", "Paper-style" = "#4C78A8")) +
  labs(
    title = "Hallmark NES correlation with primary female SAT model",
    x = "Matching ratio",
    y = "Spearman correlation of NES",
    fill = NULL
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(replicate_dirs$figures, "07_hallmark_nes_correlation.png"), p_corr, width = 8.2, height = 5, dpi = 300)

p_seed1 <- seed1_nes %>%
  select(pathway, primary_NES, exact11_NES, exact12_NES, exact13_NES, paper13_NES) %>%
  pivot_longer(-pathway, names_to = "model", values_to = "NES") %>%
  mutate(
    model = recode(
      model,
      primary_NES = "Primary",
      exact11_NES = "Exact 1:1",
      exact12_NES = "Exact 1:2",
      exact13_NES = "Exact 1:3",
      paper13_NES = "Paper 1:3"
    ),
    pathway = forcats::fct_reorder(pathway, NES, .fun = mean, .na_rm = TRUE)
  ) %>%
  ggplot(aes(x = NES, y = pathway, color = model)) +
  geom_point(size = 1.8, alpha = 0.9) +
  labs(
    title = "Seed-1 Hallmark NES across primary and matched analyses",
    x = "Normalized enrichment score",
    y = NULL,
    color = NULL
  ) +
  theme_bw(base_size = 9)

ggsave(file.path(replicate_dirs$figures, "07_seed1_hallmark_nes.png"), p_seed1, width = 8.5, height = 8, dpi = 300)

saveRDS(
  list(
    primary_gsea = primary_gsea,
    iterations = iter_results,
    summary = summary_tbl,
    seed1_nes = seed1_nes
  ),
  file.path(replicate_dirs$rds, "07_hallmark_nes_sensitivity.rds")
)

write_session_info(file.path(replicate_dirs$logs, "07_hallmark_nes_sensitivity_sessionInfo.txt"))

cat("\nHallmark NES sensitivity summary:\n")
print(summary_tbl)
