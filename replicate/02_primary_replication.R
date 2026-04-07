#!/usr/bin/env Rscript

source("kidney/manuscript/replicate/00_common.R")

suppressPackageStartupMessages({
  library(qvalue)
})

count_degs <- function(tt, fdr = 0.05) {
  tibble(
    n_deg = sum(tt$adj.P.Val < fdr),
    n_up = sum(tt$adj.P.Val < fdr & tt$logFC > 0),
    n_down = sum(tt$adj.P.Val < fdr & tt$logFC < 0),
    pi0 = tryCatch(qvalue::pi0est(tt$P.Value)$pi0, error = function(e) NA_real_)
  )
}

compare_effects <- function(tt_a, tt_b, label_a, label_b) {
  inner_join(
    tt_a %>% select(gene_id, logFC_a = logFC, t_a = t, padj_a = adj.P.Val),
    tt_b %>% select(gene_id, logFC_b = logFC, t_b = t, padj_b = adj.P.Val),
    by = "gene_id"
  ) %>%
    mutate(
      label_a = label_a,
      label_b = label_b
    )
}

get_interaction_coef <- function(design) {
  hits <- colnames(design)[grepl("SEX.*:MHRNLFLR|MHRNLFLR:SEX", colnames(design))]
  stopifnot(length(hits) == 1)
  hits[[1]]
}

primary_formula <- as.formula(paste(
  "~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
  paste(joint_factors, collapse = " + ")
))

interaction_formula <- ~ AGE + SEX * MHRNLFLR + RACE + BMI + DTHHRDY +
  ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D

marker_genes <- c("ADIPOQ", "IL6", "LEP", "PPARG", "FASN", "LPL", "CCL2", "ACTA2", "EGFR")

cat("02_primary_replication.R\n")
cat(strrep("=", 80), "\n")

inputs <- read_inputs()
sat <- prepare_tissue(inputs$pheno, inputs$counts_all, "Adipose - Subcutaneous")
vat <- prepare_tissue(inputs$pheno, inputs$counts_all, "Adipose - Visceral (Omentum)")
gene_annot <- readRDS("data/annotation/gencode_v47_genes.rds") %>%
  as_tibble() %>%
  transmute(gene_id = sub("[.].*", "", gene_id), gene_symbol = symbol) %>%
  distinct(gene_id, .keep_all = TRUE)

sat_f <- sat$tb %>% filter(SEX == "Female")
sat_m <- sat$tb %>% filter(SEX == "Male")
vat_f <- vat$tb %>% filter(SEX == "Female")
vat_m <- vat$tb %>% filter(SEX == "Male")

stopifnot(
  nrow(sat_f) == 218,
  sum(sat_f$MHRNLFLR == 1) == 27,
  nrow(sat_m) == 444,
  sum(sat_m$MHRNLFLR == 1) == 56
)

fit_sat_f <- fit_limma(sat_f, sat$counts[, sat_f$SAMPID, drop = FALSE], primary_formula, "MHRNLFLR")
fit_sat_m <- fit_limma(sat_m, sat$counts[, sat_m$SAMPID, drop = FALSE], primary_formula, "MHRNLFLR")
fit_vat_f <- fit_limma(vat_f, vat$counts[, vat_f$SAMPID, drop = FALSE], primary_formula, "MHRNLFLR")
fit_vat_m <- fit_limma(vat_m, vat$counts[, vat_m$SAMPID, drop = FALSE], primary_formula, "MHRNLFLR")

sat_design <- model.matrix(interaction_formula, data = sat$tb)
sat_design_qr <- qr(sat_design)
if (sat_design_qr$rank < ncol(sat_design)) {
  sat_design <- sat_design[, sat_design_qr$pivot[seq_len(sat_design_qr$rank)], drop = FALSE]
}
interaction_coef <- get_interaction_coef(sat_design)
fit_interaction <- fit_limma(sat$tb, sat$counts, interaction_formula, interaction_coef)

sat_f_stats <- count_degs(fit_sat_f$tt) %>% mutate(tissue = "SAT", sex = "Female")
sat_m_stats <- count_degs(fit_sat_m$tt) %>% mutate(tissue = "SAT", sex = "Male")
vat_f_stats <- count_degs(fit_vat_f$tt) %>% mutate(tissue = "VAT", sex = "Female")
vat_m_stats <- count_degs(fit_vat_m$tt) %>% mutate(tissue = "VAT", sex = "Male")
interaction_stats <- count_degs(fit_interaction$tt) %>% mutate(tissue = "SAT", sex = "Interaction")

deg_summary <- bind_rows(sat_f_stats, sat_m_stats, vat_f_stats, vat_m_stats, interaction_stats) %>%
  select(tissue, sex, everything())

write_csv(deg_summary, file.path(replicate_dirs$tables, "02_deg_summary.csv"))

sat_sex_cmp <- compare_effects(fit_sat_f$tt, fit_sat_m$tt, "Female SAT", "Male SAT")
vat_cmp <- compare_effects(fit_sat_f$tt, fit_vat_f$tt, "Female SAT", "Female VAT")

sex_corr <- tibble(
  comparison = "Female SAT vs Male SAT",
  logFC_pearson = cor(sat_sex_cmp$logFC_a, sat_sex_cmp$logFC_b, method = "pearson"),
  t_pearson = cor(sat_sex_cmp$t_a, sat_sex_cmp$t_b, method = "pearson"),
  female_only = sum(sat_sex_cmp$padj_a < 0.05 & sat_sex_cmp$padj_b >= 0.05),
  male_only = sum(sat_sex_cmp$padj_b < 0.05 & sat_sex_cmp$padj_a >= 0.05),
  both = sum(sat_sex_cmp$padj_a < 0.05 & sat_sex_cmp$padj_b < 0.05)
)

depot_corr <- tibble(
  comparison = "Female SAT vs Female VAT",
  logFC_spearman = cor(vat_cmp$logFC_a, vat_cmp$logFC_b, method = "spearman"),
  t_spearman = cor(vat_cmp$t_a, vat_cmp$t_b, method = "spearman"),
  sat_only = sum(vat_cmp$padj_a < 0.05 & vat_cmp$padj_b >= 0.05),
  vat_only = sum(vat_cmp$padj_b < 0.05 & vat_cmp$padj_a >= 0.05),
  both = sum(vat_cmp$padj_a < 0.05 & vat_cmp$padj_b < 0.05)
)

write_csv(sex_corr, file.path(replicate_dirs$tables, "02_sex_correlation_summary.csv"))
write_csv(depot_corr, file.path(replicate_dirs$tables, "02_depot_correlation_summary.csv"))

write_csv(fit_sat_f$tt, file.path(replicate_dirs$tables, "02_sat_female_limma.csv"))
write_csv(fit_sat_m$tt, file.path(replicate_dirs$tables, "02_sat_male_limma.csv"))
write_csv(fit_vat_f$tt, file.path(replicate_dirs$tables, "02_vat_female_limma.csv"))
write_csv(fit_vat_m$tt, file.path(replicate_dirs$tables, "02_vat_male_limma.csv"))
write_csv(fit_interaction$tt, file.path(replicate_dirs$tables, "02_sat_sex_ckd_interaction_limma.csv"))

marker_hits <- bind_rows(
  fit_sat_f$tt %>% mutate(model = "Female SAT"),
  fit_sat_m$tt %>% mutate(model = "Male SAT"),
  fit_vat_f$tt %>% mutate(model = "Female VAT"),
  fit_interaction$tt %>% mutate(model = "Sex x CKD interaction")
) %>%
  left_join(gene_annot, by = "gene_id") %>%
  filter(gene_symbol %in% marker_genes)

write_csv(marker_hits, file.path(replicate_dirs$tables, "02_marker_hits.csv"))

volcano_df <- bind_rows(
  fit_sat_f$tt %>% mutate(panel = "Female SAT"),
  fit_sat_m$tt %>% mutate(panel = "Male SAT")
) %>%
  mutate(sig = case_when(
    adj.P.Val < 0.05 & logFC > 0 ~ "Up",
    adj.P.Val < 0.05 & logFC < 0 ~ "Down",
    TRUE ~ "NS"
  ))

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(size = 0.35, alpha = 0.45) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  facet_wrap(~ panel) +
  scale_color_manual(values = c("Down" = "#4C78A8", "NS" = "grey75", "Up" = "#E15759"), name = NULL) +
  labs(
    title = "Primary CKD differential expression replication",
    x = "log2 fold change (CKD vs control)",
    y = expression(-log[10](FDR))
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(replicate_dirs$figures, "02_primary_volcano.png"), p_volcano, width = 8.5, height = 4.8, dpi = 300)

p_sex <- sat_sex_cmp %>%
  mutate(sig_group = case_when(
    padj_a < 0.05 & padj_b < 0.05 ~ "Both",
    padj_a < 0.05 ~ "Female only",
    padj_b < 0.05 ~ "Male only",
    TRUE ~ "Neither"
  )) %>%
  ggplot(aes(x = logFC_b, y = logFC_a, color = sig_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 0.45, alpha = 0.55) +
  scale_color_manual(
    values = c("Both" = "#B279A2", "Female only" = "#E15759", "Male only" = "#4C78A8", "Neither" = "grey80"),
    name = NULL
  ) +
  labs(
    title = "Female vs male SAT CKD effect sizes",
    x = "Male SAT log2 fold change",
    y = "Female SAT log2 fold change"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(replicate_dirs$figures, "02_sex_effect_scatter.png"), p_sex, width = 5.4, height = 5.4, dpi = 300)

p_depot <- vat_cmp %>%
  mutate(sig_group = case_when(
    padj_a < 0.05 & padj_b < 0.05 ~ "Both",
    padj_a < 0.05 ~ "SAT only",
    padj_b < 0.05 ~ "VAT only",
    TRUE ~ "Neither"
  )) %>%
  ggplot(aes(x = logFC_b, y = logFC_a, color = sig_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 0.45, alpha = 0.55) +
  scale_color_manual(
    values = c("Both" = "#B279A2", "SAT only" = "#E15759", "VAT only" = "#4C78A8", "Neither" = "grey80"),
    name = NULL
  ) +
  labs(
    title = "Female SAT vs female VAT CKD effect sizes",
    x = "Female VAT log2 fold change",
    y = "Female SAT log2 fold change"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(replicate_dirs$figures, "02_depot_effect_scatter.png"), p_depot, width = 5.4, height = 5.4, dpi = 300)

saveRDS(
  list(
    sat_f = fit_sat_f,
    sat_m = fit_sat_m,
    vat_f = fit_vat_f,
    vat_m = fit_vat_m,
    interaction = fit_interaction,
    deg_summary = deg_summary,
    sex_corr = sex_corr,
    depot_corr = depot_corr
  ),
  file.path(replicate_dirs$rds, "02_primary_replication.rds")
)

write_session_info(file.path(replicate_dirs$logs, "02_primary_replication_sessionInfo.txt"))

cat("\nPrimary replication summary:\n")
print(deg_summary)
cat("\nFemale vs male SAT correlation summary:\n")
print(sex_corr)
cat("\nFemale SAT vs VAT correlation summary:\n")
print(depot_corr)
cat("\nInteraction coefficient tested:", interaction_coef, "\n")
