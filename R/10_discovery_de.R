#!/usr/bin/env Rscript
# =============================================================================
# 10: Discovery DE — Sex-Stratified CKD Response in Subcutaneous Adipose
#
# Sections:
#   1 — Setup & Data Prep
#   2 — Sex-Stratified DE
#   3 — Depot Specificity (VAT shows no CKD signal)
#   4 — Formal Interaction Test
#   5 — Table 1
#   6 — Supplementary Analyses
#
# Sex coding: SEX = Female (218 samples, 27 CKD), Male (444, 56)
# =============================================================================
library(tidyverse)
library(edgeR)
library(limma)
library(qvalue)
library(parallel)
library(msigdbr)
library(fgsea)

source("kidney/R/functions.R")

set.seed(42)

res_dir <- "kidney/results"
fig_dir <- "kidney/figures"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

SEX_LABELS <- c("Male" = "Male", "Female" = "Female")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1: SETUP & DATA PREP                                           ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 1: Setup & Data Prep\n")
cat(strrep("=", 70), "\n")

s1 <- readRDS(file.path(res_dir, "stage1_data_adipose_subcutaneous.rds"))
tb <- s1$tissue_df
counts <- s1$counts

# Joint factors (≥50 DEGs in stage1 screening)
joint_factors <- c("MHABNWBC", "MHLVRDIS", "MHRNLFLR", "MHT2D")
cat("Joint factors:", paste(joint_factors, collapse = ", "), "\n")

base_covariates <- c("AGE", "RACE", "BMI", "DTHHRDY", "ischemic_hrs", "SMRIN", "SMCENTER", "SMEXNCRT")


# Complete cases on all model variables
all_vars <- c(base_covariates, joint_factors)
tb <- tb[complete.cases(tb[, all_vars]), ]
counts <- counts[, tb$SAMPID]

cat("Samples:", nrow(tb), "\n")
cat("Genes:", nrow(counts), "\n")
cat(sprintf("Female: %d samples, %d CKD cases\n",
            sum(tb$SEX == "Female"), sum(tb$SEX == "Female" & tb$MHRNLFLR == 1)))
cat(sprintf("Male:   %d samples, %d CKD cases\n",
            sum(tb$SEX == "Male"), sum(tb$SEX == "Male" & tb$MHRNLFLR == 1)))

stopifnot(sum(tb$SEX == "Female") >= 200, sum(tb$SEX == "Male") >= 400)

dat <- list(tb = tb, counts = counts, joint_factors = joint_factors)
saveRDS(dat, file.path(res_dir, "prep_adipose_subcutaneous.rds"))
cat("Saved:", file.path(res_dir, "prep_adipose_subcutaneous.rds"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2: SEX-STRATIFIED DE                                           ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 2: Sex-Stratified DE\n")
cat(strrep("=", 70), "\n")

# Run both sexes
res_female <- fit_sex_stratum(tb, counts, joint_factors, "Female", SEX_LABELS)
res_male   <- fit_sex_stratum(tb, counts, joint_factors, "Male", SEX_LABELS)

# Merge by gene ID
tt_f <- res_female$tt %>%
  dplyr::select(gene_id, logFC_F = logFC, t_F = t, P_F = P.Value, padj_F = adj.P.Val)
tt_m <- res_male$tt %>%
  dplyr::select(gene_id, logFC_M = logFC, t_M = t, P_M = P.Value, padj_M = adj.P.Val)

merged <- inner_join(tt_f, tt_m, by = "gene_id") %>%
  mutate(delta_logFC = logFC_F - logFC_M,
         abs_delta = abs(delta_logFC),
         sig_any = padj_F < 0.05 | padj_M < 0.05)

rho_lfc <- cor(merged$logFC_F, merged$logFC_M, method = "pearson")
rho_t <- cor(merged$t_F, merged$t_M, method = "pearson")
cat(sprintf("\nGenome-wide correlation: logFC r=%.3f, t-stat r=%.3f\n",
            rho_lfc, rho_t))

sig_f_only <- merged %>% filter(padj_F < 0.05 & padj_M >= 0.05)
sig_m_only <- merged %>% filter(padj_M < 0.05 & padj_F >= 0.05)
sig_both   <- merged %>% filter(padj_F < 0.05 & padj_M < 0.05)
cat(sprintf("DEGs: Female-only=%d, Male-only=%d, Both=%d\n",
            nrow(sig_f_only), nrow(sig_m_only), nrow(sig_both)))

# Save
saveRDS(res_female, file.path(res_dir, "de_female.rds"))
saveRDS(res_male, file.path(res_dir, "de_male.rds"))
write_csv(merged %>% arrange(desc(abs_delta)),
          file.path(res_dir, "de_merged.csv"))

# ── Figure: Volcano plots ────────────────────────────────────────────────────
cat("\nGenerating volcano plots...\n")

volcano_data <- bind_rows(
  res_female$tt %>%
    mutate(sex = "Female (SEX=1)",
           sig = ifelse(adj.P.Val < 0.05, ifelse(logFC > 0, "Up", "Down"), "NS")),
  res_male$tt %>%
    mutate(sex = "Male (SEX=0)",
           sig = ifelse(adj.P.Val < 0.05, ifelse(logFC > 0, "Up", "Down"), "NS"))
)

p_volcano <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(size = 0.4, alpha = 0.5) +
  scale_color_manual(values = c(Down = "blue", NS = "grey70", Up = "red"),
                     name = "FDR < 0.05") +
  facet_wrap(~ sex) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(x = "log2 Fold Change (CKD vs Control)",
       y = expression(-log[10](FDR)),
       title = "CKD Differential Expression by Sex") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"))

ggsave(file.path(fig_dir, "10_volcano.png"), p_volcano,
       width = 10, height = 5, dpi = 300)
cat("Saved:", file.path(fig_dir, "10_volcano.png"), "\n")

# ── Figure: logFC scatter ─────────────────────────────────────────────────────
cat("Generating logFC scatter...\n")

merged_plot <- merged %>%
  mutate(sig_cat = case_when(
    padj_F < 0.05 & padj_M < 0.05 ~ "Both",
    padj_F < 0.05 ~ "Female only",
    padj_M < 0.05 ~ "Male only",
    TRUE ~ "Neither"
  ))

p_scatter <- make_point_scatter(
  data = merged_plot,
  x_var = "logFC_M",
  y_var = "logFC_F",
  color_var = "sig_cat",
  palette = c(Both = "purple", `Female only` = "red",
              `Male only` = "blue", Neither = "grey80"),
  x_label = "logFC Male (CKD vs Control)",
  y_label = "logFC Female (CKD vs Control)",
  title = sprintf("Sex-Stratified Effect Sizes (rho = %.3f)", rho_lfc),
  legend_name = "Significant in"
)

ggsave(file.path(fig_dir, "10_logfc_scatter.png"), p_scatter,
       width = 6, height = 6, dpi = 300)
cat("Saved:", file.path(fig_dir, "10_logfc_scatter.png"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 3: DEPOT SPECIFICITY — VAT SHOWS NO CKD SIGNAL                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 3: Depot Specificity — VAT\n")
cat(strrep("=", 70), "\n")

# Load VAT data
s1_vat <- readRDS("kidney/results/stage1_data_adipose_visceral_omentum.rds")
tb_vat <- s1_vat$tissue_df
counts_vat <- s1_vat$counts

# Complete cases on same variables
tb_vat <- tb_vat[complete.cases(tb_vat[, all_vars]), ]
counts_vat <- counts_vat[, tb_vat$SAMPID]

cat(sprintf("VAT samples: %d\n", nrow(tb_vat)))
cat(sprintf("  Female: %d (%d CKD) | Male: %d (%d CKD)\n",
            sum(tb_vat$SEX == "Female"), sum(tb_vat$SEX == "Female" & tb_vat$MHRNLFLR == 1),
            sum(tb_vat$SEX == "Male"), sum(tb_vat$SEX == "Male" & tb_vat$MHRNLFLR == 1)))

# Fit sex-stratified models in VAT (reuse same helper)
res_vat_female <- fit_sex_stratum(tb_vat, counts_vat, joint_factors, "Female", SEX_LABELS)
res_vat_male   <- fit_sex_stratum(tb_vat, counts_vat, joint_factors, "Male", SEX_LABELS)

# Merge VAT results
tt_vf <- res_vat_female$tt %>%
  dplyr::select(gene_id, logFC_F = logFC, t_F = t, P_F = P.Value, padj_F = adj.P.Val)
tt_vm <- res_vat_male$tt %>%
  dplyr::select(gene_id, logFC_M = logFC, t_M = t, P_M = P.Value, padj_M = adj.P.Val)

merged_vat <- inner_join(tt_vf, tt_vm, by = "gene_id") %>%
  mutate(sig_any = padj_F < 0.05 | padj_M < 0.05)

cat(sprintf("\nVAT DEGs: Female=%d, Male=%d\n",
            sum(merged_vat$padj_F < 0.05), sum(merged_vat$padj_M < 0.05)))

# Compare SAT vs VAT female effect sizes (shared genes)
comp <- inner_join(
  merged %>% dplyr::select(gene_id, logFC_F_SAT = logFC_F, padj_F_SAT = padj_F),
  merged_vat %>% dplyr::select(gene_id, logFC_F_VAT = logFC_F, padj_F_VAT = padj_F),
  by = "gene_id"
)

rho_depot <- cor(comp$logFC_F_SAT, comp$logFC_F_VAT, method = "spearman")
cat(sprintf("SAT vs VAT female logFC correlation: rho=%.3f (%d shared genes)\n",
            rho_depot, nrow(comp)))
cat(sprintf("SAT female DEGs: %d | VAT female DEGs: %d\n",
            sum(comp$padj_F_SAT < 0.05), sum(comp$padj_F_VAT < 0.05)))

# Save VAT results
saveRDS(res_vat_female, file.path(res_dir, "de_vat_female.rds"))
saveRDS(res_vat_male, file.path(res_dir, "de_vat_male.rds"))

# ── Figure: SAT vs VAT depot comparison ─────────────────────────────────────
cat("Generating depot comparison figure...\n")

comp_plot <- comp %>%
  mutate(sig_cat = case_when(
    padj_F_SAT < 0.05 & padj_F_VAT < 0.05 ~ "Both depots",
    padj_F_SAT < 0.05 ~ "SAT only",
    padj_F_VAT < 0.05 ~ "VAT only",
    TRUE ~ "Neither"
  ))

p_depot <- make_point_scatter(
  data = comp_plot,
  x_var = "logFC_F_SAT",
  y_var = "logFC_F_VAT",
  color_var = "sig_cat",
  palette = c(`Both depots` = "purple", `SAT only` = "red",
              `VAT only` = "blue", Neither = "grey80"),
  x_label = "logFC Female CKD — Subcutaneous (SAT)",
  y_label = "logFC Female CKD — Visceral (VAT)",
  title = sprintf("Depot Specificity: SAT vs VAT (rho = %.2f)", rho_depot),
  legend_name = "FDR < 0.05 in"
)

ggsave(file.path(fig_dir, "10_depot_sat_vs_vat.png"), p_depot,
       width = 6, height = 6, dpi = 300)
cat("Saved:", file.path(fig_dir, "10_depot_sat_vs_vat.png"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 4: FORMAL INTERACTION TEST                                     ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 4: Formal Interaction Test\n")
cat(strrep("=", 70), "\n")

# Main effects + interaction formula
form_ix <- as.formula(paste(
  "~ AGE + SEX * MHRNLFLR + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
  paste(setdiff(joint_factors, "MHRNLFLR"), collapse = " + ")))

design_ix <- build_design(form_ix, tb)

ix_coef <- get_required_interaction_coef(design_ix, "MHRNLFLR",
                                         "Section 4 interaction model")
cat(sprintf("Design: %d samples x %d columns\n", nrow(design_ix), ncol(design_ix)))
cat(sprintf("Interaction column: %s\n", ix_coef))

# Fit
counts_all <- counts[, tb$SAMPID]
ix_res <- run_voom_de(counts_all, design_ix, ix_coef)
fit_ix <- ix_res$fit
tt_ix <- ix_res$tt

n_ix <- sum(tt_ix$adj.P.Val < 0.05)
n_ix_up <- sum(tt_ix$adj.P.Val < 0.05 & tt_ix$logFC > 0)
n_ix_down <- sum(tt_ix$adj.P.Val < 0.05 & tt_ix$logFC < 0)
pi0_ix <- tryCatch(pi0est(tt_ix$P.Value)$pi0, error = function(e) NA)

cat(sprintf("MHRNLFLR x SEX interaction DEGs: %d (%d up, %d down)\n",
            n_ix, n_ix_up, n_ix_down))
cat(sprintf("Pi0: %.3f\n", pi0_ix))

# Extract main MHRNLFLR effect (= male CKD effect with Male as ref)
main_coef <- get_required_coef(design_ix, "MHRNLFLR",
                               "Section 4 interaction model main effect")
tt_main <- topTable(fit_ix, coef = main_coef, number = Inf, sort.by = "none")
n_main <- sum(tt_main$adj.P.Val < 0.05)
cat(sprintf("MHRNLFLR main effect DEGs (male CKD): %d\n", n_main))

# ── Compare other phenotype x SEX interactions ────────────────────────────────
cat("\nComparing phenotype x SEX interactions...\n")

ix_comparison <- list()
for (fac in joint_factors) {
  form_fac <- as.formula(paste(
    "~ AGE + SEX *", fac, "+ RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
    paste(setdiff(joint_factors, fac), collapse = " + ")))

  design_fac <- tryCatch(build_design(form_fac, tb), error = function(e) NULL)
  if (is.null(design_fac)) next

  fac_ix_coef <- get_required_interaction_coef(
    design_fac, fac, sprintf("Section 4 interaction comparison for %s", fac)
  )

  fac_res <- run_voom_de(counts_all, design_fac, fac_ix_coef)
  tt_f <- fac_res$tt
  n_f <- sum(tt_f$adj.P.Val < 0.05)
  pi0_f <- tryCatch(pi0est(tt_f$P.Value)$pi0, error = function(e) NA)

  ix_comparison[[fac]] <- data.frame(
    factor = fac, n_interaction_degs = n_f,
    pi1 = 1 - pi0_f,
    n_cases_m = sum(tb$SEX == "Male" & tb[[fac]] == 1, na.rm = TRUE),
    n_cases_f = sum(tb$SEX == "Female" & tb[[fac]] == 1, na.rm = TRUE),
    stringsAsFactors = FALSE)
  cat(sprintf("  %s x SEX: %d interaction DEGs, pi1=%.3f\n", fac, n_f, 1 - pi0_f))
}

ix_df <- bind_rows(ix_comparison) %>% arrange(desc(n_interaction_degs))
cat(sprintf("\nCKD x SEX rank: %d of %d phenotypes tested\n",
            which(ix_df$factor == "MHRNLFLR"), nrow(ix_df)))

# Save interaction results
interaction_results <- list(
  tt_interaction = tt_ix,
  n_interaction_deg = n_ix,
  n_ix_up = n_ix_up, n_ix_down = n_ix_down,
  pi0 = pi0_ix,
  fit = fit_ix, design = design_ix,
  interaction_comparison = ix_df
)
saveRDS(interaction_results, file.path(res_dir, "interaction_results.rds"))

# ── Figure: Interaction comparison bar plot ───────────────────────────────────
cat("Generating interaction comparison bar plot...\n")

ix_plot <- ix_df %>%
  mutate(is_ckd = factor == "MHRNLFLR",
         factor = fct_reorder(factor, n_interaction_degs))

p_bar <- ggplot(ix_plot, aes(x = factor, y = n_interaction_degs, fill = is_ckd)) +
  geom_col() +
  scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey60"),
                    guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Interaction DEGs (FDR < 0.05)",
       title = "Phenotype x Sex Interaction Strength") +
  theme_bw(base_size = 12)

ggsave(file.path(fig_dir, "10_interaction_bar.png"), p_bar,
       width = 7, height = 5, dpi = 300)
cat("Saved:", file.path(fig_dir, "10_interaction_bar.png"), "\n")



# ── Female subset objects (needed for supplementary analyses) ────────────────
tb_f <- tb %>% filter(SEX == "Female")
counts_f <- counts[, tb_f$SAMPID]

form_f <- as.formula(paste(
  "~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
  paste(joint_factors, collapse = " + ")))

design_f <- build_design(form_f, tb_f)

n_obs_f <- res_female$n_deg

# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 5: TABLE 1                                                     ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 5: Table 1\n")
cat(strrep("=", 70), "\n")

tb <- tb %>%
  mutate(group = case_when(
    SEX == "Female" & MHRNLFLR == 1 ~ "Female_CKD",
    SEX == "Female" & MHRNLFLR == 0 ~ "Female_Control",
    SEX == "Male" & MHRNLFLR == 1 ~ "Male_CKD",
    SEX == "Male" & MHRNLFLR == 0 ~ "Male_Control"
  ),
  group = factor(group, levels = c("Female_CKD", "Female_Control",
                                    "Male_CKD", "Male_Control")))

groups <- levels(tb$group)

rows <- list()

# n
row_n <- c(Variable = "n")
for (g in groups) row_n[g] <- as.character(sum(tb$group == g))
rows[[length(rows) + 1]] <- row_n

# Age
row_age <- c(Variable = "Age, median [IQR]")
for (g in groups) row_age[g] <- fmt_median_iqr(tb$AGE[tb$group == g])
rows[[length(rows) + 1]] <- row_age

# BMI
row_bmi <- c(Variable = "BMI, median [IQR]")
for (g in groups) row_bmi[g] <- fmt_median_iqr(tb$BMI[tb$group == g])
rows[[length(rows) + 1]] <- row_bmi

# Race (White = 0, Black = 1)
row_race <- c(Variable = "Race (White), n (%)")
for (g in groups) {
  n_white <- sum(tb$RACE[tb$group == g] == "White", na.rm = TRUE)
  n_total <- sum(tb$group == g)
  row_race[g] <- sprintf("%d (%.1f%%)", n_white, 100 * n_white / n_total)
}
rows[[length(rows) + 1]] <- row_race

# Hardy distribution
row_hardy <- c(Variable = "Hardy Score")
rows[[length(rows) + 1]] <- row_hardy
for (h in sort(unique(tb$DTHHRDY))) {
  row_h <- c(Variable = paste0("  ", h))
  for (g in groups) {
    n_h <- sum(tb$DTHHRDY[tb$group == g] == h, na.rm = TRUE)
    n_total <- sum(tb$group == g)
    row_h[g] <- sprintf("%d (%.1f%%)", n_h, 100 * n_h / n_total)
  }
  rows[[length(rows) + 1]] <- row_h
}

# Ischemic time
row_isch <- c(Variable = "Ischemic time (hrs), median [IQR]")
for (g in groups) row_isch[g] <- fmt_median_iqr(tb$ischemic_hrs[tb$group == g])
rows[[length(rows) + 1]] <- row_isch

# RIN
row_rin <- c(Variable = "RIN, median [IQR]")
for (g in groups) row_rin[g] <- fmt_median_iqr(tb$SMRIN[tb$group == g])
rows[[length(rows) + 1]] <- row_rin

# Comorbidities
comorbid_labels <- c(MHHTN = "Hypertension", MHT2D = "Type 2 Diabetes",
                     MHHRTDIS = "Heart Disease", MHCOPD = "COPD",
                     MHHRTATT = "Heart Attack", MHCVD = "Cerebrovascular Disease",
                     MHLVRDIS = "Liver Disease", MHABNWBC = "Abnormal WBC")

row_comorbid_hdr <- c(Variable = "Comorbidities, n (%)")
rows[[length(rows) + 1]] <- row_comorbid_hdr

for (v in names(comorbid_labels)) {
  row_c <- c(Variable = paste0("  ", comorbid_labels[v]))
  for (g in groups) {
    vals <- as.numeric(tb[[v]][tb$group == g])
    vals <- vals[!is.na(vals)]
    row_c[g] <- sprintf("%d (%.1f%%)", sum(vals == 1), 100 * mean(vals == 1))
  }
  rows[[length(rows) + 1]] <- row_c
}

# Smoking
row_smoke <- c(Variable = "Smoker, n (%)")
for (g in groups) {
  vals <- as.numeric(tb$smoker[tb$group == g])
  vals <- vals[!is.na(vals)]
  row_smoke[g] <- sprintf("%d (%.1f%%)", sum(vals == 1), 100 * mean(vals == 1))
}
rows[[length(rows) + 1]] <- row_smoke

# CKD severity (within CKD groups only)
row_sev <- c(Variable = "CKD Severity (within CKD)")
rows[[length(rows) + 1]] <- row_sev

ckd_groups <- c("Female_CKD", "Male_CKD")
ctrl_groups <- c("Female_Control", "Male_Control")

for (varname in c("MHORGNTP", "MHDLYSIS")) {
  if (!(varname %in% names(tb))) next
  label <- if (varname == "MHORGNTP") "  Organ Transplant" else "  Dialysis"
  row_v <- c(Variable = label)
  for (g in groups) {
    if (g %in% ctrl_groups) {
      row_v[g] <- "\u2014"
    } else {
      vals <- as.numeric(tb[[varname]][tb$group == g])
      n_known <- sum(!is.na(vals))
      n_pos <- sum(vals == 1, na.rm = TRUE)
      if (n_known > 0) {
        row_v[g] <- sprintf("%d/%d (%.1f%%)", n_pos, n_known, 100 * n_pos / n_known)
      } else {
        row_v[g] <- "0/0"
      }
    }
  }
  rows[[length(rows) + 1]] <- row_v
}

# Build table
table1 <- bind_rows(lapply(rows, function(r) as.data.frame(t(r), stringsAsFactors = FALSE)))
colnames(table1) <- c("Variable", groups)

# Print
cat("\n")
cat(sprintf("%-35s %15s %15s %15s %15s\n",
            "Variable", "Female CKD", "Female Ctrl", "Male CKD", "Male Ctrl"))
cat(strrep("-", 95), "\n")
for (i in seq_len(nrow(table1))) {
  cat(sprintf("%-35s %15s %15s %15s %15s\n",
              table1$Variable[i],
              table1$Female_CKD[i] %||% "",
              table1$Female_Control[i] %||% "",
              table1$Male_CKD[i] %||% "",
              table1$Male_Control[i] %||% ""))
}

# Statistical tests
cat("\n--- Statistical tests ---\n")

# Female CKD vs Female Control
cat("\nFemale CKD vs Female Control:\n")
tb_females <- tb %>% filter(SEX == "Female")
cat(sprintf("  Age: p=%.3f (Wilcoxon)\n",
            wilcox.test(AGE ~ MHRNLFLR, data = tb_females)$p.value))
cat(sprintf("  BMI: p=%.3f (Wilcoxon)\n",
            wilcox.test(BMI ~ MHRNLFLR, data = tb_females)$p.value))
for (v in names(comorbid_labels)) {
  vals <- as.numeric(tb_females[[v]])
  if (sum(!is.na(vals)) > 0 && length(unique(vals[!is.na(vals)])) > 1) {
    ft <- tryCatch(fisher.test(table(tb_females$MHRNLFLR, vals)),
                   error = function(e) NULL)
    if (!is.null(ft)) {
      sig <- if (ft$p.value < 0.05) " *" else ""
      cat(sprintf("  %s: p=%.3f (Fisher)%s\n", comorbid_labels[v], ft$p.value, sig))
    }
  }
}

# Male CKD vs Male Control
cat("\nMale CKD vs Male Control:\n")
tb_males <- tb %>% filter(SEX == "Male")
cat(sprintf("  Age: p=%.3f (Wilcoxon)\n",
            wilcox.test(AGE ~ MHRNLFLR, data = tb_males)$p.value))
cat(sprintf("  BMI: p=%.3f (Wilcoxon)\n",
            wilcox.test(BMI ~ MHRNLFLR, data = tb_males)$p.value))
for (v in names(comorbid_labels)) {
  vals <- as.numeric(tb_males[[v]])
  if (sum(!is.na(vals)) > 0 && length(unique(vals[!is.na(vals)])) > 1) {
    ft <- tryCatch(fisher.test(table(tb_males$MHRNLFLR, vals)),
                   error = function(e) NULL)
    if (!is.null(ft)) {
      sig <- if (ft$p.value < 0.05) " *" else ""
      cat(sprintf("  %s: p=%.3f (Fisher)%s\n", comorbid_labels[v], ft$p.value, sig))
    }
  }
}

write_csv(table1, file.path(res_dir, "table1.csv"))
cat("\nSaved:", file.path(res_dir, "table1.csv"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 6: SUPPLEMENTARY ANALYSES                                      ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 6: Supplementary Analyses\n")
cat(strrep("=", 70), "\n")

tt_tmm <- res_female$tt

# ── Hardy/agonal subset (Ventilator only) ────────────────────────────────────
cat("\n--- Hardy 0 (Ventilator) subset (female stratum) ---\n")

tb_f_h0 <- tb_f %>% filter(DTHHRDY == "Ventilator")
n_cases_h0 <- sum(tb_f_h0$MHRNLFLR == 1)
cat(sprintf("Hardy 0 females: %d samples, %d CKD cases\n",
            nrow(tb_f_h0), n_cases_h0))

n_deg_hardy <- NA
rho_hardy <- NA
if (n_cases_h0 >= 5 && nrow(tb_f_h0) >= 30) {
  counts_h0 <- counts[, tb_f_h0$SAMPID]
  form_h0 <- as.formula(paste(
    "~ AGE + RACE + BMI + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
    paste(joint_factors, collapse = " + ")))

  design_h0 <- tryCatch(build_design(form_h0, tb_f_h0), error = function(e) NULL)

  if (!is.null(design_h0) && nrow(tb_f_h0) > ncol(design_h0) + 5) {
    hf <- get_required_coef(design_h0, "MHRNLFLR", "Hardy 0 sensitivity model")
    hardy_res <- run_voom_de(counts_h0, design_h0, hf)
    tt_hardy <- hardy_res$tt
    n_deg_hardy <- sum(tt_hardy$adj.P.Val < 0.05)
    cat(sprintf("Hardy 0 DEGs: %d (%.0f%% retained)\n",
                n_deg_hardy, 100 * n_deg_hardy / n_obs_f))

    tt_hardy$gene_id <- rownames(tt_hardy)
    common_h <- intersect(tt_tmm$gene_id, tt_hardy$gene_id)
    rho_hardy <- cor(tt_tmm$t[match(common_h, tt_tmm$gene_id)],
                     tt_hardy$t[match(common_h, tt_hardy$gene_id)],
                     method = "spearman", use = "complete.obs")
    cat(sprintf("Full vs Hardy 0 t-stat rho: %.3f\n", rho_hardy))
  } else {
    cat("Hardy 0 model not estimable (insufficient samples vs columns)\n")
  }
} else {
  cat("Too few female CKD cases in Hardy 0 for stratified analysis\n")
}

# ── Center/batch confounding (female stratum) ────────────────────────────────
cat("\n--- Center/batch confounding (female stratum) ---\n")

center_tab_f <- table(tb_f$SMCENTER, tb_f$MHRNLFLR)
ft_center_f <- fisher.test(center_tab_f, simulate.p.value = TRUE, B = 10000)
cat(sprintf("Center x MHRNLFLR Fisher p = %.4f\n", ft_center_f$p.value))

batch_tab_f <- table(tb_f$SMNABTCHD_Y, tb_f$MHRNLFLR)
ft_batch_f <- fisher.test(batch_tab_f, simulate.p.value = TRUE, B = 10000)
cat(sprintf("Batch year x MHRNLFLR Fisher p = %.4f (not in model)\n", ft_batch_f$p.value))

# ── Sample quality confounding (female stratum) ──────────────────────────────
cat("\n--- Sample quality (female CKD vs control) ---\n")

quality_vars <- c("SMRIN", "SMEXNCRT", "SMMAPRT", "SMRRNART", "ischemic_hrs",
                   "BMI")
quality_tests <- list()
for (qv in quality_vars) {
  wt <- wilcox.test(tb_f[[qv]][tb_f$MHRNLFLR == 1],
                    tb_f[[qv]][tb_f$MHRNLFLR == 0])
  cat(sprintf("  %-12s CKD=%.2f Ctrl=%.2f p=%.3f\n",
              qv,
              median(tb_f[[qv]][tb_f$MHRNLFLR == 1], na.rm = TRUE),
              median(tb_f[[qv]][tb_f$MHRNLFLR == 0], na.rm = TRUE),
              wt$p.value))
  quality_tests[[qv]] <- wt$p.value
}

# ── CKD severity: mild-only sensitivity (female stratum) ─────────────────────
cat("\n--- CKD severity: mild-only vs full model ---\n")

# Classify female CKD cases into severity buckets
# Severe = transplant and/or dialysis; Mild = neither; drop unknowns
tb_f <- tb_f %>%
  mutate(ckd_severity = case_when(
    MHRNLFLR == 0 ~ "control",
    MHRNLFLR == 1 & (MHORGNTP == 1 | MHDLYSIS == 1) ~ "severe",
    MHRNLFLR == 1 & MHORGNTP == 0 & MHDLYSIS == 0 ~ "mild",
    TRUE ~ NA_character_
  ))

n_mild <- sum(tb_f$ckd_severity == "mild", na.rm = TRUE)
n_severe <- sum(tb_f$ckd_severity == "severe", na.rm = TRUE)
n_ctrl_f <- sum(tb_f$ckd_severity == "control", na.rm = TRUE)
n_unk <- sum(is.na(tb_f$ckd_severity))
cat(sprintf(
  "  Mild CKD: %d | Severe (tx/dialysis): %d | Control: %d | Unknown: %d\n",
  n_mild, n_severe, n_ctrl_f, n_unk))

# Mild-only model: exclude severe + unknown CKD cases
tb_f_mild <- tb_f %>%
  filter(ckd_severity %in% c("control", "mild"))
counts_f_mild <- counts[, tb_f_mild$SAMPID]

d_mild <- build_design(form_f, tb_f_mild)

fc_mild <- get_required_coef(d_mild, "MHRNLFLR", "Female mild-only model")
mild_res <- run_voom_de(counts_f_mild, d_mild, fc_mild)
tt_mild <- mild_res$tt

n_deg_mild <- sum(tt_mild$adj.P.Val < 0.05)
pi0_mild <- tryCatch(pi0est(tt_mild$P.Value)$pi0,
                     error = function(e) NA)

cat(sprintf(
  "  Mild-only DEGs: %d (pi0=%.3f)\n", n_deg_mild, pi0_mild
))

# t-stat correlation: full vs mild-only
common_mild <- intersect(tt_tmm$gene_id, tt_mild$gene_id)
rho_mild <- cor(
  tt_tmm$t[match(common_mild, tt_tmm$gene_id)],
  tt_mild$t[match(common_mild, tt_mild$gene_id)],
  method = "spearman", use = "complete.obs"
)
cat(sprintf("  Full vs mild-only t-stat rho: %.3f\n", rho_mild))

# ── GSEA: mild-only vs full ─────────────────────────────────────────────────
cat("\n--- GSEA: mild-only vs full (Hallmark) ---\n")

hallmark_list <- load_hallmark_sets(filter_empty = FALSE)

gsea_full <- run_gsea(tt_tmm, hallmark_list, sprintf("Full (%d cases)", res_female$n_cases))
gsea_mild <- run_gsea(tt_mild, hallmark_list, sprintf("Mild (%d cases)", n_mild))

# Merge and compare NES
gsea_cmp_res <- compare_gsea_results(
  gsea_full, gsea_mild,
  label1 = "full", label2 = "mild",
  core_paths = CORE_PATHWAYS[1:7]
)
gsea_cmp <- gsea_cmp_res$merged
rho_nes <- gsea_cmp_res$rho_nes
cat(sprintf("  Pathway NES rho (full vs mild): %.3f\n", rho_nes))

# ── Figure: NES full vs mild dot plot ────────────────────────────────────────
cat("\nGenerating severity NES comparison figure...\n")

gsea_long <- bind_rows(
  gsea_full %>% filter(padj < 0.05) %>%
    dplyr::select(pathway, NES, padj) %>%
    mutate(model = sprintf("Full (%d CKD)", res_female$n_cases)),
  gsea_mild %>%
    filter(pathway %in% (gsea_full %>%
                           filter(padj < 0.05) %>%
                           pull(pathway))) %>%
    dplyr::select(pathway, NES, padj) %>%
    mutate(model = sprintf("Mild only (%d CKD)", n_mild))
) %>%
  mutate(
    pathway = gsub("HALLMARK_", "", pathway),
    pathway = gsub("_", " ", pathway)
  )

# Order by full-model NES
full_label <- sprintf("Full (%d CKD)", res_female$n_cases)
mild_label <- sprintf("Mild only (%d CKD)", n_mild)

path_order <- gsea_long %>%
  filter(model == full_label) %>%
  arrange(NES) %>%
  pull(pathway)
gsea_long$pathway <- factor(gsea_long$pathway, levels = path_order)

p_sev <- ggplot(gsea_long,
                aes(x = NES, y = pathway, color = model,
                    shape = padj < 0.05)) +
  geom_point(size = 3, alpha = 0.8,
             position = position_dodge(width = 0.5)) +
  scale_color_manual(
    values = setNames(c("firebrick", "steelblue"), c(full_label, mild_label)),
    name = "Model") +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    name = "FDR < 0.05") +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40") +
  labs(x = "Normalized Enrichment Score",
       y = NULL,
       title = "CKD Severity Sensitivity",
       subtitle = paste0(
         "Mild CKD (no transplant/dialysis) vs full cohort | ",
         sprintf("NES rho = %.2f", rho_nes))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "10_severity_nes.png"), p_sev,
       width = 8, height = 7, dpi = 300)
cat("Saved:", file.path(fig_dir, "10_severity_nes.png"), "\n")

# ── Save severity results ───────────────────────────────────────────────────
severity_results <- list(
  n_mild = n_mild, n_severe = n_severe,
  n_deg_mild = n_deg_mild, pi0_mild = pi0_mild,
  rho_tstat_full_vs_mild = rho_mild,
  rho_nes_full_vs_mild = rho_nes,
  gsea_full = gsea_full, gsea_mild = gsea_mild,
  gsea_comparison = gsea_cmp,
  tt_mild = tt_mild
)
write_csv(gsea_cmp, file.path(res_dir, "severity_nes_comparison.csv"))

# ── Male mild-only comparison ────────────────────────────────────────────────
cat("\n--- Male mild-only CKD (sex-specificity control) ---\n")

tb_m <- tb %>% filter(SEX == "Male") %>%
  mutate(ckd_severity = case_when(
    MHRNLFLR == 0 ~ "control",
    MHRNLFLR == 1 & (MHORGNTP == 1 | MHDLYSIS == 1) ~ "severe",
    MHRNLFLR == 1 & MHORGNTP == 0 & MHDLYSIS == 0 ~ "mild",
    TRUE ~ NA_character_
  ))

n_mild_m <- sum(tb_m$ckd_severity == "mild", na.rm = TRUE)
n_severe_m <- sum(tb_m$ckd_severity == "severe", na.rm = TRUE)
cat(sprintf("  Male: %d mild, %d severe CKD cases\n",
            n_mild_m, n_severe_m))

tb_m_mild <- tb_m %>%
  filter(ckd_severity %in% c("control", "mild"))
counts_m_mild <- counts[, tb_m_mild$SAMPID]

form_m <- as.formula(paste(
  "~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
  paste(joint_factors, collapse = " + ")))

d_m_mild <- build_design(form_m, tb_m_mild)

fc_mm <- get_required_coef(d_m_mild, "MHRNLFLR", "Male mild-only model")
male_mild_res <- run_voom_de(counts_m_mild, d_m_mild, fc_mm)
tt_m_mild <- male_mild_res$tt
n_deg_m_mild <- sum(tt_m_mild$adj.P.Val < 0.05)
pi0_m_mild <- tryCatch(pi0est(tt_m_mild$P.Value)$pi0,
                       error = function(e) NA)

cat(sprintf("  Male mild-only DEGs: %d (pi0=%.3f)\n",
            n_deg_m_mild, pi0_m_mild))
cat(sprintf(
  "  Compare: Female mild %d DEGs (%d cases) vs Male mild %d DEGs (%d cases)\n",
  n_deg_mild, n_mild, n_deg_m_mild, n_mild_m))

severity_results$n_mild_male <- n_mild_m
severity_results$n_severe_male <- n_severe_m
severity_results$n_deg_mild_male <- n_deg_m_mild
severity_results$pi0_mild_male <- pi0_m_mild

saveRDS(severity_results,
        file.path(res_dir, "severity_results.rds"))

# ── Save supplementary results ────────────────────────────────────────────────
supplementary_results <- list(
  n_deg_hardy = n_deg_hardy,
  rho_hardy = rho_hardy,
  center_p = ft_center_f$p.value,
  batch_p = ft_batch_f$p.value,
  quality_tests = quality_tests,
  severity = severity_results
)
saveRDS(supplementary_results, file.path(res_dir, "supplementary_stage1.rds"))
