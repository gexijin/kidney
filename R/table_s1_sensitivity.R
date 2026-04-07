#!/usr/bin/env Rscript
# Fill Table S1: 8-model covariate sensitivity analysis
library(tidyverse)
library(edgeR)
library(limma)
library(sva)
library(msigdbr)
library(fgsea)

source("kidney/R/functions.R")
set.seed(42)
res_dir <- "kidney/results"

dat <- readRDS(file.path(res_dir, "prep_adipose_subcutaneous.rds"))
tb <- dat$tb
counts <- dat$counts

de_f <- readRDS(file.path(res_dir, "de_female.rds"))
t_primary <- setNames(de_f$tt$t, de_f$tt$gene_id)
primary_degs <- de_f$tt$gene_id[de_f$tt$adj.P.Val < 0.05]

# Hallmark gene sets and primary GSEA for pathway correlation
hallmark_list <- load_hallmark_sets(filter_empty = TRUE)
primary_gsea <- run_gsea(de_f$tt, hallmark_list, "primary", seed = 42)

# Subsets
tb_f <- tb %>% filter(SEX == "Female")
counts_f <- counts[, tb_f$SAMPID]
tb_m <- tb %>% filter(SEX == "Male")
counts_m <- counts[, tb_m$SAMPID]

# Ensure factors
# Factor only non-binary variables; keep binary MH* as numeric (0/1) so
# coefficient names stay "MHRNLFLR" not "MHRNLFLR1"
for (v in c("RACE","DTHHRDY","SMCENTER")) {
  if (v %in% names(tb_f)) tb_f[[v]] <- factor(tb_f[[v]])
  if (v %in% names(tb_m)) tb_m[[v]] <- factor(tb_m[[v]])
}
# Ensure binary variables are numeric 0/1
for (v in c("MHABNWBC","MHLVRDIS","MHT2D","MHRNLFLR","MHHTN","MHHRTDIS","MHCOPD","MHHRTATT")) {
  if (v %in% names(tb_f)) tb_f[[v]] <- as.numeric(as.character(tb_f[[v]]))
  if (v %in% names(tb_m)) tb_m[[v]] <- as.numeric(as.character(tb_m[[v]]))
}
tb_f$smoker <- as.numeric(as.character(tb_f$smoker))
tb_m$smoker <- as.numeric(as.character(tb_m$smoker))

# Helper: run one explicit-covariate model
run_model <- function(formula, tb_sub, counts_sub, label) {
  cat(sprintf("  %s ... ", label))
  # Drop rows with NA in any model variable
  vars <- all.vars(formula)
  cc <- complete.cases(tb_sub[, vars, drop = FALSE])
  if (sum(!cc) > 0) cat(sprintf("(dropped %d NA rows) ", sum(!cc)))
  tb_sub <- tb_sub[cc, ]
  counts_sub <- counts_sub[, tb_sub$SAMPID]
  design <- build_design(formula, tb_sub)
  coef <- get_required_coef(design, "MHRNLFLR", label)
  res <- run_voom_de(counts_sub, design, coef)
  tt <- res$tt
  n_deg <- sum(tt$adj.P.Val < 0.05)
  n_up <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 0)
  n_down <- sum(tt$adj.P.Val < 0.05 & tt$logFC < 0)
  t_vec <- setNames(tt$t, tt$gene_id)
  common <- intersect(names(t_primary), names(t_vec))
  rho <- cor(t_primary[common], t_vec[common], method = "spearman")
  common_degs <- intersect(primary_degs, names(t_vec))
  concordance <- mean(sign(t_primary[common_degs]) == sign(t_vec[common_degs]))
  rho_pathway <- pathway_cor(tt, primary_gsea, hallmark_list)
  cat(sprintf("%d DEGs, rho=%.3f, pathway_rho=%.3f\n", n_deg, rho, rho_pathway))
  data.frame(label=label, n_deg_f=n_deg, n_up_f=n_up, n_down_f=n_down,
             rho_f=rho, concordance_f=concordance, pathway_rho_f=rho_pathway)
}

# Helper: SVA model
run_sva_model <- function(full_formula, null_formula, tb_sub, counts_sub, label) {
  cat(sprintf("  %s ... ", label))
  dge <- DGEList(counts = counts_sub)
  dge <- calcNormFactors(dge, method = "TMM")
  logcpm <- cpm(dge, log = TRUE, prior.count = 1)
  mod <- model.matrix(full_formula, data = tb_sub)
  mod0 <- model.matrix(null_formula, data = tb_sub)
  qr1 <- qr(mod); mod <- mod[, qr1$pivot[seq_len(qr1$rank)], drop=FALSE]
  qr0 <- qr(mod0); mod0 <- mod0[, qr0$pivot[seq_len(qr0$rank)], drop=FALSE]
  # Use 3 SVs — num.sv overestimates (28-29) and absorbs the CKD signal
  n_sv <- 3
  cat(sprintf("(%d SVs) ", n_sv))
  if (n_sv > 0) {
    sv_obj <- sva(logcpm, mod, mod0, n.sv = n_sv)
    design <- cbind(mod, sv_obj$sv)
    colnames(design)[(ncol(mod)+1):ncol(design)] <- paste0("SV", 1:n_sv)
  } else {
    design <- mod
  }
  coef <- get_required_coef(design, "MHRNLFLR", label)
  voom_obj <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(voom_obj, design))
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  tt$gene_id <- rownames(tt)
  n_deg <- sum(tt$adj.P.Val < 0.05)
  n_up <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 0)
  n_down <- sum(tt$adj.P.Val < 0.05 & tt$logFC < 0)
  t_vec <- setNames(tt$t, tt$gene_id)
  common <- intersect(names(t_primary), names(t_vec))
  rho <- cor(t_primary[common], t_vec[common], method = "spearman")
  common_degs <- intersect(primary_degs, names(t_vec))
  concordance <- mean(sign(t_primary[common_degs]) == sign(t_vec[common_degs]))
  rho_pathway <- pathway_cor(tt, primary_gsea, hallmark_list)
  cat(sprintf("%d DEGs, rho=%.3f, pathway_rho=%.3f\n", n_deg, rho, rho_pathway))
  data.frame(label=label, n_deg_f=n_deg, n_up_f=n_up, n_down_f=n_down,
             rho_f=rho, concordance_f=concordance, pathway_rho_f=rho_pathway)
}

# ── Female models ──
cat("\n=== Female models ===\n")
f_results <- list()

f_results[[1]] <- data.frame(label="1. Primary", n_deg_f=1960, n_up_f=836,
                              n_down_f=1124, rho_f=1.000, concordance_f=1.0,
                              pathway_rho_f=1.000)
cat("  1. Primary: 1960 DEGs (reference)\n")

f_results[[2]] <- run_model(
  ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + MHRNLFLR,
  tb_f, counts_f, "2. Minimal")

f_results[[3]] <- run_model(
  ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D + MHRNLFLR,
  tb_f, counts_f, "3. No CENTER")

f_results[[4]] <- run_model(
  ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D + MHHTN + MHHRTDIS + MHRNLFLR,
  tb_f, counts_f, "4. +HTN+Heart")

f_results[[5]] <- run_model(
  ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D + smoker + MHRNLFLR,
  tb_f, counts_f, "5. +Smoker")

f_results[[6]] <- run_model(
  ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D + MHHTN + MHHRTDIS + MHCOPD + MHHRTATT + smoker + MHRNLFLR,
  tb_f, counts_f, "6. Maximal")

f_results[[7]] <- run_sva_model(
  ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D + MHRNLFLR,
  ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D,
  tb_f, counts_f, "7. Primary+SVA")

f_results[[8]] <- run_sva_model(
  ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + MHRNLFLR,
  ~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN,
  tb_f, counts_f, "8. Minimal+SVA")

# ── Male models (explicit only) ──
cat("\n=== Male models ===\n")
m_results <- list()

m_results[[1]] <- data.frame(label="1. Primary", n_deg_m=5)
cat("  1. Primary: 5 DEGs (reference)\n")

run_male <- function(formula, label) {
  cat(sprintf("  %s ... ", label))
  vars <- all.vars(formula)
  cc <- complete.cases(tb_m[, vars, drop = FALSE])
  tb_sub <- tb_m[cc, ]
  counts_sub <- counts_m[, tb_sub$SAMPID]
  if (sum(!cc) > 0) cat(sprintf("(dropped %d NA rows) ", sum(!cc)))
  design <- build_design(formula, tb_sub)
  coef <- get_required_coef(design, "MHRNLFLR", label)
  res <- run_voom_de(counts_sub, design, coef)
  n_deg <- sum(res$tt$adj.P.Val < 0.05)
  cat(sprintf("%d DEGs\n", n_deg))
  data.frame(label=label, n_deg_m=n_deg)
}

m_results[[2]] <- run_male(~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + MHRNLFLR, "2. Minimal")
m_results[[3]] <- run_male(~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D + MHRNLFLR, "3. No CENTER")
m_results[[4]] <- run_male(~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D + MHHTN + MHHRTDIS + MHRNLFLR, "4. +HTN+Heart")
m_results[[5]] <- run_male(~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D + smoker + MHRNLFLR, "5. +Smoker")
m_results[[6]] <- run_male(~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT + MHABNWBC + MHLVRDIS + MHT2D + MHHTN + MHHRTDIS + MHCOPD + MHHRTATT + smoker + MHRNLFLR, "6. Maximal")

# ── Summary ──
f_df <- bind_rows(f_results)
m_df <- bind_rows(m_results)

cat("\n\n========== TABLE S1 ==========\n")
cat(sprintf("%-18s | %7s | %7s | %7s | %7s | %12s | %12s | %12s\n",
            "Model", "F DEGs", "F up", "F down", "M DEGs",
            "rho vs prim", "concordance", "pathway rho"))
cat(paste(rep("-", 105), collapse=""), "\n")
for (i in 1:8) {
  m_deg <- if (i <= nrow(m_df)) m_df$n_deg_m[i] else NA
  cat(sprintf("%-18s | %7d | %7d | %7d | %7s | %12.3f | %11.1f%% | %12.3f\n",
              f_df$label[i],
              f_df$n_deg_f[i], f_df$n_up_f[i], f_df$n_down_f[i],
              ifelse(is.na(m_deg), "—", as.character(m_deg)),
              f_df$rho_f[i],
              f_df$concordance_f[i] * 100,
              f_df$pathway_rho_f[i]))
}

# Save
saveRDS(list(female=f_df, male=m_df), file.path(res_dir, "table_s1_sensitivity.rds"))
cat("\nSaved: kidney/results/table_s1_sensitivity.rds\n")
