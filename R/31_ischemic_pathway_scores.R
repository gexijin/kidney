#!/usr/bin/env Rscript
# =============================================================================
# Stage 1c: Per-Donor Pathway Scores — CKD Cases + Ischemic-Matched Controls
#
# For each female CKD case, match 2 controls from the same ischemic time bin
# (nearest on ischemic_hrs + AGE). Compute sample-level pathway z-scores for
# core Hallmark pathways. Visualize as a heatmap to assess whether the CKD
# signature is present even at short ischemic times.
#
# Key question: Does GTEX-13OVI (2.0 hrs ischemic, the only Short-bin CKD case)
# show the CKD pathway profile despite minimal post-mortem degradation?
# =============================================================================

library(tidyverse)
library(edgeR)
library(limma)
library(msigdbr)
library(pheatmap)
library(RColorBrewer)

source("kidney/R/functions.R")

set.seed(42)

res_dir <- "kidney/results"
fig_dir <- "kidney/figures"


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 0: LOAD & PREP                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 0: Load & Prep\n")
cat(strrep("=", 70), "\n")

dat <- readRDS(file.path(res_dir, "prep_adipose_subcutaneous.rds"))
tb <- dat$tb
counts <- dat$counts
joint_factors <- dat$joint_factors

# Female only
tb_f <- tb %>% filter(SEX == "Female") %>%
  add_ischemic_bin()

counts_f <- counts[, tb_f$SAMPID]

cat(sprintf("Female samples: %d (CKD=%d, Ctrl=%d)\n",
            nrow(tb_f), sum(tb_f$MHRNLFLR == 1), sum(tb_f$MHRNLFLR == 0)))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1: MATCH 2 CONTROLS PER CKD CASE                               ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 1: Match 2 Controls per CKD Case\n")
cat(strrep("=", 70), "\n")

tb_ckd <- tb_f %>% filter(MHRNLFLR == 1) %>% arrange(ischemic_hrs)
tb_ctrl <- tb_f %>% filter(MHRNLFLR == 0)

cat(sprintf("CKD cases: %d | Controls available: %d\n",
            nrow(tb_ckd), nrow(tb_ctrl)))

# For each CKD case, find 2 nearest controls in same isch_bin
# Distance = |ischemic_hrs_diff| + 0.2 * |AGE_diff|  (weighted)
matched_ctrl_ids <- c()
match_log <- list()

for (i in seq_len(nrow(tb_ckd))) {
  case <- tb_ckd[i, ]

  # Pool = same isch_bin controls, not already used
  pool <- tb_ctrl %>%
    filter(isch_bin == case$isch_bin,
           !(SAMPID %in% matched_ctrl_ids))

  if (nrow(pool) < 2) {
    # Fallback: expand to all controls not yet used
    pool <- tb_ctrl %>%
      filter(!(SAMPID %in% matched_ctrl_ids))
    cat(sprintf("  WARNING: %s (bin=%s) — insufficient same-bin controls, using nearest overall\n",
                case$SAMPID, case$isch_bin))
  }

  pool$dist <- abs(pool$ischemic_hrs - case$ischemic_hrs) +
    0.2 * abs(pool$AGE - case$AGE)
  pool <- pool %>% arrange(dist)

  selected <- pool$SAMPID[1:2]
  matched_ctrl_ids <- c(matched_ctrl_ids, selected)

  match_log[[i]] <- data.frame(
    case_id = case$SAMPID,
    case_isch = case$ischemic_hrs,
    case_age = case$AGE,
    case_bin = as.character(case$isch_bin),
    ctrl1_id = selected[1],
    ctrl1_isch = pool$ischemic_hrs[1],
    ctrl1_age = pool$AGE[1],
    ctrl2_id = selected[2],
    ctrl2_isch = pool$ischemic_hrs[2],
    ctrl2_age = pool$AGE[2],
    stringsAsFactors = FALSE
  )
}

match_df <- bind_rows(match_log)

cat("\nMatch summary:\n")
cat(sprintf("  CKD cases: %d | Matched controls: %d (unique: %d)\n",
            nrow(match_df), length(matched_ctrl_ids),
            length(unique(matched_ctrl_ids))))

# Print matching quality for first few
cat("\nFirst 5 matches:\n")
cat(sprintf("  %-30s %6s %4s %6s | %-30s %6s | %-30s %6s\n",
            "Case", "isch", "age", "bin", "Ctrl1", "isch", "Ctrl2", "isch"))
for (j in 1:min(5, nrow(match_df))) {
  m <- match_df[j, ]
  cat(sprintf("  %-30s %6.1f %4.0f %6s | %-30s %6.1f | %-30s %6.1f\n",
              m$case_id, m$case_isch, m$case_age, m$case_bin,
              m$ctrl1_id, m$ctrl1_isch,
              m$ctrl2_id, m$ctrl2_isch))
}

write_csv(match_df, file.path(res_dir, "ischemic_matched_pairs.csv"))


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2: COVARIATE-ADJUSTED PATHWAY Z-SCORES                         ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 2: Covariate-Adjusted Pathway Z-Scores\n")
cat(strrep("=", 70), "\n")

all_sample_ids <- unique(c(tb_ckd$SAMPID, matched_ctrl_ids))

cat(sprintf("Samples in heatmap: %d (%d CKD + %d ctrl)\n",
            length(all_sample_ids), nrow(tb_ckd), length(matched_ctrl_ids)))

# ── 2a. Fit covariate model on ALL female samples, extract residuals ────────
# Regress out everything EXCEPT MHRNLFLR (the treatment of interest)
# This leaves only CKD-specific + residual variation in each gene

cat("\nFitting covariate model on all 218 female samples...\n")

isch_covars <- setdiff(joint_factors, "MHRNLFLR")
form_covonly <- as.formula(paste(
  "~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
  paste(isch_covars, collapse = " + ")))

design_cov <- build_design(form_covonly, tb_f)

cat(sprintf("  Covariate design: %d samples x %d columns\n",
            nrow(design_cov), ncol(design_cov)))
cat(sprintf("  Covariates: %s\n",
            paste(colnames(design_cov)[-1], collapse = ", ")))

dge_all_f <- DGEList(counts = counts_f)
dge_all_f <- calcNormFactors(dge_all_f, method = "TMM")
v_all_f <- voom(dge_all_f, design_cov, plot = FALSE)

# Fit and extract residuals: expression after removing covariate effects
fit_cov <- lmFit(v_all_f, design_cov)
resid_expr <- residuals(fit_cov, v_all_f)

cat(sprintf("  Residual matrix: %d genes x %d samples\n",
            nrow(resid_expr), ncol(resid_expr)))

# ── 2b. Z-score residuals across all females ────────────────────────────────
resid_means <- rowMeans(resid_expr)
resid_sds <- apply(resid_expr, 1, sd)
resid_sds[resid_sds == 0] <- 1

resid_z <- (resid_expr - resid_means) / resid_sds

# Subset to heatmap samples
resid_z_sub <- resid_z[, all_sample_ids]

# ── 2c. Define pathways and compute scores ──────────────────────────────────
core_path_names <- c(CORE_PATHWAYS, "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_APOPTOSIS")
path_labels <- PATHWAY_LABELS[core_path_names]

hallmark_genes <- load_hallmark_sets(pathways = core_path_names)

# Pathway scores = mean z-score of covariate-adjusted residuals
pathway_scores <- matrix(NA, nrow = length(core_path_names),
                         ncol = length(all_sample_ids),
                         dimnames = list(core_path_names, all_sample_ids))

for (p in core_path_names) {
  genes <- intersect(hallmark_genes[[p]], rownames(resid_z_sub))
  if (length(genes) >= 10) {
    pathway_scores[p, ] <- colMeans(resid_z_sub[genes, , drop = FALSE])
  }
  cat(sprintf("  %s: %d/%d genes matched\n",
              path_labels[p], length(genes), length(hallmark_genes[[p]])))
}

# Drop any pathways with all NA
pathway_scores <- pathway_scores[complete.cases(pathway_scores), , drop = FALSE]
rownames(pathway_scores) <- path_labels[rownames(pathway_scores)]

cat(sprintf("\nPathway score matrix: %d pathways x %d samples\n",
            nrow(pathway_scores), ncol(pathway_scores)))
cat("Note: scores are covariate-adjusted (AGE, RACE, DTHHRDY, ischemic_hrs, SMRIN, SMCENTER, comorbidities regressed out)\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 3: HEATMAP                                                      ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 3: Heatmap\n")
cat(strrep("=", 70), "\n")

# Build annotation df
sample_info <- tb_f %>%
  filter(SAMPID %in% all_sample_ids) %>%
  dplyr::select(SAMPID, MHRNLFLR, ischemic_hrs, AGE, isch_bin, DTHHRDY) %>%
  mutate(
    Status = ifelse(MHRNLFLR == 1, "CKD", "Control"),
    `Ischemic bin` = as.character(isch_bin)
  ) %>%
  arrange(Status, ischemic_hrs)

# Order samples: CKD cases by ischemic time, then their matched controls
# Build ordered sample list: for each CKD case (by ischemic time), case then ctrls
ordered_ids <- c()
for (i in seq_len(nrow(match_df))) {
  ordered_ids <- c(ordered_ids,
                   match_df$case_id[i],
                   match_df$ctrl1_id[i],
                   match_df$ctrl2_id[i])
}
# Remove duplicates (controls matched to multiple cases), keep first occurrence
ordered_ids <- unique(ordered_ids)
# Make sure all samples are included
ordered_ids <- c(ordered_ids, setdiff(all_sample_ids, ordered_ids))

pathway_scores_ordered <- pathway_scores[, ordered_ids]

# Annotation
ann_df <- sample_info %>%
  column_to_rownames("SAMPID") %>%
  dplyr::select(Status, `Ischemic bin`)
ann_df <- ann_df[ordered_ids, , drop = FALSE]

# Highlight GTEX-13OVI
is_13ovi <- grepl("13OVI", ordered_ids)
cat(sprintf("GTEX-13OVI index: %d\n", which(is_13ovi)))

# Color scales
ann_colors <- list(
  Status = c(CKD = "firebrick", Control = "steelblue"),
  `Ischemic bin` = c(Short = "#4DAF4A", Medium = "#FF7F00", Long = "#E41A1C")
)

# Truncate extreme values for better color range
score_cap <- 1.5
pathway_scores_capped <- pmin(pmax(pathway_scores_ordered, -score_cap), score_cap)

# Create column labels: short donor ID + ischemic time
col_labels <- sapply(ordered_ids, function(sid) {
  donor <- sub("^(GTEX-[A-Z0-9]+)-.*", "\\1", sid)
  isch <- sample_info$ischemic_hrs[sample_info$SAMPID == sid]
  sprintf("%s (%.0fh)", donor, isch)
})

# Gaps between CKD triplets (case + 2 controls)
# Mark gaps after every 3rd sample in the CKD-case-ordered section
gaps <- c()
pos <- 0
for (i in seq_len(nrow(match_df))) {
  # Count how many of this triplet's samples appear in ordered_ids
  triplet <- c(match_df$case_id[i], match_df$ctrl1_id[i], match_df$ctrl2_id[i])
  triplet_in <- triplet[triplet %in% ordered_ids[1:length(ordered_ids)]]
  pos <- pos + length(unique(triplet_in))
  # Subtract any that were already counted
}
# Simpler: just put gaps at every 3rd column
gaps_cols <- seq(3, ncol(pathway_scores_capped) - 1, by = 3)

cat("Generating heatmap...\n")

png(file.path(fig_dir, "31_donor_pathway_heatmap.png"),
    width = 22, height = 8, units = "in", res = 300)

pheatmap(
  pathway_scores_capped,
  cluster_rows = TRUE,
  cluster_cols = FALSE,  # keep manual ordering
  annotation_col = ann_df,
  annotation_colors = ann_colors,
  labels_col = col_labels,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-score_cap, score_cap, length.out = 101),
  fontsize_col = 5.5,
  fontsize_row = 10,
  angle_col = 90,
  gaps_col = gaps_cols,
  main = "Core Pathway Scores: Female CKD Cases + 2x Matched Controls\n(covariate-adjusted residuals, z-scored across all 218 females)",
  border_color = NA
)

dev.off()
cat("Saved:", file.path(fig_dir, "31_donor_pathway_heatmap.png"), "\n")


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 4: GTEX-13OVI SPOTLIGHT                                         ║
# ╚════════════════════════════════════════════════════════════════════════════╝

cat("\n", strrep("=", 70), "\n")
cat("SECTION 4: GTEX-13OVI Spotlight (2.0 hrs ischemic)\n")
cat(strrep("=", 70), "\n")

ovi_id <- grep("13OVI", all_sample_ids, value = TRUE)
vi_id <- grep("148VI", all_sample_ids, value = TRUE)
ovi_scores <- pathway_scores[, ovi_id]

# Compare to ALL Short-bin controls using covariate-adjusted residuals
short_ctrl_ids <- tb_f %>%
  filter(MHRNLFLR == 0, isch_bin == "Short") %>%
  pull(SAMPID)
short_ctrl_ids <- intersect(short_ctrl_ids, colnames(resid_z))

# Compute pathway scores for all Short controls from adjusted residuals
short_scores <- matrix(NA, nrow = nrow(pathway_scores),
                       ncol = length(short_ctrl_ids),
                       dimnames = list(rownames(pathway_scores), short_ctrl_ids))
for (p in names(path_labels)[names(path_labels) %in% names(hallmark_genes)]) {
  genes <- intersect(hallmark_genes[[p]], rownames(resid_z))
  if (length(genes) >= 10) {
    pl <- path_labels[p]
    if (pl %in% rownames(short_scores)) {
      short_scores[pl, ] <- colMeans(resid_z[genes, short_ctrl_ids, drop = FALSE])
    }
  }
}

cat(sprintf("\nGTEX-13OVI (2.0h) pathway scores vs Short-bin controls (n=%d):\n",
            length(short_ctrl_ids)))
cat("(covariate-adjusted: AGE, RACE, DTHHRDY, ischemic_hrs, SMRIN, SMCENTER, comorbidities removed)\n")
cat(sprintf("%-20s %8s %10s %10s %8s\n",
            "Pathway", "13OVI", "Ctrl mean", "Ctrl SD", "Pctl"))
cat(strrep("-", 60), "\n")

for (pw in rownames(pathway_scores)) {
  ovi_val <- ovi_scores[pw]
  ctrl_vals <- short_scores[pw, ]
  ctrl_vals <- ctrl_vals[!is.na(ctrl_vals)]
  pctl <- 100 * mean(ctrl_vals <= ovi_val)
  cat(sprintf("%-20s %+8.3f %+10.3f %10.3f %7.0f%%\n",
              pw, ovi_val, mean(ctrl_vals), sd(ctrl_vals), pctl))
}

# GTEX-148VI (3.8h) vs Medium-bin controls
cat(sprintf("\nGTEX-148VI (3.8h) pathway scores vs Medium-bin controls:\n"))
med_ctrl_ids <- tb_f %>%
  filter(MHRNLFLR == 0, isch_bin == "Medium") %>%
  pull(SAMPID)
med_ctrl_ids <- intersect(med_ctrl_ids, colnames(resid_z))

med_scores <- matrix(NA, nrow = nrow(pathway_scores),
                     ncol = length(med_ctrl_ids),
                     dimnames = list(rownames(pathway_scores), med_ctrl_ids))
for (p in names(path_labels)[names(path_labels) %in% names(hallmark_genes)]) {
  genes <- intersect(hallmark_genes[[p]], rownames(resid_z))
  if (length(genes) >= 10) {
    pl <- path_labels[p]
    if (pl %in% rownames(med_scores)) {
      med_scores[pl, ] <- colMeans(resid_z[genes, med_ctrl_ids, drop = FALSE])
    }
  }
}

vi_scores <- pathway_scores[, vi_id]
cat(sprintf("%-20s %8s %10s %10s %8s\n",
            "Pathway", "148VI", "Ctrl mean", "Ctrl SD", "Pctl"))
cat(strrep("-", 60), "\n")
for (pw in rownames(pathway_scores)) {
  vi_val <- vi_scores[pw]
  ctrl_vals <- med_scores[pw, ]
  ctrl_vals <- ctrl_vals[!is.na(ctrl_vals)]
  pctl <- 100 * mean(ctrl_vals <= vi_val)
  cat(sprintf("%-20s %+8.3f %+10.3f %10.3f %7.0f%%\n",
              pw, vi_val, mean(ctrl_vals), sd(ctrl_vals), pctl))
}

# Also compare to Long-bin CKD cases
long_ckd_ids <- tb_f %>%
  filter(MHRNLFLR == 1, isch_bin == "Long") %>%
  pull(SAMPID)

cat(sprintf("\nGTEX-13OVI vs Long-bin CKD cases (n=%d):\n", length(long_ckd_ids)))
cat(sprintf("%-20s %8s %10s %10s\n",
            "Pathway", "13OVI", "Long CKD", "Direction"))
cat(strrep("-", 55), "\n")

for (pw in rownames(pathway_scores)) {
  ovi_val <- ovi_scores[pw]
  long_vals <- pathway_scores[pw, long_ckd_ids]
  long_mean <- mean(long_vals, na.rm = TRUE)
  direction <- ifelse(sign(ovi_val) == sign(long_mean), "same", "OPPOSITE")
  cat(sprintf("%-20s %+8.3f %+10.3f %10s\n",
              pw, ovi_val, long_mean, direction))
}


# ╔════════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 5: SAVE                                                         ║
# ╚════════════════════════════════════════════════════════════════════════════╝

donor_pathway_results <- list(
  pathway_scores = pathway_scores,
  match_df = match_df,
  ordered_ids = ordered_ids,
  sample_info = sample_info,
  ovi_scores = ovi_scores
)
saveRDS(donor_pathway_results, file.path(res_dir, "donor_pathway_scores.rds"))

cat("\nDone.\n")
