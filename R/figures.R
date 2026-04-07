#!/usr/bin/env Rscript
# =============================================================================
# Manuscript Figures — Self-contained
#
# Loads stage1_data_adipose_subcutaneous.rds and performs all sample selection
# and matching internally. DEGs loaded from de_female.rds (stage1_01.R).
# =============================================================================

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(edgeR)
library(MatchIt)
library(ggrepel)

source("kidney/R/functions.R")

res_dir <- "kidney/results"
fig_dir <- "kidney/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── Gene annotation (one-time load) ──────────────────────────
gene_anno <- readRDS("data/annotation/gencode_v47_genes.rds")
gene_map  <- setNames(gene_anno$symbol, rownames(gene_anno))

# =============================================================================
# Data loading & prep
# =============================================================================

cat("=== Data loading & prep ===\n")

dat <- load_sat_data(res_dir = res_dir)
tb <- dat$tb
counts <- dat$counts
tb$MHT2D[is.na(tb$MHT2D)] <- 0L  # rescue 1 female CKD donor (GTEX-1CB4H)

all_vars <- c("AGE", "RACE", "BMI", "DTHHRDY", "ischemic_hrs", "SMRIN",
              "SMCENTER", "SMEXNCRT", "MHABNWBC", "MHLVRDIS", "MHRNLFLR", "MHT2D")
tb <- tb[complete.cases(tb[, all_vars]), ]
counts <- counts[, tb$SAMPID]

cat(sprintf("Samples: %d | Genes: %d\n", nrow(tb), nrow(counts)))
cat(sprintf("Female: %d (CKD=%d) | Male: %d (CKD=%d)\n",
            sum(tb$SEX == "Female"),
            sum(tb$SEX == "Female" & tb$MHRNLFLR == 1),
            sum(tb$SEX == "Male"),
            sum(tb$SEX == "Male" & tb$MHRNLFLR == 1)))

# =============================================================================
# 8-group sample selection
# =============================================================================

cat("\n=== Sample selection (8 groups) ===\n")
set.seed(2)

# ── G1: Female CKD (all available) ───────────────────────────
f_ckd_df  <- tb %>% filter(SEX == "Female", MHRNLFLR == 1)
f_ckd_ids <- f_ckd_df$SAMPID
cat(sprintf("G1 Female CKD: n=%d\n", length(f_ckd_ids)))

# ── G2: Female Control (PSM, ratio=1) ────────────────────────
fem <- tb %>% filter(SEX == "Female", !is.na(MHRNLFLR))
fem <- add_ischemic_bin(fem)

set.seed(42)
psm <- matchit(
  MHRNLFLR ~ AGE + RACE + DTHHRDY + ischemic_hrs + SMRIN,
  data = fem, method = "nearest", ratio = 1,
  exact = ~isch_bin
)
f_ctrl_ids <- match.data(psm) %>%
  filter(MHRNLFLR == 0) %>% pull(SAMPID)
cat(sprintf("G2 Female Control (PSM): n=%d\n", length(f_ctrl_ids)))

# ── G5: Male CKD (Mahalanobis matched to female CKD) ─────────
pool_mckd <- bind_rows(
  f_ckd_df %>% mutate(treat = 1),
  tb %>% filter(SEX == "Male", MHRNLFLR == 1) %>% mutate(treat = 0)
) %>%
  mutate(DTHHRDY = factor(DTHHRDY), RACE = factor(RACE))

set.seed(42)
m_ckd_match <- matchit(
  treat ~ AGE + BMI + ischemic_hrs + SMRIN,
  data = pool_mckd, method = "nearest",
  distance = "mahalanobis",
  exact = ~ DTHHRDY + RACE, ratio = 1, replace = FALSE
)
m_ckd_ids <- match.data(m_ckd_match) %>%
  filter(treat == 0) %>% pull(SAMPID)
cat(sprintf("G5 Male CKD: n=%d\n", length(m_ckd_ids)))

# ── G6: Male Control (Mahalanobis matched to female CKD) ─────
pool_mc <- bind_rows(
  f_ckd_df %>% mutate(treat = 1),
  tb %>% filter(SEX == "Male", MHRNLFLR == 0) %>% mutate(treat = 0)
) %>%
  mutate(DTHHRDY = factor(DTHHRDY), RACE = factor(RACE))

set.seed(42)
m_ctrl_match <- matchit(
  treat ~ AGE + BMI + ischemic_hrs + SMRIN,
  data = pool_mc, method = "nearest",
  distance = "mahalanobis",
  exact = ~ DTHHRDY + RACE, ratio = 1, replace = FALSE
)
m_ctrl_ids <- match.data(m_ctrl_match) %>%
  filter(treat == 0) %>% pull(SAMPID)
cat(sprintf("G6 Male Control: n=%d\n", length(m_ctrl_ids)))

# ── Low-ischemia groups (Mahalanobis on BMI+SMRIN) ───────────
used_ids <- c(f_ckd_ids, f_ctrl_ids, m_ckd_ids, m_ctrl_ids)

match_lo_isch <- function(pool_df, label) {
  pool_df <- pool_df %>%
    filter(!SAMPID %in% used_ids) %>%
    mutate(RACE = factor(RACE))
  combined <- bind_rows(
    f_ckd_df %>% mutate(treat = 1, RACE = factor(RACE)),
    pool_df  %>% mutate(treat = 0)
  )
  set.seed(42)
  m <- matchit(
    treat ~ BMI + SMRIN,
    data = combined, method = "nearest",
    distance = "mahalanobis",
    exact = ~ RACE, ratio = 1, replace = FALSE
  )
  ids <- match.data(m) %>% filter(treat == 0) %>% pull(SAMPID)
  used_ids <<- c(used_ids, ids)
  cat(sprintf("%s: n=%d\n", label, length(ids)))
  ids
}

lo_fem <- tb %>% filter(SEX == "Female", MHRNLFLR == 0, ischemic_hrs <= 3)
lo_mal <- tb %>% filter(SEX == "Male",   MHRNLFLR == 0, ischemic_hrs <= 3)

f_lo_old   <- match_lo_isch(lo_fem %>% filter(AGE >= 50),
                             "G3 Female Lo-Isch Old")
f_lo_young <- match_lo_isch(lo_fem %>% filter(AGE < 50),
                             "G4 Female Lo-Isch Young")
m_lo_old   <- match_lo_isch(lo_mal %>% filter(AGE >= 50),
                             "G7 Male Lo-Isch Old")
m_lo_young <- match_lo_isch(lo_mal %>% filter(AGE < 50),
                             "G8 Male Lo-Isch Young")

# ── Assemble all samples ─────────────────────────────────────
samp_ids <- c(
  f_ckd_ids, f_ctrl_ids, f_lo_old,
  m_ckd_ids, m_ctrl_ids, m_lo_old
)

sex <- ifelse(
  samp_ids %in% c(f_ckd_ids, f_ctrl_ids, f_lo_old),
  "Female", "Male"
)

condition <- case_when(
  samp_ids %in% f_ckd_ids  ~ "CKD",
  samp_ids %in% m_ckd_ids  ~ "CKD",
  samp_ids %in% f_ctrl_ids ~ "Control",
  samp_ids %in% m_ctrl_ids ~ "Control",
  samp_ids %in% c(f_lo_old, m_lo_old) ~ "Lo-Isch"
)

group <- paste(sex, condition)

cat(sprintf("\nTotal: %d samples\n", length(samp_ids)))
print(as.data.frame(table(Group = group)))

# ── Table 1: Group characteristics ────────────────────────────
cat("\n=== Table 1: Group characteristics ===\n")
tab1_df <- tb %>%
  filter(SAMPID %in% samp_ids) %>%
  mutate(Group = factor(group[match(SAMPID, samp_ids)], levels = c(
    "Female CKD", "Female Control", "Female Lo-Isch",
    "Male CKD", "Male Control", "Male Lo-Isch"
  )))

tab1 <- tab1_df %>%
  group_by(Group) %>%
  summarise(
    n           = n(),
    Age         = sprintf("%.0f [%.0f-%.0f]", median(AGE), min(AGE), max(AGE)),
    BMI         = sprintf("%.1f (%.1f)", mean(BMI), sd(BMI)),
    Ischemic_hr = sprintf("%.1f (%.1f)", mean(ischemic_hrs), sd(ischemic_hrs)),
    RIN         = sprintf("%.1f (%.1f)", mean(SMRIN), sd(SMRIN)),
    Hardy_Vent  = sprintf("%d (%.0f%%)",
                          sum(DTHHRDY == "Ventilator"),
                          100 * mean(DTHHRDY == "Ventilator")),
    Hardy_Fast  = sprintf("%d (%.0f%%)",
                          sum(DTHHRDY == "Fast"),
                          100 * mean(DTHHRDY == "Fast")),
    Hardy_Slow  = sprintf("%d (%.0f%%)",
                          sum(DTHHRDY == "Slow"),
                          100 * mean(DTHHRDY == "Slow")),
    Race_White  = sprintf("%d (%.0f%%)",
                          sum(RACE == "White"),
                          100 * mean(RACE == "White")),
    CKD         = sum(MHRNLFLR == 1),
    HTN         = sprintf("%d (%.0f%%)",
                          sum(MHHTN == 1, na.rm = TRUE),
                          100 * mean(MHHTN == 1, na.rm = TRUE)),
    T2D         = sprintf("%d (%.0f%%)",
                          sum(MHT2D == 1, na.rm = TRUE),
                          100 * mean(MHT2D == 1, na.rm = TRUE)),
    .groups = "drop"
  )

print(as.data.frame(tab1), right = FALSE, row.names = FALSE)

# ── Save matched raw counts ────────────────────────────────────
samp_ids_clean <- gsub("-", "_", samp_ids)
samp_ids_short <- paste0("G", sub("^GTEX_([^_]+)_.*", "\\1", samp_ids_clean))
stopifnot(!anyDuplicated(samp_ids_short))
group_clean    <- gsub(" ", "", group)

counts_out <- counts[, samp_ids]
colnames(counts_out) <- samp_ids_short
write.csv(counts_out, file.path(res_dir, "matched_raw_counts.csv"))
cat("Saved:", file.path(res_dir, "matched_raw_counts.csv"), "\n")

# ── iDEP design file ──────────────────────────────────────────
idep_design <- data.frame(row.names = c("ID", "group"),
                          matrix(c(samp_ids_short, group_clean),
                                 nrow = 2, byrow = TRUE))
colnames(idep_design) <- samp_ids_short
write.table(idep_design, file.path(res_dir, "matched_design_idep.csv"),
            sep = ",", col.names = FALSE)
cat("Saved:", file.path(res_dir, "matched_design_idep.csv"), "\n")

# ── Normalize ─────────────────────────────────────────────────
dge <- DGEList(counts = counts[, samp_ids])
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 4)

# ── DEGs from stage1_01.R ─────────────────────────────────────
de_fem <- readRDS(file.path(res_dir, "de_female.rds"))
degs <- de_fem$tt %>%
  filter(adj.P.Val < 0.05) %>%
  pull(gene_id)
cat(sprintf("DEGs: %d\n", length(degs)))

# ── DEG matrix (median-centered) ─────────────────────────────
mat <- logcpm[degs, ]
mat_centered <- mat - apply(mat, 1, median)
cat(sprintf("DEG matrix: %d genes x %d samples\n",
            nrow(mat_centered), ncol(mat_centered)))

# ── Gene modules (used for heatmap row split and Fig 3) ─────
set.seed(42)
gene_dist    <- dist(mat_centered)
gene_hc      <- hclust(gene_dist, method = "complete")
gene_modules <- cutree(gene_hc, k = 2)

mod_means <- sapply(1:2, function(m) {
  mean(mat_centered[
    gene_modules == m,
    colnames(mat_centered) %in% f_ckd_ids
  ])
})
up_mod <- which.max(mod_means)
dn_mod <- which.min(mod_means)
cat(sprintf("Up module: %d genes | Down module: %d genes\n",
            sum(gene_modules == up_mod),
            sum(gene_modules == dn_mod)))

# Gene counts per cluster
n_up <- sum(gene_modules == up_mod)
n_dn <- sum(gene_modules == dn_mod)

row_module <- factor(
  ifelse(gene_modules == up_mod,
         sprintf("Up (%d): Inflammation/EMT", n_up),
         sprintf("Down (%d): Adipogenesis/Metabolic", n_dn)),
  levels = c(sprintf("Up (%d): Inflammation/EMT", n_up),
             sprintf("Down (%d): Adipogenesis/Metabolic", n_dn))
)

# ── Shared aesthetics ────────────────────────────────────────
sex_colors <- c("Female" = "#CC79A7", "Male" = "#009E73")
cond_shapes <- c(
  "CKD" = 16, "Control" = 17, "Lo-Isch" = 0
)
cond_sizes <- c(
  "CKD" = 6, "Control" = 6, "Lo-Isch" = 5
)
cond_colors <- c(
  "CKD" = "#D6604D", "Control" = "#4393C3", "Lo-Isch" = "#33A02C"
)

group_levels <- c(
  "Female CKD", "Female Control", "Female Lo-Isch",
  "Male CKD", "Male Control", "Male Lo-Isch"
)

# =============================================================================
# Figure 0: Top 2000 most variable genes heatmap
# =============================================================================

cat("\n=== Figure 0: Most Variable Genes Heatmap ===\n")

gene_sd <- apply(logcpm[, samp_ids], 1, sd)
top2k   <- names(sort(gene_sd, decreasing = TRUE))[1:2000]
cat(sprintf("Top 2000 HVG SD range: %.2f – %.2f\n",
            min(gene_sd[top2k]), max(gene_sd[top2k])))

mat_hvg <- logcpm[top2k, samp_ids]
mat_hvg_c <- mat_hvg - apply(mat_hvg, 1, median)

ckd_status_hvg  <- ifelse(condition == "CKD", "Yes", "No")
isch_status_hvg <- ifelse(condition %in% c("CKD", "Control"), "Long", "Short")

ha_col_hvg <- HeatmapAnnotation(
  Sex = sex,
  CKD = ckd_status_hvg,
  Ischemia = isch_status_hvg,
  col = list(
    Sex = sex_colors,
    CKD = c("Yes" = "#333333", "No" = "#BBBBBB"),
    Ischemia = c("Long" = "#E69F00", "Short" = "#56B4E9")
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 18),
  annotation_legend_param = list(
    Sex      = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 24)),
    CKD      = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 24)),
    Ischemia = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 24))
  )
)

col_fun_hvg <- colorRamp2(c(-1.5, 0, 1.5), c("#0066FF", "white", "#FF0000"))

set.seed(42)
ht_hvg <- Heatmap(
  mat_hvg_c,
  name = " ",
  col = col_fun_hvg,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 18),
    labels_gp = gpar(fontsize = 16),
    legend_direction = "vertical",
    grid_height = unit(5, "mm")
  ),
  top_annotation = ha_col_hvg,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = TRUE,
  raster_quality = 3
)

out_file_hvg <- file.path(fig_dir, "figures_heatmap_hvg.png")
png(out_file_hvg, width = 16, height = 10, units = "in", res = 300)
draw(ht_hvg, heatmap_legend_side = "right", annotation_legend_side = "right",
     merge_legend = TRUE, adjust_annotation_extension = TRUE,
     legend_gap = unit(12, "mm"))
dev.off()
cat("Saved:", out_file_hvg, "\n")

# =============================================================================
# Figure 1: DEG heatmap
# =============================================================================

cat("\n=== Figure 1: DEG Heatmap ===\n")

ckd_status <- ifelse(condition == "CKD", "Yes", "No")
isch_status <- ifelse(condition %in% c("CKD", "Control"), "Long", "Short")

# Group labels for column_split (single line)
group_labels <- c(
  "Female CKD"     = "Female CKD",
  "Female Control"  = "Female Control",
  "Female Lo-Isch" = "Female Lo-Isch",
  "Male CKD"       = "Male CKD",
  "Male Control"    = "Male Control",
  "Male Lo-Isch"   = "Male Lo-Isch"
)

col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#0066FF", "white", "#FF0000"))

# Marker gene labels
marker_symbols <- c(
  "LEP", "ADIPOQ", "PPARG", "PLIN1", "FASN", "LPL",
  "IL6", "CCL2", "CXCL10", "EGR1"
)
row_symbols <- gene_map[rownames(mat_centered)]
marker_idx  <- which(row_symbols %in% marker_symbols)

ha_row <- rowAnnotation(
  Markers = anno_mark(
    at = marker_idx,
    labels = row_symbols[marker_idx],
    labels_gp = gpar(fontsize = 16, fontface = "italic")
  )
)

set.seed(42)
ht <- Heatmap(
  mat_centered,
  name = "logCPM",
  col = col_fun,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 18),
    labels_gp = gpar(fontsize = 16),
    legend_direction = "vertical",
    grid_height = unit(5, "mm"),
    title_gap = unit(8, "mm")
  ),
  right_annotation = ha_row,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  row_split = row_module,
  column_split = factor(group, levels = group_levels),
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_title_rot = 90,
  row_title_gp = gpar(fontsize = 18, fontface = "bold"),
  column_title = group_labels[group_levels],
  column_title_gp = gpar(fontsize = 18, fontface = "bold"),
  use_raster = TRUE,
  raster_quality = 3
)

out_file <- file.path(fig_dir, "figures_1_heatmap.png")

library(grid)
lgd <- Legend(
  col_fun = col_fun,
  title = "logCPM",
  title_gp = gpar(fontsize = 18),
  labels_gp = gpar(fontsize = 16),
  grid_height = unit(5, "mm"),
  title_gap = unit(5, "mm"),
  direction = "vertical"
)

png(out_file, width = 16, height = 10, units = "in", res = 300)
draw(ht, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
draw(lgd, x = unit(0.99, "npc"), y = unit(0.22, "npc"), just = c("right", "center"))
dev.off()
cat("Saved:", out_file, "\n")

# ── Shared scatter plot function ────────────────────────────
# df must have columns: X, Y, Sex, Condition, Group
# cents must have columns: X_mean, Y_mean, Sex, Condition
plot_group_scatter <- function(df, cents, x_lab = "X", y_lab = "Y",
                               trend = TRUE) {
  p <- ggplot()
  if (trend) {
    p <- p + geom_smooth(
      data = cents,
      aes(x = X_mean, y = Y_mean),
      method = "lm", color = "grey40", linewidth = 0.6,
      linetype = "dashed", se = FALSE, fullrange = TRUE
    )
  }
  p <- p +
    geom_point(
      data = df,
      aes(x = X, y = Y, color = Sex, shape = Condition),
      size = 1.5, alpha = 0.4
    ) +
    geom_point(
      data = cents,
      aes(x = X_mean, y = Y_mean,
          color = Sex, shape = Condition,
          size = Condition),
      stroke = 1.5
    ) +
    geom_text_repel(
      data = cents,
      aes(x = X_mean, y = Y_mean,
          label = paste(Sex, Condition),
          color = Sex),
      size = 4.2, fontface = "plain",
      show.legend = FALSE,
      nudge_x = 0.4,
      direction = "y",
      hjust = 0,
      segment.color = NA,
      max.overlaps = Inf,
      seed = 42
    ) +
    scale_color_manual(values = sex_colors) +
    scale_shape_manual(values = cond_shapes) +
    scale_size_manual(values = cond_sizes) +
    guides(shape = guide_legend(override.aes = list(size = 3)),
           size = "none") +
    labs(x = x_lab, y = y_lab, title = NULL) +
    theme_classic(base_size = 14) +
    theme(panel.border = element_rect(
      color = "black", fill = NA, linewidth = 0.5
    ))
  p
}

# =============================================================================
# Figure 3: Two-module scatter
# =============================================================================

cat("\n=== Figure 3: Two-module scatter ===\n")

score_up <- colMeans(mat_centered[gene_modules == up_mod, ])
score_dn <- colMeans(mat_centered[gene_modules == dn_mod, ])

scatter_df <- data.frame(
  X = score_up,
  Y = score_dn,
  Sex = sex, Condition = condition,
  Group = group,
  SAMPID = samp_ids
)

centroids3 <- scatter_df %>%
  group_by(Sex, Condition, Group) %>%
  summarise(
    X_mean = mean(X), Y_mean = mean(Y),
    X_se = sd(X) / sqrt(n()),
    Y_se = sd(Y) / sqrt(n()),
    n = n(), .groups = "drop"
  )

# Two female CKD donors with shortest ischemic time
short_isch_ids <- tb %>%
  filter(SAMPID %in% f_ckd_ids) %>%
  arrange(ischemic_hrs) %>%
  slice_head(n = 2) %>%
  pull(SAMPID)
highlight_df <- scatter_df %>% filter(SAMPID %in% short_isch_ids)

p_fig3 <- plot_group_scatter(
  scatter_df, centroids3,
  x_lab = sprintf("Up (%d): Inflammation/EMT score", n_up),
  y_lab = sprintf("Down (%d): Adipogenesis/Metabolic score", n_dn)
) +
  geom_point(data = highlight_df, aes(x = X, y = Y),
             shape = 1, size = 3.5, stroke = 0.7, color = "grey50")


out_file3 <- file.path(fig_dir, "figures_module_scatter.png")
ggsave(out_file3, p_fig3, width = 8, height = 6, dpi = 300)
cat("Saved:", out_file3, "\n")

# =============================================================================
# Figure 4: PCA using ALL filtered genes
# =============================================================================

cat("\n=== Figure 4: PCA (all filtered genes) ===\n")

mat_all   <- logcpm
mat_all_c <- mat_all - apply(mat_all, 1, median)
cat(sprintf("PCA matrix: %d genes x %d samples\n",
            nrow(mat_all_c), ncol(mat_all_c)))

pca_all <- prcomp(t(mat_all_c), center = TRUE, scale. = FALSE)
vpct <- round(
  100 * summary(pca_all)$importance[2, 1:2], 1
)
cat(sprintf("PC1: %.1f%% | PC2: %.1f%%\n", vpct[1], vpct[2]))

pca_all_df <- data.frame(
  X = pca_all$x[, 1], Y = pca_all$x[, 2],
  Sex = sex, Condition = condition,
  Group = group
)

centroids4 <- pca_all_df %>%
  group_by(Sex, Condition, Group) %>%
  summarise(
    X_mean = mean(X), Y_mean = mean(Y),
    X_se = sd(X) / sqrt(n()),
    Y_se = sd(Y) / sqrt(n()),
    n = n(), .groups = "drop"
  )

p_fig4 <- ggplot() +
  geom_point(
    data = pca_all_df,
    aes(x = X, y = Y, color = Sex, shape = Condition),
    size = 1.5, alpha = 0.4
  ) +
  geom_point(
    data = centroids4,
    aes(x = X_mean, y = Y_mean,
        color = Sex, shape = Condition,
        size = Condition),
    stroke = 1.5
  ) +
  scale_color_manual(values = sex_colors) +
  scale_shape_manual(values = cond_shapes) +
  scale_size_manual(values = cond_sizes) +
  guides(shape = guide_legend(override.aes = list(size = 3)),
         size = "none") +
  labs(
    x = sprintf("PC1 (%.1f%%)", vpct[1]),
    y = sprintf("PC2 (%.1f%%)", vpct[2]),
    title = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(
    color = "black", fill = NA, linewidth = 0.5
  ))

out_file4 <- file.path(fig_dir, "figures_pca_all_genes.png")
ggsave(out_file4, p_fig4, width = 8, height = 6, dpi = 300)
cat("Saved:", out_file4, "\n")

# =============================================================================
# Figure 4b: PCA loading GSEA (Hallmark) — PC1 vs PC2
# =============================================================================

cat("\n=== Figure 4b: PCA loading GSEA (DEGs) ===\n")

library(msigdbr)
library(fgsea)
library(ggrepel)

hallmark_list <- load_hallmark_sets()

# PCA on DEG matrix
pca_deg <- prcomp(t(mat_centered), center = TRUE, scale. = FALSE)
vpct_deg <- round(100 * summary(pca_deg)$importance[2, 1:2], 1)
cat(sprintf("DEG PCA — PC1: %.1f%% | PC2: %.1f%%\n", vpct_deg[1], vpct_deg[2]))

# Extract and sort loadings from DEG PCA
load_pc1 <- sort(pca_deg$rotation[, 1], decreasing = TRUE)
load_pc2 <- sort(pca_deg$rotation[, 2], decreasing = TRUE)

set.seed(42)
gsea_pc1 <- fgsea(pathways = hallmark_list, stats = load_pc1,
                   minSize = 15, maxSize = 500, nPermSimple = 10000)
set.seed(42)
gsea_pc2 <- fgsea(pathways = hallmark_list, stats = load_pc2,
                   minSize = 15, maxSize = 500, nPermSimple = 10000)

cat(sprintf("PC1 GSEA: %d pathways FDR < 0.05\n", sum(gsea_pc1$padj < 0.05)))
cat(sprintf("PC2 GSEA: %d pathways FDR < 0.05\n", sum(gsea_pc2$padj < 0.05)))

# Merge and filter to significant in either
gsea_merged <- merge(
  gsea_pc1[, .(pathway, NES_PC1 = NES, padj_PC1 = padj)],
  gsea_pc2[, .(pathway, NES_PC2 = NES, padj_PC2 = padj)],
  by = "pathway"
)
gsea_merged$sig <- case_when(
  gsea_merged$padj_PC1 < 0.001 & gsea_merged$padj_PC2 < 0.001 ~ "Both",
  gsea_merged$padj_PC1 < 0.001 ~ "PC1 only",
  gsea_merged$padj_PC2 < 0.001 ~ "PC2 only",
  TRUE ~ "NS"
)
gsea_sig <- gsea_merged %>%
  filter(sig != "NS", abs(NES_PC1) > 2 | abs(NES_PC2) > 2)

cat(sprintf("Significant in either: %d pathways\n", nrow(gsea_sig)))

gsea_sig$label <- clean_pathway(gsea_sig$pathway)

sig_colors <- c("Both" = "#7570B3", "PC1 only" = "#D95F02", "PC2 only" = "#1B9E77")

p_fig4b <- ggplot(gsea_sig, aes(x = NES_PC1, y = NES_PC2, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3) +
  geom_text_repel(aes(label = label), size = 3.2, max.overlaps = 30,
                  segment.color = "grey50", segment.size = 0.3) +
  scale_color_manual(values = sig_colors, name = "FDR < 0.001") +
  labs(x = sprintf("NES on PC1 (%.1f%%)", vpct_deg[1]),
       y = sprintf("NES on PC2 (%.1f%%)", vpct_deg[2])) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "top"
  )

out_file4b <- file.path(fig_dir, "figures_pca_loading_gsea.png")
ggsave(out_file4b, p_fig4b, width = 10, height = 9, dpi = 300)
cat("Saved:", out_file4b, "\n")

# Save results table
write.csv(gsea_sig, file.path(res_dir, "pca_loading_gsea.csv"), row.names = FALSE)

# =============================================================================
# Figure 5+: Marker gene scatter plots
# =============================================================================

plot_marker_scatter <- function(x_sym, y_sym, out_path,
                                logcpm, gene_anno, sex, condition) {
  x_id <- rownames(gene_anno)[gene_anno$symbol == x_sym]
  y_id <- rownames(gene_anno)[gene_anno$symbol == y_sym]
  stopifnot(length(x_id) == 1, length(y_id) == 1)

  df <- data.frame(
    X = logcpm[x_id, ], Y = logcpm[y_id, ],
    Sex = sex, Condition = condition,
    Group = paste(sex, condition)
  )

  # Remove 2 most extreme samples from each end of X and Y
  x_rank <- rank(df$X)
  y_rank <- rank(df$Y)
  keep <- x_rank > 2 & x_rank <= (nrow(df) - 2) &
          y_rank > 2 & y_rank <= (nrow(df) - 2)
  df <- df[keep, ]

  cents <- df %>%
    group_by(Sex, Condition, Group) %>%
    summarise(X_mean = mean(X), Y_mean = mean(Y),
              n = n(), .groups = "drop")

  r_val <- cor(df$X, df$Y, method = "spearman")
  cat(sprintf("  %s vs %s: rho = %.3f (n=%d after trimming)\n",
              x_sym, y_sym, r_val, nrow(df)))

  p <- plot_group_scatter(
    df, cents,
    x_lab = bquote(italic(.(x_sym)) ~ "(logCPM)"),
    y_lab = bquote(italic(.(y_sym)) ~ "(logCPM)")
  )

  ggsave(out_path, p, width = 8, height = 6, dpi = 300)
  cat("  Saved:", out_path, "\n")
  invisible(p)
}

cat("\n=== Figure 5: Marker scatter plots ===\n")

marker_pairs <- list(
  c("PTGS2",  "ADIPOQ"),
  c("ACTA2",  "LEP"),
  c("INHBA",  "PLIN1"),
  c("CCL2",   "LPL"),
  c("IL6",    "FASN")
)

for (i in seq_along(marker_pairs)) {
  x_sym <- marker_pairs[[i]][1]
  y_sym <- marker_pairs[[i]][2]
  fname <- sprintf("figures_%s_vs_%s.png",
                   tolower(x_sym), tolower(y_sym))
  plot_marker_scatter(
    x_sym, y_sym, file.path(fig_dir, fname),
    logcpm, gene_anno, sex, condition
  )
}

# =============================================================================
# Figure 6: Marker gene boxplots
# =============================================================================

library(ggsignif)
library(patchwork)

make_gene_boxplot <- function(sym, logcpm, gene_df, ylim_max = NULL) {
  eid <- names(gene_map)[gene_map == sym]
  eid <- intersect(eid, rownames(logcpm))
  if (length(eid) == 0) return(NULL)
  eid <- eid[1]

  df <- gene_df %>%
    mutate(logCPM = logcpm[eid, SAMPID]) %>%
    mutate(Condition = factor(Condition, levels = c("CKD", "Control", "Lo-Isch")))

  y_max <- max(df$logCPM, na.rm = TRUE)
  y_rng <- diff(range(df$logCPM, na.rm = TRUE))

  p <- ggplot(df, aes(x = Condition, y = logCPM, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(alpha = 0.5, width = 0.2, size = 0.8, shape = 16) +
    facet_wrap(~Sex) +
    scale_fill_manual(values = cond_colors) +
    geom_signif(
      comparisons = list(c("CKD", "Control"), c("Control", "Lo-Isch")),
      test = "wilcox.test", map_signif_level = TRUE,
      textsize = 3.5, step_increase = 0.12,
      y_position = y_max + y_rng * 0.04
    ) +
    labs(
      y = bquote(italic(.(sym))),
      x = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 13),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  if (!is.null(ylim_max)) {
    p <- p + coord_cartesian(ylim = c(NA, ylim_max))
  }
  p
}

cat("\n=== Figure 6: Marker gene boxplots ===\n")

gene_df <- data.frame(
  SAMPID = samp_ids, Sex = sex, Condition = condition,
  Group = group, stringsAsFactors = FALSE
)

# Combined ADIPOQ + IL6 panel
p_adipoq <- make_gene_boxplot("ADIPOQ", logcpm, gene_df, ylim_max = 15)
p_il6    <- make_gene_boxplot("IL6",    logcpm, gene_df, ylim_max = 12)

p_adipoq <- p_adipoq + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  plot.margin = margin(2, 5, 0, 5))
p_il6 <- p_il6 +
  theme(strip.text = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold"),
        plot.margin = margin(0, 5, 2, 5))
p_left <- p_adipoq / p_il6
p_fig3_noleg <- p_fig3 + theme(legend.position = "none")
p_combined <- (p_left | p_fig3_noleg) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(file.path(fig_dir, "figures_2_combined.png"), p_combined,
       width = 10.4, height = 5.6, dpi = 300)
cat("  Saved: figures_adipoq_il6.png\n")

# Supplementary: remaining markers in one panel
supp_genes <- c("LEP", "PPARG", "PLIN1", "FASN", "LPL", "CCL2", "CXCL10", "EGR1",
                 "ACTA2")
supp_plots <- lapply(supp_genes, make_gene_boxplot, logcpm = logcpm, gene_df = gene_df)
names(supp_plots) <- supp_genes

p_supp <- wrap_plots(supp_plots, ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave(file.path(fig_dir, "figures_boxplots_all_markers.png"), p_supp,
       width = 14, height = 13, dpi = 300)
cat("  Saved: figures_boxplots_all_markers.png\n")

# =============================================================================
# Staircase pattern analysis
# =============================================================================

cat("\n=== Staircase pattern analysis ===\n")

staircase_df <- do.call(rbind, lapply(degs, function(g) {
  x_ckd  <- logcpm[g, f_ckd_ids]
  x_ctrl <- logcpm[g, f_ctrl_ids]
  x_lo   <- logcpm[g, f_lo_old]

  m_ckd  <- mean(x_ckd)
  m_ctrl <- mean(x_ctrl)
  m_lo   <- mean(x_lo)

  p1 <- wilcox.test(x_ckd, x_ctrl)$p.value
  p2 <- wilcox.test(x_ctrl, x_lo)$p.value

  data.frame(
    gene_id = g,
    symbol  = gene_map[g],
    cluster = ifelse(gene_modules[g] == up_mod, "Up", "Down"),
    mean_CKD = m_ckd, mean_Ctrl = m_ctrl, mean_LoIsch = m_lo,
    p_CKD_vs_Ctrl = p1, p_Ctrl_vs_LoIsch = p2,
    stringsAsFactors = FALSE
  )
}))

# Staircase: both steps significant AND monotonic in expected direction
staircase_df$staircase <- with(staircase_df,
  ifelse(cluster == "Up",
    p_CKD_vs_Ctrl < 0.05 & p_Ctrl_vs_LoIsch < 0.05 &
      mean_CKD > mean_Ctrl & mean_Ctrl > mean_LoIsch,
    p_CKD_vs_Ctrl < 0.05 & p_Ctrl_vs_LoIsch < 0.05 &
      mean_CKD < mean_Ctrl & mean_Ctrl < mean_LoIsch
  )
)

# Summary
up_genes <- staircase_df[staircase_df$cluster == "Up", ]
dn_genes <- staircase_df[staircase_df$cluster == "Down", ]

cat(sprintf("Up cluster:   %d / %d staircase (%.1f%%)\n",
            sum(up_genes$staircase), nrow(up_genes),
            100 * mean(up_genes$staircase)))
cat(sprintf("Down cluster: %d / %d staircase (%.1f%%)\n",
            sum(dn_genes$staircase), nrow(dn_genes),
            100 * mean(dn_genes$staircase)))

write.csv(staircase_df, file.path(res_dir, "staircase_pattern.csv"),
          row.names = FALSE)
cat("Saved:", file.path(res_dir, "staircase_pattern.csv"), "\n")
