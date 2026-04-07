# =============================================================================
# Shared constants
# =============================================================================

#' Core Hallmark Pathways
#'
#' The canonical set of 10 core Hallmark pathways used across analyses.
CORE_PATHWAYS <- c(
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_HYPOXIA",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"
)

#' Display Labels for Hallmark Pathways
#'
#' Superset lookup: covers CORE_PATHWAYS plus additional pathways used in
#' specific analyses (e.g., ischemic pathway scores).
PATHWAY_LABELS <- c(
  HALLMARK_ADIPOGENESIS = "Adipogenesis",
  HALLMARK_OXIDATIVE_PHOSPHORYLATION = "Oxphos",
  HALLMARK_FATTY_ACID_METABOLISM = "FA Metabolism",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFa/NF-kB",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
  HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory",
  HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-gamma",
  HALLMARK_HYPOXIA = "Hypoxia",
  HALLMARK_COMPLEMENT = "Complement",
  HALLMARK_ESTROGEN_RESPONSE_LATE = "Estrogen (late)",
  HALLMARK_IL6_JAK_STAT3_SIGNALING = "IL6/JAK/STAT3",
  HALLMARK_APOPTOSIS = "Apoptosis"
)

# =============================================================================
# Shared helpers
# =============================================================================

#' Print a Section Header
#'
#' @param title Section title string.
#' @param char Separator character.
#' @param width Total width of separator line.
section_header <- function(title, char = "=", width = 70) {
  cat("\n", strrep(char, width), "\n")
  cat(title, "\n")
  cat(strrep(char, width), "\n")
}

#' Load and Prepare SAT Data
#'
#' Loads the stage1 adipose subcutaneous data, enforces factor types, converts
#' binary MH flags to integer, and optionally subsets to complete cases.
#'
#' @param res_dir Directory containing the RDS file.
#' @param mh_flags Character vector of MH flag columns to convert to integer.
#' @param complete_vars Optional character vector of columns for complete-case
#'   subsetting. If NULL, no subsetting is performed.
#'
#' @return A list with `tb` (sample metadata), `counts` (gene count matrix),
#'   and `joint_factors` (the standard 4-factor vector).
load_sat_data <- function(res_dir = "kidney/results",
                          mh_flags = c("MHABNWBC", "MHLVRDIS", "MHRNLFLR", "MHT2D"),
                          complete_vars = NULL) {
  s1 <- readRDS(file.path(res_dir, "stage1_data_adipose_subcutaneous.rds"))
  tb <- s1$tissue_df
  counts <- s1$counts

  # Binary MH flags as integer

  for (f in mh_flags) {
    if (f %in% names(tb)) tb[[f]] <- as.integer(tb[[f]])
  }

  joint_factors <- c("MHABNWBC", "MHLVRDIS", "MHRNLFLR", "MHT2D")

  # Complete cases
  if (!is.null(complete_vars)) {
    tb <- tb[complete.cases(tb[, complete_vars]), ]
    counts <- counts[, tb$SAMPID]
  }

  list(tb = tb, counts = counts, joint_factors = joint_factors)
}

#' Summarize DEG Counts
#'
#' Computes number of DEGs (up/down) and pi0 estimate from a topTable result.
#'
#' @param tt A topTable data frame.
#' @param fdr FDR threshold.
#'
#' @return A list with n_deg, n_up, n_down, and pi0.
summarize_degs <- function(tt, fdr = 0.05) {
  list(
    n_deg  = sum(tt$adj.P.Val < fdr),
    n_up   = sum(tt$adj.P.Val < fdr & tt$logFC > 0),
    n_down = sum(tt$adj.P.Val < fdr & tt$logFC < 0),
    pi0    = tryCatch(qvalue::pi0est(tt$P.Value)$pi0, error = function(e) NA_real_)
  )
}

#' Compare Two GSEA Results by NES Correlation
#'
#' Joins two GSEA result tables by pathway and computes NES Spearman
#' correlation. Optionally prints a core pathway comparison table.
#'
#' @param gsea1 First GSEA result (from `fgsea`).
#' @param gsea2 Second GSEA result (from `fgsea`).
#' @param label1 Display label for first result.
#' @param label2 Display label for second result.
#' @param core_paths Optional character vector of pathways for a side-by-side
#'   comparison table. If NULL, no table is printed.
#' @param print_table Whether to print the core pathway table.
#'
#' @return A list with `merged` (joined data frame), `rho_nes` (Spearman rho),
#'   and `core_table` (filtered to core_paths, if provided).
compare_gsea_results <- function(gsea1, gsea2,
                                  label1 = "A", label2 = "B",
                                  core_paths = NULL,
                                  print_table = TRUE) {
  merged <- dplyr::inner_join(
    gsea1 %>% dplyr::select(pathway, NES1 = NES, padj1 = padj),
    gsea2 %>% dplyr::select(pathway, NES2 = NES, padj2 = padj),
    by = "pathway"
  )

  rho_nes <- cor(merged$NES1, merged$NES2, method = "spearman")

  # Rename columns to labels
  names(merged)[names(merged) == "NES1"]  <- paste0("NES_", label1)
  names(merged)[names(merged) == "padj1"] <- paste0("padj_", label1)
  names(merged)[names(merged) == "NES2"]  <- paste0("NES_", label2)
  names(merged)[names(merged) == "padj2"] <- paste0("padj_", label2)

  core_table <- NULL
  if (!is.null(core_paths) && print_table) {
    nes1_col <- paste0("NES_", label1)
    nes2_col <- paste0("NES_", label2)
    padj1_col <- paste0("padj_", label1)
    padj2_col <- paste0("padj_", label2)

    core_table <- merged %>%
      dplyr::filter(pathway %in% core_paths) %>%
      dplyr::mutate(pathway = gsub("HALLMARK_", "", pathway)) %>%
      dplyr::arrange(.data[[nes1_col]])

    cat(sprintf("\n  %-40s %8s %8s %9s %9s\n",
                "Pathway", label1, label2,
                paste0("p_", label1), paste0("p_", label2)))
    cat(sprintf("  %s\n", strrep("-", 76)))
    for (i in seq_len(nrow(core_table))) {
      cat(sprintf("  %-40s %+8.2f %+8.2f %9.1e %9.1e\n",
                  core_table$pathway[i],
                  core_table[[nes1_col]][i], core_table[[nes2_col]][i],
                  core_table[[padj1_col]][i], core_table[[padj2_col]][i]))
    }
  }

  list(merged = merged, rho_nes = rho_nes, core_table = core_table)
}

# =============================================================================
# Coefficient resolution helpers
# =============================================================================

#' Resolve a Required Coefficient by Exact Name
#'
#' Validates that a design matrix contains an expected coefficient column and
#' returns its exact name.
#'
#' @param design A model matrix.
#' @param coef_name A single coefficient name expected in `colnames(design)`.
#' @param context A short label used in error messages.
#'
#' @return A length-1 character vector containing `coef_name`.
get_required_coef <- function(design, coef_name, context) {
  if (!(coef_name %in% colnames(design))) {
    stop(sprintf(
      "%s: expected coefficient '%s' not found. Available columns: %s",
      context, coef_name, paste(colnames(design), collapse = ", ")
    ))
  }
  coef_name
}

#' Resolve a Required Coefficient by Prefix
#'
#' Finds the unique coefficient whose name starts with a requested prefix.
#'
#' @param design A model matrix.
#' @param coef_prefix A single coefficient prefix expected in `colnames(design)`.
#' @param context A short label used in error messages.
#'
#' @return A length-1 character vector containing the matched coefficient name.
get_required_coef_prefix <- function(design, coef_prefix, context) {
  matches <- grep(paste0("^", coef_prefix), colnames(design), value = TRUE)
  if (length(matches) != 1) {
    stop(sprintf(
      "%s: expected exactly 1 coefficient starting with '%s'; found %d. Candidates: %s",
      context, coef_prefix, length(matches),
      if (length(matches) > 0) paste(matches, collapse = ", ") else "<none>"
    ))
  }
  matches[[1]]
}

#' Resolve a Required Sex Interaction Coefficient
#'
#' Finds the unique `SEX` interaction term for a requested phenotype term in a
#' design matrix. The interaction must be represented as a two-way `:` term.
#'
#' @param design A model matrix.
#' @param term_name The non-sex term expected to participate in the interaction.
#' @param context A short label used in error messages.
#'
#' @return A length-1 character vector containing the interaction coefficient
#'   name.
get_required_interaction_coef <- function(design, term_name, context) {
  interaction_cols <- colnames(design)[grepl(":", colnames(design), fixed = TRUE)]
  matches <- interaction_cols[vapply(strsplit(interaction_cols, ":", fixed = TRUE), function(parts) {
    length(parts) == 2 && term_name %in% parts && any(startsWith(parts, "SEX"))
  }, logical(1))]

  if (length(matches) != 1) {
    stop(sprintf(
      "%s: expected exactly 1 SEX interaction coefficient for '%s'; found %d. Candidates: %s",
      context, term_name, length(matches),
      if (length(matches) > 0) paste(matches, collapse = ", ") else "<none>"
    ))
  }

  matches[[1]]
}

#' Build a Full-Rank Design Matrix
#'
#' Builds a model matrix from a formula and data frame, then drops aliased
#' columns when the design is rank deficient.
#'
#' @param formula A model formula.
#' @param data A data frame used to evaluate `formula`.
#'
#' @return A numeric design matrix.
build_design <- function(formula, data) {
  design <- model.matrix(formula, data = data)
  qr_obj <- qr(design)
  if (qr_obj$rank < ncol(design)) {
    design <- design[, qr_obj$pivot[seq_len(qr_obj$rank)], drop = FALSE]
  }
  design
}

#' Run a limma-voom Differential Expression Fit
#'
#' Runs TMM normalization, `voom`, `lmFit`, `eBayes`, and `topTable` for a
#' single coefficient.
#'
#' @param counts A genes x samples count matrix.
#' @param design A model matrix.
#' @param coef The coefficient name or index passed to `limma::topTable()`.
#'
#' @return A list containing the `topTable` output and fitted analysis objects.
run_voom_de <- function(counts, design, coef) {
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  voom_obj <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(voom_obj, design))
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  tt$gene_id <- rownames(tt)

  list(tt = tt, fit = fit, voom_obj = voom_obj, dge = dge)
}

#' Run a limma-voom Differential Expression Analysis from Sample Metadata
#'
#' Subsets a counts matrix to the requested samples, builds the design matrix,
#' resolves a coefficient, and returns both summary metrics and fitted objects.
#'
#' @param tb_sub Sample metadata table for the analysis subset.
#' @param counts_all A genes x samples count matrix containing `tb_sub$SAMPID`.
#' @param formula A model formula.
#' @param coef_name Coefficient name or prefix to test.
#' @param coef_match Either `"exact"` or `"prefix"`.
#' @param context A short label used in error messages.
#'
#' @return A list with top-table results, DEG summary metrics, and fitted
#'   analysis objects.
run_limma_analysis <- function(tb_sub, counts_all, formula, coef_name,
                               coef_match = c("exact", "prefix"),
                               context = deparse(formula)) {
  coef_match <- match.arg(coef_match)
  counts_sub <- counts_all[, tb_sub$SAMPID, drop = FALSE]
  design <- build_design(formula, tb_sub)

  coef <- switch(
    coef_match,
    exact = get_required_coef(design, coef_name, context),
    prefix = get_required_coef_prefix(design, coef_name, context)
  )

  de_res <- run_voom_de(counts_sub, design, coef)
  tt <- de_res$tt

  list(
    tt = tt,
    fit = de_res$fit,
    voom_obj = de_res$voom_obj,
    dge = de_res$dge,
    design = design,
    n_deg = sum(tt$adj.P.Val < 0.05),
    n_up = sum(tt$adj.P.Val < 0.05 & tt$logFC > 0),
    n_down = sum(tt$adj.P.Val < 0.05 & tt$logFC < 0),
    n_samples = nrow(tb_sub)
  )
}

#' Create a Highlighted Point Scatter Plot
#'
#' Draws a scatter plot with non-significant points in the background and
#' highlighted categories overlaid on top. The x and y axes are forced to share
#' the same numeric limits.
#'
#' @param data A data frame containing plot columns.
#' @param x_var Name of the x-axis column.
#' @param y_var Name of the y-axis column.
#' @param color_var Name of the color/grouping column.
#' @param palette A named character vector of colors for `color_var`.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param title Plot title.
#' @param legend_name Legend title.
#'
#' @return A ggplot object.
make_point_scatter <- function(data, x_var, y_var, color_var, palette,
                               x_label, y_label, title, legend_name) {
  vals <- c(data[[x_var]], data[[y_var]])
  vals <- vals[is.finite(vals)]
  lims <- range(vals)
  pad <- 0.05 * diff(lims)
  if (!is.finite(pad) || pad == 0) pad <- 0.1
  lims <- lims + c(-pad, pad)

  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[color_var]])) +
    geom_point(data = data %>% filter(.data[[color_var]] == "Neither"),
               size = 0.3, alpha = 0.3) +
    geom_point(data = data %>% filter(.data[[color_var]] != "Neither"),
               size = 0.8, alpha = 0.7) +
    scale_color_manual(values = palette, name = legend_name) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = x_label, y = y_label, title = title) +
    coord_fixed(xlim = lims, ylim = lims) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}

#' Fit a Sex-Stratified Differential Expression Model
#'
#' Subsets samples to one sex stratum, builds the design matrix, runs the
#' differential expression model, and returns both summary metrics and fitted
#' objects.
#'
#' @param tb Sample metadata table.
#' @param counts A genes x samples count matrix.
#' @param joint_factors Character vector of phenotype covariates in the model.
#' @param sex_value The sex stratum to analyze.
#' @param sex_labels Named character vector mapping `sex_value` to display
#'   labels.
#'
#' @return A list with top-table results, DEG summary metrics, and fitted
#'   objects.
fit_sex_stratum <- function(tb, counts, joint_factors, sex_value, sex_labels) {
  sex_label <- sex_labels[as.character(sex_value)]
  tb_s <- tb %>% filter(SEX == sex_value)
  counts_s <- counts[, tb_s$SAMPID]

  cat(sprintf("\n--- %s (SEX=%s) ---\n", sex_label, sex_value))
  cat(sprintf("  Samples: %d | CKD cases: %d | Controls: %d\n",
              nrow(tb_s), sum(tb_s$MHRNLFLR == 1), sum(tb_s$MHRNLFLR == 0)))

  form <- as.formula(paste(
    "~ AGE + RACE + BMI + DTHHRDY + ischemic_hrs + SMRIN + SMCENTER + SMEXNCRT +",
    paste(joint_factors, collapse = " + ")))

  design <- build_design(form, tb_s)

  cat(sprintf("  Design: %d cols, %d rows\n", ncol(design), nrow(tb_s)))

  f_coef <- get_required_coef(design, "MHRNLFLR",
                              sprintf("fit_sex_stratum(%s)", sex_label))

  de_res <- run_voom_de(counts_s, design, f_coef)
  tt <- de_res$tt

  n_deg <- sum(tt$adj.P.Val < 0.05)
  n_up <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 0)
  n_down <- sum(tt$adj.P.Val < 0.05 & tt$logFC < 0)
  pi0 <- tryCatch(pi0est(tt$P.Value)$pi0, error = function(e) NA)

  cat(sprintf("  DEGs (FDR<0.05): %d (%d up, %d down)\n", n_deg, n_up, n_down))
  cat(sprintf("  Pi0: %.3f (%.1f%% of genes affected)\n", pi0, 100 * (1 - pi0)))

  list(tt = tt, n_deg = n_deg, n_up = n_up, n_down = n_down,
       pi0 = pi0, sex_value = sex_value, sex_label = sex_label,
       n_samples = nrow(tb_s), n_cases = sum(tb_s$MHRNLFLR == 1),
       fit = de_res$fit, voom_obj = de_res$voom_obj, design = design,
       dge = de_res$dge)
}

#' Format a Numeric Vector as Median and IQR
#'
#' Formats a numeric vector as `median [Q1-Q3]` after dropping missing values.
#'
#' @param x A numeric vector.
#'
#' @return A length-1 character string.
fmt_median_iqr <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.1f [%.1f-%.1f]", median(x), quantile(x, 0.25), quantile(x, 0.75))
}

#' Load Hallmark Gene Sets
#'
#' Loads Hallmark gene sets from `msigdbr` and returns them as a named list.
#'
#' @param species Species name passed to `msigdbr`.
#' @param collection Collection code passed to `msigdbr`.
#' @param gene_col Column containing gene identifiers.
#' @param filter_empty Whether to drop empty identifiers.
#' @param pathways Optional character vector of pathway names to retain.
#'
#' @return A named list mapping pathway names to gene vectors.
load_hallmark_sets <- function(species = "Homo sapiens", collection = "H",
                               gene_col = "ensembl_gene",
                               filter_empty = TRUE, pathways = NULL) {
  sets <- msigdbr(species = species, collection = collection) %>%
    dplyr::select(gs_name, !!rlang::sym(gene_col))

  if (filter_empty) {
    sets <- sets %>% filter(.data[[gene_col]] != "")
  }
  if (!is.null(pathways)) {
    sets <- sets %>% filter(gs_name %in% pathways)
  }

  split(sets[[gene_col]], sets$gs_name)
}

#' Run GSEA on Ranked Test Statistics
#'
#' Uses a specified statistic column and identifier column to run pathway
#' enrichment with `fgsea`.
#'
#' @param tt A result table containing at least `stat_col` and `id_col`.
#' @param pathways A named list of pathways.
#' @param label Optional model label to attach to the enrichment results.
#' @param id_col Name of the identifier column.
#' @param stat_col Name of the statistic column.
#' @param seed Optional seed for reproducible fgsea execution.
#'
#' @return A data frame returned by `fgsea`, with an added `label` column when
#'   requested.
run_gsea <- function(tt, pathways, label = NULL, id_col = "gene_id",
                     stat_col = "t", seed = NULL) {
  ranks <- setNames(tt[[stat_col]], tt[[id_col]])
  ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()))
    set.seed(seed)
  }
  res <- fgsea(pathways = pathways, stats = ranks,
               minSize = 15, maxSize = 500, nPermSimple = 10000)
  if (!is.null(label)) {
    res$label <- label
  }
  res
}

#' Compute Pathway NES Correlation Against a Reference GSEA
#'
#' Runs GSEA on a topTable and returns the Spearman correlation of NES values
#' versus a reference GSEA result. Convenience wrapper around `run_gsea` and
#' `compare_gsea_results`.
#'
#' @param tt A topTable data frame (must contain `id_col` and `stat_col`).
#' @param ref_gsea Reference GSEA result (from `fgsea`).
#' @param pathways A named list of pathway gene sets.
#' @param id_col Name of the gene identifier column in `tt`.
#' @param stat_col Name of the statistic column in `tt`.
#' @param seed Optional seed for reproducible fgsea execution.
#'
#' @return A single numeric value: Spearman rho of NES (query vs reference).
pathway_cor <- function(tt, ref_gsea, pathways,
                        id_col = "gene_id", stat_col = "t",
                        seed = 42) {
  gsea_new <- run_gsea(tt, pathways, id_col = id_col,
                       stat_col = stat_col, seed = seed)
  merged <- dplyr::inner_join(
    ref_gsea %>% dplyr::select(pathway, NES_ref = NES),
    gsea_new %>% dplyr::select(pathway, NES_new = NES),
    by = "pathway"
  )
  cor(merged$NES_ref, merged$NES_new, method = "spearman")
}

#' Compare Two Differential Expression Result Tables
#'
#' Aligns two top-table outputs and computes gene-level concordance metrics.
#'
#' @param tt_ref Reference differential expression table.
#' @param tt_new Comparison differential expression table.
#' @param method Correlation method passed to `cor()`.
#'
#' @return A list containing the merged table and concordance summary metrics.
compare_de_results <- function(tt_ref, tt_new, method = "spearman") {
  merged <- inner_join(
    tt_ref %>% dplyr::select(gene_id, logFC_ref = logFC, t_ref = t, padj_ref = adj.P.Val),
    tt_new %>% dplyr::select(gene_id, logFC_new = logFC, t_new = t, padj_new = adj.P.Val),
    by = "gene_id"
  )

  sig_ref <- merged %>% filter(padj_ref < 0.05) %>% pull(gene_id)
  sig_new <- merged %>% filter(padj_new < 0.05) %>% pull(gene_id)
  overlap <- intersect(sig_ref, sig_new)

  list(
    merged = merged,
    rho_t = cor(merged$t_ref, merged$t_new, method = method),
    rho_lfc = cor(merged$logFC_ref, merged$logFC_new, method = method),
    n_deg_ref = length(sig_ref),
    n_deg_new = length(sig_new),
    overlap = length(overlap),
    overlap_genes = overlap,
    direction_concordance = if (length(overlap) > 0) {
      shared <- merged %>% filter(gene_id %in% overlap)
      mean(sign(shared$logFC_ref) == sign(shared$logFC_new))
    } else {
      NA_real_
    },
    jaccard = length(overlap) / max(length(union(sig_ref, sig_new)), 1)
  )
}

#' Add Ischemic Time Bins
#'
#' Adds an ordered `isch_bin` factor using the cutpoints shared across the
#' robustness analyses.
#'
#' @param tb A sample metadata table containing `ischemic_hrs`.
#'
#' @return `tb` with an added ordered `isch_bin` column.
add_ischemic_bin <- function(tb) {
  tb %>%
    mutate(isch_bin = factor(case_when(
      ischemic_hrs <= 3  ~ "Short",
      ischemic_hrs <= 10 ~ "Medium",
      TRUE               ~ "Long"
    ), levels = c("Short", "Medium", "Long")))
}

#' Clean Hallmark Pathway Names for Display
#'
#' Converts machine-readable Hallmark pathway identifiers into display labels.
#'
#' @param x A character vector of pathway identifiers.
#'
#' @return A character vector of cleaned pathway labels.
clean_pathway <- function(x) {
  x %>%
    str_remove("^HALLMARK_") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_replace("Tnfa", "TNFa") %>%
    str_replace("Nfkb", "NF-kB") %>%
    str_replace("Ifn", "IFN") %>%
    str_replace("Uv ", "UV ") %>%
    str_replace("Dna ", "DNA ") %>%
    str_replace("Kras ", "KRAS ") %>%
    str_replace("Il2 ", "IL2 ") %>%
    str_replace("Il6 ", "IL6 ") %>%
    str_replace("Myc ", "MYC ") %>%
    str_replace("E2f ", "E2F ") %>%
    str_replace("P53 ", "p53 ") %>%
    str_replace("G2m ", "G2M ") %>%
    str_replace("Mtorc1", "mTORC1") %>%
    str_replace("Tgf ", "TGF ") %>%
    str_replace("Wnt ", "WNT ") %>%
    str_replace("Pi3k ", "PI3K ")
}
