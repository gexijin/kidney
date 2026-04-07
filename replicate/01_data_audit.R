#!/usr/bin/env Rscript

source("kidney/manuscript/replicate/00_common.R")

cat("01_data_audit.R\n")
cat(strrep("=", 80), "\n")

inputs <- read_inputs()
sat <- prepare_tissue(inputs$pheno, inputs$counts_all, "Adipose - Subcutaneous")
vat <- prepare_tissue(inputs$pheno, inputs$counts_all, "Adipose - Visceral (Omentum)")

cat("SAT samples:", nrow(sat$tb), "genes:", nrow(sat$counts), "\n")
cat("VAT samples:", nrow(vat$tb), "genes:", nrow(vat$counts), "\n")

sat_counts <- sex_ckd_counts(sat$tb)
vat_counts <- sex_ckd_counts(vat$tb)

write_csv(sat_counts, file.path(replicate_dirs$tables, "01_sat_sex_ckd_counts.csv"))
write_csv(vat_counts, file.path(replicate_dirs$tables, "01_vat_sex_ckd_counts.csv"))

female_sat <- sat$tb %>% filter(SEX == "Female") %>% mutate(CKD = as.integer(MHRNLFLR == 1))
male_sat <- sat$tb %>% filter(SEX == "Male") %>% mutate(CKD = as.integer(MHRNLFLR == 1))

stopifnot(sum(female_sat$CKD) == 27, sum(male_sat$CKD) == 56)
stopifnot(nrow(female_sat) == 218, nrow(male_sat) == 444)

female_cont_support <- female_sat %>% filter(CKD == 0)
female_case_support <- female_sat %>% filter(CKD == 1)

range_overlap <- tibble(
  variable = c("AGE", "BMI", "ischemic_hrs", "SMRIN", "SMEXNCRT"),
  case_min = c(
    min(female_case_support$AGE), min(female_case_support$BMI),
    min(female_case_support$ischemic_hrs), min(female_case_support$SMRIN),
    min(female_case_support$SMEXNCRT)
  ),
  case_max = c(
    max(female_case_support$AGE), max(female_case_support$BMI),
    max(female_case_support$ischemic_hrs), max(female_case_support$SMRIN),
    max(female_case_support$SMEXNCRT)
  ),
  control_min = c(
    min(female_cont_support$AGE), min(female_cont_support$BMI),
    min(female_cont_support$ischemic_hrs), min(female_cont_support$SMRIN),
    min(female_cont_support$SMEXNCRT)
  ),
  control_max = c(
    max(female_cont_support$AGE), max(female_cont_support$BMI),
    max(female_cont_support$ischemic_hrs), max(female_cont_support$SMRIN),
    max(female_cont_support$SMEXNCRT)
  )
) %>%
  mutate(
    overlap_min = pmax(case_min, control_min),
    overlap_max = pmin(case_max, control_max),
    overlap_exists = overlap_min <= overlap_max
  )

write_csv(range_overlap, file.path(replicate_dirs$tables, "01_female_overlap_ranges.csv"))

female_balance <- bind_rows(
  tibble(
    variable = c("AGE", "BMI", "ischemic_hrs", "SMRIN", "SMEXNCRT"),
    smd = c(
      numeric_smd(female_sat$AGE, female_sat$CKD),
      numeric_smd(female_sat$BMI, female_sat$CKD),
      numeric_smd(female_sat$ischemic_hrs, female_sat$CKD),
      numeric_smd(female_sat$SMRIN, female_sat$CKD),
      numeric_smd(female_sat$SMEXNCRT, female_sat$CKD)
    )
  ),
  factor_level_smds(female_sat$DTHHRDY, female_sat$CKD, "DTHHRDY"),
  factor_level_smds(female_sat$RACE, female_sat$CKD, "RACE"),
  tibble(
    variable = c("MHHTN", "MHT2D", "MHLVRDIS", "MHABNWBC", "MHHRTDIS", "MHHRTATT", "MHCOPD", "MHORGNTP", "MHDLYSIS"),
    smd = c(
      binary_smd(female_sat$MHHTN, female_sat$CKD),
      binary_smd(female_sat$MHT2D, female_sat$CKD),
      binary_smd(female_sat$MHLVRDIS, female_sat$CKD),
      binary_smd(female_sat$MHABNWBC, female_sat$CKD),
      binary_smd(female_sat$MHHRTDIS, female_sat$CKD),
      binary_smd(female_sat$MHHRTATT, female_sat$CKD),
      binary_smd(female_sat$MHCOPD, female_sat$CKD),
      binary_smd(replace_na(female_sat$MHORGNTP, 0), female_sat$CKD),
      binary_smd(replace_na(female_sat$MHDLYSIS, 0), female_sat$CKD)
    )
  )
) %>%
  mutate(abs_smd = abs(smd)) %>%
  arrange(desc(abs_smd))

write_csv(female_balance, file.path(replicate_dirs$tables, "01_female_balance_smd.csv"))

female_demographics <- tibble(
  group = c("Female CKD", "Female control"),
  n = c(sum(female_sat$CKD == 1), sum(female_sat$CKD == 0)),
  age = c(fmt_median_iqr(female_case_support$AGE), fmt_median_iqr(female_cont_support$AGE)),
  bmi = c(fmt_median_iqr(female_case_support$BMI), fmt_median_iqr(female_cont_support$BMI)),
  ischemic_hrs = c(fmt_median_iqr(female_case_support$ischemic_hrs), fmt_median_iqr(female_cont_support$ischemic_hrs)),
  rin = c(fmt_median_iqr(female_case_support$SMRIN), fmt_median_iqr(female_cont_support$SMRIN)),
  exonic_rate = c(fmt_median_iqr(female_case_support$SMEXNCRT), fmt_median_iqr(female_cont_support$SMEXNCRT))
)

write_csv(female_demographics, file.path(replicate_dirs$tables, "01_female_demographics.csv"))

hardy_table <- female_sat %>%
  count(CKD, DTHHRDY, name = "n") %>%
  mutate(group = if_else(CKD == 1, "Female CKD", "Female control")) %>%
  select(group, DTHHRDY, n)

write_csv(hardy_table, file.path(replicate_dirs$tables, "01_female_hardy_counts.csv"))

female_long_overlap <- female_sat %>%
  mutate(
    within_case_range = between(ischemic_hrs, min(female_case_support$ischemic_hrs), max(female_case_support$ischemic_hrs)),
    within_control_range = between(ischemic_hrs, min(female_cont_support$ischemic_hrs), max(female_cont_support$ischemic_hrs))
  )

overlap_summary <- tibble(
  metric = c(
    "Female CKD cases",
    "Female controls",
    "Female CKD with ischemic_hrs <= control max",
    "Female controls with ischemic_hrs >= case min",
    "Female CKD slow deaths",
    "Female controls slow deaths"
  ),
  value = c(
    nrow(female_case_support),
    nrow(female_cont_support),
    sum(female_case_support$ischemic_hrs <= max(female_cont_support$ischemic_hrs)),
    sum(female_cont_support$ischemic_hrs >= min(female_case_support$ischemic_hrs)),
    sum(female_case_support$DTHHRDY == "Slow"),
    sum(female_cont_support$DTHHRDY == "Slow")
  )
)

write_csv(overlap_summary, file.path(replicate_dirs$tables, "01_female_overlap_summary.csv"))

p_ischemia <- ggplot(female_sat, aes(x = ischemic_hrs, fill = factor(CKD))) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 25) +
  scale_fill_manual(values = c("0" = "#4C78A8", "1" = "#E15759"), labels = c("Control", "CKD"), name = NULL) +
  labs(
    title = "Female SAT ischemic time overlap",
    x = "Ischemic time (hours)",
    y = "Donor count"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(replicate_dirs$figures, "01_female_ischemic_overlap.png"),
  p_ischemia, width = 7, height = 4.5, dpi = 300
)

p_hardy <- female_sat %>%
  count(CKD, DTHHRDY, name = "n") %>%
  mutate(group = if_else(CKD == 1, "CKD", "Control")) %>%
  ggplot(aes(x = group, y = n, fill = DTHHRDY)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "Set2", name = "Hardy") +
  labs(
    title = "Female SAT Hardy death-mode composition",
    x = NULL,
    y = "Composition"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(replicate_dirs$figures, "01_female_hardy_composition.png"),
  p_hardy, width = 5.5, height = 4.5, dpi = 300
)

p_smd <- female_balance %>%
  mutate(variable = forcats::fct_reorder(variable, abs_smd)) %>%
  ggplot(aes(x = smd, y = variable)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dotted", color = "grey60") +
  geom_point(color = "#E15759", size = 2) +
  labs(
    title = "Female SAT CKD-control imbalance",
    x = "Standardized mean difference",
    y = NULL
  ) +
  theme_bw(base_size = 10)

ggsave(
  file.path(replicate_dirs$figures, "01_female_balance_smd.png"),
  p_smd, width = 7.5, height = 6.5, dpi = 300
)

saveRDS(
  list(
    sat = sat,
    vat = vat,
    sat_counts = sat_counts,
    vat_counts = vat_counts,
    female_balance = female_balance,
    range_overlap = range_overlap,
    overlap_summary = overlap_summary
  ),
  file.path(replicate_dirs$rds, "01_data_audit.rds")
)

write_session_info(file.path(replicate_dirs$logs, "01_data_audit_sessionInfo.txt"))

cat("\nChecks completed.\n")
cat("SAT sex x CKD counts:\n")
print(sat_counts)
cat("\nVAT sex x CKD counts:\n")
print(vat_counts)
cat("\nWorst female SAT imbalances by absolute SMD:\n")
print(head(female_balance, 10))
cat("\nFemale SAT overlap summary:\n")
print(overlap_summary)
