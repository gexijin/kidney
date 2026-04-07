# Replication Findings

## Scope

This folder contains a fresh replication and challenge analysis of the manuscript claim that female subcutaneous adipose tissue in CKD shows a large ischemia-like transcriptomic signature.

I did **not** modify any code in `kidney/R/`. All new work is in:

- [00_common.R](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/00_common.R)
- [01_data_audit.R](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/01_data_audit.R)
- [02_primary_replication.R](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/02_primary_replication.R)
- [03_confounding_diagnostics.R](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/03_confounding_diagnostics.R)
- [04_matching_sensitivity.R](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/04_matching_sensitivity.R)
- [05_matching_variations.R](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/05_matching_variations.R)
- [06_exact_matching_variants.R](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/06_exact_matching_variants.R)
- [07_hallmark_nes_sensitivity.R](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/07_hallmark_nes_sensitivity.R)

## What reproduced exactly

Using only:

- `data/phenotype/merged/merged_phenotype.rds`
- `data/reads/gene_reads_v11.rds`

and reconstructing the cohort exactly as described in the manuscript:

- tissue subset by `SMTSD`
- sample intersection by `SAMPID`
- gene filter `>=10` counts in `>=20%` of samples
- complete cases on the manuscript model variables

I reproduced the manuscript's primary female SAT results exactly:

- SAT cohort: `218` females (`27` CKD) and `444` males (`56` CKD)
- Female SAT: `1,960` DEGs
- Male SAT: `5` DEGs
- Female VAT: `1` DEG
- SAT sex-by-CKD interaction: `917` DEGs
- Female vs male SAT logFC correlation: Pearson `r = 0.198`

Key tables and plots:

- [02_deg_summary.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/02_deg_summary.csv)
- [02_sex_correlation_summary.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/02_sex_correlation_summary.csv)
- [02_primary_volcano.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/02_primary_volcano.png)
- [02_sex_effect_scatter.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/02_sex_effect_scatter.png)
- [02_depot_effect_scatter.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/02_depot_effect_scatter.png)

## Why confounding is a serious issue

Before matching or restriction, the female CKD and female control groups are extremely imbalanced on agonal-state variables.

Largest female SAT standardized mean differences:

- `DTHHRDY=Slow`: `0.949`
- `ischemic_hrs`: `0.880`
- `AGE`: `0.661`
- `SMRIN`: `-0.646`
- `SMEXNCRT`: `-0.608`

That is documented in:

- [01_female_balance_smd.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/01_female_balance_smd.csv)
- [01_female_balance_smd.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/01_female_balance_smd.png)
- [01_female_ischemic_overlap.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/01_female_ischemic_overlap.png)
- [01_female_hardy_composition.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/01_female_hardy_composition.png)

The raw female SAT summary matches the manuscript:

- female CKD ischemic time median `12.8 h`
- female control ischemic time median `3.8 h`
- female CKD slow deaths `16/27`
- female controls slow deaths `26/191`

See:

- [01_female_demographics.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/01_female_demographics.csv)
- [01_female_hardy_counts.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/01_female_hardy_counts.csv)

## Gene-level attenuation under stricter control

I tested several stricter female SAT analyses in [03_confounding_diagnostics.R](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/03_confounding_diagnostics.R).

Summary:

- Primary female SAT model: `1,960` DEGs
- Spline ischemia adjustment: `884` DEGs
- Ventilator-only females: `687` DEGs
- Empirical ischemia-overlap subset: `17` DEGs
- Exact `Hardy + ischemia-bin` matching, `1:1`, with continuous ischemic time still in the model:
  - `21` CKD + `21` controls: `0` DEGs
  - alternate exact variant `17` CKD + `17` controls: `0` DEGs

Key files:

- [03_model_restriction_summary.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/03_model_restriction_summary.csv)
- [03_exact_matching_summary.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/03_exact_matching_summary.csv)
- [03_propensity_summary.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/03_propensity_summary.csv)
- [03_propensity_overlap.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/03_propensity_overlap.png)
- [03_deg_attenuation.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/03_deg_attenuation.png)
- [03_exact_matching_balance.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/03_exact_matching_balance.png)

At the gene level, this is strong evidence that the exact number of DEGs is highly sensitive to overlap restrictions and agonal-state control.

## Matching-ratio sensitivity

The first exact `1:1` match was too easy to dismiss as a power problem, so I expanded this in [05_matching_variations.R](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/05_matching_variations.R).

I compared:

- `paper_style`
  - exact on ischemia bin only
  - nearest-neighbor on age, race, Hardy, ischemic time, RIN

- `exact_hardy`
  - exact on Hardy death mode and ischemia bin
  - nearest-neighbor on age, BMI, ischemic time, RIN, exonic rate

Repeated over `10` seeds for ratios `1:1`, `1:2`, `1:3`.

Median DEG counts:

- `exact_hardy 1:1`: `0`
- `exact_hardy 1:2`: `686`
- `exact_hardy 1:3`: `264`
- `paper_style 1:1`: `182.5`
- `paper_style 1:2`: `585.5`
- `paper_style 1:3`: `742`

Ranges:

- `exact_hardy 1:1`: `0` to `1`
- `exact_hardy 1:2`: `352` to `779`
- `exact_hardy 1:3`: `0` to `291`
- `paper_style 1:1`: `0` to `381`
- `paper_style 1:2`: `18` to `886`
- `paper_style 1:3`: `67` to `1080`

Key point: the stricter exact match does **not** produce a monotone “more controls = more stable biology” pattern. The signal is unstable across ratio and control allocation inside tight strata.

See:

- [05_matching_variations_summary.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/05_matching_variations_summary.csv)
- [05_matching_variations_iterations.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/05_matching_variations_iterations.csv)
- [05_matching_variations_deg.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/05_matching_variations_deg.png)
- [05_matching_variations_balance.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/05_matching_variations_balance.png)

## Targeted exact-matching variants

I then tested whether the rebound in exact `1:2` matching was tied to residual imbalance in other variables, especially BMI, in [06_exact_matching_variants.R](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/06_exact_matching_variants.R).

Seed-1 targeted variants:

- Exact `Hardy+isch 1:1`
  - `21` CKD / `21` controls
  - BMI SMD `0.683`
  - DEGs `0`

- Exact `Hardy+isch 1:2`
  - `21` / `34`
  - BMI SMD `0.461`
  - DEGs `470`

- Exact `Hardy+isch 1:2 Mahalanobis`
  - `27` / `47`
  - BMI SMD `0.298`
  - DEGs `1422`

- Exact `Hardy+isch+bmi_bin 1:2`
  - `15` / `20`
  - BMI SMD `-0.066`
  - DEGs `0`

- Paper-style `1:2`
  - `22` / `36`
  - BMI SMD `0.336`
  - DEGs `1003`

This does not prove BMI is the whole story, but it shows the recovered signal is sensitive to which residual imbalances are allowed even after exact `Hardy + ischemia` control.

See:

- [06_exact_matching_variants_summary.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/06_exact_matching_variants_summary.csv)
- [06_exact_matching_variants_deg.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/06_exact_matching_variants_deg.png)
- [06_exact_matching_variants_balance.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/06_exact_matching_variants_balance.png)

## Hallmark GSEA changes the interpretation

The strongest update is pathway-level.

Because DEG counts are very cutoff-sensitive, I ran Hallmark GSEA on the primary female SAT model and on all matched variants in [07_hallmark_nes_sensitivity.R](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/07_hallmark_nes_sensitivity.R).

For each matched analysis, I ranked genes by limma moderated `t` and computed Hallmark NES. I then correlated matched-model NES with primary-model NES.

Median Hallmark NES Spearman correlations:

- `exact_hardy 1:1`: `0.867`
- `exact_hardy 1:2`: `0.905`
- `exact_hardy 1:3`: `0.918`
- `paper_style 1:1`: `0.888`
- `paper_style 1:2`: `0.946`
- `paper_style 1:3`: `0.963`

Ranges across 10 seeds:

- `exact_hardy 1:1`: `0.858` to `0.874`
- `exact_hardy 1:2`: `0.880` to `0.932`
- `exact_hardy 1:3`: `0.905` to `0.940`
- `paper_style 1:1`: `0.861` to `0.917`
- `paper_style 1:2`: `0.920` to `0.954`
- `paper_style 1:3`: `0.952` to `0.979`

This means that even when gene-level DEGs collapse under strict matching, the overall Hallmark pathway architecture remains highly similar to the primary female SAT result.

Key files:

- [07_hallmark_nes_summary.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/07_hallmark_nes_summary.csv)
- [07_hallmark_nes_iterations.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/07_hallmark_nes_iterations.csv)
- [07_hallmark_nes_correlation.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/07_hallmark_nes_correlation.png)

## Core Hallmark pathways are preserved

From the seed-1 comparison in [07_seed1_hallmark_nes.csv](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/tables/07_seed1_hallmark_nes.csv):

- Adipogenesis
  - primary NES `-3.18`
  - exact `1:1` `-2.21`
  - exact `1:2` `-2.84`
  - exact `1:3` `-3.02`

- Oxidative phosphorylation
  - primary `-2.55`
  - exact `1:1` `-1.84`
  - exact `1:2` `-1.96`
  - exact `1:3` `-2.34`

- Fatty acid metabolism
  - primary `-2.51`
  - exact `1:1` `-1.98`
  - exact `1:2` `-2.34`
  - exact `1:3` `-2.34`

- Hypoxia
  - primary `2.43`
  - exact `1:1` `2.69`
  - exact `1:2` `2.42`
  - exact `1:3` `2.31`

- Inflammatory response
  - primary `2.68`
  - exact `1:1` `2.75`
  - exact `1:2` `2.46`
  - exact `1:3` `2.48`

- EMT
  - primary `2.84`
  - exact `1:1` `2.45`
  - exact `1:2` `2.71`
  - exact `1:3` `2.68`

- TNFa via NFkB
  - primary `3.17`
  - exact `1:1` `3.19`
  - exact `1:2` `2.92`
  - exact `1:3` `2.86`

These pathways remained strongly significant even in the exact `1:1` run.

Useful plot:

- [07_seed1_hallmark_nes.png](/Users/ge/Documents/research/gtex/kidney/manuscript/replicate/figures/07_seed1_hallmark_nes.png)

## Current interpretation

The most defensible interpretation after the full follow-up is:

- The manuscript's headline female SAT result is reproducible.
- The exact number of gene-level DEGs is highly sensitive to overlap, matching ratio, and residual imbalance.
- There is strong evidence of substantial confounding by post-mortem state.
- But the stricter matching analyses do **not** erase the pathway-level biology.
- The core metabolic-suppression / inflammatory-activation Hallmark program is preserved across strict matching variants.

So the challenge to the manuscript is now narrower and more precise:

- The manuscript likely **overstates gene-level robustness**.
- The manuscript's matching-based robustness claims should be presented with much more nuance.
- A strong statement like “the whole signal is just confounding” is **not supported** by the pathway-level analysis.
- A stronger statement like “the exact DEG count is not robust, but the pathway architecture is substantially preserved” **is** supported.

## What I have not found

- I have **not** found a fatal coding bug in `kidney/R/`.
- I have **not** shown that the entire female SAT CKD signature is an artifact.
- I have **not** invalidated the manuscript's core Hallmark pathway pattern.

## Best current conclusion

The current evidence supports this conclusion:

The female SAT CKD signal is **partly confounded and gene-level fragile**, but the underlying Hallmark pathway program is **surprisingly stable** under stricter matching. The manuscript should probably soften claims about gene-level robustness and stronger causal interpretation, but the pathway-level claim remains substantially credible.
