# =============================================================================
# LIPID MELANOMA IMMUNOTHERAPY RESEARCH  
# Script 08: Lipid-Immune Interaction Survival Analysis
# Purpose: Test if Lipid-High + CD8-Low combination predicts poor survival
# Author: madhavellinilab
# Date: Week 2 - Final Analysis
# =============================================================================

library(SummarizedExperiment)
library(tidyverse)
library(survival)
library(survminer)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("manuscript/figures", showWarnings = FALSE, recursive = TRUE)

cat("\n========================================================")
cat("\nLIPID-IMMUNE INTERACTION SURVIVAL ANALYSIS")
cat("\n========================================================\n")

# =============================================================================
# STEP 1: LOAD INTEGRATED DATA
# =============================================================================

cat("\n--- Step 1: Loading integrated dataset ---\n")

integrated_data <- read.csv("results/integrated_lipid_immune.csv")
analysis_data <- read.csv("results/analysis_dataset.csv")

# Merge with clinical data
final_data <- analysis_data %>%
  left_join(integrated_data, by = c("barcode" = "sample_id"))

# Filter complete cases
final_data <- final_data %>%
  filter(!is.na(combined_group), !is.na(OS), !is.na(OS.time), OS.time > 0)

cat(sprintf("Final dataset: %d samples\n", nrow(final_data)))
cat(sprintf("Events: %d deaths\n", sum(final_data$OS)))

# =============================================================================
# STEP 2: 4-GROUP SURVIVAL ANALYSIS
# =============================================================================

cat("\n--- Step 2: Four-group survival analysis ---\n")

# Create survival object
surv_obj <- Surv(time = final_data$OS.time / 365.25, event = final_data$OS)

# Fit Kaplan-Meier by 4 groups
fit_4groups <- survfit(surv_obj ~ combined_group, data = final_data)

# Log-rank test
surv_diff <- survdiff(surv_obj ~ combined_group, data = final_data)
pval_4group <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

cat(sprintf("\nLog-rank test (4 groups): p = %.4f\n", pval_4group))

# Median survival by group
median_surv <- summary(fit_4groups)$table
cat("\nMedian survival (years) by group:\n")
for (i in 1:nrow(median_surv)) {
  cat(sprintf("  %s: %.2f [95%% CI: %.2f-%.2f]\n",
              rownames(median_surv)[i],
              median_surv[i, "median"] / 365.25,
              median_surv[i, "0.95LCL"] / 365.25,
              median_surv[i, "0.95UCL"] / 365.25))
}

# =============================================================================
# STEP 3: PAIRWISE COMPARISONS
# =============================================================================

cat("\n--- Step 3: Pairwise comparisons ---\n")

# Best vs Worst: Lipid-Low_CD8-High vs Lipid-High_CD8-Low
best_worst <- final_data %>%
  filter(combined_group %in% c("Lipid-Low_CD8-High", "Lipid-High_CD8-Low"))

if (nrow(best_worst) > 0) {
  surv_bw <- Surv(time = best_worst$OS.time / 365.25, event = best_worst$OS)
  fit_bw <- survfit(surv_bw ~ combined_group, data = best_worst)
  diff_bw <- survdiff(surv_bw ~ combined_group, data = best_worst)
  pval_bw <- 1 - pchisq(diff_bw$chisq, 1)
  
  cat(sprintf("\nBest vs Worst group comparison: p = %.4f\n", pval_bw))
}

# =============================================================================
# STEP 4: COX REGRESSION WITH INTERACTION TERM
# =============================================================================

cat("\n--- Step 4: Cox regression with interaction ---\n")

# Test interaction effect
cox_interaction <- coxph(surv_obj ~ lipid_score * cd8_score, data = final_data)

cat("\nCox model with interaction:\n")
print(summary(cox_interaction))

# Extract interaction term p-value
interaction_p <- summary(cox_interaction)$coefficients["lipid_score:cd8_score", "Pr(>|z|)"]
cat(sprintf("\nInteraction term p-value: %.4f\n", interaction_p))

# =============================================================================
# STEP 5: GENERATE VISUALIZATIONS
# =============================================================================

cat("\n--- Step 5: Generating plots ---\n")

# Four-group Kaplan-Meier plot
pdf("figures/04_survival_4groups.pdf", width = 10, height = 7)
gg <- ggsurvplot(
  fit_4groups,
  data = final_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  legend.title = "Lipid-Immune Group",
  legend.labs = levels(final_data$combined_group),
  palette = c("#2E7D32", "#1976D2", "#F57C00", "#C62828"),
  xlab = "Time (years)",
  ylab = "Overall Survival",
  title = "Survival by Lipid-Immune Stratification\nTCGA-SKCM Melanoma (n=335)",
  ggtheme = theme_bw(base_size = 12)
)
print(gg)
dev.off()

cat("Saved: figures/04_survival_4groups.pdf\n")

# Copy to manuscript folder
file.copy("figures/04_survival_4groups.pdf",
          "manuscript/figures/Figure4_interaction_survival.pdf",
          overwrite = TRUE)

# Scatter plot: Lipid vs CD8 scores
pdf("figures/05_lipid_cd8_correlation.pdf", width = 7, height = 6)
ggplot(final_data, aes(x = lipid_score, y = cd8_score, color = combined_group)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#2E7D32", "#1976D2", "#F57C00", "#C62828")) +
  labs(
    title = "Lipid Metabolism vs CD8+ T-Cell Infiltration",
    subtitle = sprintf("r = -0.244, p < 0.0001"),
    x = "Lipid Metabolism Score",
    y = "CD8+ T-Cell Score",
    color = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")
dev.off()

cat("Saved: figures/05_lipid_cd8_correlation.pdf\n")

# =============================================================================
# STEP 6: SAVE RESULTS
# =============================================================================

cat("\n--- Step 6: Saving final results ---\n")

# Create results summary
interaction_results <- data.frame(
  Analysis = c("Four-group log-rank", "Best vs Worst", "Interaction term"),
  P_value = c(pval_4group, pval_bw, interaction_p),
  N_samples = nrow(final_data),
  N_events = sum(final_data$OS)
)

write.csv(interaction_results, "results/interaction_results.csv", row.names = FALSE)
write.csv(final_data, "results/final_analysis_dataset.csv", row.names = FALSE)

cat("\n========================================================\n")
cat("ANALYSIS COMPLETE!")
cat("\n========================================================\n")
cat("\nKey Findings:\n")
if (pval_4group < 0.05) {
  cat(sprintf("  ✓ SIGNIFICANT 4-group difference: p=%.4f\n", pval_4group))
} else {
  cat(sprintf("  ✗ No significant 4-group difference: p=%.4f\n", pval_4group))
}
if (pval_bw < 0.05) {
  cat(sprintf("  ✓ SIGNIFICANT best-vs-worst: p=%.4f\n", pval_bw))
} else {
  cat(sprintf("  ✗ No significant best-vs-worst: p=%.4f\n", pval_bw))
}
if (interaction_p < 0.05) {
  cat(sprintf("  ✓ SIGNIFICANT interaction: p=%.4f\n", interaction_p))
} else {
  cat(sprintf("  ✗ No significant interaction: p=%.4f\n", interaction_p))
}
cat("\nFiles generated:\n")
cat("  - figures/04_survival_4groups.pdf\n")
cat("  - figures/05_lipid_cd8_correlation.pdf\n")
cat("  - results/interaction_results.csv\n")
cat("  - results/final_analysis_dataset.csv\n")
cat("\n")
