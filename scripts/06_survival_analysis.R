# =============================================================================
# LIPID MELANOMA IMMUNOTHERAPY RESEARCH
# Script 06: Survival Analysis and Cox Regression
# Purpose: Analyze association between lipid signature and clinical outcomes
# Author: madhavellinilab
# Date: Week 2 - Phase 3 Integration
# =============================================================================

library(SummarizedExperiment)
library(tidyverse)
library(survival)
library(survminer)

# Set working directory
setwd("/path/to/project")

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

cat("\n================================================")
cat("\nSURVIVAL ANALYSIS: Lipid Signature & Outcomes")
cat("\n================================================\n")

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

cat("\n--- Step 1: Loading data ---\n")

# Load QC-filtered TCGA data
if (!file.exists("data/tcga_skcm_qc.rds")) {
  stop("Error: QC-filtered data not found.")
}

tcga_skcm_qc <- readRDS("data/tcga_skcm_qc.rds")

# Load lipid signature scores
if (!file.exists("results/lipid_score.rds")) {
  stop("Error: Lipid scores not found. Please run 05_lipid_signature.R first.")
}

lipid_scores <- readRDS("results/lipid_score.rds")

cat(sprintf("Loaded TCGA data: %d samples\n", ncol(tcga_skcm_qc)))
cat(sprintf("Loaded lipid scores: %d samples\n", nrow(lipid_scores)))

# =============================================================================
# STEP 2: PREPARE ANALYSIS DATASET
# =============================================================================

cat("\n--- Step 2: Preparing analysis dataset ---\n")

# Extract clinical data
clinical <- as.data.frame(colData(tcga_skcm_qc))

# Merge with lipid scores
analysis_df <- clinical %>%
  left_join(lipid_scores, by = c("barcode" = "sample_id"))

# Check for survival data availability
cat("\nChecking available clinical variables:\n")
cat(sprintf("  - Overall survival time: %s\n", 
            ifelse("OS.time" %in% colnames(analysis_df), "Available", "Missing")))
cat(sprintf("  - Overall survival event: %s\n", 
            ifelse("OS" %in% colnames(analysis_df), "Available", "Missing")))
cat(sprintf("  - Disease-free survival: %s\n", 
            ifelse("DSS.time" %in% colnames(analysis_df), "Available", "Missing")))

# Create survival variables if available in standard TCGA format
if (!"OS.time" %in% colnames(analysis_df)) {
  # Try alternative column names
  if ("days_to_death" %in% colnames(analysis_df) || 
      "days_to_last_followup" %in% colnames(analysis_df)) {
    
    analysis_df <- analysis_df %>%
      mutate(
        OS.time = coalesce(days_to_death, days_to_last_follow_up),
        OS = ifelse(!is.na(days_to_death), 1, 0)
      )
    cat("\nCreated OS variables from TCGA clinical data\n")
  }
}

# Filter to samples with complete data
analysis_complete <- analysis_df %>%
  filter(
    !is.na(lipid_score),
    !is.na(OS.time),
    !is.na(OS),
    OS.time > 0
  )

cat(sprintf("\nSamples with complete data: %d\n", nrow(analysis_complete)))
cat(sprintf("Events (deaths): %d\n", sum(analysis_complete$OS)))
cat(sprintf("Censored: %d\n", sum(1 - analysis_complete$OS)))

# Save analysis dataset
write_csv(analysis_complete, "results/analysis_dataset.csv")
saveRDS(analysis_complete, "results/analysis_dataset.rds")
cat("\nSaved: results/analysis_dataset.{csv,rds}\n")

# =============================================================================
# STEP 3: KAPLAN-MEIER SURVIVAL ANALYSIS
# =============================================================================

cat("\n--- Step 3: Kaplan-Meier survival analysis ---\n")

# Create survival object
surv_obj <- Surv(time = analysis_complete$OS.time / 365.25, 
                 event = analysis_complete$OS)

# Fit Kaplan-Meier curves by lipid group
fit_km <- survfit(surv_obj ~ lipid_group, data = analysis_complete)

# Log-rank test
surv_diff <- survdiff(surv_obj ~ lipid_group, data = analysis_complete)
pval_logrank <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

cat(sprintf("\nLog-rank test p-value: %.4f\n", pval_logrank))

# Get median survival times
median_surv <- summary(fit_km)$table
cat("\nMedian survival times:\n")
print(median_surv[, c("median", "0.95LCL", "0.95UCL")])

# Generate Kaplan-Meier plot
cat("\nGenerating Kaplan-Meier plot...\n")

pdf("figures/03_survival_km_lipid.pdf", width = 8, height = 6)
ggsurvplot(
  fit_km,
  data = analysis_complete,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  legend.title = "Lipid Metabolism",
  legend.labs = c("Low", "High"),
  palette = c("#377EB8", "#E41A1C"),
  xlab = "Time (years)",
  ylab = "Overall Survival Probability",
  title = "Survival by Lipid Metabolism Signature\nTCGA-SKCM Melanoma",
  ggtheme = theme_bw(),
  font.main = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.legend = c(11)
)
dev.off()

cat("Saved: figures/03_survival_km_lipid.pdf\n")

# =============================================================================
# STEP 4: UNIVARIATE COX REGRESSION
# =============================================================================

cat("\n--- Step 4: Univariate Cox regression ---\n")

# Continuous lipid score
cox_univ_cont <- coxph(surv_obj ~ lipid_score, 
                       data = analysis_complete)

cat("\nCox model (continuous lipid score):\n")
print(summary(cox_univ_cont))

# Categorical lipid group
cox_univ_cat <- coxph(surv_obj ~ lipid_group, 
                      data = analysis_complete)

cat("\nCox model (categorical lipid group):\n")
print(summary(cox_univ_cat))

# Extract hazard ratios
hr_cont <- exp(coef(cox_univ_cont))
hr_cont_ci <- exp(confint(cox_univ_cont))
p_cont <- summary(cox_univ_cont)$coefficients["lipid_score", "Pr(>|z|)"]

hr_cat <- exp(coef(cox_univ_cat))
hr_cat_ci <- exp(confint(cox_univ_cat))
p_cat <- summary(cox_univ_cat)$coefficients["lipid_groupLipid-High", "Pr(>|z|)"]

cat(sprintf("\nHazard Ratio (continuous): %.3f (95%% CI: %.3f-%.3f), p=%.4f\n",
            hr_cont, hr_cont_ci[1], hr_cont_ci[2], p_cont))
cat(sprintf("Hazard Ratio (High vs Low): %.3f (95%% CI: %.3f-%.3f), p=%.4f\n",
            hr_cat, hr_cat_ci[1], hr_cat_ci[2], p_cat))

# =============================================================================
# STEP 5: MULTIVARIATE COX REGRESSION
# =============================================================================

cat("\n--- Step 5: Multivariate Cox regression ---\n")

# Check for available covariates
available_covars <- c()

if ("age_at_diagnosis" %in% colnames(analysis_complete)) {
  available_covars <- c(available_covars, "age_at_diagnosis")
}
if ("gender" %in% colnames(analysis_complete)) {
  available_covars <- c(available_covars, "gender")
}
if ("ajcc_pathologic_stage" %in% colnames(analysis_complete)) {
  available_covars <- c(available_covars, "ajcc_pathologic_stage")
}

if (length(available_covars) > 0) {
  cat("\nAvailable covariates for adjustment:\n")
  cat(paste("  -", available_covars, collapse = "\n"))
  
  # Build multivariate formula
  formula_str <- paste("surv_obj ~ lipid_score +", 
                       paste(available_covars, collapse = " + "))
  
  # Fit multivariate model
  cox_multiv <- coxph(as.formula(formula_str), 
                      data = analysis_complete)
  
  cat("\n\nMultivariate Cox model:\n")
  print(summary(cox_multiv))
  
  # Extract adjusted HR for lipid score
  hr_adj <- exp(coef(cox_multiv)["lipid_score"])
  hr_adj_ci <- exp(confint(cox_multiv)["lipid_score", ])
  p_adj <- summary(cox_multiv)$coefficients["lipid_score", "Pr(>|z|)"]
  
  cat(sprintf("\nAdjusted Hazard Ratio: %.3f (95%% CI: %.3f-%.3f), p=%.4f\n",
              hr_adj, hr_adj_ci[1], hr_adj_ci[2], p_adj))
} else {
  cat("\nInsufficient covariates available for multivariate analysis\n")
}

# =============================================================================
# STEP 6: SAVE RESULTS SUMMARY
# =============================================================================

cat("\n--- Step 6: Saving results summary ---\n")

# Create results summary
results_summary <- data.frame(
  analysis = c("Univariate (continuous)", "Univariate (categorical)"),
  HR = c(hr_cont, hr_cat),
  CI_lower = c(hr_cont_ci[1], hr_cat_ci[1]),
  CI_upper = c(hr_cont_ci[2], hr_cat_ci[2]),
  p_value = c(p_cont, p_cat),
  n_samples = rep(nrow(analysis_complete), 2),
  n_events = rep(sum(analysis_complete$OS), 2)
)

if (exists("hr_adj")) {
  results_summary <- rbind(
    results_summary,
    data.frame(
      analysis = "Multivariate (adjusted)",
      HR = hr_adj,
      CI_lower = hr_adj_ci[1],
      CI_upper = hr_adj_ci[2],
      p_value = p_adj,
      n_samples = nrow(analysis_complete),
      n_events = sum(analysis_complete$OS)
    )
  )
}

write_csv(results_summary, "results/survival_results.csv")
cat("Saved: results/survival_results.csv\n")

cat("\n================================================\n")
cat("SURVIVAL ANALYSIS COMPLETE")
cat("\n================================================\n")
cat("\nKey findings:\n")
if (p_cat < 0.05) {
  cat(sprintf("  ✓ Significant association: HR=%.2f, p=%.4f\n", hr_cat, p_cat))
  cat("  ✓ High lipid metabolism associated with", 
      ifelse(hr_cat > 1, "worse", "better"), "survival\n")
} else {
  cat(sprintf("  ✗ No significant association: p=%.4f\n", p_cat))
}

cat("\nFiles generated:\n")
cat("  - figures/03_survival_km_lipid.pdf\n")
cat("  - results/analysis_dataset.csv\n")
cat("  - results/survival_results.csv\n")
cat("\n")
