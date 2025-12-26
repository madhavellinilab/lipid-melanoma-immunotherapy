# STUDY PROTOCOL: Lipid Metabolism & Immunotherapy Resistance

## Pre-Specified Analysis Plan
(Written BEFORE analyzing data—prevents p-hacking)

### Primary Research Question
Does elevated lipid metabolism gene expression predict poor immunotherapy response through CD8+ T cell exclusion in melanoma?

### Hypothesis
Lipid-high tumors (high expression of ACLY, FASN, CPT1A, SCD1, FABP5, ACSL5, HMGCS2, SREBP1) have:
1. Lower CD8+ T-cell infiltration (Spearman ρ < -0.2, p < 0.05)
2. Worse overall survival (log-rank p < 0.05, HR > 1.2)

### Primary Outcomes
1. Correlation: Lipid score ↔ CD8+ fraction (Spearman ρ, 95% CI)
2. Survival: Lipid-high vs low (HR, 95% CI, log-rank p-value)

### Secondary Outcomes
1. Effect adjustment (lipid effect independent of stage/age?)
2. Combination analysis (lipid + CD8 together)
3. Pathway analysis (which lipid genes drive effect?)

### Sample Size
N = 369 primary melanomas (TCGA-SKCM), survival events ~150

### Statistical Approach
- Spearman correlation (non-parametric)
- Kaplan-Meier curves (log-rank test)
- Cox proportional hazards regression
- Multiple testing: FDR correction if >5 tests
- Sensitivity analyses: exclude stage IV, exclude short follow-up (<6 mo)

### Pre-Specified Figures
1. Lipid signature definition + gene expression heatmap
2. Lipid vs CD8 scatter plot + boxplot
3. Kaplan-Meier survival curves
4. Cox model forest plot
5. Combination analysis (lipid + CD8)

### Exclusion/Inclusion Criteria
- Include: Primary melanomas with complete survival follow-up
- Exclude: Metastatic tumors, samples with missing RNA-seq

### Success Criteria
✓ Lipid-CD8 correlation: p < 0.01
✓ Survival difference: log-rank p < 0.05
✓ Publication: Accepted at Nature-tier journal

---
Protocol Date: December 26, 2025
PI: [Your Name]
Analyst: [Your Name]
