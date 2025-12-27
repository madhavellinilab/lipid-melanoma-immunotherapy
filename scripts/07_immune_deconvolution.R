# =============================================================================
# LIPID MELANOMA IMMUNOTHERAPY RESEARCH
# Script 07: CD8+ T-Cell Immune Deconvolution
# Purpose: Estimate CD8+ T-cell infiltration using gene expression signatures
# Author: madhavellinilab
# Date: Week 2 - Phase 3.2
# =============================================================================

library(SummarizedExperiment)
library(tidyverse)

# Create output directories
dir.create("results", showWarnings = FALSE)

cat("\n================================================")
cat("\nIMMUNE DECONVOLUTION: CD8+ T-Cell Estimation")
cat("\n================================================\n")

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

cat("\n--- Step 1: Loading data ---\n")

# Load QC-filtered TCGA data
tcga_skcm_qc <- readRDS("data/processed/tcga_skcm_qc.rds")

# Load lipid scores
lipid_scores <- readRDS("results/lipid_score.rds")

cat(sprintf("Loaded %d samples\n", ncol(tcga_skcm_qc)))

# =============================================================================
# STEP 2: DEFINE CD8+ T-CELL SIGNATURE GENES
# =============================================================================

cat("\n--- Step 2: Defining CD8+ T-cell signature ---\n")

# CD8+ T-cell marker genes (established from literature)
cd8_genes <- c(
  "CD8A",    # CD8 alpha chain (primary marker)
  "CD8B",    # CD8 beta chain
  "CD3D",    # T-cell receptor
  "CD3E",    # T-cell receptor
  "GZMA",    # Granzyme A (cytotoxicity)
  "GZMB",    # Granzyme B (cytotoxicity)
  "PRF1",    # Perforin (cytotoxicity)
  "IFNG"     # Interferon gamma (effector function)
)

cat("\nCD8+ T-cell signature genes:\n")
for (gene in cd8_genes) {
  cat(sprintf("  - %s\n", gene))
}

# =============================================================================
# STEP 3: EXTRACT AND NORMALIZE EXPRESSION
# =============================================================================

cat("\n--- Step 3: Extracting CD8 gene expression ---\n")

# Extract expression matrix
expr_mat <- assay(tcga_skcm_qc, "unstranded")
gene_symbols <- rowData(tcga_skcm_qc)$gene_name

# Verify gene presence
cd8_found <- cd8_genes %in% gene_symbols
cat("\nGene verification:\n")
for (i in 1:length(cd8_genes)) {
  status <- if (cd8_found[i]) "✓ FOUND" else "✗ MISSING"
  cat(sprintf("  %s %s\n", status, cd8_genes[i]))
}

cat(sprintf("\nTotal: %d/%d genes found\n", sum(cd8_found), length(cd8_genes)))

# Extract available CD8 genes
available_cd8 <- cd8_genes[cd8_found]
cd8_expr_idx <- which(gene_symbols %in% available_cd8)
cd8_expr <- expr_mat[cd8_expr_idx, , drop = FALSE]
rownames(cd8_expr) <- gene_symbols[cd8_expr_idx]

# Log2 transform if needed
if (max(cd8_expr) > 100) {
  cd8_expr <- log2(cd8_expr + 1)
  cat("Applied log2 transformation\n")
}

# =============================================================================
# STEP 4: CALCULATE CD8+ T-CELL SCORE
# =============================================================================

cat("\n--- Step 4: Computing CD8+ T-cell infiltration score ---\n")

# Compute CD8 score (mean expression of signature genes)
cd8_score <- colMeans(cd8_expr, na.rm = TRUE)

cat(sprintf("Computed CD8 scores for %d samples\n", length(cd8_score)))

# Create categorical groups (median split)
median_cd8 <- median(cd8_score, na.rm = TRUE)
cd8_group <- ifelse(cd8_score > median_cd8, "CD8-High", "CD8-Low")
cd8_group <- factor(cd8_group, levels = c("CD8-Low", "CD8-High"))

cat(sprintf("\nMedian CD8 score: %.3f\n", median_cd8))
cat("\nGroup distribution:\n")
print(table(cd8_group))

# Summary statistics
cat("\nCD8 score statistics:\n")
cat(sprintf("  Mean: %.3f\n", mean(cd8_score)))
cat(sprintf("  SD: %.3f\n", sd(cd8_score)))
cat(sprintf("  Range: [%.3f, %.3f]\n", min(cd8_score), max(cd8_score)))

# =============================================================================
# STEP 5: CORRELATE CD8 WITH LIPID SIGNATURE
# =============================================================================

cat("\n--- Step 5: Analyzing lipid-immune correlation ---\n")

# Merge CD8 and lipid scores
integrated_data <- data.frame(
  sample_id = names(cd8_score),
  cd8_score = cd8_score,
  cd8_group = cd8_group,
  lipid_score = lipid_scores$lipid_score[match(names(cd8_score), lipid_scores$sample_id)],
  lipid_group = lipid_scores$lipid_group[match(names(cd8_score), lipid_scores$sample_id)]
)

# Test correlation
cor_test <- cor.test(integrated_data$lipid_score, integrated_data$cd8_score)

cat(sprintf("\nCorrelation between lipid and CD8 scores:\n"))
cat(sprintf("  Pearson r = %.3f\n", cor_test$estimate))
cat(sprintf("  p-value = %.4f\n", cor_test$p.value))

# =============================================================================
# STEP 6: CREATE 2x2 STRATIFICATION
# =============================================================================

cat("\n--- Step 6: Creating lipid-immune stratification ---\n")

# Create 2x2 groups
integrated_data$combined_group <- paste0(
  integrated_data$lipid_group, "_", integrated_data$cd8_group
)

integrated_data$combined_group <- factor(
  integrated_data$combined_group,
  levels = c("Lipid-Low_CD8-High", "Lipid-Low_CD8-Low",
             "Lipid-High_CD8-High", "Lipid-High_CD8-Low")
)

cat("\n2x2 Stratification:\n")
print(table(integrated_data$combined_group))

cat("\nContingency table:\n")
print(table(integrated_data$lipid_group, integrated_data$cd8_group))

# Chi-square test for independence
chi_test <- chisq.test(table(integrated_data$lipid_group, integrated_data$cd8_group))
cat(sprintf("\nChi-square test p-value: %.4f\n", chi_test$p.value))

# =============================================================================
# STEP 7: SAVE INTEGRATED RESULTS
# =============================================================================

cat("\n--- Step 7: Saving integrated dataset ---\n")

# Save RDS
saveRDS(integrated_data, "results/integrated_lipid_immune.rds")
cat("Saved: results/integrated_lipid_immune.rds\n")

# Save CSV
write_csv(integrated_data, "results/integrated_lipid_immune.csv")
cat("Saved: results/integrated_lipid_immune.csv\n")

# Save CD8 expression matrix
saveRDS(cd8_expr, "results/cd8_expression_log2.rds")
cat("Saved: results/cd8_expression_log2.rds\n")

cat("\n================================================\n")
cat("IMMUNE DECONVOLUTION COMPLETE")
cat("\n================================================\n")
cat("\nNext steps:\n")
cat("  1. Run lipid-immune interaction survival analysis\n")
cat("  2. Generate integrated visualizations\n")
cat("  3. Test hypothesis: High-Lipid + Low-CD8 = poor outcome\n")
cat("\nFiles generated:\n")
cat("  - results/integrated_lipid_immune.csv\n")
cat("  - results/cd8_expression_log2.rds\n")
cat("\n")
