# =============================================================================
# LIPID MELANOMA IMMUNOTHERAPY RESEARCH
# Script 04: Batch Processing Pipeline
# Purpose: Orchestrate all analysis steps in sequence
# =============================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(limma)
library(survival)

# Set up environment
setwd("/path/to/project")
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# =============================================================================
# STEP 1: DATA DOWNLOAD AND PREPARATION
# =============================================================================

cat("\n=== STEP 1: Downloading TCGA-SKCM data ===\n")
source("scripts/01_download_tcga.R")

# =============================================================================
# STEP 2: DATA EXPLORATION
# =============================================================================

cat("\n=== STEP 2: Exploring TCGA-SKCM data ===\n")
source("scripts/02_explore_tcga.R")

# =============================================================================
# STEP 3: QUALITY CONTROL
# =============================================================================

cat("\n=== STEP 3: Quality control and filtering ===\n")
source("scripts/03_quality_control.R")

# =============================================================================
# STEP 4: LIPID SIGNATURE SCORING (PHASE 3)
# =============================================================================

cat("\n=== STEP 4: Creating lipid metabolism signature ===\n")

# Define 8-gene lipid panel
lipid_genes <- c(
  "ACLY",    # de novo lipogenesis
  "FASN",    # fatty acid synthesis
  "SCD1",    # desaturation
  "SREBP1",  # transcriptional regulation (SREBF1)
  "CPT1A",   # fatty acid oxidation
  "FABP5",   # lipid transport
  "ACSL5",   # fatty acid activation
  "HMGCS2"   # ketogenesis
)

# Map gene symbols to proper identifiers
lipid_genes_map <- c(
  "ACLY" = "ACLY",
  "FASN" = "FASN",
  "SCD1" = "SCD",     # Note: Gene symbol is SCD
  "SREBP1" = "SREBF1", # Note: Gene symbol is SREBF1
  "CPT1A" = "CPT1A",
  "FABP5" = "FABP5",
  "ACSL5" = "ACSL5",
  "HMGCS2" = "HMGCS2"
)

# Load QC-filtered data
load("data/tcga_skcm_qc.rds")

# Extract expression matrix
expr_mat <- assay(tcga_skcm_qc)
gene_symbols <- rowData(tcga_skcm_qc)$gene_name

# Verify gene presence
cat("\nChecking lipid gene panel availability:\n")
for (gene in names(lipid_genes_map)) {
  actual_gene <- lipid_genes_map[gene]
  if (actual_gene %in% gene_symbols) {
    cat(sprintf("  ✓ %s (%s): FOUND\n", gene, actual_gene))
  } else {
    cat(sprintf("  ✗ %s (%s): MISSING\n", gene, actual_gene))
  }
}

# Extract lipid gene expression
lipid_expr <- expr_mat[gene_symbols %in% lipid_genes_map, ]
rownames(lipid_expr) <- gene_symbols[gene_symbols %in% lipid_genes_map]

# Log2 transform if needed (TCGA data is usually log2(x+1) already)
if (max(lipid_expr) > 100) {
  lipid_expr <- log2(lipid_expr + 1)
  cat("\nApplied log2 transformation\n")
}

# Z-score normalization within each gene
lipid_expr_z <- t(scale(t(lipid_expr)))

# Calculate lipid signature score (mean z-score)
lipid_signature <- colMeans(lipid_expr_z, na.rm = TRUE)

# Save lipid signature scores
lipid_scores <- data.frame(
  sample_id = names(lipid_signature),
  lipid_score = lipid_signature,
  lipid_score_percentile = percent_rank(lipid_signature)
)

write_csv(lipid_scores, "results/lipid_signature_scores.csv")
cat("\nLipid signature scores saved to results/lipid_signature_scores.csv\n")

# =============================================================================
# STEP 5: CLINICAL CORRELATION ANALYSIS
# =============================================================================

cat("\n=== STEP 5: Correlating lipid signature with clinical outcomes ===\n")

# Extract clinical data
clinical <- colData(tcga_skcm_qc)

# Merge with lipid scores
clinical_lipid <- cbind(clinical, lipid_scores)

# Define high vs low lipid metabolism groups (median split)
median_score <- median(lipid_signature)
clinical_lipid$lipid_group <- ifelse(
  clinical_lipid$lipid_score > median_score,
  "High", "Low"
)

# Survival analysis
if ("OS.time" %in% colnames(clinical_lipid) && "OS" %in% colnames(clinical_lipid)) {
  surv_obj <- Surv(clinical_lipid$OS.time, clinical_lipid$OS)
  fit <- survfit(surv_obj ~ lipid_group, data = clinical_lipid)
  
  # Save survival plot
  pdf("figures/lipid_signature_survival.pdf", width = 8, height = 6)
  plot(fit, col = c("blue", "red"), lwd = 2,
       xlab = "Time (days)", ylab = "Overall Survival",
       main = "Lipid Metabolism Signature and Survival")
  legend("topright", legend = c("High", "Low"), 
         col = c("red", "blue"), lwd = 2)
  dev.off()
  
  # Log-rank test
  surv_test <- survdiff(surv_obj ~ lipid_group, data = clinical_lipid)
  cat(sprintf("\nLog-rank test p-value: %.4f\n", 
              1 - pchisq(surv_test$chisq, 1)))
}

# =============================================================================
# STEP 6: IMMUNOTHERAPY RESPONSE CORRELATION
# =============================================================================

cat("\n=== STEP 6: Analyzing immunotherapy biomarkers ===\n")

# If immune checkpoint gene expression available
immune_genes <- c("PDCD1", "CD274", "CTLA4", "LAG3")
available_immune <- immune_genes[immune_genes %in% gene_symbols]

if (length(available_immune) > 0) {
  immune_expr <- expr_mat[gene_symbols %in% available_immune, ]
  
  # Correlate with lipid signature
  for (gene in available_immune) {
    cor_test <- cor.test(lipid_signature, 
                         immune_expr[gene_symbols == gene, ])
    cat(sprintf("Correlation with %s: r=%.3f, p=%.4f\n",
                gene, cor_test$estimate, cor_test$p.value))
  }
}

# =============================================================================
# PIPELINE COMPLETION
# =============================================================================

cat("\n=== BATCH PROCESSING COMPLETE ===\n")
cat("Results saved to:\n")
cat("  - results/lipid_signature_scores.csv\n")
cat("  - figures/lipid_signature_survival.pdf\n")
cat("\nNext steps:\n")
cat("  1. Review QC metrics in results/\n")
cat("  2. Examine lipid signature distribution\n")
cat("  3. Validate findings with external datasets\n")
cat("  4. Proceed to mechanistic validation (Week 3)\n")
