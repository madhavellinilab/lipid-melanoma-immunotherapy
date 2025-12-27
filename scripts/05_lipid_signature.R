# =============================================================================
# LIPID MELANOMA IMMUNOTHERAPY RESEARCH
# Script 05: Lipid Metabolism Signature Scoring (Phase 3)
# Purpose: Define 8-gene lipid panel and compute per-sample signature scores
# Author: madhavellinilab
# Date: Week 2 - Phase 3.1
# =============================================================================

library(SummarizedExperiment)
library(tidyverse)
library(pheatmap)

# Set working directory
setwd("/path/to/project")

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

cat("\n=======================================")
cat("\nPHASE 3: LIPID SIGNATURE SCORING")
cat("\n=======================================\n")

# =============================================================================
# STEP 1: DEFINE LIPID GENE PANEL
# =============================================================================

cat("\n--- Step 1: Defining 8-gene lipid metabolism panel ---\n")

# Fixed 8-gene lipid panel with metabolic roles
lipid_panel <- data.frame(
  display_name = c("ACLY", "FASN", "SCD1", "SREBP1", "CPT1A", "FABP5", "ACSL5", "HMGCS2"),
  gene_symbol = c("ACLY", "FASN", "SCD", "SREBF1", "CPT1A", "FABP5", "ACSL5", "HMGCS2"),
  metabolic_role = c(
    "de novo lipogenesis",
    "fatty acid synthesis",
    "fatty acid desaturation",
    "lipid transcription factor",
    "fatty acid oxidation",
    "lipid transport",
    "fatty acid activation",
    "ketogenesis"
  ),
  stringsAsFactors = FALSE
)

cat("\nLipid metabolism gene panel:\n")
print(lipid_panel, row.names = FALSE)

# =============================================================================
# STEP 2: LOAD QC-FILTERED DATA AND VERIFY GENE PRESENCE
# =============================================================================

cat("\n--- Step 2: Loading QC-filtered TCGA-SKCM data ---\n")

# Load the QC-filtered SummarizedExperiment
if (!file.exists("data/tcga_skcm_qc.rds")) {
  stop("Error: QC-filtered data not found. Please run 03_quality_control.R first.")
}

tcga_skcm_qc <- readRDS("data/tcga_skcm_qc.rds")

cat(sprintf("Loaded data: %d genes x %d samples\n", 
            nrow(tcga_skcm_qc), ncol(tcga_skcm_qc)))

# Extract expression matrix and gene symbols
expr_mat <- assay(tcga_skcm_qc, "unstranded")
gene_symbols <- rowData(tcga_skcm_qc)$gene_name

# Verify gene presence in TCGA-SKCM data
cat("\n--- Verifying lipid gene presence in TCGA-SKCM ---\n")

lipid_panel$found <- lipid_panel$gene_symbol %in% gene_symbols
lipid_panel$count <- sapply(lipid_panel$gene_symbol, function(g) {
  sum(gene_symbols == g)
})

cat("\nGene verification results:\n")
for (i in 1:nrow(lipid_panel)) {
  status <- if (lipid_panel$found[i]) "✓ FOUND" else "✗ MISSING"
  cat(sprintf("  %s %s (%s): %s [n=%d]\n",
              status,
              lipid_panel$display_name[i],
              lipid_panel$gene_symbol[i],
              lipid_panel$metabolic_role[i],
              lipid_panel$count[i]))
}

# Check if all genes are present
genes_found <- sum(lipid_panel$found)
cat(sprintf("\nSummary: %d/%d genes found in dataset\n", 
            genes_found, nrow(lipid_panel)))

if (genes_found < nrow(lipid_panel)) {
  warning("Not all lipid genes found in dataset. Proceeding with available genes.")
}

# =============================================================================
# STEP 3: EXTRACT AND NORMALIZE EXPRESSION
# =============================================================================

cat("\n--- Step 3: Extracting and normalizing lipid gene expression ---\n")

# Get available lipid genes
available_genes <- lipid_panel$gene_symbol[lipid_panel$found]

# Extract lipid gene expression matrix
lipid_expr_idx <- which(gene_symbols %in% available_genes)
lipid_expr <- expr_mat[lipid_expr_idx, , drop = FALSE]

# Add gene symbols as rownames
rownames(lipid_expr) <- gene_symbols[lipid_expr_idx]

cat(sprintf("Extracted expression for %d lipid genes\n", nrow(lipid_expr)))

# Apply log2(x+1) transformation for normalization
cat("\nApplying log2(x+1) transformation...\n")
lipid_expr_log <- log2(lipid_expr + 1)

cat(sprintf("Expression range before: [%.2f, %.2f]\n", 
            min(lipid_expr), max(lipid_expr)))
cat(sprintf("Expression range after: [%.2f, %.2f]\n", 
            min(lipid_expr_log), max(lipid_expr_log)))

# =============================================================================
# STEP 4: COMPUTE LIPID SCORE AND GROUPS
# =============================================================================

cat("\n--- Step 4: Computing per-sample lipid metabolism score ---\n")

# Compute per-sample lipid score (arithmetic mean of log2 expression)
lipid_score <- colMeans(lipid_expr_log, na.rm = TRUE)

cat(sprintf("Computed lipid scores for %d samples\n", length(lipid_score)))

# Derive categorical variable by median split
median_score <- median(lipid_score, na.rm = TRUE)
lipid_group <- ifelse(lipid_score > median_score, "Lipid-High", "Lipid-Low")
lipid_group <- factor(lipid_group, levels = c("Lipid-Low", "Lipid-High"))

cat(sprintf("\nMedian lipid score: %.3f\n", median_score))
cat("\nGroup distribution:\n")
print(table(lipid_group))

# =============================================================================
# STEP 5: QUALITY CONTROL & VISUALIZATION OF SIGNATURE
# =============================================================================

cat("\n--- Step 5: Quality control and visualization ---\n")

# Summary statistics
cat("\nLipid score summary statistics:\n")
cat(sprintf("  Mean: %.3f\n", mean(lipid_score, na.rm = TRUE)))
cat(sprintf("  SD: %.3f\n", sd(lipid_score, na.rm = TRUE)))
cat(sprintf("  Median: %.3f\n", median_score))
cat(sprintf("  Range: [%.3f, %.3f]\n", 
            min(lipid_score, na.rm = TRUE), 
            max(lipid_score, na.rm = TRUE)))
cat(sprintf("  IQR: %.3f\n", IQR(lipid_score, na.rm = TRUE)))

cat("\nCounts per group:\n")
cat(sprintf("  Lipid-Low: %d\n", sum(lipid_group == "Lipid-Low")))
cat(sprintf("  Lipid-High: %d\n", sum(lipid_group == "Lipid-High")))

# Plot histogram of lipid score
cat("\nGenerating lipid score distribution plot...\n")

pdf("figures/02_lipid_distribution.pdf", width = 8, height = 6)
par(mar = c(5, 5, 4, 2))
hist(lipid_score, 
     breaks = 30, 
     col = "lightblue", 
     border = "white",
     xlab = "Lipid Metabolism Score (mean log2 expression)",
     ylab = "Frequency",
     main = "Distribution of Lipid Metabolism Signature Score\nTCGA-SKCM (n=samples)",
     cex.lab = 1.2,
     cex.main = 1.3)
abline(v = median_score, col = "red", lwd = 2, lty = 2)
legend("topright", 
       legend = sprintf("Median = %.3f", median_score),
       col = "red", lwd = 2, lty = 2,
       bty = "n", cex = 1.1)
dev.off()

cat("Saved: figures/02_lipid_distribution.pdf\n")

# =============================================================================
# STEP 6: GENE-LEVEL SIGNATURE VISUALIZATION
# =============================================================================

cat("\n--- Step 6: Generating gene-level heatmap ---\n")

# Z-score scale each gene across samples
lipid_expr_z <- t(scale(t(lipid_expr_log)))

# Order samples by lipid score
sample_order <- order(lipid_score)

# Prepare annotation for heatmap
annotation_col <- data.frame(
  Lipid_Group = lipid_group,
  row.names = colnames(lipid_expr_z)
)

annotation_colors <- list(
  Lipid_Group = c("Lipid-Low" = "#377EB8", "Lipid-High" = "#E41A1C")
)

# Generate heatmap without column clustering
cat("Generating heatmap...\n")

pdf("figures/lipid_genes_heatmap.pdf", width = 12, height = 6)
pheatmap(lipid_expr_z[, sample_order],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-3, 3, length.out = 101),
         main = "Lipid Metabolism Gene Expression Signature\n(Z-scores ordered by lipid score)",
         fontsize = 10,
         fontsize_row = 9)
dev.off()

cat("Saved: figures/lipid_genes_heatmap.pdf\n")

# =============================================================================
# STEP 7: PERSIST SIGNATURE FOR DOWNSTREAM ANALYSIS
# =============================================================================

cat("\n--- Step 7: Saving lipid signature results ---\n")

# Create results dataframe
lipid_signature_results <- data.frame(
  sample_id = colnames(tcga_skcm_qc),
  lipid_score = lipid_score,
  lipid_group = lipid_group,
  stringsAsFactors = FALSE
)

# Save as RDS
saveRDS(lipid_signature_results, "results/lipid_score.rds")
cat("Saved: results/lipid_score.rds\n")

# Save as CSV
write_csv(lipid_signature_results, "results/lipid_score.csv")
cat("Saved: results/lipid_score.csv\n")

# Save lipid expression matrix for further analysis
saveRDS(lipid_expr_log, "results/lipid_expression_log2.rds")
cat("Saved: results/lipid_expression_log2.rds\n")

cat("\n=======================================\n")
cat("PHASE 3 COMPLETE: Lipid signature scoring finished")
cat("\n=======================================\n")
cat("\nNext steps:\n")
cat("  1. Run CD8+ T-cell deconvolution (Phase 3.2)\n")
cat("  2. Merge with clinical data for survival analysis\n")
cat("  3. Perform correlation analyses with immune markers\n")
cat("\nFiles generated:\n")
cat("  - figures/02_lipid_distribution.pdf\n")
cat("  - figures/lipid_genes_heatmap.pdf\n")
cat("  - results/lipid_score.rds\n")
cat("  - results/lipid_score.csv\n")
cat("  - results/lipid_expression_log2.rds\n")
cat("\n")
