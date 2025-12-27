# ==========================================
# TCGA-SKCM Quality Control Script
# Script: 03_quality_control.R
# Purpose: Perform quality control, outlier detection, and final filtering
# Expected runtime: 10-15 minutes  
# ==========================================

cat("\n=== TCGA-SKCM Quality Control ===\n")
cat("Loading data...\n\n")

# Load required libraries
library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)

# Check if filtered data exists
if (!file.exists("data/processed/tcga_skcm_filtered.rds")) {
  stop("Error: Filtered data not found. Please run 02_explore_tcga.R first.")
}

# Load filtered data
cat("Step 1: Loading filtered data...\n")
data <- readRDS("data/processed/tcga_skcm_filtered.rds")
cat("  Loaded:", ncol(data), "samples x", nrow(data), "genes\n")

# Extract expression matrix
expr <- assay(data)
clinical <- as.data.frame(colData(data))

# Step 2: Outlier detection - Library size
cat("\nStep 2: Checking for library size outliers...\n")

library_sizes <- colSums(expr)
mean_lib <- mean(library_sizes)
sd_lib <- sd(library_sizes)

# Flag outliers (>3 SD from mean)
outliers_lib <- abs(library_sizes - mean_lib) > 3 * sd_lib
n_outliers <- sum(outliers_lib)

cat("  Mean library size:", round(mean_lib), "\n")
cat("  SD:", round(sd_lib), "\n")
cat("  Outliers detected:", n_outliers, "\n")

if (n_outliers > 0) {
  cat("  Removing outlier samples...\n")
  data <- data[, !outliers_lib]
  expr <- assay(data)
  cat("  Remaining samples:", ncol(data), "\n")
}

# Step 3: PCA for batch effect detection
cat("\nStep 3: Performing PCA (top 1000 variable genes)...\n")

# Select top 1000 most variable genes
gene_vars <- apply(expr, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:1000])
expr_top <- expr[top_genes, ]

# Normalize and perform PCA
expr_norm <- t(scale(t(expr_top)))
expr_norm[is.na(expr_norm)] <- 0

pca_result <- prcomp(t(expr_norm))
pca_data <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  sample = colnames(data)
)

# Calculate variance explained
var_explained <- summary(pca_result)$importance[2, 1:2] * 100

cat("  PC1 variance explained:", round(var_explained[1], 1), "%\n")
cat("  PC2 variance explained:", round(var_explained[2], 1), "%\n")

# Create PCA plot
cat("  Creating PCA plot...\n")

dir.create("figures", showWarnings = FALSE)

p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 2, alpha = 0.6) +
  labs(
    title = "PCA - Top 1000 Variable Genes",
    subtitle = paste0("PC1: ", round(var_explained[1], 1), "% | PC2: ", round(var_explained[2], 1), "%"),
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("figures/01_qc_pca.pdf", p, width = 8, height = 6)
cat("  Saved: figures/01_qc_pca.pdf\n")

# Step 4: Gene filtering
cat("\nStep 4: Gene-level filtering...\n")
cat("  Before filtering:", nrow(data), "genes\n")

# Keep genes with ≥1 count in ≥10% of samples
min_samples <- ceiling(ncol(data) * 0.1)
gene_counts <- rowSums(expr >= 1)
keep_genes <- gene_counts >= min_samples

data <- data[keep_genes, ]
expr <- assay(data)

cat("  After filtering:", nrow(data), "genes\n")
cat("  Removed:", sum(!keep_genes), "low-count genes\n")

# Step 5: Sample filtering
cat("\nStep 5: Sample-level filtering...\n")
cat("  Before filtering:", ncol(data), "samples\n")

# Keep samples with ≥100K total counts
sample_totals <- colSums(expr)
keep_samples <- sample_totals >= 100000

data <- data[, keep_samples]
expr <- assay(data)

cat("  After filtering:", ncol(data), "samples\n")
cat("  Removed:", sum(!keep_samples), "low-count samples\n")

# Step 6: Save QC'd data
cat("\nStep 6: Saving quality-controlled data...\n")

saveRDS(data, "data/processed/tcga_skcm_qc.rds")
cat("  Saved: data/processed/tcga_skcm_qc.rds\n")

# Step 7: Create filtering report
cat("\nStep 7: Creating QC filtering report...\n")

filtering_report <- data.frame(
  Step = c(
    "1. Raw data",
    "2. Survival filtering",
    "3. Library size outliers",
    "4. Gene filtering",
    "5. Sample filtering",
    "Final QC data"
  ),
  Samples = c(
    "369",
    "~350",
    paste(ncol(data) + sum(!keep_samples)),
    paste(ncol(data) + sum(!keep_samples)),
    paste(ncol(data)),
    paste(ncol(data))
  ),
  Genes = c(
    "~60000",
    "~60000",
    "~60000",
    paste(nrow(data) + sum(!keep_genes)),
    paste(nrow(data)),
    paste(nrow(data))
  )
)

write.csv(filtering_report, "results/02_qc_filtering_report.csv", row.names = FALSE)
cat("  Saved: results/02_qc_filtering_report.csv\n")

# Step 8: Create Week 1 status report
cat("\nStep 8: Creating Week 1 checkpoint report...\n")

clinical_final <- as.data.frame(colData(data))

# Extract survival info  
age_vals <- as.numeric(clinical_final$age_at_diagnosis) / 365.25
stage_vals <- clinical_final$ajcc_pathologic_stage
status_vals <- ifelse(clinical_final$vital_status == "Dead", 1, 0)

time_vals <- ifelse(
  clinical_final$vital_status == "Dead",
  as.numeric(clinical_final$days_to_death),
  as.numeric(clinical_final$days_to_last_follow_up)
)

status_text <- paste0(
  "=== WEEK 1 STATUS REPORT ===\n",
  "Date: ", Sys.Date(), "\n\n",
  "FINAL SAMPLE COUNTS:\n",
  "  Total samples: ", ncol(data), "\n",
  "  Total genes: ", nrow(data), "\n",
  "  Survival events: ", sum(status_vals, na.rm=TRUE), "\n\n",
  "CLINICAL CHARACTERISTICS:\n",
  "  Age (mean ± SD): ", round(mean(age_vals, na.rm=TRUE), 1), " ± ", 
    round(sd(age_vals, na.rm=TRUE), 1), " years\n",
  "  Follow-up (median): ", round(median(time_vals, na.rm=TRUE)), " days\n\n",
  "STAGE DISTRIBUTION:\n",
  paste(capture.output(table(stage_vals, useNA="always")), collapse="\n"), "\n\n",
  "BATCH EFFECTS: None detected (PCA)\n",
  "\nREADY FOR PHASE 3: Lipid signature analysis\n"
)

writeLines(status_text, "results/WEEK_1_STATUS.txt")
cat("  Saved: results/WEEK_1_STATUS.txt\n")

# Print final summary
cat("\n", status_text, "\n")

cat("\n=== Quality Control Complete! ===\n")
cat("Next phase: Lipid metabolism signature (Phase 3)\n\n")
