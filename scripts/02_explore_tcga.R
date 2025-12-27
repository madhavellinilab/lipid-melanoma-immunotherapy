# ==========================================
# TCGA-SKCM Data Exploration Script
# Script: 02_explore_tcga.R
# Purpose: Extract clinical data and filter to survival-complete samples
# Expected runtime: 5-10 minutes
# ==========================================

cat("\n=== TCGA-SKCM Data Exploration ===\n")
cat("Loading data...\n\n")

# Load required libraries
library(SummarizedExperiment)
library(tidyverse)

# Check if raw data exists
if (!file.exists("data/TCGA/tcga_skcm_raw.rds")) {
  stop("Error: Raw data not found. Please run 01_download_tcga.R first.")
}

# Load raw data
cat("Step 1: Loading raw TCGA data...\n")
data <- readRDS("data/TCGA/tcga_skcm_raw.rds")
cat("  Loaded:", ncol(data), "samples x", nrow(data), "genes\n")

# Extract clinical data
cat("\nStep 2: Extracting clinical and survival data...\n")

clinical <- as.data.frame(colData(data))

# Key variables for survival analysis
cat("  Extracting key variables:\n")
cat("    - Age\n")
cat("    - Stage\n")
cat("    - Vital status\n")
cat("    - Days to death/last follow-up\n")

# Create survival variables
clinical_clean <- clinical %>%
  mutate(
    # Age
    age = as.numeric(age_at_diagnosis) / 365.25,
    
    # Stage (simplified)
    stage = case_when(
      grepl("Stage I", ajcc_pathologic_stage, ignore.case = TRUE) ~ "I",
      grepl("Stage II", ajcc_pathologic_stage, ignore.case = TRUE) ~ "II",
      grepl("Stage III", ajcc_pathologic_stage, ignore.case = TRUE) ~ "III",
      grepl("Stage IV", ajcc_pathologic_stage, ignore.case = TRUE) ~ "IV",
      TRUE ~ NA_character_
    ),
    
    # Vital status (1 = dead, 0 = alive)
    status = ifelse(vital_status == "Dead", 1, 0),
    
    # Time to event (days)
    time = ifelse(
      vital_status == "Dead",
      as.numeric(days_to_death),
      as.numeric(days_to_last_follow_up)
    )
  )

# Filter to samples with complete survival data
cat("\nStep 3: Filtering to survival-complete samples...\n")
cat("  Before filtering:", nrow(clinical_clean), "samples\n")

clinical_filtered <- clinical_clean %>%
  filter(
    !is.na(age),
    !is.na(time),
    time > 0,
    !is.na(status)
  )

cat("  After filtering:", nrow(clinical_filtered), "samples\n")
cat("  Survival events (deaths):", sum(clinical_filtered$status), "\n")

# Filter expression data to match
cat("\nStep 4: Filtering expression data...\n")

keep_samples <- rownames(clinical_filtered)
data_filtered <- data[, colnames(data) %in% keep_samples]

cat("  Filtered data:", ncol(data_filtered), "samples x", nrow(data_filtered), "genes\n")

# Remove zero-variance genes
cat("\nStep 5: Removing zero-variance genes...\n")

expr_matrix <- assay(data_filtered)
gene_vars <- apply(expr_matrix, 1, var)
keep_genes <- gene_vars > 0

data_filtered <- data_filtered[keep_genes, ]

cat("  Removed", sum(!keep_genes), "zero-variance genes\n")
cat("  Remaining:", nrow(data_filtered), "genes\n")

# Save filtered data
cat("\nStep 6: Saving filtered data...\n")

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
saveRDS(data_filtered, "data/processed/tcga_skcm_filtered.rds")
cat("  Saved: data/processed/tcga_skcm_filtered.rds\n")

# Save expression matrix (compressed)
expr_df <- as.data.frame(assay(data_filtered))
write.csv(expr_df, gzfile("data/processed/expression_matrix.csv.gz"), row.names = TRUE)
cat("  Saved: data/processed/expression_matrix.csv.gz\n")

# Create and save summary report
cat("\nStep 7: Creating summary report...\n")

dir.create("results", showWarnings = FALSE)

summary_text <- paste0(
  "=== TCGA-SKCM Data Summary ===\n",
  "Date: ", Sys.Date(), "\n\n",
  "Sample Statistics:\n",
  "  Total samples: ", ncol(data_filtered), "\n",
  "  Total genes: ", nrow(data_filtered), "\n",
  "  Survival events: ", sum(clinical_filtered$status), "\n\n",
  "Clinical Variables:\n",
  "  Age (mean ± SD): ", round(mean(clinical_filtered$age, na.rm=TRUE), 1), 
  " ± ", round(sd(clinical_filtered$age, na.rm=TRUE), 1), " years\n",
  "  Follow-up (median): ", round(median(clinical_filtered$time, na.rm=TRUE)), " days\n\n",
  "Stage Distribution:\n",
  paste(capture.output(table(clinical_filtered$stage, useNA="always")), collapse="\n")
)

writeLines(summary_text, "results/01_data_summary.txt")
cat("  Saved: results/01_data_summary.txt\n")

# Print summary
cat("\n", summary_text, "\n")

cat("\n=== Exploration Complete! ===\n")
cat("Next step: Run 03_quality_control.R\n\n")
