# ==========================================
# TCGA-SKCM Data Download Script
# Script: 01_download_tcga.R  
# Purpose: Download TCGA-SKCM RNA-seq data
# Expected runtime: 30-45 minutes
# ==========================================

cat("\n=== TCGA-SKCM Data Download ===")
cat("\nStarting download...\n\n")

# Load required libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

# Create output directory
dir.create("data/TCGA", showWarnings = FALSE, recursive = TRUE)

# Query TCGA-SKCM data
cat("Step 1: Querying TCGA database for SKCM data...\n")

query <- GDCquery(
  project = "TCGA-SKCM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

cat("  Found", length(query$results[[1]]$cases), "samples\n")

# Download data (this takes ~30-45 minutes)
cat("\nStep 2: Downloading data (this may take 30-45 minutes)...\n")
cat("  Progress will be shown below:\n\n")

GDCdownload(
  query,
  method = "api",
  files.per.chunk = 10
)

cat("\n  Download complete!\n")

# Prepare data
cat("\nStep 3: Preparing data...\n")

data <- GDCprepare(query)

cat("  Data dimensions:", dim(data), "\n")
cat("  Samples:", ncol(data), "\n")
cat("  Genes:", nrow(data), "\n")

# Save raw data
cat("\nStep 4: Saving raw data...\n")

saveRDS(data, "data/TCGA/tcga_skcm_raw.rds")
cat("  Saved: data/TCGA/tcga_skcm_raw.rds\n")

# Extract and save metadata
metadata <- as.data.frame(colData(data))
write.csv(metadata, "data/metadata.csv", row.names = FALSE)
cat("  Saved: data/metadata.csv\n")

# Print summary
cat("\n=== Download Summary ===\n")
cat("Total samples:", ncol(data), "\n")
cat("Total genes:", nrow(data), "\n")
cat("File size:", round(file.size("data/TCGA/tcga_skcm_raw.rds") / 1024^2, 1), "MB\n")

cat("\n=== Download Complete! ===\n")
cat("Next step: Run 02_explore_tcga.R\n\n")
