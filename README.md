# Lipid Metabolic Reprogramming & Immunotherapy Resistance in Melanoma

## Project Overview
Analysis of TCGA-SKCM (n=369) exploring lipid metabolism as a biomarker for immunotherapy resistance in melanoma.

## Research Question
Does elevated lipid metabolism gene expression predict poor anti-PD-1 response through CD8+ T cell exclusion?

## Analysis Pipeline
1. Download TCGA-SKCM RNA-seq data (60K genes, 369 samples)
2. Score lipid metabolism signature (8-gene panel)
3. Estimate CD8+ T cell infiltration (immune deconvolution)
4. Test association with survival (Cox regression, Kaplan-Meier)

## Key Files
- `scripts/01_download_tcga.R` → Download TCGA data
- `scripts/06_survival_analysis.R` → Cox models + figures

## Status
[Update weekly]
