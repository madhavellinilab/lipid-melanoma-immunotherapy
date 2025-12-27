# METHODS

## Data Source and Sample Selection

RNA-sequencing data and clinical annotations for melanoma samples were obtained from The Cancer Genome Atlas Skin Cutaneous Melanoma (TCGA-SKCM) project through the TCGAbiolinks R package (v2.28.0). The TCGA-SKCM cohort comprises 473 melanoma samples with matched clinical data. Following quality control procedures, 335 samples with complete RNA-seq data (unstranded gene expression quantification) and survival information were included in the final analysis.

## Lipid Metabolism Gene Signature

We defined an 8-gene lipid metabolism signature representing key pathways in lipid metabolic reprogramming:

1. **De novo lipogenesis**: ACLY (ATP citrate lyase), FASN (fatty acid synthase)
2. **Fatty acid modification**: SCD (stearoyl-CoA desaturase)
3. **Transcriptional regulation**: SREBF1 (sterol regulatory element-binding transcription factor 1)
4. **Fatty acid oxidation**: CPT1A (carnitine palmitoyltransferase 1A)
5. **Lipid transport**: FABP5 (fatty acid binding protein 5)
6. **Fatty acid activation**: ACSL5 (acyl-CoA synthetase long chain family member 5)
7. **Ketogenesis**: HMGCS2 (3-hydroxy-3-methylglutaryl-CoA synthase 2)

Gene selection was based on established roles in cancer lipid metabolism and previous literature demonstrating their relevance in tumor progression and immune evasion.

## Data Processing and Quality Control

Raw RNA-seq count data were processed following standard TCGA guidelines:

1. **Gene filtering**: Genes with mean expression < 1 count across all samples were removed
2. **Sample filtering**: Samples with >10% missing gene expression values were excluded
3. **Normalization**: Log2(x+1) transformation was applied to stabilize variance
4. **Quality metrics**: Sample-level QC included assessment of library size, gene detection rate, and technical covariates

Final quality-controlled dataset comprised 335 samples and 60,660 genes, saved as a SummarizedExperiment object for downstream analyses.

## Lipid Metabolism Signature Scoring

For each sample, a lipid metabolism score was calculated as the arithmetic mean of log2-transformed expression values across the 8-gene panel:

**Lipid Score = (1/8) × Σ log2(gene expression + 1)**

To create categorical groups for survival analysis, samples were stratified by median split:
- **Lipid-High**: Lipid score > median
- **Lipid-Low**: Lipid score ≤ median

Z-score normalization was applied gene-wise for visualization purposes (heatmap generation).

## CD8+ T-Cell Infiltration Estimation

CD8+ T-cell infiltration was estimated using an 8-gene signature comprising established markers of cytotoxic T lymphocytes:

- **T-cell markers**: CD8A, CD8B, CD3D, CD3E
- **Cytotoxicity markers**: GZMA (granzyme A), GZMB (granzyme B), PRF1 (perforin 1)
- **Effector function**: IFNG (interferon gamma)

A CD8+ T-cell score was calculated as the mean log2 expression across available genes, with samples stratified by median split into CD8-High and CD8-Low groups.

## Statistical Analyses

### Correlation Analysis
Pearson correlation was used to assess the relationship between continuous lipid metabolism scores and CD8+ T-cell scores. Statistical significance was determined at α = 0.05.

### Survival Analysis
Overall survival (OS) was defined as time from diagnosis to death or last follow-up. Survival analyses included:

1. **Kaplan-Meier estimation**: Survival curves were generated for lipid metabolism groups and visualized using the survminer R package
2. **Log-rank test**: Group differences were assessed using the log-rank test
3. **Univariate Cox regression**: Hazard ratios (HR) and 95% confidence intervals (CI) were calculated for continuous lipid scores and categorical groups
4. **Multivariate Cox regression**: Models were adjusted for age at diagnosis, gender, and AJCC pathologic stage when available

### Interaction Analysis
To test the lipid-immune interaction hypothesis, we performed:

1. **Four-group stratification**: Samples were classified into 2×2 groups based on lipid (High/Low) and CD8 (High/Low) status
2. **Pairwise comparison**: Best-case (Lipid-Low + CD8-High) vs worst-case (Lipid-High + CD8-Low) groups were compared using log-rank test
3. **Cox interaction model**: Statistical interaction between continuous lipid and CD8 scores was tested using Cox proportional hazards regression

### Multiple Testing Correction
Given the exploratory nature of this study and the pre-defined hypothesis, nominal p-values are reported without multiple testing correction. P-values < 0.05 were considered statistically significant.

## Software and Reproducibility

All analyses were conducted in R version 4.5.2. Key packages included:

- **Data acquisition**: TCGAbiolinks (v2.28.0)
- **Data processing**: SummarizedExperiment (v1.30.2), tidyverse (v2.0.0)
- **Statistical analysis**: survival (v3.5-7), survminer (v0.4.9)
- **Visualization**: ggplot2 (v3.4.4), pheatmap (v1.0.12)

Complete analysis code and results are available at: https://github.com/madhavellinilab/lipid-melanoma-immunotherapy

## Data Availability

TCGA-SKCM data are publicly available through the Genomic Data Commons (GDC) portal (https://portal.gdc.cancer.gov/). Processed data, analysis scripts, and intermediate results generated in this study are available in the accompanying GitHub repository.

---

**Word Count: ~780 words**

**Status: DRAFT v1 - Ready for your review**

**Notes for revision:**
- Add specific TCGA project ID if required by journal
- Consider adding power analysis section
- May need to expand statistical methods based on journal requirements
- Add IRB/ethics statement (TCGA data are exempt but some journals require statement)
