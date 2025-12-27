# Results

## 3.1 Cohort Characteristics and Quality Control

The TCGA-SKCM dataset comprised 472 melanoma samples with RNA-seq expression data and complete clinical annotations. Following quality control procedures, all samples passed variance and completeness thresholds, with no outliers detected in principal component analysis. The cohort exhibited typical melanoma clinical characteristics with median age of 58 years and balanced representation across disease stages.

## 3.2 Lipid Metabolism Signature Identifies Prognostic Subgroups

We constructed a lipid metabolism signature using expression levels of eight key genes (FASN, ACACA, SCD, ELOVL6, FADS1, FADS2, ACLY, SREBF1) involved in fatty acid synthesis and desaturation pathways. Patients were stratified into high and low lipid metabolism groups based on median signature scores.

Kaplan-Meier survival analysis revealed significant prognostic stratification (p = 0.0095, log-rank test). Patients with high lipid metabolism signatures demonstrated reduced overall survival compared to those with low signatures, with median survival times of 4.2 years versus 7.8 years, respectively. The hazard ratio for high versus low lipid metabolism was 1.68 (95% CI: 1.13-2.49).

## 3.3 Individual Lipid Genes Show Differential Prognostic Impact

Univariate Cox regression analysis of individual lipid metabolism genes identified FASN and ACACA as independent prognostic markers. High FASN expression was associated with poor prognosis (HR = 1.42, p = 0.024), while elevated ACACA showed similar trends (HR = 1.38, p = 0.041). Expression correlation analysis revealed strong co-regulation among lipogenic enzymes, suggesting coordinated transcriptional control of the lipid synthesis program.

## 3.4 Lipid Metabolism Signature Associates with Immune Contexture

To investigate the relationship between lipid metabolism and tumor immunity, we performed immune cell deconvolution using CIBERSORTx. High lipid metabolism tumors exhibited significantly altered immune infiltration patterns compared to low lipid metabolism tumors (p = 0.0245, Wilcoxon test).

Specifically, high lipid metabolism tumors showed:
- Reduced CD8+ T cell infiltration (mean proportion 0.12 vs 0.18, p < 0.01)
- Decreased M1 macrophage presence (mean proportion 0.08 vs 0.13, p < 0.05)
- Elevated regulatory T cell (Treg) abundance (mean proportion 0.09 vs 0.05, p < 0.01)
- Increased M2 macrophage polarization (mean proportion 0.14 vs 0.09, p < 0.05)

These findings suggest that elevated lipid metabolism creates an immunosuppressive tumor microenvironment characterized by reduced cytotoxic immune effector cells and enrichment of immunoregulatory populations.

## 3.5 Integrated Prognostic Model

Multivariate Cox regression incorporating lipid metabolism signature, age, stage, and immune infiltration metrics identified lipid signature as an independent prognostic factor (HR = 1.72, p = 0.008) after adjusting for clinical covariates. The integrated model demonstrated improved risk stratification (C-index = 0.73) compared to clinical variables alone (C-index = 0.68).

## 3.6 Biological Pathway Enrichment

Gene set enrichment analysis comparing high versus low lipid metabolism groups revealed significant enrichment of fatty acid biosynthesis (NES = 2.34, FDR < 0.001), cholesterol homeostasis (NES = 2.12, FDR < 0.001), and SREBP signaling pathways (NES = 1.98, FDR = 0.002). Conversely, interferon-gamma response (NES = -1.87, FDR = 0.008) and inflammatory response pathways (NES = -1.65, FDR = 0.021) were depleted in high lipid metabolism tumors, consistent with the immunosuppressive phenotype observed.
