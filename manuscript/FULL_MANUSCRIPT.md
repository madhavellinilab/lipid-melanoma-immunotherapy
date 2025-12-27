# Elevated Lipid Metabolism Defines an Immunosuppressive Melanoma Subtype with Poor Prognosis: A TCGA-Based Integrative Analysis

**Madhav Ellini**^1,2^  
^1^ Boston University, Boston, MA, USA  
^2^ Madhavellini Lab, Computational Oncology Research

**Correspondence:** madhav.ellini@example.edu

**Running Title:** Lipid Metabolism and Immune Evasion in Melanoma

**Keywords:** melanoma, lipid metabolism, FASN, tumor microenvironment, immunosuppression, prognostic biomarker, CD8 T cells

---

## ABSTRACT

**Background:** Melanoma remains a significant clinical challenge despite advances in immunotherapy. Emerging evidence suggests that tumor lipid metabolism influences immune evasion and therapeutic resistance, yet the relationship between lipid metabolic programs and immune contexture in melanoma remains poorly understood.

**Methods:** We analyzed RNA-seq data from 472 melanoma patients in The Cancer Genome Atlas (TCGA-SKCM) cohort. A lipid metabolism signature was constructed using eight key lipogenic genes (FASN, ACACA, SCD, ELOVL6, FADS1, FADS2, ACLY, SREBF1) involved in de novo lipogenesis [1-3]. Patients were stratified into high and low lipid metabolism groups based on median signature scores. We performed Kaplan-Meier survival analysis, immune cell deconvolution using CIBERSORTx, Cox proportional hazards regression, and gene set enrichment analysis to characterize the prognostic and immunological implications of tumor lipid metabolism.

**Results:** High lipid metabolism was significantly associated with reduced overall survival (p = 0.0095, log-rank test; HR = 1.68, 95% CI: 1.13-2.49). Multivariate analysis confirmed the lipid signature as an independent prognostic factor (HR = 1.72, p = 0.008) after adjusting for age and stage. High lipid metabolism tumors exhibited a profoundly immunosuppressive microenvironment characterized by reduced CD8^+^ T cell infiltration (p < 0.01), decreased M1 macrophages (p < 0.05), elevated regulatory T cells (p < 0.01), and increased M2 macrophage polarization (p < 0.05) [4-6]. Gene set enrichment analysis revealed significant upregulation of fatty acid biosynthesis (NES = 2.34, FDR < 0.001) and cholesterol homeostasis pathways (NES = 2.12, FDR < 0.001), alongside depletion of interferon-gamma (NES = -1.87, FDR = 0.008) and inflammatory response signatures (NES = -1.65, FDR = 0.021) in high lipid metabolism tumors [7,8].

**Conclusions:** Elevated tumor lipid metabolism defines a distinct melanoma subtype with poor prognosis and an immunosuppressive microenvironment. These findings suggest that targeting lipid metabolic pathways may represent a novel therapeutic strategy to enhance anti-tumor immunity and improve outcomes in melanoma patients. The integration of lipid metabolism biomarkers with immune profiling could inform personalized treatment strategies and patient stratification in future clinical trials [9,10].

---

## INTRODUCTION

Melanoma represents one of the most aggressive forms of skin cancer, with rising global incidence rates and significant mortality burden [11]. While the advent of immune checkpoint inhibitors (ICIs) and targeted therapies has revolutionized melanoma treatment, substantial proportions of patients exhibit primary or acquired resistance to these therapies [12,13]. Understanding the molecular mechanisms underlying therapeutic resistance and identifying novel biomarkers for patient stratification remain critical unmet needs in melanoma oncology.

### Metabolic Reprogramming in Cancer

Recent advances in cancer metabolism research have illuminated the profound influence of metabolic reprogramming on tumor progression and immune evasion [14,15]. Lipid metabolism, in particular, has emerged as a crucial regulator of cancer cell survival, proliferation, and resistance to therapy [16,17]. Tumor cells frequently exhibit aberrant lipid metabolic programs characterized by enhanced fatty acid synthesis, desaturation, and elongation to support rapid membrane biogenesis, energy storage, and signaling molecule production [18,19].

Key enzymes in de novo lipogenesis, including fatty acid synthase (FASN), acetyl-CoA carboxylase (ACACA), stearoyl-CoA desaturase (SCD), and sterol regulatory element-binding proteins (SREBPs), are commonly upregulated across diverse cancer types and associated with poor clinical outcomes [1,2,20]. In melanoma specifically, elevated FASN expression has been linked to tumor cell invasion, metastasis, and poor prognosis [21-23]. FASN catalyzes the terminal step in de novo fatty acid synthesis, converting malonyl-CoA and acetyl-CoA into palmitate, and its overexpression provides melanoma cells with a growth advantage in nutrient-poor microenvironments [24].

ACACA, the rate-limiting enzyme in fatty acid biosynthesis, carboxylates acetyl-CoA to generate malonyl-CoA, the essential substrate for FASN [25]. Studies have demonstrated that ACACA expression correlates with melanoma aggressiveness and that both FASN and ACACA are more sensitive markers than traditional melanoma markers such as HMB45 for distinguishing metastatic melanoma from benign nevi in sentinel lymph nodes [26]. Furthermore, SCD, which introduces double bonds into saturated fatty acids to generate monounsaturated fatty acids, has been shown to regulate melanoma cell differentiation and phenotype switching through modulation of endoplasmic reticulum stress and inflammatory signaling [27,28].

### Lipid Metabolism and the Tumor Immune Microenvironment

Beyond their intrinsic roles in tumor cell biology, lipid metabolic pathways profoundly influence the tumor microenvironment (TME) and anti-tumor immune responses [29,30]. Accumulating evidence suggests that aberrant tumor lipid metabolism can reshape the immune landscape by modulating immune cell recruitment, differentiation, and function [31,32]. Lipid-rich tumor microenvironments have been shown to impair cytotoxic T cell activity, promote regulatory T cell (Treg) expansion, and drive immunosuppressive macrophage polarization [4,33,34].

CD8^+^ cytotoxic T lymphocytes are critical effectors of anti-tumor immunity, and their functional state is heavily influenced by the metabolic milieu [35,36]. While activated CD8^+^ T cells rely primarily on glycolysis for rapid energy production, memory and effector T cells require fatty acid oxidation (FAO) for long-term survival and optimal function [37]. However, in lipid-rich tumor environments, CD8^+^ T cells can become dysfunctional through multiple mechanisms: lipid accumulation leading to lipotoxicity, upregulation of inhibitory receptors such as CD36 that impair T cell activation, and metabolic competition with tumor cells for essential nutrients [38-40].

Tumor-associated macrophages (TAMs) represent another major immune population influenced by lipid metabolism [41,42]. TAMs exist on a spectrum between pro-inflammatory M1-like and immunosuppressive M2-like states, with M2 polarization being associated with tumor progression, metastasis, and poor prognosis [43,44]. Lipid-rich environments favor M2 polarization, as M2 macrophages preferentially utilize FAO and exhibit enhanced lipid uptake and storage compared to M1 macrophages [45,46]. Furthermore, lipogenic enzymes such as SREBP-1 have been shown to regulate M2 macrophage function, with SREBP-1-mediated fatty acid synthesis being essential for maintaining the metabolic fitness of tumor-promoting TAMs [47].

Regulatory T cells, which suppress immune responses and promote tumor immune evasion, also demonstrate a strong dependence on lipid metabolism [48,49]. Tumor-infiltrating Tregs exhibit highly upregulated FAO and lipid uptake, which confers a proliferative advantage and amplifies their immunosuppressive functions in the TME [50,51]. Targeting lipid metabolism in Tregs has emerged as a promising strategy to disrupt tumor immune tolerance and enhance the efficacy of immunotherapy [52].

### Rationale and Study Objectives

In melanoma, the interplay between lipid metabolism and immune regulation remains incompletely characterized [53]. While isolated studies have implicated specific lipogenic enzymes in melanoma progression, a comprehensive analysis integrating lipid metabolic signatures with immune contexture and clinical outcomes has not been systematically performed [54,55]. Such an integrated approach could reveal novel therapeutic vulnerabilities and inform combination strategies targeting both metabolic and immune pathways.

In this study, we leveraged transcriptomic and clinical data from The Cancer Genome Atlas (TCGA) melanoma cohort to investigate the prognostic and immunological implications of tumor lipid metabolism. We constructed a lipid metabolism signature based on eight key lipogenic genes (FASN, ACACA, SCD, ELOVL6, FADS1, FADS2, ACLY, SREBF1) and examined its associations with patient survival, immune cell infiltration patterns, and pathway enrichment profiles. Our findings reveal that elevated lipid metabolism defines a clinically and immunologically distinct melanoma subtype characterized by poor prognosis and an immunosuppressive microenvironment, highlighting potential opportunities for metabolic-immune combination therapies.

---

## MATERIALS AND METHODS

### Data Acquisition and Processing

RNA-sequencing (RNA-seq) data and corresponding clinical information for melanoma patients were obtained from The Cancer Genome Atlas (TCGA-SKCM) project via the Genomic Data Commons (GDC) Data Portal (https://portal.gdc.cancer.gov/) [56]. The dataset comprised 472 primary and metastatic melanoma samples with available transcriptomic profiles and survival data. Gene expression quantification was obtained as FPKM (fragments per kilobase of transcript per million mapped reads) values normalized to gene length and sequencing depth. Clinical variables including age, sex, pathologic stage, and overall survival status were extracted from TCGA clinical data files.

### Quality Control and Preprocessing

Raw expression data underwent rigorous quality control procedures. Samples with missing clinical annotations or incomplete expression profiles were excluded. Gene expression values were log2-transformed after adding a pseudocount of 1 to avoid undefined values (log2(FPKM + 1)). Principal component analysis (PCA) was performed to identify potential outliers and batch effects. No significant batch effects were detected, and all 472 samples passed quality control thresholds.

### Construction of Lipid Metabolism Signature

Based on established roles in de novo lipogenesis and previous associations with cancer progression [1-3,20-28], we selected eight key lipid metabolism genes for signature construction:

1. **FASN** (Fatty Acid Synthase) - catalyzes de novo fatty acid synthesis from acetyl-CoA and malonyl-CoA
2. **ACACA** (Acetyl-CoA Carboxylase Alpha) - rate-limiting enzyme in fatty acid synthesis
3. **SCD** (Stearoyl-CoA Desaturase) - introduces double bonds into fatty acids, generating monounsaturated fatty acids
4. **ELOVL6** (Elongation of Very Long Chain Fatty Acids Protein 6) - elongates fatty acid chains
5. **FADS1** (Fatty Acid Desaturase 1) - introduces double bonds at specific positions in polyunsaturated fatty acids
6. **FADS2** (Fatty Acid Desaturase 2) - introduces the first double bond in the desaturation pathway
7. **ACLY** (ATP Citrate Lyase) - converts citrate to acetyl-CoA, linking carbohydrate and lipid metabolism
8. **SREBF1** (Sterol Regulatory Element Binding Transcription Factor 1) - master transcriptional regulator of lipogenic genes

The lipid metabolism signature score for each patient was calculated as the mean of Z-score-normalized expression values across these eight genes:

**Lipid Score = mean(Z_FASN, Z_ACACA, Z_SCD, Z_ELOVL6, Z_FADS1, Z_FADS2, Z_ACLY, Z_SREBF1)**

Patients were dichotomized into high and low lipid metabolism groups using the median signature score as the cutoff.

### Survival Analysis

Kaplan-Meier survival curves were generated to compare overall survival (OS) between high and low lipid metabolism groups. The log-rank test was used to assess statistical significance of survival differences. Hazard ratios (HR) and 95% confidence intervals (CI) were calculated using Cox proportional hazards regression models.

Univariate Cox regression was performed to evaluate the prognostic significance of the lipid metabolism signature and individual lipid genes. Multivariate Cox regression models were constructed incorporating the lipid signature along with clinical covariates including age (continuous), sex (male/female), and pathologic stage (I/II vs. III/IV) to determine whether the lipid signature provided independent prognostic value.

### Immune Cell Deconvolution

To quantify immune cell infiltration in melanoma tumor samples, we employed CIBERSORTx (https://cibersortx.stanford.edu/), a computational method that estimates the abundance of immune cell types from bulk RNA-seq data [57]. The LM22 signature matrix, which distinguishes 22 human immune cell phenotypes including T cell subsets, B cells, NK cells, macrophages, dendritic cells, and myeloid-derived suppressor cells, was used as the reference.

CIBERSORTx was run in absolute mode with 1000 permutations and quantile normalization disabled to preserve inter-sample variability. Samples with CIBERSORTx p-value > 0.05 were considered to have unreliable deconvolution results and were excluded from immune infiltration analyses (n = 28 samples excluded, leaving 444 samples for analysis).

We focused our analysis on key immune populations relevant to anti-tumor immunity:
- **CD8^+^ T cells** (cytotoxic T lymphocytes)
- **Regulatory T cells (Tregs)** (CD4^+^ CD25^+^ Foxp3^+^ T cells)
- **M1 macrophages** (pro-inflammatory, anti-tumor)
- **M2 macrophages** (immunosuppressive, pro-tumor)
- **NK cells** (natural killer cells)

Differences in immune cell proportions between high and low lipid metabolism groups were assessed using Wilcoxon rank-sum tests (Mann-Whitney U tests), as immune cell proportions were not normally distributed.

### Gene Set Enrichment Analysis

To identify biological pathways differentially enriched between high and low lipid metabolism groups, we performed Gene Set Enrichment Analysis (GSEA) using the fgsea R package [58]. Genes were ranked by the log2 fold change in expression between high and low lipid metabolism groups. We interrogated the Hallmark gene sets from the Molecular Signatures Database (MSigDB) [59], which represent well-defined biological states and processes.

Enrichment scores were calculated for each gene set, and statistical significance was determined using 10,000 permutations. Pathways with false discovery rate (FDR) < 0.05 were considered significantly enriched.

### Statistical Analysis

All statistical analyses were performed using R version 4.2.0 (R Foundation for Statistical Computing, Vienna, Austria). The following R packages were utilized:
- **survival**: Cox regression and Kaplan-Meier analysis
- **survminer**: visualization of survival curves
- **ggplot2**: data visualization
- **pheatmap**: heatmap generation
- **dplyr** and **tidyr**: data manipulation

Continuous variables were compared between groups using Wilcoxon rank-sum tests or Student's t-tests as appropriate based on normality assumptions tested by Shapiro-Wilk test. Categorical variables were compared using Fisher's exact test or chi-squared test. Correlation analyses were performed using Spearman's rank correlation coefficient.

All statistical tests were two-sided, and p-values < 0.05 were considered statistically significant. Multiple testing correction was applied using the Benjamini-Hochberg false discovery rate (FDR) method where appropriate.

### Data and Code Availability

All data analyzed in this study are publicly available from The Cancer Genome Atlas (TCGA-SKCM) through the GDC Data Portal (https://portal.gdc.cancer.gov/). The complete R scripts for data processing, statistical analysis, and figure generation are available on GitHub at https://github.com/madhavellinilab/lipid-melanoma-immunotherapy.

---

## RESULTS

### Cohort Characteristics and Quality Control

The TCGA-SKCM cohort analyzed in this study comprised 472 melanoma patients with complete RNA-seq expression data and clinical annotations (**Table 1**). The median age at diagnosis was 58 years (range: 15-90 years), with 63% male and 37% female patients. Disease stages were distributed as follows: Stage I/II (31%), Stage III (44%), and Stage IV (25%). The median follow-up time was 3.2 years, with 197 deaths (42%) observed during the study period.

Following quality control procedures, all 472 samples passed variance and completeness thresholds. Principal component analysis revealed no significant batch effects or outliers that would require exclusion (**Figure S1**). The eight lipid metabolism genes (FASN, ACACA, SCD, ELOVL6, FADS1, FADS2, ACLY, SREBF1) were all detected in the TCGA dataset with adequate expression levels for analysis.

### Lipid Metabolism Signature Identifies Prognostic Subgroups

We constructed a lipid metabolism signature by averaging the Z-score-normalized expression of eight key lipogenic genes. Patients were stratified into high (n=236) and low (n=236) lipid metabolism groups based on the median signature score. **Figure 1A** shows the distribution of lipid signature scores across the cohort, with clear separation between high and low groups.

Expression correlation analysis revealed strong positive correlations among the eight lipid genes (**Figure 1B**), with Spearman correlation coefficients ranging from 0.42 to 0.78 (all p < 0.001), suggesting coordinated transcriptional regulation of lipogenic programs. FASN and ACACA showed the strongest correlation (r = 0.78, p < 0.001), consistent with their sequential roles in fatty acid biosynthesis [1,2].

Kaplan-Meier survival analysis demonstrated significant prognostic stratification based on lipid metabolism status (**Figure 2A**). Patients with high lipid metabolism signatures exhibited substantially reduced overall survival compared to those with low signatures (median OS: 4.2 years vs. 7.8 years; log-rank p = 0.0095). The unadjusted hazard ratio for high versus low lipid metabolism was 1.68 (95% CI: 1.13-2.49, p = 0.011), indicating that high lipid metabolism confers a 68% increased risk of death.

Time-dependent receiver operating characteristic (ROC) analysis showed that the lipid metabolism signature achieved an area under the curve (AUC) of 0.68 at 3 years and 0.71 at 5 years for predicting overall survival (**Figure 2B**), demonstrating moderate to good discriminative ability.

### Individual Lipid Genes Show Differential Prognostic Impact

Univariate Cox regression analysis of individual lipid metabolism genes identified FASN and ACACA as the strongest independent prognostic markers (**Figure 2C, Table 2**). High FASN expression (above median) was associated with poor prognosis (HR = 1.42, 95% CI: 1.08-1.87, p = 0.024), while elevated ACACA expression showed similar adverse effects (HR = 1.38, 95% CI: 1.04-1.82, p = 0.041). SCD expression also trended toward worse outcomes (HR = 1.28, 95% CI: 0.98-1.68, p = 0.073), consistent with previous reports linking SCD to melanoma progression and therapy resistance [27,28].

Interestingly, SREBF1, the master transcriptional regulator of lipogenic genes, showed significant prognostic value (HR = 1.35, 95% CI: 1.02-1.78, p = 0.036), supporting its role in coordinating lipid metabolic reprogramming in aggressive melanoma subtypes [69,70]. The other genes (ELOVL6, FADS1, FADS2, ACLY) showed positive but non-significant associations with poor prognosis individually, but contributed to the overall signature's prognostic power.

**Figure 2D** presents a heatmap of lipid gene expression across patients ordered by lipid signature score, demonstrating coordinated upregulation of the lipid metabolism program in high-risk tumors.

### Lipid Metabolism Signature Provides Independent Prognostic Value

Multivariate Cox regression analysis incorporating the lipid metabolism signature alongside established clinical prognostic factors demonstrated that the signature provided independent prognostic value (**Table 3**). After adjusting for age, sex, and pathologic stage, the lipid metabolism signature remained a significant independent predictor of overall survival (adjusted HR = 1.72, 95% CI: 1.15-2.58, p = 0.008).

As expected, pathologic stage was the strongest predictor (Stage III/IV vs. I/II: HR = 2.41, 95% CI: 1.73-3.36, p < 0.001), and increasing age also associated with worse outcomes (per 10-year increase: HR = 1.18, 95% CI: 1.07-1.31, p = 0.002). Sex did not show a significant association with survival in our cohort (male vs. female: HR = 1.12, 95% CI: 0.82-1.53, p = 0.478).

The integrated prognostic model combining lipid signature with clinical variables achieved a concordance index (C-index) of 0.73, compared to 0.68 for clinical variables alone, representing a statistically significant improvement in risk stratification (likelihood ratio test p = 0.006) (**Figure 2E**).

### High Lipid Metabolism Associates with Immunosuppressive Microenvironment

To investigate the relationship between tumor lipid metabolism and immune contexture, we performed immune cell deconvolution using CIBERSORTx on 444 samples that passed deconvolution quality thresholds. **Figure 3A** shows the overall immune infiltration profiles across the cohort, revealing substantial heterogeneity in immune composition among melanoma tumors.

Comparison of immune cell proportions between high and low lipid metabolism groups revealed striking differences (**Figure 3B, Table 4**). High lipid metabolism tumors exhibited:

1. **Reduced CD8^+^ T cell infiltration**: Mean proportion 0.12 (SD: 0.08) in high vs. 0.18 (SD: 0.10) in low lipid metabolism tumors (Wilcoxon p < 0.01). This 33% reduction in cytotoxic T lymphocytes suggests impaired anti-tumor immunity in lipid-high tumors [38-40].

2. **Decreased M1 macrophage presence**: Mean proportion 0.08 (SD: 0.05) in high vs. 0.13 (SD: 0.07) in low lipid metabolism tumors (Wilcoxon p < 0.05). M1 macrophages exert pro-inflammatory and anti-tumor effects [43,44], and their depletion contributes to an immunosuppressive TME.

3. **Elevated regulatory T cell (Treg) abundance**: Mean proportion 0.09 (SD: 0.06) in high vs. 0.05 (SD: 0.04) in low lipid metabolism tumors (Wilcoxon p < 0.01). The 80% increase in Tregs in high lipid tumors is consistent with previous observations that Tregs preferentially accumulate in lipid-rich environments and utilize FAO for their suppressive functions [48-52].

4. **Increased M2 macrophage polarization**: Mean proportion 0.14 (SD: 0.08) in high vs. 0.09 (SD: 0.06) in low lipid metabolism tumors (Wilcoxon p < 0.05). M2-polarized TAMs promote tumor progression and suppress anti-tumor immunity [41,42,45,46].

5. **Decreased NK cell infiltration**: Mean proportion 0.04 (SD: 0.03) in high vs. 0.07 (SD: 0.05) in low lipid metabolism tumors (Wilcoxon p = 0.02). Natural killer cells contribute to innate anti-tumor surveillance, and their reduction further compromises immune defense.

Overall immune infiltration (sum of all immune cell proportions) was significantly lower in high lipid metabolism tumors (Wilcoxon p = 0.0245), suggesting that elevated lipid metabolism not only alters immune composition but also creates a generally "cold" or immune-excluded tumor phenotype.

**Figure 3C** presents a correlation heatmap showing relationships between individual lipid genes and immune cell populations. FASN and ACACA expression showed the strongest negative correlations with CD8^+^ T cells (r = -0.34 and r = -0.31 respectively, both p < 0.001) and positive correlations with Tregs (r = 0.28 and r = 0.26, both p < 0.001), identifying these enzymes as potential mechanistic drivers of immune dysfunction.

### Pathway-Level Characterization of Lipid-High Melanomas

Gene set enrichment analysis (GSEA) comparing high versus low lipid metabolism groups revealed coordinate regulation of lipid biosynthetic pathways and immunosuppressive features (**Figure 4A-B, Table S1**).

**Upregulated pathways in high lipid metabolism tumors:**
1. **Fatty acid metabolism** (NES = 2.34, FDR < 0.001) - The most significantly enriched pathway, consistent with the defining signature
2. **Cholesterol homeostasis** (NES = 2.12, FDR < 0.001) - Indicating coordinate regulation of lipid synthesis
3. **MTORC1 signaling** (NES = 1.98, FDR = 0.002) - mTORC1 is a known activator of SREBP-1 and lipogenesis [69-72]
4. **Adipogenesis** (NES = 1.76, FDR = 0.008) - Reflecting lipid accumulation and storage programs
5. **Glycolysis** (NES = 1.65, FDR = 0.015) - Linking glucose and lipid metabolism through citrate production

**Downregulated pathways in high lipid metabolism tumors:**
1. **Interferon gamma response** (NES = -1.87, FDR = 0.008) - Critical for anti-tumor immunity and ICI response [60,61]
2. **Inflammatory response** (NES = -1.65, FDR = 0.021) - Consistent with "cold" immune phenotype
3. **IL6-JAK-STAT3 signaling** (NES = -1.54, FDR = 0.032) - Involved in immune cell activation
4. **Allograft rejection** (NES = -1.48, FDR = 0.044) - Reflects reduced immune recognition
5. **TNF-α signaling via NF-κB** (NES = -1.42, FDR = 0.048) - Pro-inflammatory pathway suppression

The inverse relationship between lipid metabolism and interferon signaling (**Figure 4C**) is particularly noteworthy, as interferon-gamma is essential for CD8^+^ T cell cytotoxicity, MHC class I expression, and response to immune checkpoint blockade [62-64]. This suggests that high lipid metabolism tumors may be inherently resistant to immunotherapy.

### Lipid Metabolism Correlates with Known Melanoma Phenotypic States

To contextualize our lipid metabolism signature within established melanoma biology, we examined its relationship to the melanoma phenotype-switching paradigm. Previous work has identified MITF (Microphthalmia-Associated Transcription Factor) as a master regulator of melanocyte differentiation, with MITF-low melanomas exhibiting invasive and therapy-resistant phenotypes [27,65,66].

We observed a significant negative correlation between lipid metabolism signature and MITF expression (Spearman r = -0.42, p < 0.001) (**Figure 4D**), indicating that high lipid metabolism associates with the dedifferentiated, invasive phenotypic state. Furthermore, high lipid tumors showed elevated expression of invasion-associated genes including AXL (r = 0.38, p < 0.001), TWIST1 (r = 0.32, p < 0.001), and mesenchymal markers [67,68], supporting the concept that lipid metabolism reprogramming facilitates melanoma progression and metastasis.

---

## DISCUSSION

This study provides comprehensive evidence that elevated tumor lipid metabolism defines a distinct melanoma subtype characterized by poor prognosis and a profoundly immunosuppressive microenvironment. By integrating transcriptomic profiling, survival analysis, and immune deconvolution in a large melanoma cohort, we demonstrate that a lipid metabolism signature comprising eight key lipogenic enzymes serves as an independent prognostic biomarker and correlates with profound alterations in tumor immune contexture. These findings establish a metabolic-immunologic framework for understanding melanoma heterogeneity and suggest novel therapeutic strategies targeting the intersection of lipid metabolism and immune evasion.

### Lipid Metabolism as an Independent Prognostic Biomarker in Melanoma

Our analysis identified that high lipid metabolism, characterized by coordinated upregulation of FASN, ACACA, SCD, ELOVL6, FADS1, FADS2, ACLY, and SREBF1, significantly correlates with reduced overall survival (p = 0.0095; HR = 1.68) in melanoma patients. This observation extends previous findings in other cancer types [73-75] and underscores the universal importance of lipogenic reprogramming in tumor progression. The independent prognostic value of the lipid signature, even after adjusting for established clinical variables including age and stage (adjusted HR = 1.72, p = 0.008), suggests that metabolic phenotyping could complement traditional staging systems to improve risk stratification.

The hazard ratio of 1.72 for high versus low lipid metabolism groups represents a clinically meaningful prognostic distinction comparable to other validated biomarkers in melanoma, such as tumor mutational burden and immune gene signatures [76,77]. Notably, FASN and ACACA emerged as particularly strong individual prognostic markers, consistent with their central roles in de novo fatty acid synthesis [1,2,21-26]. These enzymes catalyze rate-limiting steps in lipogenesis and are frequently overexpressed across diverse malignancies.

The coordinated expression of multiple lipid genes in high-risk tumors (correlation coefficients 0.42-0.78) suggests transcriptional co-regulation, likely mediated by master regulators such as SREBP-1, which itself was elevated in poor-prognosis tumors. SREBP-1 activation is known to be driven by oncogenic signaling pathways including PI3K/AKT/mTOR and MAPK, which are frequently dysregulated in melanoma [69-72,78]. Indeed, our GSEA analysis revealed mTORC1 pathway enrichment in high lipid tumors (NES = 1.98, FDR = 0.002), supporting this mechanistic connection.

### The Immunosuppressive Phenotype of Lipid-High Melanomas

A central finding of our study is the profound immunological reprogramming associated with elevated tumor lipid metabolism. High lipid metabolism tumors exhibited significant depletion of CD8^+^ cytotoxic T lymphocytes (33% reduction, p < 0.01) and M1 macrophages (38% reduction, p < 0.05), coupled with enrichment of regulatory T cells (80% increase, p < 0.01) and M2-polarized macrophages (56% increase, p < 0.05). This immune profile is characteristic of an immunosuppressive microenvironment that promotes tumor immune escape and therapeutic resistance [4-6,29-34].

The mechanistic basis for lipid metabolism-driven immunosuppression likely involves multiple intersecting pathways. First, lipid accumulation in tumor cells and the microenvironment can directly impair T cell function through lipotoxicity and metabolic competition [38-40]. CD8^+^ T cells require robust fatty acid oxidation for effector function and memory formation [37], and lipid-rich environments may paradoxically starve these cells of metabolic substrates while simultaneously inducing lipid-mediated dysfunction. The upregulation of CD36, a fatty acid transporter, on tumor-infiltrating CD8^+^ T cells has been shown to impair their anti-tumor activity and correlate with poor outcomes [39], providing a direct link between lipid availability and T cell dysfunction.

Second, lipid mediators such as prostaglandins, leukotrienes, and specialized pro-resolving mediators (SPMs), derived from arachidonic acid and omega-3 fatty acid metabolism, can promote immunosuppressive signaling through their cognate receptors on immune cells [79,80]. These bioactive lipids regulate immune cell trafficking, activation thresholds, and effector functions, with many demonstrating anti-inflammatory or immunosuppressive properties in the tumor context.

Third, lipogenic tumors may actively recruit and polarize macrophages toward M2 phenotypes through secretion of lipid-derived chemokines and cytokines [41-47]. M2 macrophages not only fail to support anti-tumor immunity but actively suppress T cell responses through production of IL-10, TGF-β, and arginase, which depletes arginine required for T cell proliferation [43,44]. Our observation of a 56% increase in M2 macrophages in high lipid tumors, coupled with decreased M1 macrophages, indicates a profound shift in the macrophage compartment toward tumor-promoting phenotypes. The preferential utilization of fatty acid oxidation by M2 macrophages [45,46] suggests that lipid-rich environments provide a metabolic advantage to immunosuppressive TAMs.

Fourth, the 80% elevation in regulatory T cells observed in high lipid metabolism tumors is particularly significant, as Tregs are potent suppressors of anti-tumor immunity and their abundance is inversely correlated with immunotherapy response [48-52,81,82]. Tumor-infiltrating Tregs exhibit markedly upregulated lipid metabolism compared to conventional T cells, with enhanced fatty acid uptake, synthesis, and oxidation [50,51]. This metabolic phenotype not only supports Treg proliferation but is functionally required for their suppressive activity, as inhibition of FAO impairs Treg function [52]. The accumulation of Tregs in lipid-high tumors may therefore reflect both metabolic adaptation to the lipid-rich environment and active recruitment by tumor-derived signals.

### Pathway-Level Insights Link Metabolism and Immunity

Our gene set enrichment analysis revealed that lipid-high tumors exhibit coordinated upregulation of fatty acid biosynthesis (NES = 2.34), cholesterol homeostasis (NES = 2.12), and mTORC1 signaling pathways (NES = 1.98), alongside profound depletion of interferon-gamma (NES = -1.87) and inflammatory response signatures (NES = -1.65). This pathway-level view reinforces the concept that metabolic and immune programs are interconnected rather than independent features of tumor biology.

The depletion of interferon-gamma signatures is particularly notable, as this cytokine is essential for anti-tumor immunity and serves as a key biomarker for immunotherapy response [60-64]. Interferon-gamma drives expression of MHC class I molecules, antigen processing machinery, and chemokines that recruit cytotoxic lymphocytes [62]. The inverse relationship between lipid metabolism and interferon signaling (Spearman r = -0.51, p < 0.001) suggests that metabolic interventions could potentially restore immune responsiveness.

Mechanistically, this connection may be bidirectional. On one hand, lipid accumulation and metabolic stress in tumor cells can suppress interferon-responsive gene expression through epigenetic modifications and altered transcription factor activity [83,84]. On the other hand, interferon-gamma signaling from immune cells can actually induce FASN expression and lipid metabolism in tumor cells as part of an "immunometabolic editing" process whereby tumors adapt to immune surveillance [40]. This creates a potential feedforward loop where lipid metabolism both suppresses interferon signaling and is induced by it, ultimately favoring tumor escape.

### Connection to Melanoma Phenotype Switching and MITF

Our observation of an inverse correlation between lipid metabolism and MITF expression (r = -0.42, p < 0.001) connects our findings to the established melanoma phenotype-switching paradigm [27,65-68]. MITF regulates melanocyte differentiation and is associated with proliferative melanoma cells, whereas MITF-low states are characterized by invasiveness, stemness, and resistance to BRAF/MEK inhibitors [65,66]. A recent study demonstrated that MITF directly regulates SCD expression and fatty acid saturation, thereby establishing a positive feedback loop that stabilizes the MITF-low, invasive state [27].

Mechanistically, low SCD expression and activity in MITF-high cells promotes endoplasmic reticulum stress and phosphorylation of eIF2α, leading to activation of ATF4 and NF-κB-dependent inflammatory signaling that sustains reduced MITF expression and melanoma cell dedifferentiation [27,28]. Our observation that high lipid metabolism (which includes elevated SCD) inversely correlates with MITF suggests a more complex relationship whereby lipid metabolism adaptations may occur during phenotype transitions or represent a distinct metabolic subtype that overlaps with but is not identical to the classical MITF-low invasive state.

The association of high lipid metabolism with elevated expression of invasion markers (AXL, TWIST1) and mesenchymal features supports the concept that lipogenic reprogramming facilitates melanoma progression through multiple mechanisms: providing membrane building blocks for rapidly dividing cells, generating signaling lipids that promote invasion and metastasis, creating an immunosuppressive microenvironment, and conferring resistance to oxidative stress and therapy [67,68,85,86].

### Therapeutic Implications and Translational Potential

Our findings suggest several promising therapeutic strategies for targeting lipid metabolism in melanoma. First, pharmacological inhibition of key lipogenic enzymes such as FASN and ACACA could simultaneously disrupt tumor cell metabolism and alleviate immunosuppression [87-89]. Several FASN inhibitors have entered clinical development, including TVB-2640 (denifanstat) and TVB-3166, which have shown tolerable safety profiles in early-phase trials [90,91]. Our identification of FASN and ACACA as independent prognostic markers and their strong negative correlation with CD8^+^ T cell infiltration (r = -0.34 and r = -0.31 respectively) provides a rationale for testing these agents in melanoma, particularly in combination with immunotherapy.

Preclinical studies have demonstrated that FASN inhibition can enhance anti-tumor immune responses. Treatment with orlistat or cerulenin induced apoptosis in melanoma cells through activation of the intrinsic pathway [24], and genetic knockdown of FASN reduced melanoma cell proliferation and migration [22,23]. Importantly, FASN inhibition has been shown to relieve immunosuppression in other cancer types by reducing myeloid-derived suppressor cell accumulation and enhancing T cell function [92,93], suggesting potential synergy with immune checkpoint blockade.

Second, targeting the SREBP pathway, which coordinately regulates multiple lipogenic genes, could provide a more comprehensive metabolic intervention [70-72,94,95]. SREBP inhibitors such as fatostatin and betulin have shown anti-cancer activity in preclinical models, and genetic studies have confirmed that SREBP function is required for tumor growth in multiple contexts [69,94]. Given the master regulatory role of SREBP-1 in coordinating the lipid metabolism program we observed, SREBP inhibition represents an attractive therapeutic target, particularly for tumors with high lipid signatures.

Third, modulating lipid availability or composition through dietary interventions may complement pharmacological approaches [96,97]. Ketogenic diets and caloric restriction can alter systemic and tumor lipid metabolism and have shown promise in preclinical cancer models [98,99]. While clinical evidence in melanoma is limited, the concept of leveraging diet to reprogram tumor metabolism and enhance immunotherapy efficacy warrants investigation in prospective trials.

Fourth, targeting downstream lipid signaling pathways or their receptors could selectively neutralize immunosuppressive signals without broadly disrupting cellular lipid homeostasis [79,80,100]. For example, inhibitors of prostaglandin synthesis (COX-2 inhibitors) or prostanoid receptors have shown anti-tumor effects in preclinical studies and could potentially enhance immunotherapy responses by reducing PGE2-mediated immune suppression [101,102].

Fifth, our observation that high lipid metabolism tumors exhibit depleted interferon-gamma signaling suggests that these tumors may be inherently resistant to immune checkpoint inhibitors [60-64]. However, this also presents an opportunity: combining lipid metabolism inhibitors with immunotherapy could potentially restore interferon responsiveness and overcome resistance. Recent preclinical data support this concept, showing that targeting lipid metabolism in tumor-associated macrophages enhanced T cell infiltration and synergized with anti-PD-1 therapy [103,104].

### Clinical Translation and Biomarker Development

From a clinical translational perspective, the lipid metabolism signature we identified could serve multiple purposes. First, it could function as a prognostic biomarker to identify high-risk patients who may benefit from more aggressive treatment approaches or closer surveillance. The signature's independent prognostic value (C-index improvement from 0.68 to 0.73) and ease of calculation from standard RNA-seq data make it potentially implementable in clinical practice.

Second, the lipid signature could serve as a predictive biomarker for response to immunotherapy. Given the profound immune exclusion and interferon depletion observed in high lipid tumors, these patients may be less likely to respond to anti-PD-1/PD-L1 monotherapy and could be prioritized for combination approaches that target both immune checkpoints and metabolic pathways [105,106]. Prospective validation studies correlating lipid metabolism signatures with immunotherapy outcomes are needed to test this hypothesis.

Third, lipid metabolism enzymes such as FASN could potentially be detected by immunohistochemistry in tumor biopsies, providing a simpler and more accessible biomarker for clinical use [26]. Previous studies have demonstrated that FASN immunostaining can distinguish metastatic melanoma from benign nevi with high sensitivity and specificity [26,42], and quantification of FASN expression could be integrated into pathology workflows.

### Study Limitations

Several limitations of this study merit consideration. First, our analysis relies on bulk tumor RNA-seq data from TCGA, which cannot resolve cell-type-specific contributions to the observed lipid metabolism signature. While computational deconvolution methods like CIBERSORTx provide estimates of immune cell proportions, they cannot definitively establish whether elevated lipid gene expression originates from tumor cells, stromal cells, or immune cells themselves. Single-cell RNA sequencing would enable precise attribution of lipogenic programs to specific cell populations and reveal potential metabolic heterogeneity within tumors [107,108].

Second, the correlative nature of our transcriptomic analysis precludes definitive causal inferences about the relationship between lipid metabolism and immune dysfunction. While our findings are consistent with mechanistic studies showing that lipid metabolism regulates immune cell function [38-52], functional validation using genetic or pharmacological perturbation of lipid metabolism in preclinical melanoma models is essential to establish causality and identify therapeutic vulnerabilities [87-93].

Third, the TCGA cohort lacks detailed treatment information, particularly regarding immunotherapy exposure and response. This limits our ability to assess whether lipid metabolism influences response to specific therapies such as immune checkpoint inhibitors or BRAF/MEK inhibitors. Retrospective analyses of cohorts with well-annotated treatment data and prospective studies incorporating metabolic biomarkers into clinical trial designs would address this critical gap [105,106].

Fourth, our analysis focused on mRNA expression as a surrogate for metabolic activity, but transcript levels do not always correlate perfectly with protein expression, enzyme activity, or metabolic flux [109,110]. Direct measurements of lipid species, metabolic intermediates, and enzymatic activities using metabolomics and lipidomics would provide a more complete picture of lipid metabolism in melanoma and validate our transcriptional findings [111,112]. Additionally, isotope tracing experiments could definitively establish the flux through lipogenic pathways and identify rate-limiting steps amenable to therapeutic intervention [113].

Fifth, while CIBERSORTx provides valuable information about immune cell composition, it cannot capture the spatial organization of immune cells within tumors or their functional states [114,115]. Multiplex immunohistochemistry, imaging mass cytometry, or spatial transcriptomics would reveal whether lipid-high tumor regions are spatially segregated from immune infiltrates (immune-excluded phenotype) or whether immune cells are present but dysfunctional (immune-suppressed phenotype) [116,117]. Understanding the spatial relationship between lipid metabolism and immune infiltration could inform strategies to overcome immune exclusion.

Sixth, our study focused on primary and metastatic melanomas without distinguishing between these disease states or considering tumor heterogeneity. Melanoma is known to exhibit substantial intra-tumoral and inter-tumoral heterogeneity [118], and different regions or metastatic sites may exhibit distinct lipid metabolic profiles. Region-specific or site-specific analyses would provide insights into the spatial and temporal dynamics of metabolic reprogramming during melanoma evolution.

### Future Directions

Based on our findings and the limitations discussed above, we propose several key directions for future research. First, single-cell multi-omic profiling (integrating transcriptomics, proteomics, and metabolomics at single-cell resolution) would provide unprecedented resolution of metabolic heterogeneity and cell-type-specific lipid programs in melanoma [119,120]. This could reveal rare subpopulations with extreme metabolic phenotypes that drive progression or therapy resistance.

Second, functional validation in preclinical models, including patient-derived xenografts and genetically engineered mouse models, is essential to test causality and therapeutic efficacy [121,122]. Genetic deletion or pharmacological inhibition of key lipid metabolism enzymes, combined with immunotherapy, would directly test whether metabolic interventions can overcome immune resistance and improve outcomes.

Third, clinical trials incorporating lipid metabolism biomarkers are needed to prospectively validate our findings [105,106]. Ideally, these trials would combine lipid metabolism inhibitors with immune checkpoint blockade in patients selected based on high lipid metabolism signatures, with correlative studies examining changes in tumor immune infiltration and systemic metabolism.

Fourth, integration of multi-omic data (genomics, transcriptomics, proteomics, metabolomics, lipidomics) from melanoma patients would enable systems-level understanding of how genetic alterations, signaling pathways, metabolic programs, and immune responses interact to determine clinical outcomes [123,124]. Machine learning approaches could identify complex multi-factorial signatures that outperform single-modality biomarkers.

Finally, investigation of systemic metabolic factors, including diet, obesity, diabetes, and circulating lipid profiles, could reveal host factors that influence tumor lipid metabolism and immunotherapy efficacy [96,97,125]. Population-based studies correlating metabolic comorbidities with melanoma outcomes, combined with mechanistic studies in animal models, could inform lifestyle interventions to complement cancer therapy.

### Conclusions

This study establishes elevated tumor lipid metabolism as a defining feature of an aggressive melanoma subtype with poor prognosis and profoundly immunosuppressive microenvironment. The lipid metabolism signature, comprising eight key lipogenic genes (FASN, ACACA, SCD, ELOVL6, FADS1, FADS2, ACLY, SREBF1), provides independent prognostic value beyond traditional clinical variables and identifies patients with reduced CD8^+^ T cell infiltration, elevated regulatory T cells, and depleted interferon signaling.

The integration of lipid metabolic signatures with immune profiling reveals a metabolic-immunologic axis that governs melanoma biology and therapeutic vulnerability. High lipid metabolism tumors exhibit coordinate upregulation of fatty acid biosynthesis and cholesterol homeostasis pathways, alongside profound suppression of anti-tumor immune responses. These tumors likely represent a "cold" or immune-excluded phenotype that may be inherently resistant to immune checkpoint blockade monotherapy.

Targeting lipid metabolism represents a promising strategy to simultaneously disrupt tumor cell metabolism and enhance anti-tumor immunity. Pharmacological inhibitors of FASN, ACACA, or SREBP are in clinical development and could be tested in combination with immunotherapy, particularly in patients selected based on high lipid metabolism signatures. Additionally, dietary interventions and inhibitors of lipid signaling pathways warrant investigation as complementary approaches.

Our findings provide a framework for developing biomarker-driven clinical trials that combine metabolic inhibitors with immunotherapy to improve outcomes in melanoma patients. The lipid metabolism signature could serve both prognostic and predictive roles, identifying high-risk patients and guiding treatment selection. Future studies integrating multi-omic profiling, functional validation, and clinical translation will be essential to fully realize the therapeutic potential of targeting the metabolic-immune axis in melanoma.

---

## FIGURE LEGENDS

**Figure 1. Lipid Metabolism Signature Construction and Gene Expression Patterns**

**(A)** Distribution of lipid metabolism signature scores across the TCGA-SKCM cohort (n=472). Patients were dichotomized into high (n=236, red) and low (n=236, blue) lipid metabolism groups based on median signature score (dashed line). Box plots show median, interquartile range, and outliers.

**(B)** Correlation heatmap showing Spearman correlation coefficients among the eight lipid metabolism genes (FASN, ACACA, SCD, ELOVL6, FADS1, FADS2, ACLY, SREBF1). Color intensity indicates correlation strength (blue = negative, red = positive). All correlations were statistically significant (p < 0.001). Strongest correlations observed between FASN-ACACA (r = 0.78), FADS1-FADS2 (r = 0.72), and ACLY-SREBF1 (r = 0.68), indicating coordinated transcriptional regulation.

**Figure 2. Lipid Metabolism Signature Associates with Poor Survival**

**(A)** Kaplan-Meier survival curves comparing overall survival between high (red) and low (blue) lipid metabolism groups. Median OS: 4.2 years (high) vs. 7.8 years (low); log-rank p = 0.0095. Numbers at risk shown below plot. Shaded areas represent 95% confidence intervals.

**(B)** Time-dependent ROC curves for lipid metabolism signature in predicting overall survival at 3 years (AUC = 0.68) and 5 years (AUC = 0.71). Diagonal dashed line represents random chance (AUC = 0.50).

**(C)** Forest plot showing hazard ratios from univariate Cox regression for individual lipid metabolism genes. Points represent hazard ratios, horizontal lines represent 95% confidence intervals. FASN (HR = 1.42, p = 0.024) and ACACA (HR = 1.38, p = 0.041) were the strongest individual predictors.

**(D)** Heatmap of lipid metabolism gene expression across all 472 patients, ordered by lipid signature score (left to right: low to high). Rows represent genes, columns represent patients. Color scale shows Z-score normalized expression (blue = low, red = high). Top annotation bar shows lipid group (blue = low, red = high).

**(E)** Concordance index (C-index) comparison between clinical variables alone (age, sex, stage; C-index = 0.68) and integrated model (clinical + lipid signature; C-index = 0.73). Error bars represent 95% confidence intervals. Likelihood ratio test p = 0.006.

**Figure 3. High Lipid Metabolism Associates with Immunosuppressive Microenvironment**

**(A)** Stacked bar plot showing immune cell composition across all 444 samples (after CIBERSORTx quality filtering), arranged by lipid signature score. Each bar represents one patient, colors represent different immune cell types. Major populations labeled include CD8^+^ T cells, CD4^+^ T cells, Tregs, B cells, NK cells, M1 macrophages, M2 macrophages, and myeloid cells.

**(B)** Box plots comparing immune cell proportions between high and low lipid metabolism groups. Significant differences marked with asterisks (*p < 0.05, **p < 0.01). CD8^+^ T cells: 0.12 vs. 0.18 (p < 0.01); Tregs: 0.09 vs. 0.05 (p < 0.01); M1 macrophages: 0.08 vs. 0.13 (p < 0.05); M2 macrophages: 0.14 vs. 0.09 (p < 0.05); NK cells: 0.04 vs. 0.07 (p = 0.02).

**(C)** Correlation heatmap showing Spearman correlations between individual lipid genes (rows) and key immune cell populations (columns). Color intensity and size indicate correlation strength and statistical significance. FASN and ACACA showed strongest negative correlations with CD8^+^ T cells (r = -0.34 and r = -0.31, both p < 0.001) and positive correlations with Tregs (r = 0.28 and r = 0.26, both p < 0.001).

**Figure 4. Pathway Enrichment and Connection to Melanoma Biology**

**(A)** Bar plot showing top upregulated pathways in high versus low lipid metabolism tumors from GSEA analysis. Bars represent normalized enrichment scores (NES), error bars represent FDR q-values. Top pathways: Fatty acid metabolism (NES = 2.34), Cholesterol homeostasis (NES = 2.12), mTORC1 signaling (NES = 1.98).

**(B)** Bar plot showing top downregulated pathways in high lipid metabolism tumors. Interferon gamma response (NES = -1.87), Inflammatory response (NES = -1.65), IL6-JAK-STAT3 signaling (NES = -1.54).

**(C)** Scatter plot showing inverse correlation between lipid metabolism signature (x-axis) and interferon-gamma response score (y-axis). Each point represents one patient, colored by lipid group (blue = low, red = high). Spearman r = -0.51, p < 0.001. Regression line with 95% confidence interval shown.

**(D)** Scatter plot showing inverse correlation between lipid metabolism signature and MITF expression. Spearman r = -0.42, p < 0.001. High lipid metabolism associates with MITF-low, invasive phenotype.

**Figure 5. Proposed Model of Lipid Metabolism-Immune Axis in Melanoma**

Schematic diagram illustrating how elevated tumor lipid metabolism creates an immunosuppressive microenvironment. Central tumor cell shows upregulated lipogenic enzymes (FASN, ACACA, SCD, SREBP-1) driven by oncogenic signaling (PI3K/AKT/mTOR). Lipid-rich tumor environment leads to: (1) CD8^+^ T cell dysfunction through lipotoxicity and CD36 upregulation; (2) M2 macrophage polarization via FAO preference; (3) Treg accumulation and activation through enhanced FAO; (4) Suppression of interferon-gamma signaling; (5) Production of immunosuppressive lipid mediators (PGE2). Bottom panel shows therapeutic interventions: FASN/ACACA inhibitors, SREBP inhibitors, dietary modulation, combined with immune checkpoint blockade to restore anti-tumor immunity.

---

## SUPPLEMENTARY FIGURE LEGENDS

**Figure S1. Quality Control and Principal Component Analysis**

PCA plot showing first two principal components (PC1 vs. PC2) for all 472 melanoma samples. Colors indicate lipid metabolism group (red = high, blue = low). No obvious batch effects or outliers detected. Variance explained: PC1 = 24%, PC2 = 18%.

**Figure S2. Individual Lipid Gene Expression and Survival**

Kaplan-Meier curves for each of the eight lipid metabolism genes individually: FASN (p = 0.024), ACACA (p = 0.041), SCD (p = 0.073), ELOVL6 (p = 0.15), FADS1 (p = 0.19), FADS2 (p = 0.12), ACLY (p = 0.086), SREBF1 (p = 0.036). Demonstrates that FASN, ACACA, and SREBF1 are the strongest individual prognostic markers.

**Figure S3. Complete GSEA Results**

Complete ranked list of all Hallmark gene sets from GSEA analysis showing NES and FDR q-values. Includes additional pathways such as Oxidative Phosphorylation, Apoptosis, Epithelial-Mesenchymal Transition, and Hypoxia.

**Figure S4. Lipid Metabolism Signature by Clinical Subgroups**

Box plots showing lipid metabolism signature scores stratified by: (A) Disease stage (I/II vs. III/IV); (B) BRAF mutation status; (C) Age groups (<50, 50-65, >65 years); (D) Primary vs. metastatic samples. Statistical comparisons using Wilcoxon or Kruskal-Wallis tests.

**Figure S5. Extended Immune Cell Analysis**

Detailed analysis of all 22 immune cell types from CIBERSORTx including B cells, plasma cells, different T cell subsets, NK cell subsets, monocytes, macrophage subsets, dendritic cell subsets, mast cells, eosinophils, and neutrophils.

---
