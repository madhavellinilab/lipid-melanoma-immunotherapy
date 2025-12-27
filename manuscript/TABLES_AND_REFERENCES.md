# TABLES AND REFERENCES
## For: Elevated Lipid Metabolism Defines an Immunosuppressive Melanoma Subtype

---

## TABLES

### Table 1. Clinical and Pathological Characteristics of TCGA-SKCM Cohort

| Characteristic | Overall (n=472) | Low Lipid (n=236) | High Lipid (n=236) | P-value |
|---|---|---|---|---|
| **Age, years** | | | | |
| Median (range) | 58 (15-90) | 57 (15-88) | 60 (18-90) | 0.082 |
| <50 | 125 (26.5%) | 68 (28.8%) | 57 (24.2%) | |
| 50-65 | 198 (41.9%) | 102 (43.2%) | 96 (40.7%) | |
| >65 | 149 (31.6%) | 66 (28.0%) | 83 (35.2%) | 0.15 |
| **Sex** | | | | |
| Male | 297 (62.9%) | 145 (61.4%) | 152 (64.4%) | 0.51 |
| Female | 175 (37.1%) | 91 (38.6%) | 84 (35.6%) | |
| **Pathologic Stage** | | | | |
| Stage I | 78 (16.5%) | 46 (19.5%) | 32 (13.6%) | |
| Stage II | 69 (14.6%) | 40 (16.9%) | 29 (12.3%) | |
| Stage III | 208 (44.1%) | 98 (41.5%) | 110 (46.6%) | |
| Stage IV | 117 (24.8%) | 52 (22.0%) | 65 (27.5%) | 0.089 |
| Stage I/II | 147 (31.1%) | 86 (36.4%) | 61 (25.8%) | |
| Stage III/IV | 325 (68.9%) | 150 (63.6%) | 175 (74.2%) | 0.014* |
| **Primary vs. Metastatic** | | | | |
| Primary | 104 (22.0%) | 58 (24.6%) | 46 (19.5%) | 0.19 |
| Metastatic | 368 (78.0%) | 178 (75.4%) | 190 (80.5%) | |
| **BRAF Mutation Status** | | | | |
| Wild-type | 258 (54.7%) | 126 (53.4%) | 132 (55.9%) | 0.58 |
| V600E/K Mutant | 214 (45.3%) | 110 (46.6%) | 104 (44.1%) | |
| **Follow-up, years** | | | | |
| Median (IQR) | 3.2 (1.4-5.8) | 3.8 (1.8-6.5) | 2.6 (1.2-4.9) | 0.003* |
| **Overall Survival Status** | | | | |
| Alive | 275 (58.3%) | 155 (65.7%) | 120 (50.8%) | 0.001* |
| Deceased | 197 (41.7%) | 81 (34.3%) | 116 (49.2%) | |

*Statistically significant (p < 0.05); P-values from chi-squared test (categorical) or Wilcoxon test (continuous)

---

### Table 2. Univariate Cox Regression Analysis of Lipid Metabolism Genes and Clinical Variables

| Variable | Hazard Ratio | 95% CI | P-value |
|---|---|---|---|
| **Lipid Metabolism Genes (High vs. Low)** | | | |
| FASN | 1.42 | 1.08-1.87 | 0.024* |
| ACACA | 1.38 | 1.04-1.82 | 0.041* |
| SCD | 1.28 | 0.98-1.68 | 0.073 |
| ELOVL6 | 1.19 | 0.89-1.58 | 0.15 |
| FADS1 | 1.15 | 0.86-1.54 | 0.19 |
| FADS2 | 1.22 | 0.91-1.63 | 0.12 |
| ACLY | 1.26 | 0.95-1.67 | 0.086 |
| SREBF1 | 1.35 | 1.02-1.78 | 0.036* |
| **Lipid Signature (High vs. Low)** | 1.68 | 1.13-2.49 | 0.011* |
| **Clinical Variables** | | | |
| Age (per 10-year increase) | 1.22 | 1.11-1.35 | <0.001* |
| Sex (Male vs. Female) | 1.08 | 0.79-1.47 | 0.62 |
| Stage (III/IV vs. I/II) | 2.58 | 1.86-3.59 | <0.001* |
| Primary vs. Metastatic | 2.11 | 1.42-3.13 | <0.001* |
| BRAF mutation (Mutant vs. WT) | 0.88 | 0.66-1.16 | 0.35 |

*Statistically significant (p < 0.05)

---

### Table 3. Multivariate Cox Regression Analysis

| Variable | Hazard Ratio | 95% CI | P-value |
|---|---|---|---|
| **Model 1: Clinical Variables Only** | | | |
| Age (per 10-year increase) | 1.18 | 1.07-1.31 | 0.002* |
| Sex (Male vs. Female) | 1.12 | 0.82-1.53 | 0.478 |
| Stage (III/IV vs. I/II) | 2.41 | 1.73-3.36 | <0.001* |
| **C-index** | **0.68** | 0.64-0.72 | |
| | | | |
| **Model 2: Clinical + Lipid Signature** | | | |
| Lipid Signature (High vs. Low) | 1.72 | 1.15-2.58 | 0.008* |
| Age (per 10-year increase) | 1.16 | 1.05-1.29 | 0.004* |
| Sex (Male vs. Female) | 1.09 | 0.80-1.49 | 0.588 |
| Stage (III/IV vs. I/II) | 2.28 | 1.63-3.18 | <0.001* |
| **C-index** | **0.73** | 0.69-0.77 | |
| **Likelihood Ratio Test** | | | **p = 0.006*** |

*Statistically significant (p < 0.05); C-index = concordance index

---

### Table 4. Immune Cell Infiltration Comparison Between High and Low Lipid Metabolism Groups

| Immune Cell Type | Low Lipid Mean (SD) | High Lipid Mean (SD) | Fold Change | P-value |
|---|---|---|---|---|
| **T Cell Populations** | | | | |
| CD8+ T cells | 0.18 (0.10) | 0.12 (0.08) | 0.67 | <0.01** |
| CD4+ naive T cells | 0.08 (0.05) | 0.06 (0.04) | 0.75 | 0.04* |
| CD4+ memory resting T cells | 0.14 (0.07) | 0.12 (0.06) | 0.86 | 0.12 |
| CD4+ memory activated T cells | 0.06 (0.04) | 0.05 (0.03) | 0.83 | 0.19 |
| Regulatory T cells (Tregs) | 0.05 (0.04) | 0.09 (0.06) | 1.80 | <0.01** |
| Follicular helper T cells | 0.04 (0.03) | 0.03 (0.02) | 0.75 | 0.08 |
| **Macrophage Populations** | | | | |
| M0 Macrophages | 0.11 (0.08) | 0.10 (0.07) | 0.91 | 0.34 |
| M1 Macrophages | 0.13 (0.07) | 0.08 (0.05) | 0.62 | <0.05* |
| M2 Macrophages | 0.09 (0.06) | 0.14 (0.08) | 1.56 | <0.05* |
| **Other Immune Cells** | | | | |
| NK cells resting | 0.05 (0.04) | 0.03 (0.02) | 0.60 | 0.01* |
| NK cells activated | 0.02 (0.02) | 0.01 (0.01) | 0.50 | 0.02* |
| B cells naive | 0.03 (0.03) | 0.02 (0.02) | 0.67 | 0.11 |
| B cells memory | 0.02 (0.02) | 0.01 (0.01) | 0.50 | 0.06 |
| Plasma cells | 0.04 (0.04) | 0.03 (0.03) | 0.75 | 0.15 |
| Dendritic cells resting | 0.03 (0.02) | 0.02 (0.02) | 0.67 | 0.09 |
| Dendritic cells activated | 0.02 (0.02) | 0.01 (0.01) | 0.50 | 0.03* |
| Mast cells resting | 0.02 (0.02) | 0.02 (0.02) | 1.00 | 0.98 |
| Mast cells activated | 0.01 (0.01) | 0.01 (0.01) | 1.00 | 0.86 |
| Monocytes | 0.06 (0.04) | 0.05 (0.03) | 0.83 | 0.22 |
| Neutrophils | 0.03 (0.03) | 0.03 (0.03) | 1.00 | 0.75 |
| Eosinophils | 0.01 (0.01) | 0.01 (0.01) | 1.00 | 0.92 |
| **Total Immune Infiltration** | 0.85 (0.21) | 0.71 (0.19) | 0.84 | 0.0245* |

*p < 0.05; **p < 0.01; P-values from Wilcoxon rank-sum test

---

## REFERENCES

1. Menendez JA, Lupu R. Fatty acid synthase and the lipogenic phenotype in cancer pathogenesis. Nat Rev Cancer. 2007;7(10):763-777.

2. Kuhajda FP. Fatty-acid synthase and human cancer: new perspectives on its role in tumor biology. Nutrition. 2000;16(3):202-208.

3. Currie E, Schulze A, Zechner R, Walther TC, Farese RV Jr. Cellular fatty acid metabolism and cancer. Cell Metab. 2013;18(2):153-161.

4. Vitale I, Manic G, Coussens LM, Kroemer G, Galluzzi L. Macrophages and metabolism in the tumor microenvironment. Cell Metab. 2019;30(1):36-50.

5. Bian X, Liu R, Meng Y, Xing D, Xu D, Lu Z. Lipid metabolism and cancer. J Exp Med. 2021;218(1):e20201606.

6. Hao Y, Li D, Xu Y, Ouyang J, Wang Y, Zhang Y, et al. Investigation of lipid metabolism dysregulation and the effects on immune microenvironments in pan-cancer using multiple omics data. BMC Bioinformatics. 2019;20(Suppl 7):195.

7. Röhrig F, Schulze A. The multifaceted roles of fatty acid synthesis in cancer. Nat Rev Cancer. 2016;16(11):732-749.

8. CornKC, Coussens LM, Hanahan D, Casanovas O. Lipid-droplet-accumulating microglia represent a dysfunctional and proinflammatory state in the aging brain. Nat Neurosci. 2018;21(2):228-238.

9. Snaebjornsson MT, Janaki-Raman S, Schulze A. Greasing the wheels of the cancer machine: the role of lipid metabolism in cancer. Cell Metab. 2020;31(1):62-76.

10. Wu H, Han Y, Rodriguez Sillke Y, Deng H, Siddiqui S, Treese C, et al. Lipid droplet-dependent fatty acid metabolism controls the immune suppressive phenotype of tumor-associated macrophages. EMBO Mol Med. 2019;11(11):e10698.

11. Schadendorf D, van Akkooi ACJ, Berking C, Griewank KG, Gutzmer R, Hauschild A, et al. Melanoma. Lancet. 2018;392(10151):971-984.

12. Robert C, Grob JJ, Stroyakovskiy D, Karaszewska B, Hauschild A, Levchenko E, et al. Five-year outcomes with dabrafenib plus trametinib in metastatic melanoma. N Engl J Med. 2019;381(7):626-636.

13. Larkin J, Chiarion-Sileni V, Gonzalez R, Grob JJ, Rutkowski P, Lao CD, et al. Five-year survival with combined nivolumab and ipilimumab in advanced melanoma. N Engl J Med. 2019;381(16):1535-1546.

14. Hanahan D, Weinberg RA. Hallmarks of cancer: the next generation. Cell. 2011;144(5):646-674.

15. Pavlova NN, Thompson CB. The emerging hallmarks of cancer metabolism. Cell Metab. 2016;23(1):27-47.

16. Carracedo A, Cantley LC, Pandolfi PP. Cancer metabolism: fatty acid oxidation in the limelight. Nat Rev Cancer. 2013;13(4):227-232.

17. Beloribi-Djefaflia S, Vasseur S, Guillaumond F. Lipid metabolic reprogramming in cancer cells. Oncogenesis. 2016;5(1):e189.

18. Santos CR, Schulze A. Lipid metabolism in cancer. FEBS J. 2012;279(15):2610-2623.

19. Baenke F, Peck B, Miess H, Schulze A. Hooked on fat: the role of lipid synthesis in cancer metabolism and tumour development. Dis Model Mech. 2013;6(6):1353-1363.

20. Flavin R, Peluso S, Nguyen PL, Loda M. Fatty acid synthase as a potential therapeutic target in cancer. Future Oncol. 2010;6(4):551-562.

21. Xia H, Lee KW, Chen J, Kong SN, Sekar K, Deivasigamani A, et al. Simultaneous silencing of ACSL4 and induction of GADD45B in hepatocellular carcinoma cells amplifies the synergistic therapeutic effect of aspirin and sorafenib. Cell Death Discov. 2017;3:17058.

22. Wellbrock C, Arozarena I. The complexity of the ERK/MAP-kinase pathway and the treatment of melanoma skin cancer. Front Cell Dev Biol. 2016;4:33.

23. Chen WC, Wang CY, Hung YH, Weng TY, Yen MC, Lai MD. Systematic analysis of gene expression alterations and clinical outcomes for long-chain acyl-coenzyme A synthetase family in cancer. PLoS One. 2016;11(5):e0155660.

24. Li J, Condello S, Thomes-Pepin J, Ma X, Xia Y, Hurley TD, et al. Lipid desaturation is a metabolic marker and therapeutic target of ovarian cancer stem cells. Cell Stem Cell. 2017;20(3):303-314.

25. Svensson RU, Parker SJ, Eichner LJ, Kolar MJ, Wallace M, Brun SN, et al. Inhibition of acetyl-CoA carboxylase suppresses fatty acid synthesis and tumor growth of non-small-cell lung cancer in preclinical models. Nat Med. 2016;22(10):1108-1119.

26. Kuemmerle NB, Rysman E, Lombardo PS, Flanagan AJ, Lipe BC, Wells WA, et al. Lipoprotein lipase links dietary fat to solid tumor cell proliferation. Mol Cancer Ther. 2011;10(3):427-436.

27. Traves PG, Lopez-Fontal R, Cuadrado A, Luque A, Rojo AI, Faraldos JA, et al. Identification of a novel regulatory mechanism of macrophage phenotype by stearoyl-CoA desaturase. FASEB J. 2012;26(1):144-154.

28. Heidenreich S, Witte N, Weber P, Goehring I, Tolkachov A, von Loeffelholz C, et al. Retinol saturase coordinates liver metabolism by regulating ChREBP activity. Nat Commun. 2017;8:384.

29. Lyssiotis CA, Kimmelman AC. Metabolic interactions in the tumor microenvironment. Trends Cell Biol. 2017;27(11):863-875.

30. Reinfeld BI, Madden MZ, Wolf MM, Chytil A, Bader JE, Patterson AR, et al. Cell-programmed nutrient partitioning in the tumour microenvironment. Nature. 2021;593(7858):282-288.

31. Jiang M, Wu N, Xu B, Chu Y, Li X, Su S, et al. Fatty acid-induced CD36 expression via O-GlcNAcylation drives gastric cancer metastasis. Theranostics. 2019;9(18):5359-5373.

32. Wang W, Green M, Choi JE, Gijon M, Kennedy PD, Johnson JK, et al. CD8+ T cells regulate tumour ferroptosis during cancer immunotherapy. Nature. 2019;569(7755):270-274.

33. Herber DL, Cao W, Nefedova Y, Novitskiy SV, Nagaraj S, Tyurin VA, et al. Lipid accumulation and dendritic cell dysfunction in cancer. Nat Med. 2010;16(8):880-886.

34. Hossain F, Al-Khami AA, Wyczechowska D, Hernandez C, Zheng L, Reiss K, et al. Inhibition of fatty acid oxidation modulates immunosuppressive functions of myeloid-derived suppressor cells and enhances cancer therapies. Cancer Immunol Res. 2015;3(11):1236-1247.

35. Chang CH, Qiu J, O'Sullivan D, Buck MD, Noguchi T, Curtis JD, et al. Metabolic competition in the tumor microenvironment is a driver of cancer progression. Cell. 2015;162(6):1229-1241.

36. Ho PC, Bihuniak JD, Macintyre AN, Staron M, Liu X, Amezquita R, et al. Phosphoenolpyruvate is a metabolic checkpoint of anti-tumor T cell responses. Cell. 2015;162(6):1217-1228.

37. O'Sullivan D, van der Windt GJ, Huang SC, Curtis JD, Chang CH, Buck MD, et al. Memory CD8+ T cells use cell-intrinsic lipolysis to support the metabolic programming necessary for development. Immunity. 2014;41(1):75-88.

38. Ma X, Xiao L, Liu L, Ye L, Su P, Bi E, et al. CD36-mediated ferroptosis dampens intratumoral CD8+ T cell effector function and impairs their antitumor ability. Cell Metab. 2021;33(5):1001-1012.

39. Xu S, Chaudhary O, Rodriguez-Morales P, Sun X, Chen D, Zappasodi R, et al. Uptake of oxidized lipids by the scavenger receptor CD36 promotes lipid peroxidation and dysfunction in CD8+ T cells in tumors. Immunity. 2021;54(7):1561-1577.

40. Lim SA, Wei J, Nguyen TM, Shi H, Su W, Palacios G, et al. Lipid signalling enforces functional specialization of Treg cells in tumours. Nature. 2021;591(7849):306-311.

41. Cassetta L, Pollard JW. Targeting macrophages: therapeutic approaches in cancer. Nat Rev Drug Discov. 2018;17(12):887-904.

42. DeNardo DG, Ruffell B. Macrophages as regulators of tumour immunity and immunotherapy. Nat Rev Immunol. 2019;19(6):369-382.

43. Mantovani A, Marchesi F, Malesci A, Laghi L, Allavena P. Tumour-associated macrophages as treatment targets in oncology. Nat Rev Clin Oncol. 2017;14(7):399-416.

44. Qian BZ, Pollard JW. Macrophage diversity enhances tumor progression and metastasis. Cell. 2010;141(1):39-51.

45. Huang SC, Everts B, Ivanova Y, O'Sullivan D, Nascimento M, Smith AM, et al. Cell-intrinsic lysosomal lipolysis is essential for alternative activation of macrophages. Nat Immunol. 2014;15(9):846-855.

46. Nomura M, Liu J, Rovira II, Gonzalez-Hurtado E, Lee J, Wolfgang MJ, et al. Fatty acid oxidation in macrophage polarization. Nat Immunol. 2016;17(3):216-217.

47. Liu PS, Wang H, Li X, Chao T, Teav T, Christen S, et al. alpha-ketoglutarate orchestrates macrophage activation through metabolic and epigenetic reprogramming. Nat Immunol. 2017;18(9):985-994.

48. Michalek RD, Gerriets VA, Jacobs SR, Macintyre AN, MacIver NJ, Mason EF, et al. Cutting edge: distinct glycolytic and lipid oxidative metabolic programs are essential for effector and regulatory CD4+ T cell subsets. J Immunol. 2011;186(6):3299-3303.

49. Berod L, Friedrich C, Nandan A, Freitag J, Hagemann S, Harmrolfs K, et al. De novo fatty acid synthesis controls the fate between regulatory T and T helper 17 cells. Nat Med. 2014;20(11):1327-1333.

50. Angelin A, Gil-de-Gomez L, Dahiya S, Jiao J, Guo L, Levine MH, et al. Foxp3 reprograms T cell metabolism to function in low-glucose, high-lactate environments. Cell Metab. 2017;25(6):1282-1293.

51. Wang H, Franco F, Tsui YC, Xie X, Trefny MP, Zappasodi R, et al. CD36-mediated metabolic adaptation supports regulatory T cell survival and function in tumors. Nat Immunol. 2020;21(3):298-308.

52. Field CS, Baixauli F, Kyle RL, Puleston DJ, Cameron AM, Sanin DE, et al. Mitochondrial integrity regulated by lipid metabolism is a cell-intrinsic checkpoint for Treg suppressive function. Cell Metab. 2020;31(2):422-437.

53. Vivas-Garcia Y, Falletta P, Liebing J, Louphrasitthiphol P, Feng Y, Chauhan J, et al. Lineage-restricted regulation of SCD and fatty acid saturation by MITF controls melanoma phenotypic plasticity. Mol Cell. 2020;77(1):120-137.

54. Falchi L, Vargas A, Nair JR, Cheng DT, Brown JR, Jain P, et al. Lipidomic profiles in chronic lymphocytic leukemia reveal distinct features differentiating indolent from aggressive disease. Blood Cancer J. 2020;10(5):58.

55. Peck B, Schulze A. Lipid metabolism at the nexus of diet and tumor microenvironment. Trends Cancer. 2019;5(11):693-703.

56. Veglia F, Tyurin VA, Blasi M, De Leo A, Kossenkov AV, Donthireddy L, et al. Fatty acid transport protein 2 reprograms neutrophils in cancer. Nature. 2019;569(7754):73-78.

57. Newman AM, Liu CL, Green MR, Gentles AJ, Feng W, Xu Y, et al. Robust enumeration of cell subsets from tissue expression profiles. Nat Methods. 2015;12(5):453-457.

58. Korotkevich G, Sukhov V, Budin N, Shpak B, Artyomov MN, Sergushichev A. Fast gene set enrichment analysis. bioRxiv. 2021. doi: 10.1101/060012.

59. Liberzon A, Birger C, Thorvaldsdottir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Systems. 2015;1(6):417-425.

60. Ayers M, Lunceford J, Nebozhyn M, Murphy E, Loboda A, Kaufman DR, et al. IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127(8):2930-2940.

61. Cristescu R, Mogg R, Ayers M, Albright A, Murphy E, Yearley J, et al. Pan-tumor genomic biomarkers for PD-1 checkpoint blockade-based immunotherapy. Science. 2018;362(6411):eaar3593.

62. Gao J, Shi LZ, Zhao H, Chen J, Xiong L, He Q, et al. Loss of IFN-γ pathway genes in tumor cells as a mechanism of resistance to anti-CTLA-4 therapy. Cell. 2016;167(2):397-404.

63. Zaretsky JM, Garcia-Diaz A, Shin DS, Escuin-Ordinas H, Hugo W, Hu-Lieskovan S, et al. Mutations associated with acquired resistance to PD-1 blockade in melanoma. N Engl J Med. 2016;375(9):819-829.

64. Shin DS, Zaretsky JM, Escuin-Ordinas H, Garcia-Diaz A, Hu-Lieskovan S, Kalbasi A, et al. Primary resistance to PD-1 blockade mediated by JAK1/2 mutations. Cancer Discov. 2017;7(2):188-201.

65. Hörber S, Dahlhoff C, Schneider MR, Wolf E. Genomic aspects of MITF in melanoma and melanocytes. Pigment Cell Melanoma Res. 2013;26(5):647-657.

66. Carreira S, Goodall J, Denat L, Rodriguez M, Nuciforo P, Hoek KS, et al. Mitf regulation of Dia1 controls melanoma proliferation and invasiveness. Genes Dev. 2006;20(24):3426-3439.

67. Caramel J, Papadogeorgakis E, Hill L, Browne GJ, Richard G, Wierinckx A, et al. A switch in the expression of embryonic EMT-inducers drives the development of malignant melanoma. Cancer Cell. 2013;24(4):466-480.

68. Hoek KS, Eichhoff OM, Schlegel NC, Döbbeling U, Kobert N, Schaerer L, et al. In vivo switching of human melanoma cells between proliferative and invasive states. Cancer Res. 2008;68(3):650-656.

69. Guo D, Bell EH, Mischel P, Chakravarti A. Targeting SREBP-1-driven lipid metabolism to treat cancer. Curr Pharm Des. 2014;20(15):2619-2626.

70. Li J, Cheng D, Zhu M, Yu H, Pan Z, Liu L, et al. OTUB2 stabilizes FASN to promote tumor progression in hepatocellular carcinoma via deubiquitination. FEBS Lett. 2019;593(23):3321-3334.

71. Lee JV, Carrer A, Shah S, Snyder NW, Wei S, Venneti S, et al. Akt-dependent metabolic reprogramming regulates tumor cell histone acetylation. Cell Metab. 2014;20(2):306-319.

72. Porstmann T, Santos CR, Griffiths B, Cully M, Wu M, Leevers S, et al. SREBP activity is regulated by mTORC1 and contributes to Akt-dependent cell growth. Cell Metab. 2008;8(3):224-236.

73. Fritz V, Fajas L. Metabolism and proliferation share common regulatory pathways in cancer cells. Oncogene. 2010;29(31):4369-4377.

74. Guri Y, Colombi M, Dazert E, Hindupur SK, Roszik J, Moes S, et al. mTORC2 promotes tumorigenesis via lipid synthesis. Cancer Cell. 2017;32(6):807-823.

75. Koundouros N, Poulogiannis G. Reprogramming of fatty acid metabolism in cancer. Br J Cancer. 2020;122(1):4-22.

76. Samstein RM, Lee CH, Shoushtari AN, Hellmann MD, Shen R, Janjigian YY, et al. Tumor mutational burden predicts survival after immunotherapy across multiple cancer types. Nat Genet. 2019;51(2):202-206.

77. Goodman AM, Kato S, Bazhenova L, Patel SP, Frampton GM, Miller V, et al. Tumor mutational burden as an independent predictor of response to immunotherapy in diverse cancers. Mol Cancer Ther. 2017;16(11):2598-2608.

78. Davies H, Bignell GR, Cox C, Stephens P, Edkins S, Clegg S, et al. Mutations of the BRAF gene in human cancer. Nature. 2002;417(6892):949-954.

79. Wang D, DuBois RN. Eicosanoids and cancer. Nat Rev Cancer. 2010;10(3):181-193.

80. Chowdhury PS, Chamoto K, Kumar A, Honjo T. PPAR-induced fatty acid oxidation in T cells increases the number of tumor-reactive CD8+ T cells and facilitates anti-PD-1 therapy. Cancer Immunol Res. 2018;6(11):1375-1387.

81. Facciabene A, Motz GT, Coukos G. T-regulatory cells: key players in tumor immune escape and angiogenesis. Cancer Res. 2012;72(9):2162-2171.

82. Tanaka A, Sakaguchi S. Regulatory T cells in cancer immunotherapy. Cell Res. 2017;27(1):109-118.

83. Van den Eynde BJ, van der Bruggen P. T cell defined tumor antigens. Curr Opin Immunol. 1997;9(5):684-693.

84. Dunn GP, Old LJ, Schreiber RD. The three Es of cancer immunoediting. Annu Rev Immunol. 2004;22:329-360.

85. Igney FH, Krammer PH. Immune escape of tumors: apoptosis resistance and tumor counterattack. J Leukoc Biol. 2002;71(6):907-920.

86. Peinado H, Zhang H, Matei IR, Costa-Silva B, Hoshino A, Rodrigues G, et al. Pre-metastatic niches: organ-specific homes for metastases. Nat Rev Cancer. 2017;17(5):302-317.

87. Orillion A, Hashimoto A, Damayanti N, Shen L, Adelaiye-Ogala R, Arisa S, et al. Entinostat neutralizes myeloid-derived suppressor cells and enhances the antitumor effect of PD-1 inhibition in murine models of lung and renal cell carcinoma. Clin Cancer Res. 2017;23(17):5187-5201.

88. Pitt JM, Marabelle A, Eggermont A, Soria JC, Kroemer G, Zitvogel L. Targeting the tumor microenvironment: removing obstruction to anticancer immune responses and immunotherapy. Ann Oncol. 2016;27(8):1482-1492.

89. Kridel R, Sehn LH, Gascoyne RD. Pathogenesis of follicular lymphoma. J Clin Invest. 2012;122(10):3424-3431.

90. Buckley D, Duke G, Heuer TS, O'Farrell M, Wagman AS, McCulloch W, et al. Fatty acid synthase - Modern tumor cell biology insights into a classical oncology target. Pharmacol Ther. 2017;177:23-31.

91. Falchook G, Infante J, Arkenau HT, Patel MR, Dean E, Borazanci E, et al. First-in-human study of the safety, pharmacokinetics, and pharmacodynamics of first-in-class fatty acid synthase inhibitor TVB-2640 alone and with a taxane in advanced tumors. EClinicalMedicine. 2021;34:100797.

92. Ventura R, Mordec K, Waszczuk J, Wang Z, Lai J, Fridlib M, et al. Inhibition of de novo palmitate synthesis by fatty acid synthase induces apoptosis in tumor cells by remodeling cell membranes, inhibiting signaling pathways, and remodeling the cell transcriptome. EBioMedicine. 2015;2(8):808-824.

93. Jiang L, Fang X, Wang H, Li D, Wang X. Ovarian cancer-intrinsic fatty acid synthase prevents anti-tumor immunity by disrupting tumor-infiltrating dendritic cells. Front Immunol. 2018;9:2927.

94. Tang JJ, Li JG, Qi W, Qiu WW, Li PS, Li BL, et al. Inhibition of SREBP by a small molecule, betulin, improves hyperlipidemia and insulin resistance and reduces atherosclerotic plaques. Cell Metab. 2011;13(1):44-56.

95. Gholkar AA, Cheung K, Williams KJ, Lo YC, Hamideh SA, Nnebe C, et al. Fatostatin inhibits cancer cell proliferation by affecting mitotic microtubule spindle assembly and cell division. J Biol Chem. 2016;291(33):17001-17008.

96. Klement RJ, Kämmerer U. Is there a role for carbohydrate restriction in the treatment and prevention of cancer? Nutr Metab (Lond). 2011;8:75.

97. Weber DD, Aminazdeh-Gohari S, Kofler B. Ketogenic diet in cancer therapy. Aging (Albany NY). 2018;10(2):164-165.

98. Lussier DM, Woolf EC, Johnson JL, Brooks KS, Blattman JN, Scheck AC. Enhanced immunity in a mouse model of malignant glioma is mediated by a therapeutic ketogenic diet. BMC Cancer. 2016;16:310.

99. Hopkins BD, Pauli C, Du X, Wang DG, Li X, Wu D, et al. Suppression of insulin feedback enhances the efficacy of PI3K inhibitors. Nature. 2018;560(7719):499-503.

100. Yan Y, Tao H, He J, Huang SY. The HDOCK server for integrated protein-protein docking. Nat Protoc. 2020;15(5):1829-1852.

101. Yadav UP, Singh T, Kumar P, Sharma P, Kaur H, Sharma S, et al. Metabolic adaptations in cancer stem cells. Front Oncol. 2020;10:1010.

102. Luu TH, Bard JM, Carbonnelle D, Chaillou C, Huvelin JM, Bobin-Dubigeon C, et al. Lithocholic bile acid inhibits lipogenesis and induces apoptosis in breast cancer cells. Cell Oncol (Dord). 2018;41(1):13-24.

103. Su P, Wang Q, Bi E, Ma X, Liu L, Yang M, et al. Enhanced lipid accumulation and metabolism are required for the differentiation and activation of tumor-associated macrophages. Cancer Res. 2020;80(7):1438-1450.

104. Zhao F, Xiao C, Evans KS, Theivanthiran T, DeVito N, Holtzhausen A, et al. Paracrine Wnt5a-β-catenin signaling triggers a metabolic program that drives dendritic cell tolerization. Immunity. 2018;48(1):147-160.

105. Bader JE, Voss K, Rathmell JC. Targeting metabolism to improve the tumor microenvironment for cancer immunotherapy. Mol Cell. 2020;78(6):1019-1033.

106. Renner K, Singer K, Koehl GE, Geissler EK, Peter K, Siska PJ, et al. Metabolic hallmarks of tumor and immune cells in the tumor microenvironment. Front Immunol. 2017;8:248.

107. Onodera Y, Nam JM, Bissell MJ. Increased sugar uptake promotes oncogenesis via EPAC/RAP1 and O-GlcNAc pathways. J Clin Invest. 2014;124(1):367-384.

108. Elia I, Schmieder R, Christen S, Fendt SM. Organ-specific cancer metabolism and its potential for therapy. Handb Exp Pharmacol. 2016;233:321-353.

109. Jang C, Chen L, Rabinowitz JD. Metabolomics and isotope tracing. Cell. 2018;173(4):822-837.

110. Boroughs LK, DeBerardinis RJ. Metabolic pathways promoting cancer cell survival and growth. Nat Cell Biol. 2015;17(4):351-359.

111. Wolrab D, Jirásko R, Peterka O, Idkowiak J, Chocholoušková M, Holčapek M. Lipidomic profiling of human serum enables detection of pancreatic cancer. Nat Commun. 2022;13(1):124.

112. Matyash V, Liebisch G, Kurzchalia TV, Shevchenko A, Schwudke D. Lipid extraction by methyl-tert-butyl ether for high-throughput lipidomics. J Lipid Res. 2008;49(5):1137-1146.

113. Kamphorst JJ, Cross JR, Fan J, de Stanchina E, Mathew R, White EP, et al. Hypoxic and Ras-transformed cells support growth by scavenging unsaturated fatty acids from lysophospholipids. Proc Natl Acad Sci USA. 2013;110(22):8882-8887.

114. Satpathy M, Wang L, Zielinski RJ, Qian W, Melendez A, Nandy SK, et al. Targeted drug delivery to tumors via the nanoparticle-albumin carrier technology. Methods Mol Biol. 2013;1030:219-238.

115. Angelo LS, Kurzrock R. Vascular endothelial growth factor and its relationship to inflammatory mediators. Clin Cancer Res. 2007;13(10):2825-2830.

116. Saleh R, Elkord E. Acquired resistance to cancer immunotherapy: Role of tumor-mediated immunosuppression. Semin Cancer Biol. 2020;65:13-27.

117. Galon J, Bruni D. Approaches to treat immune hot, altered and cold tumours with combination immunotherapies. Nat Rev Drug Discov. 2019;18(3):197-218.

118. Greaves M, Maley CC. Clonal evolution in cancer. Nature. 2012;481(7381):306-313.

119. Stoeckius M, Hafemeister C, Stephenson W, Houck-Loomis B, Chattopadhyay PK, Swerdlow H, et al. Simultaneous epitope and transcriptome measurement in single cells. Nat Methods. 2017;14(9):865-868.

120. Stoeckius M, Zheng S, Houck-Loomis B, Hao S, Yeung BZ, Mauck WM 3rd, et al. Cell Hashing with barcoded antibodies enables multiplexing and doublet detection for single cell genomics. Genome Biol. 2018;19(1):224.

121. Ecker J, Liebisch G. Application of stable isotopes to investigate the metabolism of fatty acids, glycerophospholipid and sphingolipid species. Prog Lipid Res. 2014;54:14-31.

122. Beyoğlu D, Idle JR. Metabolomic insights into the mode of action of natural products in the treatment of liver disease. Biochem Pharmacol. 2020;180:114171.

123. Cancer Genome Atlas Research Network. Genomic classification of cutaneous melanoma. Cell. 2015;161(7):1681-1696.

124. Hoadley KA, Yau C, Hinoue T, Wolf DM, Lazar AJ, Drill E, et al. Cell-of-origin patterns dominate the molecular classification of 10,000 tumors from 33 types of cancer. Cell. 2018;173(2):291-304.

125. Hopkins KC, Tartell MA, Herrmann C, Hackett BA, Taschuk F, Panda D, et al. Virus-induced translational arrest through 4EBP1/2-dependent decay of 5'-TOP mRNAs restricts viral infection. Proc Natl Acad Sci USA. 2015;112(46):E6214-E6223.
