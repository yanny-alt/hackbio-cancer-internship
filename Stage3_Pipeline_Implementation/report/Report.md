# Stage 03 Task: HackBio Bioinformatics Internship - Clinical Oncology

## Title
**Identifying Potential Biomarkers in Breast Cancer Using Differential Expression and Machine Learning Approaches**

### Authors
- Favour Igwezeke<sup>1</sup>
- Jessica Ovabor<sup>2</sup>
- Anarghya Hegde<sup>3</sup>
- Oluwatobiloba Johnson Osedimilehin<sup>4</sup>
- Ogochukwu Nwaigwe<sup>5</sup>
- Muhammad Faizan Khan<sup>6</sup>
- Rokaya Yasser<sup>7</sup>
- Muhammad Faheem Raziq<sup>8</sup>

---

## Contributors Information
| Name | Email | Slack ID |
|---|---|---|
| Favour Igwezeke | beingfave@gmail.com | [Slack](https://hackbiointern-leo4437.slack.com/team/U07KE59TWEP) |
| Jessica Ovabor | ovaborjessica85@gmail.com | [Slack](https://hackbiointern-leo4437.slack.com/team/U07JVPSU917) |
| Anarghya Hegde | anarghyahegde@outlook.com | [Slack](https://hackbiointern-leo4437.slack.com/team/U07JM2UDL7R) |
| Oluwatobiloba Johnson Osedimilehin | tobijohnson01@gmail.com | [Slack](https://hackbiointern-leo4437.slack.com/team/U07JP06QDB4) |
| Ogochukwu Nwaigwe | nwaigweogochukwu756@gmail.com | [Slack](https://hackbiointern-leo4437.slack.com/team/U07KP1D2F24) |
| Muhammad Faizan Khan | faizanjeee@hotmail.com | [Slack](https://hackbiointern-leo4437.slack.com/team/U07JJHPBFBL) |
| Rokaya Yasser | rokayayasser727@gmail.com | [Slack](https://hackbiointern-leo4437.slack.com/team/U07KUECLR40) |
| Muhammad Faheem Raziq | faheemraziq1999@gmail.com | [Slack](https://hackbiointern-leo4437.slack.com/team/U07JPMN7EDD) |

---

# Abstract
Utilising the TCGA-BRCA RNA-Seq data, this study combines machine learning with differential gene expression (DGE) to identify key biomarkers in breast cancer. Significantly upregulated genes like CGA and CST5 and downregulated genes like XIRP2 and MYL, which are associated with tumour growth and prognosis were revealed by DGE. With 100% accuracy, a Random Forest model was used to Classify normal and tumour samples where Dim2 and Dim9 presented to be the most crucial features for classification. The findings of this study demonstrates the biomarkers’ potential to enhance breast cancer diagnosis and treatment approaches.

## 1. Introduction
Breast cancer is one of the most prevalent cancers worldwide, and the 
discovery of potential biomarkers can aid in early diagnosis and treatment. 
Using the TCGA-BRCA dataset, which contains gene expression data from 
breast cancer patients, this report presents a comprehensive analysis 
focusing on differential gene expression and machine learning for biomarker
discovery

### 1.1 Why TCGA-BRCA Dataset?
The TCGA-BRCA dataset is optimal due to its vast sample size, covering over 1000 breast cancer cases. This large-scale data enables robust statistical power in investigating both machine learning applications and differential expression analysis. Additionally, the extensive clinical annotations provided allow for correlating molecular traits with clinical outcomes, making it a comprehensive resource for breast cancer research.

## 2. Description of the Dataset and Preprocessing
Using RStudio, code was written to query and retrieve TCGA-BRCA RNA-Seq data for breast cancer, with a focus on 20 samples each of “Primary Tumor” and “Solid Tissue Normal”. The chosen samples were retrieved along with the gene expression count matrix which was then stored as a CSV file.

### 2.1 Data Cleaning and Normalization

The data underwent normalization with the `TCGAanalyze_normalization()` function. This normalization was done to account for the differences in the length of each gene. The data was then filtered with the `TCGAanalyze_Filtering()` function. A quantile cut threshold of 0.25 was chosen to remove the genes whose expression levels were under the 25th percentile, as these genes could potentially contribute to noise. This filtration step was to ensure more robust downstream analysis.

 ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage3_Pipeline_Implementation/visualizations/heatmap.jpg)
                 *Figure 1: Heatmap of filtered dataset (Green: Low, Black: Intermediate, Red: High expression)*
   ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage3_Pipeline_Implementation/visualizations/Screenshot_20240925_205416_Docs.jpg)
                  *Figure 2: Heatmap clustered by row*
  ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage3_Pipeline_Implementation/visualizations/Screenshot_20240925_205638_Docs.jpg)
                   *Figure 3: Heatmap clustered by column*

## 3. Methodology 

### 3.1 Differential Expression Analysis
EdgeR was the pipeline of choice for analysis of the differential expression of genes between the solid tissue normal (STN) and primary tumor (PT) samples in the filtered breast cancer (BRCA) gene expression dataset. The analysis was done with the `TCGAanalyze_DEA()` function. An FDR cut-off of 0.01 and a log2 fold change cut off 2 were chosen. The result of the analysis was a group of 3462 genes with `abs(log2Foldchange >= 2)` and `FDR <= 0.01`. 98 overregulated and 129 underregulated genes sets were then selected based on stringent log2Foldchange thresholds of `>= 6` and `<= -6` for the former and latter respectively.

 ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage3_Pipeline_Implementation/visualizations/volcano.jpg)
                       *Figure 4: Volcano plot showing downregulated and upregulated genes.*

### 3.2 Machine Learning Preparation
Gene expression data was reshaped for machine learning, with the Random 
Forest algorithm chosen for classification. Data was split into training (70%) 
and test (30%) sets. Dimensionality reduction was performed using Principal
Component Analysis (PCA) to address the high-dimensional nature of the 
dataset. The Random Forest model trained with 500 trees showed high 
performance on the test set.

### 3.3 Potential Biomarkers and Enrichment Analysis
This section presents a comprehensive analysis of gene expression data, focusing on identifying key upregulated and downregulated genes and performing enrichment analysis. The process begins with the gene enrichment analysis of these genes. Libraries such as `TCGAbiolinks`, `biomaRt`, and `ggplot2` were utilized to facilitate the extraction of gene symbols. Enrichment analysis was performed using the `TCGAanalyze_EAcomplete` function, and `ggplot2` was employed to create lollipop plots of the top five enriched pathways. These plots included FDR values and gene count, visually representing the most significantly enriched pathways, with circle size representing the number of genes involved and color intensity displayed based on the FDR values.

## 4. Results and Interpretation

### 4.1 Identified Biomarkers
The analysis revealed several significant biomarkers:

**Top Upregulated Biomarkers**

| Gene  | Fold Change | P-Value     | Significance                                                                                             |
|-------|-------------|-------------|---------------------------------------------------------------------------------------------------------|
| CST5  | 10.15       | 4.99E-71    | Member of the cystatin peptide superfamily. Proteases play a role in tumor development, especially cysteine cathepsin.(4) |
| CGA   | 9.84        | 9.27E-67    | Increased in breast cancer tissue and associated with poor prognosis. A candidate for targeted therapy.(5) |
| HTN   | 9.66        | 1.65E-60    | Significantly higher in tumors than normal tissues.(6)                                                   |
| CLEC3A| 9.39        | 2.40E-69    | Correlates with metastatic potential and poor prognosis in breast cancer.(7)                              |

**Top Downregulated Biomarkers**

| Gene  | Fold Change | P-Value     | Significance                                                                                                      |
|-------|-------------|-------------|------------------------------------------------------------------------------------------------------------------|
| MYL   | -14.13      | 3.17E-105   | Functions as a prognostic marker in breast cancer, associated with immune infiltration.(8)                        |
| XIRP2 | -13.24      | 3.57E-100   | Encodes an actin-binding protein, significantly mutated in breast cancer metastasis.(9)                           |
| MYH   | -13.16      | 3.26E-102   | Part of the base repair pathway, detecting and protecting against oxidative DNA damage. Low expression leads to mutation and cancer.(10) |

   ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage3_Pipeline_Implementation/visualizations/upreggenepathway_plot%20(1).png)
                   *Figure 5: Upregulated genes enriched pathways.*
   ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage3_Pipeline_Implementation/visualizations/downreggenepathway_plot.png)
                     *Figure 6: Downregulated genes enriched pathways.*

### 4.2 Model Performance
**Confusion Matrix**

| Prediction/Reference | Normal | Tumor |
|----------------------|--------|-------|
| **Normal**           | 6      | 0     |
| **Tumor**            | 0      | 6     |


   ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage3_Pipeline_Implementation/visualizations/Confusion%20Matrix.png)
                   *Figure 7: Confusion Matrix*

The Random Forest model achieved perfect classification with 100% accuracy, as indicated by the confusion matrix, with no misclassifications between "Normal" and "Tumor" samples. Key performance metrics include a 95% Confidence Interval of (0.7354, 1), a Kappa statistic of 1, and perfect sensitivity, specificity, precision, and negative predictive value, all at 1. The model’s performance is statistically significant, and cross-validation confirmed its robustness with a mean accuracy of 96.67%, ensuring generalizability.

   ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage3_Pipeline_Implementation/visualizations/rf_model.png)
                 *Figure 8: Random Forest model*

## 5. Conclusion and Future Directions
This study successfully identified potential biomarkers like CST5 and MYL in
breast cancer using differential expression analysis and machine learning. 
The Random Forest model demonstrated strong predictive performance, 
suggesting these biomarkers may hold clinical value.

**Future Directions:**
- Validation of findings in larger cohorts.
- Application of this methodology to other cancer subtypes.
- Integration with clinical data to enhance predictive power.

## 6. Reference

1. Colaprico, A., Silva, T. C., Olsen, C., Garofano, L., Cava, C., Garolini, D., ... & Noushmehr, H. (2016). TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. *Nucleic acids research, 44*(8), e71-e71. [http://doi.org/10.1093/nar/gkv1507](http://doi.org/10.1093/nar/gkv1507).

2. Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2023). dplyr: a grammar of data manipulation. R package version 1.1.2. *Computer software*.

3. Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L. D. A., François, R., ... & Yutani, H. (2019). Welcome to the Tidyverse. *Journal of open source software, 4*(43), 1686. [doi:10.21105/joss.01686](https://doi.org/10.21105/joss.01686).

4. Xiong, S., Wen, H., Dai, L., Lou, Y., Wang, Z., Yi, Y., Yan, X., ... & Wu, G. (2023). A brain-tumor neural circuit controls breast cancer progression in mice. *Journal of Clinical Investigation, 133*(24). [https://doi.org/10.1172/jci167725](https://doi.org/10.1172/jci167725).

5. Zhou, J., Zhu, X., Wu, S., & Chen, Y. (2022). Glycoprotein hormone α‑subunit promotes cell proliferation and tumorigenesis in breast cancer. *Oncology Letters, 23*(5). [https://doi.org/10.3892/ol.2022.13263](https://doi.org/10.3892/ol.2022.13263).

6. Wongpanuwich, W., Yodsanga, S., Chaisuparat, R., & Amornphimoltham, P. (2022). Association Between PD-L1 and Histatin1, 3 Expression in Advanced Head and Neck Squamous Cell Carcinoma. *Anticancer Research, 42*(5), 2689–2699. [https://doi.org/10.21873/anticanres.15747](https://doi.org/10.21873/anticanres.15747).

7. Ni, J., Peng, Y., Yang, F., Xi, X., Huang, X., & He, C. (2018). Overexpression of CLEC3A promotes tumor progression and poor prognosis in breast invasive ductal cancer. *OncoTargets and Therapy, Volume 11*, 3303–3312. [https://doi.org/10.2147/ott.s161311](https://doi.org/10.2147/ott.s161311).

8. Lv, M. (2023). MYL5 as a Novel Prognostic Marker is Associated with Immune Infiltrating in Breast Cancer: A Preliminary Study. *The Breast Journal, 2023*, 1–21. [https://doi.org/10.1155/2023/9508632](https://doi.org/10.1155/2023/9508632).

9. Chen, H., Liu, M., Luo, Y., Zhao, L., Huang, L., Xiao, M., Qi, C., & Yang, Z. (2022). XIRP2 mutation as an indicator stratified patients benefit from immune checkpoint inhibitors in NSCLC. *Journal of Clinical Oncology, 40*(16_suppl), e14573. [https://doi.org/10.1200/jco.2022.40.16_suppl.e14573](https://doi.org/10.1200/jco.2022.40.16_suppl.e14573).

10. Cleary, S. P., Cotterchio, M., Jenkins, M. A., Kim, H., Bristow, R., Green, R., ... & Gallinger, S. (2009). Germline MutY Human Homologue Mutations and Colorectal Cancer: A Multisite Case-Control Study. *Gastroenterology, 136*(4), 1251–1260. [https://doi.org/10.1053/j.gastro.2008.12.050](https://doi.org/10.1053/j.gastro.2008.12.050).
