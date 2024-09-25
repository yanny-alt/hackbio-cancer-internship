# Stage 03 Task: HackBio Internship - Cancers

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

# Why TCGA-BRCA Dataset?

Breast cancer dataset is optimal due to its high prevalence and as the most frequent cancer worldwide, breast cancer is a crucial research focus. With more than 1000 samples, it provides strong statistical power for investigation of machine learning and differential expression. Molecular traits and clinical outcomes can be correlated thanks to its extensive clinical annotations.

## Dataset Download and Extraction

Using RStudio, code was written to query and retrieve TCGA-BRCA RNA-Seq data for breast cancer, with a focus on 20 samples each of “Primary Tumor” and “Solid Tissue Normal”. These samples were chosen such that the first 20 samples were solid tissue normal, and the remaining 20 samples were of the BRCA primary tumor. The chosen samples were retrieved along with the gene expression count matrix which was then stored as a CSV file.

## Data Preprocessing Workflow

The data underwent normalization with the `TCGAanalyze_normalization()` function. This normalization was done to account for the differences in the length of each gene. The data was then filtered with the `TCGAanalyze_Filtering()` function. A quantile cut threshold of 0.25 was chosen to remove the genes whose expression levels were under the 25th percentile, as these genes could potentially contribute to noise. This filtration step was to ensure more robust downstream analysis.

**Figures**:
- **Figure 1**: Heatmap of filtered dataset (Green: Low, Black: Intermediate, Red: High expression)
- **Figure 2**: Heatmap clustered by row
- **Figure 3**: Heatmap clustered by column

## 4. Differential Expression Analysis
EdgeR was the pipeline of choice for analysis of the differential expression of genes between the solid tissue normal (STN) and primary tumor (PT) samples in the filtered breast cancer (BRCA) gene expression dataset. The analysis was done with the `TCGAanalyze_DEA()` function. An FDR cut-off of 0.01 and a log2 fold change cut off 2 were chosen. The result of the analysis was a group of 3462 genes with `abs(log2Foldchange >= 2)` and `FDR <= 0.01`. 98 overregulated and 129 underregulated genes sets were then selected based on stringent log2Foldchange thresholds of `>= 6` and `<= -6` for the former and latter respectively.

**Figure 4**: Volcano plot showing downregulated and upregulated genes.

## 5. Enrichment Analysis
The gene enrichment analysis for upregulated and downregulated genes using various bioinformatics libraries. It begins by loading libraries such as `TCGAbiolinks`, `biomaRt`, and `ggplot2`. The script reads CSV files containing Ensembl gene IDs, retrieves corresponding gene symbols from the Ensembl database, and merges this information into new datasets. After preparing gene lists for enrichment analysis, it utilizes the `TCGAanalyze_EAcomplete` function to perform the analysis and save the results in CSV format. The script reshapes the enrichment data using functions from the `tidyr` package, extracting relevant details like GO terms and FDR values. Finally, it creates lollipop plots for the top five enriched pathways, employing `ggplot2` for visualization. The plots display FDR values, and the number of genes associated with each pathway, which are then saved as PNG files.

**Figure 5**: Upregulated genes enriched pathways.  
**Figure 6**: Downregulated genes enriched pathways.

## 6. Upregulated Genes
### 6.1 Platelet Activation
This pathway plays a key role in hematogenous metastasis, protecting metastatic tumor cells from immune evasion and apoptosis.

### 6.2 Hindbrain & Metencephalon Development
Neural circuits associated with cancer-related anxiety may indicate deeper neuroimmune interactions in the tumor microenvironment.

## 7. Downregulated Genes
### 7.1 Muscle Tissue Development
Downregulation is linked with muscle fiber reduction and mitochondrial dysfunction, contributing to cancer-related cachexia.

## 8. Machine Learning Analysis
### 8.1 Dataset Preparation
The dataset was split into 70% training and 30% test data. Labels were assigned as "Normal" or "Tumor" for machine learning classification. 

### 8.2 Feature Selection and Dimensionality Reduction
Using Principal Component Analysis (PCA), features were selected based on the most relevant biomarkers identified during enrichment analysis to reduce dimensionality.

## 9. Model Selection and Training
Random Forest was chosen for its ability to handle high-dimensional data. The model was trained with 500 trees and achieved high performance.

## 10. Model Performance
### 10.1 Confusion Matrix
| Prediction/Reference | Normal | Tumor |
|----------------------|--------|-------|
| **Normal**           | 6      | 0     |
| **Tumor**            | 0      | 6     |

The model achieved:
- **Accuracy**: 100%
- **95% CI**: (0.7354, 1)
- **Kappa**: 1
- **Sensitivity & Specificity**: Both 1

## 11. Feature Importance Analysis
Random Forest feature importance shows which features significantly impact the model's predictions.
