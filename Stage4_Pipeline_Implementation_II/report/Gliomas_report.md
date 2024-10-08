## Title
**Unsupervised Clustering of Gliomas Using Methylation and RNA Expression Levels to Classify IDH Status**

### Authors
- Favour Igwezeke<sup>1</sup>
- Jessica Ovabor<sup>2</sup>
- Oluwatobiloba Johnson Osedimilehin<sup>3</sup>
- Ogochukwu Nwaigwe<sup>4</sup>
- Muhammad Faizan Khan<sup>5</sup>
- Rokaya Yasser<sup>7</sup>

---

# Abstract
In this study, we aim to reproduce the unsupervised clustering of low-grade gliomas (n = 516 LGG) based on methylation levels to classify the six IDH statuses. 
Leveraging machine learning techniques, specifically the K-Nearest Neighbors (KNN) algorithm, we attempted to replicate the clustering process, as outlined in the study “Unsupervised Clustering of Gliomas Identifies Six Methylation Groups and Four RNA Expression Groups Associated with IDH Status” by Ceccarelli et al. (2016). 
The clustering was performed on methylation data obtained from The Cancer Genome Atlas (TCGA), where the gliomas were pre-labeled based on IDH mutations. 
Our results revealed that KNN alone did not provide sufficiently high classification accuracy for distinct clusters, but the analysis offered insights into potential factors impacting glioma classification and paved the way for future improvements using alternative techniques. 
We also compared our results with the findings of the referenced paper, identifying both overlaps and discrepancies in the clustering performance.

## 1. Introduction
Gliomas are the most common type of primary brain tumors, arising from glial cells. They are notoriously difficult to treat due to their highly infiltrative nature, and they are commonly associated with poor prognoses (Chen et al., 2021). Mutations in the Isocitrate Dehydrogenase (IDH) genes, particularly IDH1 and IDH2, have been found to play a significant role in gliomagenesis, especially in low-grade gliomas (LGG) (Cohen et al., 2013). 
These mutations are early events in tumor development and serve as important biomarkers for diagnosis, prognosis, and treatment strategies.These mutations often coexist with specific genetic alterations, such as TP53 mutations in astrocytic tumors and 1p/19q co-deletions in oligodendrogliomas (Cohen et al.). The presence of IDH mutations plays a critical role in treatment strategies and prognostic evaluations for low-grade gliomas, emphasizing their significance in tumorigenesis (Chen et al.).

Recent studies have demonstrated that molecular profiling of gliomas based on methylation and RNA expression data can help classify gliomas into distinct subtypes, with significant prognostic implications (Ceccarelli et al., 2016). 
This study aims to reproduce the unsupervised clustering methodology used to identify six methylation-based groups and four RNA expression-based groups, associated with IDH statuses.

## 2. Dataset and Preprocessing
This analysis utilizes gene expression data from the TCGA Low-Grade Glioma project (TCGA-LGG). After querying and downloading the transcriptome profiling data, raw count data was extracted for further analysis. Sample IDs were obtained, and IDH status was determined using subtype data.

### 2.1 Data Cleaning and Normalization
Data cleaning involved excluding samples without IDH status and factorizing the remaining samples into Wildtype and Mutant categories. The expression dataset was normalized using the 'TCGAanalyze_Normalization' function, which adjusted for gene length. Filtering was then applied to remove low-expression genes based on a 25th percentile cutoff. This ensured a high-quality dataset for differential expression analysis, enhancing the reliability of the results by minimizing noise.

## 3. Methodology
### 3.1 Differential Expression Analysis
Differential expression analysis was conducted to compare IDH Wildtype and Mutant samples using a log-fold change cutoff of 2.5 and a false discovery rate of 0.01, resulting in identification of underregulated and overregulated genes in gliomas.

### 3.2 Random Forest and KNN Classification
To evaluate the glioma dataset, both Random Forest and KNN models were employed:

**Random Forest Model:**
- Accuracy: 17.15%
- Kappa: 0.0057
- Balanced Accuracy by Class ranged from 49.95% to 50.65%.

**KNN Model (k = 3):**
- Accuracy: 17.32%
- Kappa: 0.0078
- Balanced Accuracy by Class ranged from 49.74% to 51.36%.

Both models had low overall accuracies and kappas, suggesting significant challenges in predicting glioma clusters based on gene expression data alone, possibly due to high inter-cluster similarity and insufficient discriminatory features.

### 3.3 Clustering Analysis
To further explore the structure within the glioma dataset, we applied K-means clustering with six centers. The silhouette score, which measures how well each point lies within its assigned cluster, was 0.81, indicating that the clustering generally formed well-defined groups. This was confirmed by evaluating the confusion matrices of the models and comparing the clusters to the known classes. Both Random Forest and KNN models exhibited similar misclassification patterns, with sensitivity and specificity values that suggested the difficulty in clearly distinguishing between IDH Wildtype and Mutant gliomas.

## 4.0 Enrichment Analysis (EA)
EA was carried out with the `TCGAanalyze_EAComplete()` function provided by the TCGAbiolinks package in RStudio. EA revealed significantly enriched biological processes in lists of overexpressed and underexpressed genes. Notable enriched biological processes include:

- **Overexpressed Genes:** GPCR signaling pathway, sensory perception, and neurological system processes.
- **Underexpressed Genes:** Pattern specification process, anterior/posterior pattern formation, and regionalization.

### Figure 1
Lollipop plot of the top pathways enriched with downregulated/underexpressed genes.

### Figure 2
Lollipop plot of the top pathways enriched with upregulated/overexpressed genes.

Top-ranked upregulated genes include the protein-coding gene DEFB119, which plays a role in cell/tissue defense (Martinelli et al., 2020), and long non-coding (lnc) RNA genes LINC01602, LINC03018, which may be involved in transcriptional control (Zhang et al., 2019). The top-ranked downregulated genes include TTR (transthyretin), which is often not expressed in gliomas (Albrecht et al., 1995), and HOXC12, involved in morphological processes (Kawasumi-Kita et al., 2024).
