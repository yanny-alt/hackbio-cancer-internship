## Title
**Unsupervised Clustering of Gliomas Using Methylation and RNA Expression Levels to Classify IDH Status**

### Authors
- Favour Igwezeke<sup>1</sup>
- Jessica Ovabor<sup>2</sup>
- Oluwatobiloba Johnson Osedimilehin<sup>3</sup>
- Muhammad Faizan Khan<sup>4</sup>
- Rokaya Yasser<sup>5</sup>

---

# Abstract
This study aims to reproduce the unsupervised clustering of low-grade gliomas (n = 516 LGG) based on expression levels to classify  IDH statuses. Utilizing machine learning, specifically the K-Nearest Neighbors (KNN) algorithm, we attempted to replicate the clustering from “Unsupervised Clustering of Gliomas Identifies Four RNA Expression Groups Associated with IDH Status” by Ceccarelli et al. (2016). Clustering was performed on expression data from The Cancer Genome Atlas (TCGA), where gliomas were pre-labeled by IDH mutations. Our results showed that KNN did not yield high classification accuracy, but provided insights into factors influencing glioma classification.

## 1. Introduction
Gliomas are the most common type of primary brain tumors, arising from glial cells. They are notoriously difficult to treat due to their highly infiltrative nature, and they are commonly associated with poor prognoses (Chen et al., 2021). Mutations in the Isocitrate Dehydrogenase (IDH) genes, particularly IDH1 and IDH2, have been found to play a significant role in gliomagenesis, especially in low-grade gliomas (LGG) (Cohen et al., 2013). 
These mutations are early events in tumor development and serve as important biomarkers for diagnosis, prognosis, and treatment strategies.These mutations often coexist with specific genetic alterations, such as TP53 mutations in astrocytic tumors and 1p/19q co-deletions in oligodendrogliomas (Cohen et al.). The presence of IDH mutations plays a critical role in treatment strategies and prognostic evaluations for low-grade gliomas, emphasizing their significance in tumorigenesis (Chen et al.).

## 2. Dataset and Preprocessing
This analysis utilizes gene expression data from the TCGA Low-Grade Glioma project (TCGA-LGG). After querying and downloading the transcriptome profiling data, raw count data was extracted for further analysis. Sample IDs were obtained, and IDH status was determined using subtype data.

### 2.1 Data Cleaning and Normalization
Data cleaning involved excluding samples without IDH status and factorizing the remaining samples into Wildtype and Mutant categories. The expression dataset was normalized using the 'TCGAanalyze_Normalization' function, which adjusted for gene length. Filtering was then applied to remove low-expression genes based on a 25th percentile cutoff. This ensured a high-quality dataset for differential expression analysis, enhancing the reliability of the results by minimizing noise.

   ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/heatmap%20of%20LLG%20expresion%20data.png)
                 *Figure: Heatmap of LLG expression data*

## 3. Methodology
### 3.1 Differential Expression Analysis
Differential expression analysis was conducted to compare IDH Wildtype and Mutant samples using a log-fold change cutoff of 2.5 and a false discovery rate of 0.01, resulting in identification of underregulated and overregulated genes in gliomas.

  ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/MA%20plot%203%20For%20DEA%20results.png)
             *Figure: MA plot For Differential Expression Analysis*

  ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/Density%20Plot%203.png)
               *Figure: Density Plot*

  ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/LGG%20Volcano%20Plot3.png)
               *Figure: Volcano Plot For LGG Dataset*

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


   ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/Random%20Forest%20Confusion%20Matrix.png)
                        *Figure: Random Forest Confusion Matrix*

   ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/KNN%20Confusion%20Matrix.png)
              *Figure: KNN Confusion Matrix*

  ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/Variable%20Importace.png)
              *Figure: Variable Importace Plot*

Both models had low overall accuracies and kappas, suggesting significant challenges in predicting glioma clusters based on gene expression data alone, possibly due to high inter-cluster similarity and insufficient discriminatory features.

### 3.3 Clustering Analysis
To further explore the structure within the glioma dataset, we applied K-means clustering with six centers. The silhouette score, which measures how well each point lies within its assigned cluster, was 0.81, indicating that the clustering generally formed well-defined groups. This was confirmed by evaluating the confusion matrices of the models and comparing the clusters to the known classes. Both Random Forest and KNN models exhibited similar misclassification patterns, with sensitivity and specificity values that suggested the difficulty in clearly distinguishing between IDH Wildtype and Mutant gliomas.

 ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/t-SNE%20Plot%20of%20Clustering%20Results.png)
            *Figure: t-SNE Plot of Clustering Results*

## 4. Enrichment Analysis (EA)
EA was carried out with the `TCGAanalyze_EAComplete()` function provided by the TCGAbiolinks package in RStudio. EA revealed significantly enriched biological processes in lists of overexpressed and underexpressed genes. Notable enriched biological processes include:

- **Overexpressed Genes:** GPCR signaling pathway, sensory perception, and neurological system processes.
- **Underexpressed Genes:** Pattern specification process, anterior/posterior pattern formation, and regionalization.

 ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/downreggenepathway_plot%20(1).png)
          *Figure: Lollipop plot of the top pathways enriched with downregulated/underexpressed genes*

 ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/upreggenepathway_plot%20(2).png)
          *Figure : Lollipop plot of the top pathways enriched with upregulated/overexpressed genes*

Top-ranked upregulated genes include the protein-coding gene DEFB119, which plays a role in cell/tissue defense (Martinelli et al., 2020), and long non-coding (lnc) RNA genes LINC01602, LINC03018, which may be involved in transcriptional control (Zhang et al., 2019). The top-ranked downregulated genes include TTR (transthyretin), which is often not expressed in gliomas (Albrecht et al., 1995), and HOXC12, involved in morphological processes (Kawasumi-Kita et al., 2024).

## 5. Comparison with the Reference Paper
The study by Ceccarelli et al. (2016) identified six distinct methylation groups and four RNA expression groups associated with IDH status using unsupervised clustering techniques. Our analysis, while following a similar methodology, did not achieve comparable levels of classification accuracy or clustering performance. Specifically:

- **Cluster Definitions:** Ceccarelli et al. achieved clearer separation between IDH mutant and wildtype gliomas. Our results suggest that the methylation-based features used in our analysis may not be as discriminative as those employed in the original study.
- **Dataset Quality:** The low accuracy in our models might be due to the quality and preprocessing of the data, as well as the differences in the clustering algorithms used.
- **Future Directions:** Additional steps, such as using deeper neural networks or incorporating more advanced feature selection techniques, may be necessary to improve clustering accuracy.

## 6. Discussion
Our study highlights the significant challenges in classifying gliomas into distinct subtypes based solely on gene expression and methylation levels. Although the KNN and Random Forest models provided some insights into the clustering structure of the glioma dataset, the results indicate that more sophisticated algorithms may be required to achieve the desired classification performance.

The comparison with Ceccarelli et al. (2016) reveals important discrepancies in clustering performance, yet both studies observe consistent clustering trends. Our findings underscore the necessity for further research into the molecular subtypes of gliomas. Specifically, optimizing clustering algorithms and exploring additional biological variables will be crucial for better differentiation between IDH statuses.

Furthermore, our analysis indicates that specific genes identified as upregulated or downregulated may serve as potential therapeutic targets. However, the complexity of glioma biology requires an integrative approach that combines gene expression data with other genomic, epigenomic, and clinical parameters. This comprehensive strategy will enhance our understanding of glioma subtypes and inform future therapeutic interventions.

## 7. Conclusion
In conclusion, our study aimed to classify gliomas based on gene expression data, specifically focusing on IDH status. While we achieved some insights into the gene expression profiles, the overall classification accuracy and clustering performance fell short of expectations when compared to the reference study by Ceccarelli et al. (2016). This suggests that our methodology, while valid, may require refinement and the incorporation of more advanced analytical techniques.

## 8. References

1. Ceccarelli, M., Barthel, F. P., Malta, T. M., Sabedot, T. S., Salama, S. R., Murray, B. A., et al. (2016). "Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma." *Cell*, 164(3), 550-563.

2. Chen, H., Luo, Y., Li, C., Zhan, W., Tan, Q., Xie, C., Sharma, A., Sharma, H. S., & Zhang, Z. (2021b). Multimodal imaging in the differential diagnosis of glioma recurrence from treatment-related effects: A protocol for systematic review and network meta-analysis. *Progress in Brain Research*, 377–383.

3. Cohen, A. L., Holmen, S. L., & Colman, H. (2013b). IDH1 and IDH2 Mutations in Gliomas. *Current Neurology and Neuroscience Reports*, 13(5).

4. Albrecht, S., Bayer, T. A., Kraus, J. A. & Pietsch, T. (1995). Transthyretin expression in medulloblastomas and medulloblastoma cell lines. *Neuropathology and Applied Neurobiology*, 21(5), 399-409.

5. Kawasumi-Kita, A., Lee, S. W., Ohtsuka, D., Niimi, K., Asakura, Y., Kitajima, K., Sakane, Y., Tamura, K., Ochi, H., Suzuki, K. I. T. & Morishita, Y. (2024). hoxc12/c13 as key regulators for rebooting the developmental program in Xenopus limb regeneration. *Nature Communications*, 15(1), p. 3340.

6. Martinelli, C., Gabriele, F., Manai, F., Ciccone, R., Novara, F., Sauta, E., Bellazzi, R., Patane, M., Moroni, I., Paterra, R. & Comincini, S. (2020). The search for molecular markers in a gene-orphan case study of a pediatric spinal cord pilocytic astrocytoma. *Cancer Genomics & Proteomics*, 17(2), 117-130.

7. Nilsson, F., Storm, P., Sozzi, E., Hidalgo Gil, D., Birtele, M., Sharma, Y., Parmar, M. & Fiorenzano, A. (2021). Single-cell profiling of coding and noncoding genes in human dopamine neuron differentiation. *Cells*, 10(1), p. 137.

8. Shinoda, T., Ito, H., Sudo, K., Iwamoto, I., Morishita, R. & Nagata, K. I. (2010). Septin 14 is involved in cortical neuronal migration via interaction with Septin 4. *Molecular Biology of the Cell*, 21(8), 1324-1334.

9. Zhang, X., Wang, W., Zhu, W., Dong, J., Cheng, Y., Yin, Z. & Shen, F. (2019). Mechanisms and functions of long non-coding RNAs at multiple regulatory levels. *International Journal of Molecular Sciences*, 20(22), p. 5573.

