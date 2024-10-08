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
This study aimed to reproduce unsupervised clustering of low-grade gliomas (LGG) based on expression levels to classify  IDH statuses, utilizing the K-Nearest Neighbors (KNN) algorithm. We analyzed expression data from 516 LGG samples in The Cancer Genome Atlas (TCGA) and compared our results with the findings of Ceccarelli et al. (2016). While KNN did not achieve high classification accuracy, it provided insights into factors influencing glioma classification, indicating a need for more advanced techniques.

## 1. Introduction
Gliomas, arising from glial cells, are the most common primary brain tumors and are challenging to treat. IDH mutations significantly impact gliomagenesis and serve as vital biomarkers for diagnosis, prognosis, and treatment strategies (Cohen et al., 2013). Recent molecular profiling studies suggest that methylation and RNA expression data can classify gliomas into distinct subtypes with prognostic implications (Ceccarelli et al., 2016).


## 2. Dataset and Preprocessing
We used gene expression data from the TCGA Low-Grade Glioma project. After data cleaning, IDH status was categorized into Wildtype and Mutant. Normalization was performed to enhance data quality.

   ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/heatmap%20of%20LLG%20expresion%20data.png)
                 *Figure: Heatmap of LLG expression data*

## 3. Methodology
### 3.1 Differential Expression Analysis
Differential expression analysis identified underregulated and overregulated genes using a log-fold change cutoff of 2.5 and a false discovery rate of 0.01.

  ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/MA%20plot%203%20For%20DEA%20results.png)
             *Figure: MA plot For Differential Expression Analysis*

  ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/Density%20Plot%203.png)
               *Figure: Density Plot*

  ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/LGG%20Volcano%20Plot3.png)
               *Figure: Volcano Plot For LGG Dataset*

### 3.2 Random Forest and KNN Classification
Both Random Forest and KNN models yielded low overall accuracies (17.15% for Random Forest; 17.32% for KNN). This suggests challenges in predicting glioma clusters due to high inter-cluster similarity and insufficient features.

   ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/Random%20Forest%20Confusion%20Matrix.png)
                        *Figure: Random Forest Confusion Matrix*

   ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/KNN%20Confusion%20Matrix.png)
              *Figure: KNN Confusion Matrix*

  ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/Variable%20Importace.png)
              *Figure: Variable Importace Plot*



### 3.3 Clustering Analysis
K-means clustering with six centers showed a silhouette score of 0.81, indicating well-defined groups. However, misclassification patterns were noted, emphasizing difficulties in distinguishing IDH statuses.
 ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/t-SNE%20Plot%20of%20Clustering%20Results.png)
            *Figure: t-SNE Plot of Clustering Results*

## 4. Enrichment Analysis (EA)
EA identified enriched biological processes linked to overexpressed and underexpressed genes, highlighting potential therapeutic targets.

- **Overexpressed Genes:** GPCR signaling pathway, sensory perception, and neurological system processes.
- **Underexpressed Genes:** Pattern specification process, anterior/posterior pattern formation, and regionalization.

 ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/downreggenepathway_plot%20(1).png)
          *Figure: Lollipop plot of the top pathways enriched with downregulated/underexpressed genes*

 ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage4_Pipeline_Implementation_II/visualizations/upreggenepathway_plot%20(2).png)
          *Figure : Lollipop plot of the top pathways enriched with upregulated/overexpressed genes*

Top-ranked upregulated genes include the protein-coding gene DEFB119, which plays a role in cell/tissue defense (Martinelli et al., 2020), and long non-coding (lnc) RNA genes LINC01602, LINC03018, which may be involved in transcriptional control (Zhang et al., 2019). The top-ranked downregulated genes include TTR (transthyretin), which is often not expressed in gliomas (Albrecht et al., 1995), and HOXC12, involved in morphological processes (Kawasumi-Kita et al., 2024).

## 5. Comparison with the Reference Paper
Our results did not match the classification accuracy of Ceccarelli et al. (2016), indicating possible issues with dataset quality and feature discriminability.

## 6. Discussion
This study reveals the challenges of classifying gliomas based on gene expression alone. While KNN and Random Forest provided insights, more sophisticated algorithms are necessary for improved accuracy. Future research should integrate multiple data types and advanced methods to enhance understanding of glioma subtypes.

## 7. Conclusion
Our findings underscore the need for refined methodologies in glioma classification and the exploration of molecular targets for therapeutic intervention, highlighting the complexity of glioma biology.

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

