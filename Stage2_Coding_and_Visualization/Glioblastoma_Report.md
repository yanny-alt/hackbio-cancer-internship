# Stage 02 Task

## HackBio BioInformatics Internship: Clinical Oncology

### Title
**Gene Expression Analysis and Functional Enrichment of Differentially Expressed Genes**

### Authors
Favour Igwezeke<sup>1</sup>, Jessica Ovabor<sup>2</sup>, Anarghya Hegde<sup>3</sup>, Oluwatobiloba Johnson Osedimilehin<sup>4</sup>, Ogochukwu Nwaigwe<sup>5</sup>, Muhammad Faheem Raziq<sup>6</sup>

### Contributors Information

| Authors | Email | Slack ID |
| ------- | ----- | -------- |
| 1 | beingfave@gmail.com | [Slack ID](https://hackbiointern-leo4437.slack.com/team/U07KE59TWEP) |
| 2 | ovaborjessica85@gmail.com | [Slack ID](https://hackbiointern-leo4437.slack.com/team/U07JVPSU917) |
| 3 | anarghyahegde@outlook.com | [Slack ID](https://hackbiointern-leo4437.slack.com/team/U07JM2UDL7R) |
| 4 | tobijohnson01@gmail.com | [Slack ID](https://hackbiointern-leo4437.slack.com/team/U07JP06QDB4) |
| 5 | nwaigweogochukwu756@gmail.com | [Slack ID](https://hackbiointern-leo4437.slack.com/team/U07KP1D2F24) |
| 6 | faheemraziq1999@gmail.com | [Slack ID](https://hackbiointern-leo4437.slack.com/team/U07KUECLR40) |

## 1. Abstract

The biological pathways that are enriched in a dataset of underregulated genes are 
presented in detail in this study. Important biological pathways using Gene Ontology 
(GO) enrichment analysis were found, including regulation of biosynthetic activities 
and transcription regulation by RNA polymerase II.(Ashburner et 
al., 2000; The Gene Ontology Consortium, 2021)
Understanding the regulatory mechanisms of gene expression and their effects on metabolism, development, and 
cellular differentiation depends on these pathways. This study sheds light on the 
biological relevance of these pathways by revealing the number of genes involved in 
each process and the fold enrichment values. Relevant literature in the domain of 
bioinformatics and molecular biology supports this study’s conclusion.

Dataset utilized for the study: **Glioblastoma gene expression dataset**

## 2. Data Preprocessing

Prior to generating heatmaps, necessary preprocessing steps were conducted:

- **Normalization:** Data was normalized to ensure comparable expression values across genes and samples.
- **Filtering:** Differentially expressed genes were filtered based on fold change and p-values, depending on team-decided threshold values. Genes with a fold change greater than 2 and p-values less than 0.05 were focused on for significant upregulation, and fold change less than -2 for significant downregulation.

**Tools used:** R, tidyverse, dplyr, ShinyGo, and DESeq.

## 3. Heatmap Generation

To visualize the gene expression data, the `heatmap.2()` function from the gplots package in RStudio was utilized. Heatmaps were generated using two different color palettes: diverging and sequential. The purpose of selecting these color schemes was to improve the plots’ readability and facilitate the differentiation of various expression levels across the dataset.

### 3.1 Importance of Color Selection

- **Diverging color palette:** This palette helps in clearly distinguishing between upregulated and downregulated genes. It uses two contrasting colors, making it ideal for highlighting significant changes.
        ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage2_Coding_and_Visualization/visualizations/Image%20of%20divergent%20colour%20palette.jpg)
    *Figure 1: Image of divergent color palette*


- **Sequential color palette:** This palette helps to visualize a gradient of expression changes. It is useful for observing trends where changes are more gradual rather than radically divergent.
       ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage2_Coding_and_Visualization/visualizations/Image%20of%20sequential%20colour%20palette.jpg)
    *Figure 2: Image of sequential color palette*

Both color variants of the heatmap were included to demonstrate the impact of color choices on data interpretation.

## 4. Clustering of Genes and Samples

Three different clustering strategies were explored to investigate gene expression patterns across samples:

- **Clustering genes alone (rows):** Grouping genes that show similar expression patterns across samples.
      ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage2_Coding_and_Visualization/visualizations/Clustering%20gene%20(rows).jpg)
    *Figure 3: Clustering gene (rows)*

- **Clustering samples alone (columns):** Grouping samples with similar expression patterns.
      ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage2_Coding_and_Visualization/visualizations/Clustering%20samples%20(columns).jpg)
    *Figure 4: Clustering samples (columns)*

- **Clustering both genes and samples:** A combined approach to observe how genes and samples cluster together.
       ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage2_Coding_and_Visualization/visualizations/Clustering%20of%20both%20genes%20(rows)%20and%20samples%20(columns).jpg)
    *Figure 5: Clustering of both genes (rows) and samples (columns)*

## 5. Subsetting Genes Based on Expression

To get the adjusted p-values and fold changes, DESeq was performed, and then subsets of genes were identified based on significant upregulation and downregulation using team-decided fold change and adjusted p-value cutoffs; a significance threshold of 0.05 was chosen. To get the fold change, the following formulas were used: 

**Formulars:**  
- Fold change :  Expression level in final condition / Expression level in Initial Condition
- Log2 fold change : log2(final/initial)

- **Upregulated genes:** log2Foldchange > 1 and Adjusted P value < 0.05
- **Downregulated genes:** log2Foldchange < -1 and Adjusted P value < 0.05

## 6. Functional Enrichment Analysis

Functional enrichment analysis was performed using ShinyGO to determine which pathways are overrepresented in the gene expression dataset. The analysis helped identify biological processes associated with significantly downregulated genes, but there were no biological processes that were significantly enriched with the genes in the upregulated set. (Ge et al., 2020)

### 6.1 Enrichment Results

Based on the analysis, the top 5 enriched pathways were identified. A bubble plot was generated to visualize these pathways, showing the number of genes involved and the significance of each pathway. The size of the bubbles corresponds to the negative log10 of the p-value, emphasizing pathway significance.

**Top 5 Enriched Pathways:**

1. Regulation of Transcription by RNA Polymerase II
2. Transcription by RNA Polymerase II
3. Regulation of Biosynthetic Process
4. Nervous System Development
5. Anatomical Structure Morphogenesis

## 7. Detailed Description of Top 3 Enriched Pathways

Here, a detailed description of the top 3 enriched pathways according to biological processes is provided:

### 7.1 Regulation of Transcription by RNA Polymerase II (GO:0006357)

- **Enrichment FDR:** 0.0024
- **Number of Genes Involved:** 21
- **Total Pathway Genes:** 2747
- **Fold Enrichment:** 2.34

One essential biological process that governs the transcription of DNA into messenger RNA (mRNA) is Transcription Regulation by RNA polymerase II. The transcription of genes that code for proteins is carried out by RNA polymerase II, whose regulation is necessary for cell development, differentiation, and responsiveness to external stimuli. The over-representation of this pathway in the dataset implies that transcriptional regulation is a major contributor to the observed gene under-regulation. (Linzer N et al.,  (2021) . Front. Mol. Biosci. 8:681550.)

### 7.2 Transcription by RNA Polymerase II (GO:0006366)

- **Enrichment FDR:** 0.0024
- **Number of Genes Involved:** 21
- **Total Pathway Genes:** 2859
- **Fold Enrichment:** 2.34

The expression of genes and cellular function depend on the Transcription process by RNA polymerase II. The pathway’s over-representation suggests that a large number of the genes in this dataset are crucial for transcriptional processes. The pathway’s significance in maintaining cellular homeostasis and controlling gene expression in response to both external and internal signals is highlighted by the high fold enrichment. (Acker, J et al., (2013).  Nature Reviews Molecular Cell Biology, 14(5), 283-295.) 

### 7.3 Regulation of Biosynthetic Process (GO:0009889)

- **Enrichment FDR:** 0.0024
- **Number of Genes Involved:** 25
- **Total Pathway Genes:** 4545
- **Fold Enrichment:** 2.08

In simpler terms, the process by which cells generate complex chemicals from simpler ones is known as Biosynthesis and it is essential to both metabolism and cellular growth. The regulation of biosynthetic pathways is essential for cellular growth, cellular energy-balance maintenance, and for its responsive function towards external stimuli. The pathway’s enrichment in this dataset implies that under-regulating genes related to biosynthetic activities may have important metabolic repercussions. (Mukherjee, A et al., 2014; Journal of Biological Chemistry 289(10), 6569-6579.)

## 8. Visualization of Enrichment Results

To convey the importance of these pathways, a bubble plot was constructed. Each bubble represents one of the top 5 pathways, with the size of the bubble proportional to the number of genes involved and the color representing the log10 p-value for enrichment significance.
      ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage2_Coding_and_Visualization/visualizations/Bubble%20plot%20of%20top%205%20enriched%20pathways.PNG)
    *Figure 6: Bubble plot of top 5 enriched pathways.*

## 9. Conclusion

Our investigation’s findings have demonstrated the biological significance of several enriched pathways in a dataset of underregulated genes. Our understanding of the control of gene expression, biosynthesis, and cellular differentiation is improved by focusing on biological pathways like the regulation of biosynthetic activities and transcription regulation by RNA polymerase II (Ashburner et al., 2000). These pathways are crucial to regulating metabolism, development, and cellular differentiation, as demonstrated by their high fold enrichment values and the number of genes involved. These findings advance our understanding of the molecular processes driving gene under-regulation and might lead to novel approaches to addressing biological issues.

## 10. References

1. Ashburner, M., Ball, C. A., Blake, J. A., et al. (2000). Gene ontology: Tool for the unification of biology. *Nature Genetics*, 25(1), 25-29. DOI: [10.1038/75556](https://doi.org/10.1038/75556)
2. Mukherjee, A., Rotwein, P. (2014). Transcriptional regulation of biosynthetic processes in skeletal muscle growth. *Journal of Biological Chemistry*, 289(10), 6569-6579. DOI: [10.1074/jbc.M113.523001](https://doi.org/10.1074/jbc.M113.523001)
3. Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L. D. A., François, R., ... & Yutani, H. (2019). Welcome to the Tidyverse. *Journal of Open Source Software*, 4(43), 1686. DOI: [10.21105/joss.01686](https://doi.org/10.21105/joss.01686)
4. Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2023). dplyr: A grammar of data manipulation. R package version 1.1.2. Computer software.
5. Ge, S. X., Jung, D., & Yao, R. (2020). ShinyGO: A graphical gene-set enrichment tool for animals and plants. *Bioinformatics*, 36(8), 2628-2629. DOI: [10.1093/bioinformatics/btz931](https://doi.org/10.1093/bioinformatics/btz931)
6. Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 1-21. DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
7. Acker, J., Conaway, R. C., & Conaway, J. W. (2013). The RNA polymerase II machinery: Insights into structure, function, and regulation. *Nature Reviews Molecular Cell Biology*, 14(5), 283-295. DOI: [10.1038/nrm3554](https://doi.org/10.1038/nrm3554)
8. Linzer, N., Trumbull, A., Nar, R., Gibbons, M. D., Yu, D. T., Strouboulis, J., & Bungert, J. (2021). Regulation of RNA Polymerase II Transcription Initiation and Elongation by Transcription Factor TFII-I. *Frontiers in Molecular Biosciences*, 8:681550. DOI: [10.3389/fmolb.2021.681550](https://doi.org/10.3389/fmolb.2021.681550)

