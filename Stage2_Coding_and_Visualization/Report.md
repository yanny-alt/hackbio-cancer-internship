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

## 1. Overview

In this analysis, a gene expression dataset was processed, visualized, and interpreted to gain insights into the biological importance of differential expression. To investigate the pathways associated with gene expression under specific conditions, heatmaps were generated, clustering was performed, and functional enrichment analysis was conducted. This analysis was completed by a team of 3 data scientists and 3 biomarker hunters.

Dataset utilized for the study: **Glioblastoma gene expression dataset**

The top 500+ differentially expressed genes are included in this dataset, and the goal of the study was to visualize the data, find significant gene clusters, and employ enrichment analysis to determine the biological significance of these genes.

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

To get the p-values and fold changes, DESeq was performed, and then subsets of genes were identified based on significant upregulation and downregulation using team-decided fold change and p-value cutoffs; a significance threshold of 0.05 was chosen.

- **Upregulated genes:** log2Foldchange > 1 and P value < 0.05
- **Downregulated genes:** log2Foldchange < -1 and P value < 0.05

## 6. Functional Enrichment Analysis

Functional enrichment analysis was performed using ShinyGO to determine which pathways are overrepresented in the gene expression dataset. The analysis helped identify biological processes associated with significantly downregulated genes, but there were no biological processes that were significantly enriched with the genes in the upregulated set.

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

One essential biological process that governs the transcription of DNA into messenger RNA (mRNA) is Transcription Regulation by RNA polymerase II. The transcription of genes that code for proteins is carried out by RNA polymerase II, whose regulation is necessary for cell development, differentiation, and responsiveness to external stimuli. The over-representation of this pathway in the dataset implies that transcriptional regulation is a major contributor to the observed gene under-regulation.

### 7.2 Transcription by RNA Polymerase II (GO:0006366)

- **Enrichment FDR:** 0.0024
- **Number of Genes Involved:** 21
- **Total Pathway Genes:** 2859
- **Fold Enrichment:** 2.34

The expression of genes and cellular function depend on the Transcription process by RNA polymerase II. The pathway’s over-representation suggests that a large number of the genes in this dataset are crucial for transcriptional processes. The pathway’s significance in maintaining cellular homeostasis and controlling gene expression in response to both external and internal signals is highlighted by the high fold enrichment.

### 7.3 Regulation of Biosynthetic Process (GO:0009889)

- **Enrichment FDR:** 0.0024
- **Number of Genes Involved:** 25
- **Total Pathway Genes:** 4545
- **Fold Enrichment:** 2.08

In simpler terms, the process by which cells generate complex chemicals from simpler ones is known as Biosynthesis and it is essential to both metabolism and cellular growth. The regulation of biosynthetic pathways is essential for cellular growth, cellular energy-balance maintenance, and for its responsive function towards external stimuli. The pathway’s enrichment in this dataset implies that under-regulating genes related to biosynthetic activities may have important metabolic repercussions.

## 8. Visualization of Enrichment Results

To convey the importance of these pathways, a bubble plot was constructed. Each bubble represents one of the top 5 pathways, with the size of the bubble proportional to the number of genes involved and the color representing the log10 p-value for enrichment significance.
      ![](https://github.com/yanny-alt/hackbio-cancer-internship/blob/main/Stage2_Coding_and_Visualization/visualizations/Bubble%20plot%20of%20top%205%20enriched%20pathways.PNG)
    *Figure 6: Bubble plot of top 5 enriched pathways.*

## 9. Conclusion

The enrichment analysis conducted in this study reveals the critical roles of biological processes such as transcriptional regulation and biosynthesis in the dataset of the underregulated genes. The identified pathways are heavily involved in maintaining cellular homeostasis, development, and metabolism. Their dysregulation may have broad effects in domains like metabolic diseases and neurobiology. FOXP2, HIF3A, and GBX2 are among the genes found in these pathways that are critical for functions ranging from transcription to neurogenesis. When it comes to developmental disorders and metabolic dysfunctions, more research into the functional implications of these pathways may yield important insights into disease mechanisms and causes.
