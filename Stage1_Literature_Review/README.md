# Stage1_Literature_Review

This directory contains all materials related to the literature review conducted during the Bioinformatics Internship program. The literature review focuses on the analysis of RNA-Seq transcriptome data from TCGA, investigating gene expression changes across various cancer types.

## Files and Directories

- **README.md**: Provides an overview of the literature review stage and its contents.
- **literature_review.md**: The full literature review document, including study overview, methods, results, and discussion.
- **references/**: Directory containing references and citations used in the literature review.

## Video Presentation

- **Video Title**: Large-scale RNA-Seq Transcriptome Analysis of 4043 Cancers and 548 Normal Tissue Controls across 12 TCGA Cancer Types
- **Summary**: A presentation video summarizing the findings of the literature review. [link to the video summary](https://www.loom.com/share/d1dc6145ee1a4672819c460b2b7c4900?sid=a43cac16-836f-4822-80cb-bf3b7b7bb81e)
- **Contributors**: Muhammad Faheem Raziq, Favour Igwezeke, Josiah Isong, Nnadiekwe Chigozie

## Study Overview

### Introduction
- The study utilizes large-scale TCGA RNA-seq data from 4043 cancer samples and 548 normal tissue samples across 12 cancer types.
- It aims to investigate global gene expression changes and develop predictive gene signatures to improve cancer diagnosis and therapy.

### Methods and Results

- **Data Sets**: RNA-seq data obtained from TCGA Data Portal, focusing on UNC (IlluminaHiSeq_RNASeqV2).
- **Differential Expression Analysis**: Conducted using edgeR Bioconductor package, revealing universally up and downregulated genes across cancer types.
- **Clustering and Gene Set Association Analysis**: Performed using APCluster and GSAASeqSP, identifying gene clusters and significant gene sets.
- **Pathway Enrichment and Disease Association Analysis**: Conducted using WebGestalt to explore biological relevance and disease associations.

### Results
- Identified seven cross-cancer gene signatures and a 14-gene signature with high predictive accuracy across several cancer types.

### Discussion and Conclusion
- The study highlights the potential of transcriptome analysis in identifying key gene signatures for cancer diagnosis and treatment.

## References

- Anders, S., & Huber, W. (2010). Differential expression analysis for sequence count data. Genome Biology, 11(10), R106. [https://doi.org/10.1186/gb-2010-11-10-r106](https://doi.org/10.1186/gb-2010-11-10-r106)
- Peng, L., et al. (2015). Large-scale RNA-Seq Transcriptome Analysis of 4043 Cancers and 548 Normal Tissue Controls across 12 TCGA Cancer Types. Scientific Reports, 5(1), 13413. [https://doi.org/10.1038/srep13413](https://doi.org/10.1038/srep13413)
- Wang, Z., et al. (2009). RNA-Seq: A revolutionary tool for transcriptomics. Nature Reviews Genetics, 10(1), 57-63. [https://doi.org/10.1038/nrg2484](https://doi.org/10.1038/nrg2484)
- The Cancer Genome Atlas (TCGA). (n.d.). Retrieved from [http://cancergenome.nih.gov](http://cancergenome.nih.gov)
- UCSC Cancer Genomics Hub. (n.d.). Retrieved from [https://cghub.ucsc.edu](https://cghub.ucsc.edu)

