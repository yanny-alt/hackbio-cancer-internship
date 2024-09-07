**Large-scale RNA-Seq Transcriptome Analysis of 4043 Cancers and 548 Normal **

**Tissue Controls across 12 TCGA Cancer Types**

**Muhammad Faheem Raziq, Favour Igwezeke, Josiah Isong, Nnadiekwe Chigozie **

**LinkedIn/Twitter link:**

**Study Overview2**

**Introduction**

**Through such large-scale databases such as TCGA, cancer genomics has brought about profound **

**knowledge  in the general character of cancer. In this study, we use TCGA RNA-seq data **

**containing 4043 cancer and 548 normal tissues from 12 cancer types to investigate global gene **

**expression changes and develop predictive gene signatures (Wang et al., 2009). The focus was **

**these transcriptomes, scrutinize the gene expressions, and establish signatures for improving the **

**cancer diagnosis and therapy (Peng et al., 2015).**

**Methods and Results **

**Data Sets**

**We  obtained  transcriptome  and  clinical  data  from  the  TCGA  Data  Portal  (https\://tcga- **

**data.nci.nih.gov/tcga/tcgaHome2.jsp), focusing on UNC (IlluminaHiSeq\_RNASeqV2) to avoid **

**sequencing platform variability. Twelve cancer types with both cancerous and normal tissue data **

**were selected using "primary tumor" and "solid tissue normal" samples (Table 4).3**

**Differential Expression Analysis**

**Differential  expression  analysis  was  conducted  using  the  edgeR  Bioconductor  package **

**(http\://www\.bioconductor.org/packages/release/bioc/html/edgeR.html). This can be achieved by a **

**comparison of gene expression in primary tumors and solid tissue normal for each cancer type. **

**Robust changes in gene expression were found which showed that there are indeed universally **

**down and up regulated genes across various cancer types (Anders & Huber, 2010).**

**Clustering and Gene Set Association Analysis**

**Gene clustering was performed using APCluster (http\://www\.bioinf.jku.at/software/apcluster/), **

**with  DESeq-normalized  RNA-Seq  data.  The  Pearson  correlation  coefficient  measured  gene **

**similarity, and genes in the same cluster showed highly correlated expression profiles. Gene set **

**association  analysis  was  performed  with  GSAASeqSP  (http\://gsaa.unc.edu),  identifying **

**significant gene sets with an FDR < 0.15.**

**Pathway Enrichment and Disease Association Analysis**

**To explore biological relevance, we conducted pathway enrichment and disease association **

**analyses using WebGestalt (http\://bioinfo.vanderbilt.edu/webgestalt/). GO, KEGG, and Pathway **

**Commons analyses identified over-represented GO terms and pathways, linking gene sets to **

**disease associations.**

**Results**

**Seven cross-cancer gene signatures were consistently altered across cancers. A 14-gene signature **

**showed high predictive accuracy, with rates ranging from 88.17% to 99.10% across seven cancer **

**types. A lung cancer-specific signature achieved 95.68% accuracy in TCGA data and 100% in **

**GSE5364 data.4**

**Discussion**

**Using a collection of ± 6,000 RNA-Seq experiments, 322 cancer-associated genes share conserved **

**upregulation across cancers and 127 cancer-associated genes display unique upregulation across **

**cancers. The presented gene signature provides opportunities for general diagnostic markers and **

**pointed treatments. Such high predictive accuracy of these molecular signatures indicates their **

**potential to greatly aid clinical diagnoses of cancer types and classification thereof.**

**Conclusion**

**This study underscores the power of transcriptome analysis in identifying key gene signatures that **

**differentiate cancerous from normal tissues. These findings enhance our understanding of cancer **

**biology and offer pathways for precise diagnostic and therapeutic strategies.**

**References**

**Anders, S., & Huber, W. (2010). Differential expression analysis for sequence count data. Genome **

**Biology, 11(10), R106. https\://doi.org/10.1186/gb-2010-11-10-r106**

**Peng, L., Bian, X. W., Li, D. K., Xu, C., Wang, G. M., Xia, Q. Y., & Xiong, Q. (2015). Large- **

**scale RNA-Seq Transcriptome Analysis of 4043 Cancers and 548 Normal Tissue Controls **

**12**

**Reports,**

**13413.**

**TCGA **

**Cancer **

**Types. **

**across **

**Scientific **

**5(1), **

**https\://doi.org/10.1038/srep13413**

**Wang, Z., Gerstein, M., & Snyder, M. (2009). RNA-Seq: A revolutionary tool for transcriptomics. **

**Nature Reviews. Genetics, 10(1), 57–63. https\://doi.org/10.1038/nrg2484**

**The Cancer Genome Atlas (TCGA). (n.d.). Retrieved from http\://cancergenome.nih.gov **

**UCSC Cancer Genomics Hub. (n.d.). Retrieved from https\://cghub.ucsc.edu**
