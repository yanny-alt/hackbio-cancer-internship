**Large-scale RNA-Seq Transcriptome Analysis of 4043 Cancers and 548 Normal Tissue Controls across 12 TCGA Cancer Types**

Muhammad Faheem Raziq, Favour Igwezeke, Josiah Isong, Nnadiekwe Chigozie

LinkedIn/Twitter link:

**_Study Overview_**

** **

<!--[if gte vml 1]><v:shapetype id="_x0000_t75" coordsize="21600,21600"
 o:spt="75" o:preferrelative="t" path="m@4@5l@4@11@9@11@9@5xe" filled="f"
 stroked="f">
 <v:stroke joinstyle="miter"/>
 <v:formulas>
  <v:f eqn="if lineDrawn pixelLineWidth 0"/>
  <v:f eqn="sum @0 1 0"/>
  <v:f eqn="sum 0 0 @1"/>
  <v:f eqn="prod @2 1 2"/>
  <v:f eqn="prod @3 21600 pixelWidth"/>
  <v:f eqn="prod @3 21600 pixelHeight"/>
  <v:f eqn="sum @0 0 1"/>
  <v:f eqn="prod @6 1 2"/>
  <v:f eqn="prod @7 21600 pixelWidth"/>
  <v:f eqn="sum @8 21600 0"/>
  <v:f eqn="prod @7 21600 pixelHeight"/>
  <v:f eqn="sum @10 21600 0"/>
 </v:formulas>
 <v:path o:extrusionok="f" gradientshapeok="t" o:connecttype="rect"/>
 <o:lock v:ext="edit" aspectratio="t"/>
</v:shapetype><v:shape id="Picture_x0020_1" o:spid="_x0000_i1027" type="#_x0000_t75"
 style='width:469.5pt;height:396pt;visibility:visible;mso-wrap-style:square'>
 <v:imagedata src="file:///C:/Users/Yanny/AppData/Local/Packages/oice_16_974fa576_32c1d314_122a/AC/Temp/msohtmlclip1/01/clip_image001.png"
  o:title=""/>
</v:shape><![endif]--><!--[if !vml]-->![](file:///C:/Users/Yanny/AppData/Local/Packages/oice_16_974fa576_32c1d314_122a/AC/Temp/msohtmlclip1/01/clip_image002.gif)<!--[endif]-->

 

 

 

 


# <a id="_974ajmee7jae"></a>**Introduction**

Through such large-scale databases such as TCGA, cancer genomics has brought about profound knowledge in the general character of cancer. In this study, we use TCGA RNA-seq data containing 4043 cancer and 548 normal tissues from 12 cancer types to investigate global gene expression changes and develop predictive gene signatures <!--[if supportFields]><span
lang=EN style='font-size:12.0pt;line-height:150%;font-family:"Times New Roman",serif;
background:white;mso-highlight:white'><span style='mso-element:field-begin'></span><span
style='mso-spacerun:yes'> </span>ADDIN ZOTERO_ITEM CSL_CITATION
{&quot;citationID&quot;:&quot;RJxe65at&quot;,&quot;properties&quot;:{&quot;formattedCitation&quot;:&quot;(Wang
et al., 2009)&quot;,&quot;plainCitation&quot;:&quot;(Wang et al., 2009)&quot;,&quot;noteIndex&quot;:0},&quot;citationItems&quot;:[{&quot;id&quot;:1118,&quot;uris&quot;:[&quot;http://zotero.org/users/local/GZBaIJVt/items/76BGZYXU&quot;],&quot;itemData&quot;:{&quot;id&quot;:1118,&quot;type&quot;:&quot;article-journal&quot;,&quot;abstract&quot;:&quot;RNA-Seq
is a recently developed approach to transcriptome profiling that uses
deep-sequencing technologies. Studies using this method have already altered
our view of the extent and complexity of eukaryotic transcriptomes. RNA-Seq
also provides a far more precise measurement of levels of transcripts and their
isoforms than other methods. This article describes the RNA-Seq approach, the
challenges associated with its application, and the advances made so far in
characterizing several eukaryote
transcriptomes.&quot;,&quot;container-title&quot;:&quot;Nature Reviews.
Genetics&quot;,&quot;DOI&quot;:&quot;10.1038/nrg2484&quot;,&quot;ISSN&quot;:&quot;1471-0064&quot;,&quot;issue&quot;:&quot;1&quot;,&quot;journalAbbreviation&quot;:&quot;Nat
Rev
Genet&quot;,&quot;language&quot;:&quot;eng&quot;,&quot;note&quot;:&quot;PMID:
19015660\nPMCID:
PMC2949280&quot;,&quot;page&quot;:&quot;57-63&quot;,&quot;source&quot;:&quot;PubMed&quot;,&quot;title&quot;:&quot;RNA-Seq:
a revolutionary tool for
transcriptomics&quot;,&quot;title-short&quot;:&quot;RNA-Seq&quot;,&quot;volume&quot;:&quot;10&quot;,&quot;author&quot;:[{&quot;family&quot;:&quot;Wang&quot;,&quot;given&quot;:&quot;Zhong&quot;},{&quot;family&quot;:&quot;Gerstein&quot;,&quot;given&quot;:&quot;Mark&quot;},{&quot;family&quot;:&quot;Snyder&quot;,&quot;given&quot;:&quot;Michael&quot;}],&quot;issued&quot;:{&quot;date-parts&quot;:[[&quot;2009&quot;,1]]}}}],&quot;schema&quot;:&quot;https://github.com/citation-style-language/schema/raw/master/csl-citation.json&quot;}
<span style='mso-element:field-separator'></span></span><![endif]-->(Wang et al., 2009)<!--[if supportFields]><span lang=EN style='font-size:12.0pt;
line-height:150%;font-family:"Times New Roman",serif;background:white;
mso-highlight:white'><span style='mso-element:field-end'></span></span><![endif]-->. The focus was these transcriptomes, scrutinize the gene expressions, and establish signatures for improving the cancer diagnosis and therapy <!--[if supportFields]><span lang=EN
style='font-size:12.0pt;line-height:150%;font-family:"Times New Roman",serif;
background:white;mso-highlight:white'><span style='mso-element:field-begin'></span><span
style='mso-spacerun:yes'> </span>ADDIN ZOTERO_ITEM CSL_CITATION
{&quot;citationID&quot;:&quot;WQ4QdKbu&quot;,&quot;properties&quot;:{&quot;formattedCitation&quot;:&quot;(Peng
et al., 2015)&quot;,&quot;plainCitation&quot;:&quot;(Peng et al.,
2015)&quot;,&quot;noteIndex&quot;:0},&quot;citationItems&quot;:[{&quot;id&quot;:1116,&quot;uris&quot;:[&quot;http://zotero.org/users/local/GZBaIJVt/items/AUJ8CNKW&quot;],&quot;itemData&quot;:{&quot;id&quot;:1116,&quot;type&quot;:&quot;article-journal&quot;,&quot;abstract&quot;:&quot;The
Cancer Genome Atlas (TCGA) has accrued RNA-Seq-based transcriptome data for
more than 4000 cancer tissue samples across 12 cancer types, translating these
data into biological insights remains a major challenge. We analyzed and
compared the transcriptomes of 4043 cancer and 548 normal tissue samples from
21 TCGA cancer types and created a comprehensive catalog of gene expression
alterations for each cancer type. By clustering genes into co-regulated gene
sets, we identified seven cross-cancer gene signatures altered across a diverse
panel of primary human cancer samples. A 14-gene signature extracted from these
seven cross-cancer gene signatures precisely differentiated between cancerous
and normal samples, the predictive accuracy of leave-one-out cross-validation
(LOOCV) were 92.04%, 96.23%, 91.76%, 90.05%, 88.17%, 94.29% and 99.10% for
BLCA, BRCA, COAD, HNSC, LIHC, LUAD and LUSC, respectively. A lung
cancer-specific gene signature, containing SFTPA1 and SFTPA2 genes, accurately
distinguished lung cancer from other cancer samples, the predictive accuracy of
LOOCV for TCGA and GSE5364 data were 95.68% and 100%, respectively. These gene
signatures provide rich insights into the transcriptional programs that trigger
tumorigenesis and metastasis and many genes in the signature gene panels may be
of significant value to the diagnosis and treatment of
cancer.&quot;,&quot;container-title&quot;:&quot;Scientific
Reports&quot;,&quot;DOI&quot;:&quot;10.1038/srep13413&quot;,&quot;ISSN&quot;:&quot;2045-2322&quot;,&quot;issue&quot;:&quot;1&quot;,&quot;journalAbbreviation&quot;:&quot;Sci
Rep&quot;,&quot;language&quot;:&quot;en&quot;,&quot;license&quot;:&quot;2015
The Author(s)&quot;,&quot;note&quot;:&quot;publisher: Nature Publishing
Group&quot;,&quot;page&quot;:&quot;13413&quot;,&quot;source&quot;:&quot;www.nature.com&quot;,&quot;title&quot;:&quot;Large-scale
RNA-Seq Transcriptome Analysis of 4043 Cancers and 548 Normal Tissue Controls
across 12 TCGA Cancer Types&quot;,&quot;volume&quot;:&quot;5&quot;,&quot;author&quot;:[{&quot;family&quot;:&quot;Peng&quot;,&quot;given&quot;:&quot;Li&quot;},{&quot;family&quot;:&quot;Bian&quot;,&quot;given&quot;:&quot;Xiu
Wu&quot;},{&quot;family&quot;:&quot;Li&quot;,&quot;given&quot;:&quot;Di
Kang&quot;},{&quot;family&quot;:&quot;Xu&quot;,&quot;given&quot;:&quot;Chuan&quot;},{&quot;family&quot;:&quot;Wang&quot;,&quot;given&quot;:&quot;Guang
Ming&quot;},{&quot;family&quot;:&quot;Xia&quot;,&quot;given&quot;:&quot;Qing
You&quot;},{&quot;family&quot;:&quot;Xiong&quot;,&quot;given&quot;:&quot;Qing&quot;}],&quot;issued&quot;:{&quot;date-parts&quot;:[[&quot;2015&quot;,8,21]]}}}],&quot;schema&quot;:&quot;https://github.com/citation-style-language/schema/raw/master/csl-citation.json&quot;}
<span style='mso-element:field-separator'></span></span><![endif]-->(Peng et al., 2015)<!--[if supportFields]><span lang=EN style='font-size:12.0pt;
line-height:150%;font-family:"Times New Roman",serif;background:white;
mso-highlight:white'><span style='mso-element:field-end'></span></span><![endif]-->.


# <a id="_m9n8g3t9p6ln"></a>**Methods and Results<a id="_3r7bkm5ri35f"></a>**

**Data Sets**

We obtained transcriptome and clinical data from the TCGA Data Portal (<https://tcga-data.nci.nih.gov/tcga/tcgaHome2.jsp>), focusing on UNC (IlluminaHiSeq\_RNASeqV2) to avoid sequencing platform variability. Twelve cancer types with both cancerous and normal tissue data were selected using "primary tumor" and "solid tissue normal" samples (Table 4).

**<!--[if gte vml 1]><v:shape id="Picture_x0020_2" o:spid="_x0000_i1026"
 type="#_x0000_t75" style='width:486pt;height:324pt;visibility:visible;
 mso-wrap-style:square' o:bordertopcolor="yellow pure" o:borderleftcolor="yellow pure"
 o:borderbottomcolor="yellow pure" o:borderrightcolor="yellow pure">
 <v:imagedata src="file:///C:/Users/Yanny/AppData/Local/Packages/oice_16_974fa576_32c1d314_122a/AC/Temp/msohtmlclip1/01/clip_image003.png"
  o:title="" croptop="220f" cropleft="145f" cropright="211f"/>
 <w:bordertop type="single" width="2"/>
 <w:borderleft type="single" width="2"/>
 <w:borderbottom type="single" width="2"/>
 <w:borderright type="single" width="2"/>
</v:shape><![endif]--><!--[if !vml]-->![](file:///C:/Users/Yanny/AppData/Local/Packages/oice_16_974fa576_32c1d314_122a/AC/Temp/msohtmlclip1/01/clip_image004.gif)<!--[endif]-->**

**Differential Expression Analysis**

Differential expression analysis was conducted using the edgeR Bioconductor package (<http://www.bioconductor.org/packages/release/bioc/html/edgeR.html>). This can be achieved by a comparison of gene expression in primary tumors and solid tissue normal for each cancer type. Robust changes in gene expression were found which showed that there are indeed universally down and up regulated genes across various cancer types <!--[if supportFields]><span
lang=EN style='font-size:12.0pt;line-height:150%;font-family:"Times New Roman",serif;
background:white;mso-highlight:white'><span style='mso-element:field-begin'></span><span
style='mso-spacerun:yes'> </span>ADDIN ZOTERO_ITEM CSL_CITATION
{&quot;citationID&quot;:&quot;NRNMV1m1&quot;,&quot;properties&quot;:{&quot;formattedCitation&quot;:&quot;(Anders
&amp; Huber, 2010)&quot;,&quot;plainCitation&quot;:&quot;(Anders &amp; Huber,
2010)&quot;,&quot;noteIndex&quot;:0},&quot;citationItems&quot;:[{&quot;id&quot;:1121,&quot;uris&quot;:[&quot;http://zotero.org/users/local/GZBaIJVt/items/CDR3QP26&quot;],&quot;itemData&quot;:{&quot;id&quot;:1121,&quot;type&quot;:&quot;article-journal&quot;,&quot;abstract&quot;:&quot;High-throughput
sequencing assays such as RNA-Seq, ChIP-Seq or barcode counting provide
quantitative readouts in the form of count data. To infer differential signal
in such data correctly and with good statistical power, estimation of data
variability throughout the dynamic range and a suitable error model are
required. We propose a method based on the negative binomial distribution, with
variance and mean linked by local regression and present an implementation,
DESeq, as an R/Bioconductor package.&quot;,&quot;container-title&quot;:&quot;Genome
Biology&quot;,&quot;DOI&quot;:&quot;10.1186/gb-2010-11-10-r106&quot;,&quot;ISSN&quot;:&quot;1474-760X&quot;,&quot;issue&quot;:&quot;10&quot;,&quot;journalAbbreviation&quot;:&quot;Genome
Biol&quot;,&quot;language&quot;:&quot;eng&quot;,&quot;note&quot;:&quot;PMID:
20979621\nPMCID:
PMC3218662&quot;,&quot;page&quot;:&quot;R106&quot;,&quot;source&quot;:&quot;PubMed&quot;,&quot;title&quot;:&quot;Differential
expression analysis for sequence count
data&quot;,&quot;volume&quot;:&quot;11&quot;,&quot;author&quot;:[{&quot;family&quot;:&quot;Anders&quot;,&quot;given&quot;:&quot;Simon&quot;},{&quot;family&quot;:&quot;Huber&quot;,&quot;given&quot;:&quot;Wolfgang&quot;}],&quot;issued&quot;:{&quot;date-parts&quot;:[[&quot;2010&quot;]]}}}],&quot;schema&quot;:&quot;https://github.com/citation-style-language/schema/raw/master/csl-citation.json&quot;}
<span style='mso-element:field-separator'></span></span><![endif]-->(Anders & Huber, 2010)<!--[if supportFields]><span lang=EN style='font-size:
12.0pt;line-height:150%;font-family:"Times New Roman",serif;background:white;
mso-highlight:white'><span style='mso-element:field-end'></span></span><![endif]-->.


# <a id="_72mafn40wzad"></a>**Clustering and Gene Set Association Analysis**

Gene clustering was performed using APCluster (<http://www.bioinf.jku.at/software/apcluster/>), with DESeq-normalized RNA-Seq data. The Pearson correlation coefficient measured gene similarity, and genes in the same cluster showed highly correlated expression profiles. Gene set association analysis was performed with GSAASeqSP ([http://gsaa.unc.edu](http://gsaa.unc.edu/)), identifying significant gene sets with an FDR < 0.15.


# <a id="_l3b1063cpzry"></a>**Pathway Enrichment and Disease Association Analysis**

To explore biological relevance, we conducted pathway enrichment and disease association analyses using WebGestalt (<http://bioinfo.vanderbilt.edu/webgestalt/>). GO, KEGG, and Pathway Commons analyses identified over-represented GO terms and pathways, linking gene sets to disease associations.&#x20;


# <a id="_e1zwptfpz03"></a>**Results**

Seven cross-cancer gene signatures were consistently altered across cancers. A 14-gene signature showed high predictive accuracy, with rates ranging from 88.17% to 99.10% across seven cancer types. A lung cancer-specific signature achieved 95.68% accuracy in TCGA data and 100% in GSE5364 data.

<!--[if gte vml 1]><v:shape
 id="Picture_x0020_4" o:spid="_x0000_i1025" type="#_x0000_t75" style='width:249.75pt;
 height:2in;visibility:visible;mso-wrap-style:square'>
 <v:imagedata src="file:///C:/Users/Yanny/AppData/Local/Packages/oice_16_974fa576_32c1d314_122a/AC/Temp/msohtmlclip1/01/clip_image005.png"
  o:title=""/>
</v:shape><![endif]--><!--[if !vml]-->![](file:///C:/Users/Yanny/AppData/Local/Packages/oice_16_974fa576_32c1d314_122a/AC/Temp/msohtmlclip1/01/clip_image006.gif)<!--[endif]-->


# <a id="_loe4335jqamv"></a>**Discussion**

Using a collection of ± 6,000 RNA-Seq experiments, 322 cancer-associated genes share conserved upregulation across cancers and 127 cancer-associated genes display unique upregulation across cancers. The presented gene signature provides opportunities for general diagnostic markers and pointed treatments. Such high predictive accuracy of these molecular signatures indicates their potential to greatly aid clinical diagnoses of cancer types and classification thereof.


# <a id="_cty35lyo662o"></a>**Conclusion**

This study underscores the power of transcriptome analysis in identifying key gene signatures that differentiate cancerous from normal tissues. These findings enhance our understanding of cancer biology and offer pathways for precise diagnostic and therapeutic strategies.


# <a id="_8vf3obn5bcfe"></a>**References**

<!--[if supportFields]><span
lang=EN style='font-size:12.0pt;mso-bidi-font-size:11.0pt;line-height:150%;
font-family:"Times New Roman",serif'><span style='mso-element:field-begin'></span><span
style='mso-spacerun:yes'> </span>ADDIN ZOTERO_BIBL
{&quot;uncited&quot;:[],&quot;omitted&quot;:[],&quot;custom&quot;:[]}
CSL_BIBLIOGRAPHY <span style='mso-element:field-separator'></span></span><![endif]-->Anders, S., & Huber, W. (2010). Differential expression analysis for sequence count data. _Genome Biology_, _11_(10), R106. https\://doi.org/10.1186/gb-2010-11-10-r106

Peng, L., Bian, X. W., Li, D. K., Xu, C., Wang, G. M., Xia, Q. Y., & Xiong, Q. (2015). Large-scale RNA-Seq Transcriptome Analysis of 4043 Cancers and 548 Normal Tissue Controls across 12 TCGA Cancer Types. _Scientific Reports_, _5_(1), 13413. https\://doi.org/10.1038/srep13413

Wang, Z., Gerstein, M., & Snyder, M. (2009). RNA-Seq: A revolutionary tool for transcriptomics. _Nature Reviews. Genetics_, _10_(1), 57–63. https\://doi.org/10.1038/nrg2484

<!--[if supportFields]><span
lang=EN style='font-size:12.0pt;mso-bidi-font-size:11.0pt;line-height:150%;
font-family:"Times New Roman",serif'><span style='mso-element:field-end'></span></span><![endif]-->The Cancer Genome Atlas (TCGA). (n.d.). Retrieved from[ http://cancergenome.nih.gov](http://cancergenome.nih.gov/)

UCSC Cancer Genomics Hub. (n.d.). Retrieved from[ https://cghub.ucsc.edu](https://cghub.ucsc.edu/)

 
