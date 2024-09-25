# Loading all necessary packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(limma)
library(edgeR)
library(EDASeq)
library(gplots)
library(dplyr)

# For project info:
getProjectSummary('TCGA-BRCA')
?GDCquery

# Generating a query object
brcaQ <- GDCquery(project = 'TCGA-BRCA', 
                 data.category = 'Transcriptome Profiling',
                 data.type = 'Gene Expression Quantification')

# Getting the results of the query
brcaQ_output <- getResults(brcaQ)

# Filtering for the information of the specific samples needed
sample_info <- brcaQ_output %>%
  filter(sample_type %in% c('Solid Tissue Normal', 'Primary Tumor')) %>%
  group_by(sample_type) %>%
  slice(1:20) %>%
  ungroup() %>%
  arrange(factor(sample_type, levels = c('Solid Tissue Normal', 'Primary Tumor'))) # To select the solid tissue normal (control) samples first

# Querying again to get the data on the sample cases filtered for
brcaQ2 <- GDCquery(project = 'TCGA-BRCA',
                  data.category = 'Transcriptome Profiling',
                  data.type = 'Gene Expression Quantification',
                  experimental.strategy = 'RNA-Seq',
                  workflow.type = 'STAR - Counts',
                  barcode = sample_info$cases, 
                  access = 'open'
)

# To download the above specified data
GDCdownload(brcaQ2)

# Preparing the data for analysis
brca.data <- GDCprepare(brcaQ2)
brca.assay.data <- assay(brca.data, 'unstranded')

# To obtain the names of all the samples in the brca.assay.data matrix
samples <- colnames(brca.assay.data)

# Converting the matrix into a proper gene expression dataset  
countData <- brca.assay.data[, samples]

# Saving countData
write.csv(as.data.frame(countData), file='pre_normalization_exp_data.csv')

# Normalization of the expression dataset based on the length of each gene
normalized <- TCGAanalyze_Normalization(tabDF = countData, geneInfo = geneInfoHT, method = 'geneLength')

# Filtering out noisy genes
filtered <- TCGAanalyze_Filtering(tabDF = normalized, method = 'quantile', qnt.cut = 0.25)
dim(filtered)

# Differential Expression Analysis with the EdgeR pipeline
dea <- TCGAanalyze_DEA(mat1 = filtered[, c(1:20)], 
                mat2 = filtered[, c(21:40)],
                Cond1type = 'Solid Tissue Normal',
                Cond2type = 'Primary Tumor',
                pipeline = 'edgeR',
                fdr.cut = 0.01,
                logFC.cut = 2
)

# Differential Expression Analysis with levels
dea_wl <- TCGAanalyze_LevelTab(dea, 'Solid Tissue Normal', 'Primary Tumor', 
                     filtered[, c(1:20)], 
                     filtered[, c(21:40)])

# Saving the normalized and filtered data
write.csv(as.data.frame(normalized), file = 'normalized_exp_data.csv')
write.csv(as.data.frame(filtered), file = 'filtered_exp_data.csv')

# Saving dea results
write.csv(as.data.frame(dea), file = 'stringent_edgeR_DEA_res.csv')
write.csv(as.data.frame(dea_wl), file = 'stringent_edgeR_DEAWL_res.csv')

# Separating the underregulated subset from the overregulated subset
underregulated <- subset(dea, logFC <= -6 & FDR <= 0.01)
overregulated <- subset(dea, logFC >= 6 & FDR <= 0.01)

# Saving these subsets
write.csv(as.data.frame(underregulated), file = 'EM_stringent_underregulated.csv')
write.csv(as.data.frame(overregulated), file = 'EM_stringent_overregulated.csv')

# Specifically getting the names of the genes
underregulated_genes <- rownames(underregulated)
overregulated_genes <- rownames(overregulated)

# Saving the names of the genes
write.csv(as.data.frame(underregulated_genes), file = 'EM_stringent_underregulated_genes.csv')
write.csv(as.data.frame(overregulated_genes), file = 'EM_stringent_overregulated_genes.csv')
