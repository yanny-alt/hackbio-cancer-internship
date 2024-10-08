# Loading all necessary packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(limma)
library(edgeR)
library(EDASeq)
library(gplots)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(tibble)

# Set the project code for Low-Grade Glioma
lggCode <- 'TCGA-LGG'

# Get project summary
getProjectSummary(lggCode)

# Generate a query object for gene expression data
lggQ <- GDCquery(project = lggCode,
                  data.category = 'Transcriptome Profiling',
                  data.type = 'Gene Expression Quantification',)


# Download the data
GDCdownload(lggQ)

# Prepare the data for analysis
lgg.data <- GDCprepare(lggQ)
lgg.assay.data <- assay(lgg.data, 'unstranded')

# Obtain sample IDs
samples <- colnames(lgg.assay.data)

# Convert the matrix into a proper gene expression dataset
countData <- lgg.assay.data[, samples]

# Save raw count data
write.csv(as.data.frame(countData), file = 'pre_normalization_exp_data.csv')

# Load subtype data to check IDH status
subtypeData <- TCGAquery_subtype('lgg')

# Filter samples for IDH status and merge with count data
patientIDs <- substr(samples, 1, 12)
idhStatus <- subtypeData[, c("patient", "IDH.status")]
idhStatus <- idhStatus[idhStatus$patient %in% patientIDs, ]

# Merge IDH status with count data
mergedData <- merge(idhStatus, data.frame(patient = patientIDs, sample = samples), by = "patient")

# Exclude samples without IDH status
mergedData <- mergedData[!is.na(mergedData$IDH.status), ]
countData <- countData[, mergedData$sample]

# Factorize the IDH status (Wildtype vs Mutant)
mergedData$IDH.status <- factor(mergedData$IDH.status, levels = c("WT", "Mutant"))

# Normalization of the expression dataset
normalized <- TCGAanalyze_Normalization(tabDF = countData, geneInfo = geneInfoHT, method = 'geneLength')

# Filtering out noisy genes
filtered <- TCGAanalyze_Filtering(tabDF = normalized, method = 'quantile', qnt.cut = 0.25)

# Differential Expression Analysis using TCGAanalyze_DEA
dea <- TCGAanalyze_DEA(mat1 = filtered[, mergedData$IDH.status == "WT"],
                       mat2 = filtered[, mergedData$IDH.status == "Mutant"],
                       Cond1type = 'IDH Wildtype',
                       Cond2type = 'IDH Mutant',
                       pipeline = 'edgeR',
                       fdr.cut = 0.01,
                       logFC.cut = 2.5)

# Differential Expression Analysis with levels
dea_wl <- TCGAanalyze_LevelTab(dea, 'IDH Wildtype', 'IDH Mutant', 
                               filtered[, mergedData$IDH.status == "WT"], 
                               filtered[, mergedData$IDH.status == "Mutant"])

# Saving the normalized and filtered data and merge data
write.csv(as.data.frame(normalized), file = 'normalized_exp_data.csv')
write.csv(as.data.frame(filtered), file = 'filtered_exp_data.csv')
write.csv(mergedData, file = 'merged_data.csv', row.names = FALSE)


# Saving DEA results
write.csv(as.data.frame(dea), file = 'FC_2.5_IDH_WT_vs_Mut_DE_results.csv')
write.csv(as.data.frame(dea_wl), file = 'FC_2.5_IDH_WT_vs_Mut_DEA_LevelTab.csv')

# Separating the underregulated and overregulated subsets
underregulated <- subset(dea, logFC <= -2.5 & FDR <= 0.01)
overregulated <- subset(dea, logFC >= 2.5 & FDR <= 0.01)

# Saving these subsets
write.csv(as.data.frame(underregulated), file = 'IDH_Underregulated_Genes.csv')
write.csv(as.data.frame(overregulated), file = 'IDH_Overregulated_Genes.csv')

# Specifically getting the names of the genes
underregulated_genes <- rownames(underregulated)
overregulated_genes <- rownames(overregulated)

# Saving the names of the genes
write.csv(as.data.frame(underregulated_genes), file = 'IDH_Underregulated_Gene_Names.csv')
write.csv(as.data.frame(overregulated_genes), file = 'IDH_Overregulated_Gene_Names.csv')

#Heatmap Plotting

# Read in the expression data and metadata
exp_data <- read.csv("filtered_exp_data.csv", row.names = 1)
metadata <- read.csv("merged_data.csv")

# Ensure that the column names of the expression data match the sample IDs in metadata
colnames(exp_data) <- gsub("\\.", "-", colnames(exp_data))  # Replace '.' with '-' for matching

# Order metadata based on sample order in exp_data
metadata <- metadata[match(colnames(exp_data), metadata$sample), ]

# Prepare row annotations for the mutation status
sample_annotations <- data.frame(IDH.status = metadata$IDH.status)
rownames(sample_annotations) <- metadata$sample

# Define color schemes for the annotation
annotation_colors <- list(IDH.status = c(WT = "blue", Mutant = "red"))

# Normalize the expression data using log2 transformation (adding 1 to avoid log(0))
exp_data_normalized <- log2(exp_data + 1)

# Subset the normalized data to include all genes and only the first 50 samples
exp_data_subset <- exp_data_normalized[, 1:100]

# Create a breaks list for color compression
breaksList <- seq(min(exp_data_subset), max(exp_data_subset), length.out = 100)

# Heatmap without clustering for both genes and samples
pheatmap(
  exp_data_subset, 
  annotation_col = sample_annotations[1:100, , drop = FALSE],
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,  # Do not cluster genes
  cluster_cols = FALSE,  # Do not cluster samples
  show_rownames = FALSE,  # Show gene labels
  annotation_legend = TRUE,   # Show annotation legend
  annotation_fontsize = 3,    # Decrease font size for WT and Mutant annotations
  color = colorRampPalette(c("darkgreen", "yellow"))(50), # Color gradient
  breaks = breaksList,    # Apply breaks for color scaling
  main = "Heatmap OF LLG EXPRESION",
  width = 20,             # Heatmap width
  height = 25             # Heatmap height
)

