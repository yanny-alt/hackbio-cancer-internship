counts <- read.csv('filtered_exp_data.csv', row.names = 1) # Loading the expression data for the purpose of creating heatmaps of the dataset 
res <- read.csv('stringent_edgeR_DEA_res.csv', row.names = 1) # Loading DEA results for a volcano plot

normalized <- as.matrix(log2(counts + 1)) # Log2 Transformation

library(gplots)

div_palette <- colorRampPalette(c('green', 'black', 'red'))(n=299) # To simply highlight the changes in expression of the genes
seq_palette <- colorRampPalette(c('lightblue', 'darkblue'))(n=299) # To highlight the intensity of expression of the genes

# To visualise the entire dataset with the diverging color palette
pdf('div_heatmaps.pdf', width = 10, height = 10)
heatmap.2(normalized, trace = 'none', scale = 'none', dendrogram = 'none', col = div_palette, key = T)
dev.off() 

# To visualise the entire dataset with the sequential color palette
pdf('seq_heatmaps.pdf', width = 10, height = 10)
heatmap.2(normalized, trace = 'none', scale = 'none', dendrogram = 'none', col = seq_palette, key = T)
dev.off() 

# To cluster the genes and samples in the dataset
# One has to scale the values either by rows or columns to cluster effectively

# Scaling and clustering the filtered dataset by rows
pdf('rowCluster.pdf', width=10, height=10)
heatmap.2(normalized, scale = 'row', trace = 'none', dendrogram = 'row', col=seq_palette, key = T) # To cluster the genes
dev.off()

# Scaling and clustering the filtered dataset by columns
pdf('colCluster.pdf', width=10, height=10)
heatmap.2(normalized, scale = 'col', trace = 'none', dendrogram = 'col', col=seq_palette, key = T) # To cluster the samples
dev.off()


# Loading ggplot2
library(ggplot2)

# Creating a column for the classification of the genes based on expression level
res$conditions <- with(res,
                       ifelse(FDR <= 0.01 & abs(logFC) >= 6, 
                              ifelse(logFC >= 6, 'Upregulated', 'Downregulated'), 'Insignificant'))

# Creating the volcano plot with ggplot2
volcano <- ggplot(res, aes(x=logFC, y=-log10(FDR))) + geom_point(aes(colour=conditions)) + ggtitle('Differentially Expressed Genes') + xlab('Log2 Fold Change') + ylab('-Log10 FDR')

# Saving the volcano plot as a pdf file
pdf('volcano6.pdf', width=10, height=10)
print(volcano)
dev.off()
