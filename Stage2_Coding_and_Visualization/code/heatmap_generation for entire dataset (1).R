# Load necessary libraries
library(gplots)
library(matrixStats)

# Load the dataset
glioblastoma <- read.csv(url("https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv"), row.names = 1)

# Convert to numeric matrix
data_numeric <- as.matrix(glioblastoma)
data_numeric <- apply(data_numeric, 2, as.numeric)

# Ensure row names (Ensembl IDs) are preserved
rownames(data_numeric) <- rownames(glioblastoma)

# Log transformation
data_log <- log2(data_numeric + 1)

# Z-score normalization (scale rows)
data_normalized <- t(scale(t(data_log)))

# Ensure geneIDs are preserved in the normalized data
rownames(data_normalized) <- rownames(data_log)

# Define color palettes
color_palette <- colorRampPalette(c("yellow", "white", "darkgreen"))(100)
color_palette1 <- colorRampPalette(c("blue", "white", "green"))(100)  
color_palette2 <- colorRampPalette(c("black", "white", "red"))(100)  
color_palette3 <- colorRampPalette(c("purple", "white", "brown"))(100)
color_palette4 <- colorRampPalette(c("darkred", "white", "darkblue"))(100)

# Generate heatmaps
# Heatmap1
pdf("heatmap_dataset1.pdf", width = 12, height = 12)
heatmap.2(
  data_normalized,
  col = color_palette,
  trace = "none",
  scale = "none",
  dendrogram = "none",
  Rowv = TRUE,
  Colv = TRUE,
  key = TRUE,
  margins = c(5, 10),  
  cexRow = 0.5,        
  cexCol = 0.3,        
  density.info = "none",
  key.title = NA,
  key.xlab = "Expression Level"
)
dev.off()


# Heatmap2
pdf("heatmap_dataset.pdf", width = 12, height = 12)
heatmap.2(
  data_normalized,
  col = color_palette1,
  trace = "none",
  scale = "none",
  dendrogram = "none",
  Rowv = TRUE,
  Colv = TRUE,
  key = TRUE,
  margins = c(5, 10),  
  cexRow = 0.5,        
  cexCol = 0.3,        
  density.info = "none",
  key.title = NA,
  key.xlab = "Expression Level"
)
dev.off()

# Heatmap with Clustering of Genes Alone
pdf("heatmap_cluster_genes.pdf", width = 12, height = 12)
heatmap.2(
  data_normalized,
  col = color_palette2,
  trace = "none",
  scale = "none",
  dendrogram = "row",  
  Rowv = TRUE,
  Colv = FALSE,        
  key = TRUE,
  margins = c(5, 10),  
  cexRow = 0.5,        
  cexCol = 0.3,        
  density.info = "none",
  key.title = NA,
  key.xlab = "Expression Level"
)
dev.off()

# Heatmap with Clustering of Samples Alone
pdf("heatmap_cluster_samples.pdf", width = 12, height = 12)
heatmap.2(
  data_normalized,
  col = color_palette3,
  trace = "none",
  scale = "none",
  dendrogram = "column",  
  Rowv = FALSE,           
  Colv = TRUE,
  key = TRUE,
  margins = c(5, 10),  
  cexRow = 0.5,        
  cexCol = 0.3,        
  density.info = "none",
  key.title = NA,
  key.xlab = "Expression Level"
)
dev.off()


pdf("heatmap_dataset_cluster_both.pdf", width = 12, height = 12)
heatmap.2(
  data_normalized,
  col = color_palette4,
  trace = "none",
  scale = "none",
  dendrogram = "both",
  Rowv = TRUE,
  Colv = TRUE,
  key = TRUE,
  margins = c(5, 10),  
  cexRow = 0.5,        
  cexCol = 0.3,        
  density.info = "none",
  key.title = NA,
  key.xlab = "Expression Level"
)
dev.off()

