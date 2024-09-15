# Load the CSV file
enrichment_data <- read.csv("C:/Users/Yanny/Downloads/enrichment_biological_processes1 (2).csv")

# View the first few rows
head(enrichment_data)

# Sort the data by the FDR (ascending order)
enrichment_data_sorted <- enrichment_data[order(enrichment_data$Enrichment.FDR), ]

# Select the top 5 pathways
top5_pathways <- enrichment_data_sorted[1:5, ]

# Select the top 3 pathways
top3_pathways <- enrichment_data_sorted[1:3, ]

# Extract the genes from the top 5 pathways
top5_genes <- top5_pathways$Genes

# Extract the genes from the top 3 pathways
top3_genes <- top3_pathways$Genes

# Assuming enrichment_data_sorted has the top pathways already sorted by FDR
top5_pathways <- enrichment_data_sorted[1:5, ]

# Calculate the negative log10 of the FDR for each pathway
top5_pathways$neg_log10_FDR <- -log10(top5_pathways$Enrichment.FDR)

# Load ggplot2
library(ggplot2)

# Create a bubble plot
ggplot(top5_pathways, aes(x=reorder(Pathway, neg_log10_FDR), y=neg_log10_FDR)) +
  geom_point(aes(size=nGenes, color=neg_log10_FDR)) +
  scale_color_gradient(low="blue", high="red") +
  labs(title="Top 5 Enriched Pathways",
       x="Pathway",
       y="-log10(FDR)",
       size="Number of Genes",
       color="-log10(FDR)") +
  theme_minimal() +
  coord_flip()
