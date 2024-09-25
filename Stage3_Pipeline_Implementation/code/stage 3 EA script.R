#LOAD PACKAGES
library("TCGAbiolinks")
library(SummarizedExperiment)
library(biomaRt)
library(ggplot2)
library(EDASeq)
library(dplyr)
library(tidyr)
# Read your CSV file containing Ensembl gene IDs
upreggenedata <- read.csv("upregulatedgenes.csv", stringsAsFactors = FALSE)  # Read the data without converting strings to factors
downreggenedata <- read.csv("downregulatedgenes.csv", stringsAsFactors = FALSE)

# Set up biomaRt connection to Ensembl
ensembl <- useMart("ensembl")  # Connect to the Ensembl BioMart database
dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)  # Specify the dataset for human genes

# Retrieve gene symbols corresponding to the Ensembl IDs
results <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),  # Columns to retrieve: Ensembl ID and gene symbol
  filters = "ensembl_gene_id",  # Filter based on Ensembl gene IDs
  values = upreggenedata$overregulated_genes,  # Use the Ensembl IDs from the data frame
  mart = dataset  # Query the specified dataset
)

# Merge the results with the original data frame
upreggenemergeddata <- merge(
  upreggenedata,  # Original data frame
  results,  # Data frame containing the gene symbols
  by.x = "overregulated_genes",  # Match using the Ensembl IDs in the original data
  by.y = "ensembl_gene_id",  # Match with the Ensembl IDs in the results
  all.x = TRUE  # Keep all rows from the original data, filling with NA where there are no matches
)

# Save the merged data to a new CSV file
write.csv(upreggenemergeddata, "upreggenesymbol.csv", row.names = FALSE)  # Write the results to a CSV file without row names

# Retrieve gene symbols corresponding to the Ensembl IDs
results <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),  # Columns to retrieve: Ensembl ID and gene symbol
  filters = "ensembl_gene_id",  # Filter based on Ensembl gene IDs
  values = downreggenedata$underregulated_genes,  # Use the Ensembl IDs from the data frame
  mart = dataset  # Query the specified dataset
)

# Merge the results with the original data frame
downreggenemergeddata <- merge(
  downreggenedata,  # Original data frame
  results,  # Data frame containing the gene symbols
  by.x = "underregulated_genes",  # Match using the Ensembl IDs in the original data
  by.y = "ensembl_gene_id",  # Match with the Ensembl IDs in the results
  all.x = TRUE  # Keep all rows from the original data, filling with NA where there are no matches
)

# Save the merged data to a new CSV file
write.csv(downreggenemergeddata, "downreggenesymbol.csv", row.names = FALSE)  # Write the results to a CSV file without row names


# Enrichment Analysis

# Prepare gene lists for enrichment analysis
upregulated_genes <- na.omit(upreggenemergeddata$external_gene_name)
downregulated_genes <- na.omit(downreggenemergeddata$external_gene_name)

# Perform enrichment analysis using TCGAAnalyse_EAcomplete

upreggeneEA <- TCGAanalyze_EAcomplete(TFname="UPREGULATEDGENE_EA", upregulated_genes)
downreggeneEA <- TCGAanalyze_EAcomplete(TFname="DOWNREGULATEDGENE_EA", downregulated_genes)

# Save enrichment results in csv

write.csv(upreggeneEA, "upreggeneEA.csv", row.names = FALSE)
write.csv(downreggeneEA, "downreggeneEA.csv", row.names = FALSE)

#Read Your CSV: Read the enrichment results CSV file.

upregEA <- read.csv("upreggeneEA.csv")
downregEA <- read.csv("downreggeneEA.csv")

#Use the separate function from tidyr to split the enrichment results into separate columns.

upregEAplot <- upregEA %>%
  pivot_longer(cols = everything(), 
               names_to = "Enrichment_Type", 
               values_to = "Results") %>%
  separate(Results, into = c("GO_Term", "FDR", "ng", "ncommon"), sep = "; ", extra = "merge") %>%
  mutate(
    FDR = as.numeric(gsub("FDR= ", "", FDR)),
    ng = as.numeric(gsub("\\(ng=|\\)", "", ng)),
    ncommon = as.numeric(gsub("\\(ncommon=|\\)", "", ncommon))
  )

downregEAplot <- downregEA %>%
  pivot_longer(cols = everything(), 
               names_to = "Enrichment_Type", 
               values_to = "Results") %>%
  separate(Results, into = c("GO_Term", "FDR", "ng", "ncommon"), sep = "; ", extra = "merge") %>%
  mutate(
    FDR = as.numeric(gsub("FDR= ", "", FDR)),
    ng = as.numeric(gsub("\\(ng=|\\)", "", ng)),
    ncommon = as.numeric(gsub("\\(ncommon=|\\)", "", ncommon))
  )

# Save EA results in csv

write.csv(upregEAplot, "upreggeneplot.csv", row.names = FALSE)
write.csv(downregEAplot, "downreggeneplot.csv", row.names = FALSE)

# Read the EA CSV data

upregdata <- read.csv("upreggeneplot.csv")
downregdata <- read.csv("downreggeneplot.csv")

# Select the top 5 pathways based on FDR

upregtop_5 <- upregdata[1:5, ]
downregtop_5 <- downregdata[1:5, ]

# Create a new column for -log10(FDR) to use for coloring

upregtop_5$logFDR <- -log10(upregtop_5$FDR)
downregtop_5$logFDR <- -log10(downregtop_5$FDR)

# Extract the Pathway names from the GO_Term

upregtop_5$Pathway <- gsub("GO:\\d+~", "", upregtop_5$GO_Term)
downregtop_5$Pathway <- gsub("GO:\\d+~", "", downregtop_5$GO_Term)

# Create the lollipop plot

ggplot(upregtop_5, aes(x = FDR, y = reorder(Pathway, FDR))) +
  geom_segment(aes(xend = 0, yend = Pathway), color = "grey") +
  geom_point(aes(size = ng, color = logFDR), shape = 16) +  # Use ng for size
  scale_color_gradient(low = "blue", high = "red") +  # Color gradient from blue to red
  scale_size_continuous(range = c(3, 8)) +  # Adjusting the circle size range
  theme_minimal() +
  labs(title = "Top 5 Upregulated Genes Enriched Pathways",
       x = "-log10(FDR)",
       y = "Pathways",
       color = "-log10(FDR)",
       size = "Number of Genes") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# Save the plot

ggsave("upreggenepathway_plot.png", width = 12, height = 6, dpi = 300, bg = "white")

ggplot(downregtop_5, aes(x = FDR, y = reorder(Pathway, FDR))) +
  geom_segment(aes(xend = 0, yend = Pathway), color = "grey") +
  geom_point(aes(size = ng, color = logFDR), shape = 16) +  # Use ng for size
  scale_color_gradient(low = "blue", high = "red") +  # Color gradient from blue to red
  scale_size_continuous(range = c(3, 8)) +  # Adjusting the circle size range
  theme_minimal() +
  labs(title = "Top 5 Downregulated Genes Enriched Pathways",
       x = "-log10(FDR)",
       y = "Pathways",
       color = "-log10(FDR)",
       size = "Number of Genes") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# Save the plot

ggsave("downreggenepathway_plot.png", width = 12, height = 6, dpi = 300, bg = "white")
