# Loading required datasets
upreg <- read.csv('lgg_upreg.csv', row.names = 1)
downreg <- read.csv('lgg_downreg.csv', row.names = 1)

# Adding constants to the instances of the pvalue and fdr that are 0 to avoid errors.
pv_epsilon <- 1.0e-323
fdr_epsilon <- 1.0e-321

# Replacing all zero pvalues and fdrs in the downreg dataset to avoid issues during log transformation
downreg$PValue[downreg$PValue == 0] <- pv_epsilon
downreg$FDR[downreg$FDR == 0] <- fdr_epsilon

# Creating a reliable ranking metric for the upregulated and downregulated genes 
upreg$ranking_metric <- upreg$logFC * -log10(upreg$FDR)
downreg$ranking_metric <- downreg$logFC * -log10(downreg$FDR)

# Ordering the ranking metrics for both datasets
upreg <- upreg[order(upreg$ranking_metric, decreasing = T), ] # In descending order
downreg <- downreg[order(downreg$ranking_metric), ] # In ascending order

# Saving the lists of top genes
top10UGs <- rownames(upreg)[c(1:10)]
top10DGs <- rownames(downreg)[c(1:10)]

top10UGsDF <- data.frame(top10UGs)
top10DGsDF <- data.frame(top10DGs)

mart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

# Getting gene_infos
topUGinfo <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
             filters = 'ensembl_gene_id',
             values = top10UGs, 
             mart = mart)

topDGinfo <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = top10DGs, 
                   mart = mart)

write.csv(top10UGsDF, file = 'top10UGs.csv', row.names = F)
write.csv(top10DGsDF, file = 'top10DGs.csv', row.names = F)

write.csv(topUGinfo, file = 'top10UGinfo.csv', row.names = F)
write.csv(topDGinfo, file = 'top10DGinfo.csv', row.names = F)

# Getting gene_infos
ug_ids <- rownames(upreg)
dg_ids <- rownames(downreg)

mart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

# Getting gene_infos
ugi <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
             filters = 'ensembl_gene_id',
             values = ug_ids, 
             mart = mart)

dgi <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
             filters = 'ensembl_gene_id',
             values = dg_ids, 
             mart = mart)

# Merging the gene infos with the gene csvs
upregWithSymbols <- merge(upreg, ugi, by.x= 'row.names', by.y = 'ensembl_gene_id')
downregWithSymbols <- merge(downreg, dgi, by.x= 'row.names', by.y = 'ensembl_gene_id')

# Removing genes with missing gene symbols
cleanUpregWithSymbols <- upregWithSymbols[!is.na(upregWithSymbols$hgnc_symbol) & upregWithSymbols$hgnc_symbol != '', ]
cleanDownregWithSymbols <- downregWithSymbols[!is.na(downregWithSymbols$hgnc_symbol) & downregWithSymbols$hgnc_symbol != '', ]

ugEA <- TCGAanalyze_EAcomplete(TFname = 'UpregGeneEA', cleanUpregWithSymbols$hgnc_symbol)
dgEA <- TCGAanalyze_EAcomplete(TFname = 'DownregGeneEA', cleanDownregWithSymbols$hgnc_symbol)

write.csv(ugEA, 'lggUgEA.csv', row.names=F)
write.csv(dgEA, 'lggDgEA.csv', row.names=F)

upregEA <- read.csv('lggUgEA.csv')
downregEA <- read.csv('lggDgEA.csv')

library(tidyr)

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

library(ggplot2)

ggplot(upregtop_5, aes(x = logFDR, y = reorder(Pathway, logFDR))) +
  geom_segment(aes(xend = 0, yend = Pathway), color = "grey") +
  geom_point(aes(size = ncommon, color = logFDR), shape = 16) +  # Use ng for size
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

ggplot(downregtop_5, aes(x = logFDR, y = reorder(Pathway, logFDR))) +
  geom_segment(aes(xend = 0, yend = Pathway), color = "grey") +
  geom_point(aes(size = ncommon, color = logFDR), shape = 16) +  # Use ng for size
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
