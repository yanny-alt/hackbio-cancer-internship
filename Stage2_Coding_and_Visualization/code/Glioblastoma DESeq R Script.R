setwd("C:/Users/tobij/Desktop/HB Essays and Papers/R Tasks") # To set my working directory
glioblastoma <- read.csv('C:/Users/tobij/Desktop/HB Essays and Papers/R Tasks/glioblastoma1.csv', row.names = 1) # Loading the glioblastoma dataset

View(glioblastoma) # To view the dataframe

str(glioblastoma) # To see the data types of the dataframe

glioblastomaFiltered <- glioblastoma[rowSums(glioblastoma >= 3) >= 3, ] # Filtering genes where their expression data is almost zero across all samples. AT least 3 samples for each gene must have a count of at least 3. 
# Upon filtering, the number of genes reduced from 582 to 461

sample_ids <- colnames(glioblastomaFiltered) # Getting the sample ids
tumor_types <- c(rep('recurrent', 5), rep('primary', 5)) # Creating a vector for representation of the sample types
tumor_types <- factor(tumor_types) # Factoring the tumor types to create 2 levels (recurrent and primary)
column_data <- data.frame(row.names = sample_ids, tumor_type = tumor_types) # Generating column data to host the sample ids and sample types

# install.packages("BiocManager") # Install DESeq2 if not already installed through Bioconductor
# BiocManager::install('DESeq2')

library(DESeq2) # Loading DESeq2 for analysis

ddsfm <- DESeqDataSetFromMatrix(countData = glioblastomaFiltered, colData = column_data, design = ~tumor_type) # To put the count dataset and the column dataset in one matrix for analysis
dea <- DESeq(ddsfm)
res <- results(dea)
summary(res)

write.csv(as.data.frame(res), file='glio_results.csv')
