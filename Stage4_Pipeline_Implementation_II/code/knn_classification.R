install.packages("ranger")
install.packages("Rtsne")

# Load required libraries
library(tidyverse)
library(caret)
library(ranger)
library(parallel)
library(class)
library(ggplot2)
library(cluster) # For silhouette analysis
library(Rtsne) # For t-SNE visualization

# Load the gene expression dataset
gene_expr <- read.csv("C:/Users/Yanny/Downloads/filtered_exp_data.csv")

# Inspect the data
head(gene_expr)

# Remove the Gene ID column for further processing
gene_data <- gene_expr[,-1]  # Exclude the Gene ID column

# Ensure the data is numeric (for expression levels)
gene_data <- as.data.frame(lapply(gene_data, as.numeric))

# Check for missing values and handle them (e.g., impute or remove)
gene_data <- na.omit(gene_data)

# Example: Create a mock target variable (you might have different logic here)
set.seed(123)
gene_data$Cluster <- sample(1:6, nrow(gene_data), replace = TRUE)  # Simulated clusters (6 IDH statuses)

# Split data into training and testing sets
set.seed(123)
trainIndex <- createDataPartition(gene_data$Cluster, p = .8, list = FALSE, times = 1)
geneTrain <- gene_data[trainIndex, ]
geneTest <- gene_data[-trainIndex, ]

# Dimensionality Reduction using PCA
pca <- prcomp(geneTrain[,-ncol(geneTrain)], center = TRUE, scale. = TRUE)
explained_variance <- summary(pca)$importance[3,] # Cumulative proportion of variance explained
num_components <- min(which(explained_variance >= 0.95)) # Choose number of components explaining at least 95% of variance

# Reduce dimensionality of the data
geneTrain_reduced <- data.frame(pca$x[, 1:num_components])
geneTrain_reduced$Cluster <- geneTrain$Cluster

geneTest_reduced <- predict(pca, geneTest[,-ncol(geneTest)])[, 1:num_components]
geneTest_reduced <- data.frame(geneTest_reduced)
geneTest_reduced$Cluster <- geneTest$Cluster

# Random Forest Classification using ranger (faster version)
num_cores <- detectCores() - 1 # Use all but one core
rf_model <- ranger(
  as.factor(Cluster) ~ ., 
  data = geneTrain_reduced,
  num.trees = 100,
  mtry = 5,
  importance = 'impurity',
  num.threads = num_cores
)

# Predict using the Random Forest model
rf_predictions <- predict(rf_model, data = geneTest_reduced)$predictions

# Confusion matrix to evaluate Random Forest results
rf_cm <- confusionMatrix(as.factor(rf_predictions), as.factor(geneTest_reduced$Cluster))
print(rf_cm)

# Variable importance plot
importance_scores <- rf_model$variable.importance
barplot(sort(importance_scores, decreasing = TRUE), main = "Variable Importance", col = "lightblue")

# KNN Classification using reduced features from PCA
k <- 3  # Set the number of neighbors
knn_predictions <- knn(
  train = geneTrain_reduced[, -ncol(geneTrain_reduced)], 
  test = geneTest_reduced[, -ncol(geneTest_reduced)], 
  cl = geneTrain_reduced$Cluster, 
  k = k
)

# Confusion matrix to evaluate KNN results
knn_cm <- confusionMatrix(knn_predictions, as.factor(geneTest_reduced$Cluster))
print(knn_cm)

# Function to plot confusion matrix
plot_confusion_matrix <- function(cm, title) {
  cm_table <- as.data.frame(cm$table)  # Convert confusion matrix to dataframe
  colnames(cm_table) <- c("Reference", "Prediction", "Freq")  # Ensure columns are named correctly
  
  ggplot(data = cm_table, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = Freq), color = "white") +
    scale_fill_gradient(low = "white", high = "blue") +
    geom_text(aes(label = Freq), color = "black") +
    labs(title = title, x = "Reference", y = "Prediction") +
    theme_minimal()
}

# Plot KNN confusion matrix
plot_confusion_matrix(knn_cm, "KNN Confusion Matrix")

# Plot Random Forest confusion matrix
plot_confusion_matrix(rf_cm, "Random Forest Confusion Matrix")

# K-means clustering
set.seed(123)
kmeans_result <- kmeans(geneTrain_reduced[, -ncol(geneTrain_reduced)], centers = 6, nstart = 25)
geneTrain_reduced$Cluster <- kmeans_result$cluster

# Evaluate clustering using silhouette analysis
sil <- silhouette(kmeans_result$cluster, dist(geneTrain_reduced[, -ncol(geneTrain_reduced)]))
avg_sil <- mean(sil[, 3])
print(paste("Average Silhouette Score:", round(avg_sil, 2)))

# Plot the clustering result using t-SNE for visualization
set.seed(123)
tsne_result <- Rtsne(geneTrain_reduced[, -ncol(geneTrain_reduced)], dims = 2, perplexity = 30)
plot_df <- data.frame(tsne_result$Y)
plot_df$Cluster <- as.factor(kmeans_result$cluster)

ggplot(plot_df, aes(x = X1, y = X2, color = Cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(title = "t-SNE Plot of Clustering Results", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal()
