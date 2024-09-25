# Load necessary libraries
library(tidyverse)
library(caret)
library(randomForest)
library(ggplot2)

# Load dataset
df <- read_csv("C:/Users/Yanny/Downloads/normalized_exp_data.csv")

# Inspect the dataset
head(df)
str(df)

# Restructure the data by adding labels for the classes
df_tidy <- df %>%
  pivot_longer(cols = -...1, names_to = "Sample", values_to = "Expression") %>%
  mutate(
    Class = case_when(
      Sample %in% colnames(df)[2:21] ~ "Normal",
      Sample %in% colnames(df)[22:41] ~ "Tumor"
    )
  )

# Confirm the structure of the new tidy dataset
head(df_tidy)

# Prepare the feature matrix (gene expressions) and the labels
df_ml <- df %>%
  select(-...1) %>%
  t() %>%
  as.data.frame()

colnames(df_ml) <- df$...1  # Gene names as column names
df_ml$Class <- c(rep("Normal", 20), rep("Tumor", 20))  # Assign labels

# Inspect the prepared data
head(df_ml)

set.seed(123)  # For reproducibility

# Split the data into training (70%) and test (30%) sets
trainIndex <- createDataPartition(df_ml$Class, p = 0.7, list = FALSE)
trainData <- df_ml[trainIndex, ]
testData <- df_ml[-trainIndex, ]

# Install necessary libraries if not already installed
if (!require(FactoMineR)) install.packages("FactoMineR")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(caret)) install.packages("caret")
if (!require(randomForest)) install.packages("randomForest")

# Load libraries
library(FactoMineR)
library(ggplot2)
library(caret)
library(randomForest)

# Ensure Class is a factor (for classification)
trainData$Class <- as.factor(trainData$Class)
testData$Class <- as.factor(testData$Class)

# Run PCA to reduce dimensionality (exclude the Class column)
df_pca <- PCA(trainData[, -ncol(trainData)], ncp = 20)  # Keep 20 principal components

# Extract PCA components for training and test sets
trainData_pca <- as.data.frame(df_pca$ind$coord)
testData_pca <- as.data.frame(predict(df_pca, newdata = testData[, -ncol(testData)])$coord)

# Add the class labels back to the reduced datasets
trainData_pca$Class <- trainData$Class
testData_pca$Class <- testData$Class

# Check the structure of the PCA components
str(trainData_pca)
str(testData_pca)

# Train Random Forest model on the reduced dataset
set.seed(123)
rf_model <- randomForest(Class ~ ., data = trainData_pca, importance = TRUE, ntree = 500)

# Predict on the test set
rf_predictions <- predict(rf_model, newdata = testData_pca)

# Evaluate performance using confusion matrix
confusion_matrix <- confusionMatrix(rf_predictions, testData_pca$Class)
print(confusion_matrix)

# Plot feature importance
importance <- importance(rf_model)
varImpPlot(rf_model)

# Visualizing Confusion Matrix with ggplot
cm <- as.table(confusion_matrix$table)
ggplot(as.data.frame(cm), aes(Reference, Prediction)) +
  geom_tile(aes(fill = Freq), color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = Freq)) +
  theme_minimal() +
  labs(title = "Confusion Matrix", x = "True Labels", y = "Predicted Labels")

# Apply 10-fold cross-validation on PCA-reduced dataset
train_control <- trainControl(method = "cv", number = 10)
rf_cv_model <- train(Class ~ ., data = trainData_pca, method = "rf", trControl = train_control)

# Cross-validation results
print(rf_cv_model)

