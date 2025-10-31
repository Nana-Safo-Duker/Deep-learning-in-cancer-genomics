# Deep Learning in Cancer Genomics: Comprehensive Analysis Script
#
# This R script provides a complete pipeline for analyzing cancer genomics data
# using machine learning and statistical approaches. While R is not typically used
# for deep learning (Python/PyTorch/TensorFlow are preferred), this script demonstrates
# statistical analysis, visualization, and traditional ML approaches suitable for R.
#
# Research Paper: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-024-01315-6
#
# Author: Research Analysis
# Date: 2024

# ============================================================================
# 1. LIBRARY IMPORTS AND SETUP
# ============================================================================

# Install packages if not already installed
required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "caret", "randomForest", 
  "pROC", "RColorBrewer", "corrplot", "VIM", "factoextra",
  "ggthemes", "patchwork", "matrixStats"
)

install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    install.packages(new_packages, dependencies = TRUE)
  }
  sapply(packages, require, character.only = TRUE)
}

cat("Loading required packages...\n")
suppressPackageStartupMessages({
  install_if_missing(required_packages)
})

# Set random seed for reproducibility
set.seed(42)

# Set working directory (adjust if needed)
# setwd("path/to/your/project")

cat("Libraries loaded successfully!\n")
cat("R version:", R.version.string, "\n\n")

# ============================================================================
# 2. HELPER FUNCTIONS
# ============================================================================

#' Simulate cancer genomics data for demonstration
#' 
#' @param n_samples Number of samples (patients)
#' @param n_genes Number of genes in expression matrix
#' @param n_cancer_types Number of different cancer types
#' @return List containing gene expression, mutations, and clinical data
simulate_cancer_genomics_data <- function(n_samples = 500, n_genes = 2000, n_cancer_types = 5) {
  cat("Generating simulated cancer genomics data...\n")
  
  # Simulate gene expression data (log-normal distribution, typical for RNA-seq)
  base_expression <- matrix(
    rlnorm(n_samples * n_genes, meanlog = 5, sdlog = 2),
    nrow = n_samples, ncol = n_genes
  )
  
  # Add cancer-type-specific expression patterns
  cancer_types <- paste0("Cancer_Type_", 1:n_cancer_types)
  labels <- sample(cancer_types, size = n_samples, replace = TRUE)
  
  # Introduce cancer-type-specific differentially expressed genes
  for (i in 1:n_cancer_types) {
    cancer_type <- cancer_types[i]
    mask <- labels == cancer_type
    n_de_genes <- floor(n_genes / 10)  # 10% of genes are differentially expressed
    de_indices <- sample(n_genes, size = n_de_genes, replace = FALSE)
    
    fold_change <- runif(sum(mask) * n_de_genes, min = 2, max = 5)
    base_expression[mask, de_indices] <- base_expression[mask, de_indices] * 
      matrix(fold_change, nrow = sum(mask), ncol = n_de_genes)
  }
  
  # Create gene expression data frame
  gene_expression <- as.data.frame(base_expression)
  colnames(gene_expression) <- paste0("Gene_", 1:n_genes)
  
  # Simulate mutation data (binary matrix)
  mutation_rate <- 0.05  # 5% mutation rate per gene-sample pair
  n_mut_genes <- floor(n_genes / 10)
  mutations <- as.data.frame(
    matrix(rbinom(n_samples * n_mut_genes, size = 1, prob = mutation_rate),
           nrow = n_samples, ncol = n_mut_genes)
  )
  colnames(mutations) <- paste0("Gene_", 1:n_mut_genes, "_mut")
  
  # Simulate clinical data
  ages <- rnorm(n_samples, mean = 60, sd = 15)
  ages <- pmax(pmin(ages, 90), 18)  # Clip to realistic age range
  
  # Survival times (exponential distribution with cancer-type-specific rates)
  survival_times <- numeric(n_samples)
  for (i in 1:length(labels)) {
    rate <- runif(1, min = 0.01, max = 0.05)
    survival_times[i] <- rexp(1, rate = rate)
  }
  
  clinical_data <- data.frame(
    sample_id = paste0("Sample_", 1:n_samples),
    cancer_type = labels,
    age = ages,
    survival_time = survival_times,
    stage = sample(c("I", "II", "III", "IV"), size = n_samples, 
                   replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))
  )
  
  cat("Data generation complete!\n")
  cat("  Gene Expression:", nrow(gene_expression), "samples x", ncol(gene_expression), "genes\n")
  cat("  Mutations:", nrow(mutations), "samples x", ncol(mutations), "mutations\n")
  cat("  Clinical Data:", nrow(clinical_data), "samples\n\n")
  
  return(list(
    gene_expression = gene_expression,
    mutations = mutations,
    clinical = clinical_data
  ))
}

#' Preprocess gene expression data
#' 
#' @param gene_expr Gene expression data frame
#' @param clinical Clinical data with cancer type labels
#' @param n_top_genes Number of top variable genes to select
#' @return List containing train/test splits and preprocessing objects
preprocess_data <- function(gene_expr, clinical, n_top_genes = 500) {
  cat("Preprocessing data...\n")
  
  # Normalization: Standard scaling (z-score)
  gene_expr_scaled <- as.data.frame(scale(gene_expr))
  
  # Feature selection: Select top variable genes (variance-based)
  gene_variance <- apply(gene_expr_scaled, 2, var)
  gene_variance_sorted <- sort(gene_variance, decreasing = TRUE)
  top_genes <- names(gene_variance_sorted)[1:n_top_genes]
  gene_expr_selected <- gene_expr_scaled[, top_genes]
  
  cat("Selected", n_top_genes, "most variable genes from", ncol(gene_expr), "total genes\n")
  
  # Encode labels
  cancer_types <- unique(clinical$cancer_type)
  clinical$cancer_type_encoded <- as.numeric(factor(clinical$cancer_type))
  
  # Create train-test split with stratification
  train_indices <- createDataPartition(
    clinical$cancer_type, 
    p = 0.8, 
    list = FALSE
  )
  
  X_train <- as.matrix(gene_expr_selected[train_indices, ])
  X_test <- as.matrix(gene_expr_selected[-train_indices, ])
  y_train <- clinical$cancer_type_encoded[train_indices]
  y_test <- clinical$cancer_type_encoded[-train_indices]
  
  cat("Train set size:", nrow(X_train), "\n")
  cat("Test set size:", nrow(X_test), "\n")
  cat("Number of features:", ncol(X_train), "\n")
  cat("Number of classes:", length(unique(y_train)), "\n\n")
  
  return(list(
    X_train = X_train,
    X_test = X_test,
    y_train = y_train,
    y_test = y_test,
    feature_names = top_genes,
    class_names = cancer_types
  ))
}

#' Perform statistical analysis
#' 
#' @param X_train Training features
#' @param X_test Test features
#' @param y_train Training labels
#' @param y_test Test labels
#' @param feature_names Names of features
perform_statistical_analysis <- function(X_train, X_test, y_train, y_test, feature_names) {
  cat("Performing statistical analysis...\n")
  
  # Summary statistics
  cat("\nSummary Statistics:\n")
  cat("  Mean expression:", mean(X_train), "\n")
  cat("  Median expression:", median(X_train), "\n")
  cat("  Standard deviation:", sd(X_train), "\n")
  cat("  Minimum value:", min(X_train), "\n")
  cat("  Maximum value:", max(X_train), "\n\n")
  
  # T-test between cancer types (first two types)
  cancer_types_unique <- unique(c(y_train, y_test))
  if (length(cancer_types_unique) >= 2) {
    type1_mask <- y_train == cancer_types_unique[1]
    type2_mask <- y_train == cancer_types_unique[2]
    
    # Test on first gene
    gene1_train <- X_train[, 1]
    t_test_result <- t.test(gene1_train[type1_mask], gene1_train[type2_mask])
    
    cat("T-test comparing cancer types (Gene:", feature_names[1], "):\n")
    cat("  T-statistic:", t_test_result$statistic, "\n")
    cat("  P-value:", t_test_result$p.value, "\n")
    cat("  Mean Type 1:", mean(gene1_train[type1_mask]), "\n")
    cat("  Mean Type 2:", mean(gene1_train[type2_mask]), "\n\n")
  }
  
  # Cross-validation
  cat("Performing cross-validation...\n")
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = multiClassSummary
  )
  
  # Note: This is computationally intensive, so using a subset
  # For full analysis, uncomment the following:
  # rf_model_cv <- train(
  #   x = X_train[1:min(200, nrow(X_train)), ],
  #   y = factor(y_train[1:min(200, length(y_train))]),
  #   method = "rf",
  #   trControl = ctrl,
  #   ntree = 100
  # )
  
  cat("Statistical analysis complete!\n\n")
}

# ============================================================================
# 3. DATA GENERATION AND LOADING
# ============================================================================

cat(paste0(rep("=", 60), collapse=""), "\n")
cat("Deep Learning in Cancer Genomics: Comprehensive Analysis (R)\n")
cat(paste0(rep("=", 60), collapse=""), "\n\n")

# Generate simulated data
data_list <- simulate_cancer_genomics_data(n_samples = 500, n_genes = 2000, n_cancer_types = 5)

gene_expr <- data_list$gene_expression
mutations <- data_list$mutations
clinical <- data_list$clinical

# ============================================================================
# 4. DATA EXPLORATION AND VISUALIZATION
# ============================================================================

cat("Step 1: Data Exploration\n")
cat(strrep("-", 60), "\n")

# Distribution of expression values
cat("\nCreating data exploration visualizations...\n")

# Plot 1: Distribution of expression values
p1 <- ggplot(data.frame(value = as.vector(as.matrix(gene_expr[, 1:100]))), 
             aes(x = value)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
  labs(title = "Distribution of Gene Expression Values",
       x = "Expression Level (log-scale)", y = "Frequency") +
  theme_minimal()

# Plot 2: Expression by cancer type (box plot)
gene_expr_long <- gene_expr[, 1:20] %>%
  mutate(sample = row_number()) %>%
  pivot_longer(cols = -sample, names_to = "gene", values_to = "expression") %>%
  mutate(cancer_type = clinical$cancer_type[sample])

p2 <- ggplot(gene_expr_long, aes(x = cancer_type, y = expression, fill = cancer_type)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Expression Distribution by Cancer Type (Sample Genes)",
       x = "Cancer Type", y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Plot 3: Cancer type distribution
p3 <- ggplot(clinical, aes(x = cancer_type, fill = cancer_type)) +
  geom_bar(alpha = 0.7) +
  labs(title = "Distribution of Cancer Types",
       x = "Cancer Type", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Plot 4: Survival time distribution
p4 <- ggplot(clinical, aes(x = survival_time)) +
  geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7, color = "black") +
  labs(title = "Distribution of Survival Times",
       x = "Survival Time (months)", y = "Frequency") +
  theme_minimal()

# Save plots
ggsave("r_expression_distribution.png", p1, width = 10, height = 6, dpi = 300)
ggsave("r_expression_by_type.png", p2, width = 10, height = 6, dpi = 300)
ggsave("r_cancer_type_distribution.png", p3, width = 10, height = 6, dpi = 300)
ggsave("r_survival_distribution.png", p4, width = 10, height = 6, dpi = 300)

cat("Exploratory plots saved!\n\n")

# ============================================================================
# 5. DATA PREPROCESSING
# ============================================================================

cat("Step 2: Data Preprocessing\n")
cat(strrep("-", 60), "\n")

preprocessed <- preprocess_data(gene_expr, clinical, n_top_genes = 500)

X_train <- preprocessed$X_train
X_test <- preprocessed$X_test
y_train <- preprocessed$y_train
y_test <- preprocessed$y_test
feature_names <- preprocessed$feature_names
class_names <- preprocessed$class_names

# ============================================================================
# 6. DIMENSIONALITY REDUCTION (PCA)
# ============================================================================

cat("Step 3: Dimensionality Reduction (PCA)\n")
cat(strrep("-", 60), "\n")

# Perform PCA
pca_result <- prcomp(X_train, center = TRUE, scale. = FALSE)  # Already scaled

# Explained variance
explained_variance <- summary(pca_result)$importance[2, ]
cumulative_variance <- summary(pca_result)$importance[3, ]

cat("First 10 principal components explain", 
    cumulative_variance[10] * 100, "% of variance\n")
cat("First 50 principal components explain", 
    cumulative_variance[50] * 100, "% of variance\n\n")

# Visualization: Scree plot
pca_df <- data.frame(
  PC = 1:length(explained_variance),
  Explained_Variance = explained_variance,
  Cumulative_Variance = cumulative_variance
)

p_scree <- ggplot(pca_df[1:20, ], aes(x = PC, y = Explained_Variance)) +
  geom_line(size = 1.2, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  labs(title = "PCA Scree Plot (First 20 Components)",
       x = "Principal Component", y = "Explained Variance Ratio") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

p_cumulative <- ggplot(pca_df[1:50, ], aes(x = PC, y = Cumulative_Variance)) +
  geom_line(size = 1.2, color = "darkgreen") +
  geom_point(size = 2, color = "darkgreen") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(title = "Cumulative Explained Variance",
       x = "Number of Components", y = "Cumulative Variance Ratio") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

ggsave("r_pca_scree.png", p_scree, width = 10, height = 6, dpi = 300)
ggsave("r_pca_cumulative.png", p_cumulative, width = 10, height = 6, dpi = 300)

# 2D PCA visualization
pca_2d <- prcomp(X_train, center = TRUE, scale. = FALSE, rank. = 2)
X_train_2d <- as.data.frame(pca_2d$x)
X_train_2d$cancer_type <- factor(y_train)

p_pca2d <- ggplot(X_train_2d, aes(x = PC1, y = PC2, color = cancer_type)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(title = "PCA Visualization: First Two Principal Components",
       x = paste0("PC1 (", round(explained_variance[1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(explained_variance[2] * 100, 1), "% variance)")) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

ggsave("r_pca_2d.png", p_pca2d, width = 10, height = 8, dpi = 300)

cat("PCA visualizations saved!\n\n")

# ============================================================================
# 7. MACHINE LEARNING MODEL TRAINING
# ============================================================================

cat("Step 4: Machine Learning Model Training\n")
cat(strrep("-", 60), "\n")

# Train Random Forest model
cat("Training Random Forest model...\n")

rf_model <- randomForest(
  x = X_train,
  y = factor(y_train),
  ntree = 100,
  mtry = sqrt(ncol(X_train)),
  importance = TRUE,
  proximity = TRUE
)

# Make predictions
y_pred <- predict(rf_model, newdata = X_test)
y_pred_proba <- predict(rf_model, newdata = X_test, type = "prob")

# Convert predictions to numeric for comparison
y_pred_numeric <- as.numeric(y_pred)

cat("Model training complete!\n\n")

# ============================================================================
# 8. MODEL EVALUATION
# ============================================================================

cat("Step 5: Model Evaluation\n")
cat(strrep("-", 60), "\n")

# Calculate metrics
accuracy <- mean(y_pred_numeric == y_test)
cat("Accuracy:", round(accuracy * 100, 2), "%\n")

# Confusion matrix
cm <- table(True = y_test, Predicted = y_pred_numeric)
print("Confusion Matrix:")
print(cm)
cat("\n")

# Classification report
cat("Classification Report:\n")
for (i in 1:length(class_names)) {
  tp <- sum(y_test == i & y_pred_numeric == i)
  fp <- sum(y_test != i & y_pred_numeric == i)
  fn <- sum(y_test == i & y_pred_numeric != i)
  tn <- sum(y_test != i & y_pred_numeric != i)
  
  precision <- tp / (tp + fp) if (tp + fp > 0) else 0
  recall <- tp / (tp + fn) if (tp + fn > 0) else 0
  f1 <- 2 * precision * recall / (precision + recall) if (precision + recall > 0) else 0
  
  cat(sprintf("Class %s (%s): Precision=%.3f, Recall=%.3f, F1=%.3f\n",
              i, class_names[i], precision, recall, f1))
}

# Visualization: Confusion matrix
cm_df <- as.data.frame(cm)
p_cm <- ggplot(cm_df, aes(x = Predicted, y = True, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Confusion Matrix", x = "Predicted Label", y = "True Label") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("r_confusion_matrix.png", p_cm, width = 10, height = 8, dpi = 300)

# ROC curves (multi-class)
library(pROC)

if (length(class_names) > 2) {
  # One-vs-rest approach
  roc_list <- list()
  for (i in 1:length(class_names)) {
    y_test_binary <- as.numeric(y_test == i)
    roc_obj <- roc(y_test_binary, y_pred_proba[, i])
    roc_list[[i]] <- data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Class = class_names[i],
      AUC = as.numeric(auc(roc_obj))
    )
  }
  
  roc_df <- do.call(rbind, roc_list)
  
  p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Class)) +
    geom_line(size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    labs(title = "ROC Curves for Multi-Class Classification",
         x = "False Positive Rate", y = "True Positive Rate") +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")
  
  # Add AUC to legend
  auc_labels <- paste0(roc_df$Class, " (AUC = ", round(roc_df$AUC, 3), ")")
  
  ggsave("r_roc_curves.png", p_roc, width = 10, height = 8, dpi = 300)
}

cat("\nEvaluation plots saved!\n\n")

# ============================================================================
# 9. FEATURE IMPORTANCE ANALYSIS
# ============================================================================

cat("Step 6: Feature Importance Analysis\n")
cat(strrep("-", 60), "\n")

# Extract feature importance from Random Forest
importance_scores <- importance(rf_model)[, "MeanDecreaseGini"]
top_n <- 20
top_indices <- order(importance_scores, decreasing = TRUE)[1:top_n]
top_features <- feature_names[top_indices]
top_importance <- importance_scores[top_indices]

cat("Top", top_n, "Most Important Features:\n")
for (i in 1:top_n) {
  cat(sprintf("%2d. %s: %.6f\n", i, top_features[i], top_importance[i]))
}
cat("\n")

# Visualization
importance_df <- data.frame(
  Feature = top_features,
  Importance = top_importance
)

p_importance <- ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  labs(title = paste("Top", top_n, "Most Important Features (Potential Biomarkers)"),
       x = "Feature", y = "Importance") +
  theme_minimal()

ggsave("r_feature_importance.png", p_importance, width = 12, height = 8, dpi = 300)

cat("Feature importance plot saved!\n\n")

# ============================================================================
# 10. STATISTICAL VALIDATION
# ============================================================================

cat("Step 7: Statistical Validation\n")
cat(strrep("-", 60), "\n")

# Perform statistical analysis
perform_statistical_analysis(X_train, X_test, y_train, y_test, feature_names)

# Bootstrap confidence interval for accuracy
bootstrap_accuracy <- function(y_true, y_pred, n_bootstrap = 1000) {
  n <- length(y_true)
  accuracies <- numeric(n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    indices <- sample(n, size = n, replace = TRUE)
    accuracies[i] <- mean(y_true[indices] == y_pred[indices])
  }
  
  return(accuracies)
}

boot_accs <- bootstrap_accuracy(y_test, y_pred_numeric, n_bootstrap = 1000)
ci_lower <- quantile(boot_accs, 0.025)
ci_upper <- quantile(boot_accs, 0.975)

cat("Bootstrap 95% Confidence Interval for Test Accuracy:\n")
cat(sprintf("  [%.4f, %.4f]\n\n", ci_lower, ci_upper))

# Visualization
boot_df <- data.frame(Accuracy = boot_accs)
p_bootstrap <- ggplot(boot_df, aes(x = Accuracy)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = ci_lower, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = ci_upper, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = accuracy, linetype = "solid", color = "green", size = 1.5) +
  labs(title = "Bootstrap Distribution of Accuracy",
       x = "Accuracy", y = "Frequency") +
  theme_minimal()

ggsave("r_bootstrap_accuracy.png", p_bootstrap, width = 8, height = 6, dpi = 300)

cat("Statistical validation plots saved!\n\n")

# ============================================================================
# 11. SUMMARY AND CONCLUSION
# ============================================================================

cat(paste0(rep("=", 60), collapse=""), "\n")
cat("Analysis completed successfully!\n")
cat(paste0(rep("=", 60), collapse=""), "\n\n")

cat("Results Summary:\n")
cat("  - Accuracy:", round(accuracy * 100, 2), "%\n")
cat("  - Number of samples:", length(y_test), "\n")
cat("  - Number of features:", ncol(X_test), "\n")
cat("  - Number of classes:", length(class_names), "\n")
cat("  - Top feature:", top_features[1], "\n\n")

cat("Files generated:\n")
cat("  - r_expression_distribution.png\n")
cat("  - r_expression_by_type.png\n")
cat("  - r_cancer_type_distribution.png\n")
cat("  - r_survival_distribution.png\n")
cat("  - r_pca_scree.png\n")
cat("  - r_pca_cumulative.png\n")
cat("  - r_pca_2d.png\n")
cat("  - r_confusion_matrix.png\n")
cat("  - r_roc_curves.png\n")
cat("  - r_feature_importance.png\n")
cat("  - r_bootstrap_accuracy.png\n\n")

cat("Note: For deep learning models, consider using Python with PyTorch or TensorFlow.\n")
cat("This R script demonstrates statistical analysis and traditional ML approaches.\n")

# Save session info
sink("r_session_info.txt")
cat("Session Information:\n")
cat("===================\n\n")
print(sessionInfo())
sink()

cat("\nSession info saved to r_session_info.txt\n")

