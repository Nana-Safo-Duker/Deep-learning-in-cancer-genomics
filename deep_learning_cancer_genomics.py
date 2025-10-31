"""
Deep Learning in Cancer Genomics: Comprehensive Analysis Script

This Python script provides a complete pipeline for analyzing cancer genomics data
using deep learning approaches. It includes data preprocessing, model development,
evaluation, and visualization components.

Research Paper: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-024-01315-6

Author: Research Analysis
Date: 2024
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from pathlib import Path
import pickle
import json
from datetime import datetime

warnings.filterwarnings('ignore')

# Machine learning and deep learning imports
try:
    import sklearn
    from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
    from sklearn.preprocessing import StandardScaler, LabelEncoder
    from sklearn.metrics import (
        accuracy_score, classification_report, confusion_matrix,
        roc_auc_score, roc_curve, precision_recall_curve, f1_score
    )
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    from sklearn.inspection import permutation_importance
    from sklearn.ensemble import RandomForestClassifier
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    print("Warning: scikit-learn not available. Install with: pip install scikit-learn")

# Deep learning frameworks
try:
    import torch
    import torch.nn as nn
    import torch.optim as optim
    from torch.utils.data import Dataset, DataLoader
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    print("Note: PyTorch not available. Install with: pip install torch")

try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers, models
    TF_AVAILABLE = True
except ImportError:
    TF_AVAILABLE = False
    print("Note: TensorFlow not available. Install with: pip install tensorflow")

# Statistical analysis
from scipy import stats
from scipy.stats import ttest_ind, mannwhitneyu, ttest_1samp


class CancerGenomicsNet(nn.Module):
    """
    Deep neural network for cancer type classification from gene expression data.
    
    Architecture:
    - Input layer: Gene expression features
    - Hidden layers: Multiple fully connected layers with batch normalization
    - Output layer: Softmax for multi-class classification
    """
    
    def __init__(self, input_dim, hidden_dims=[256, 128, 64], num_classes=5, dropout_rate=0.3):
        super(CancerGenomicsNet, self).__init__()
        
        layers_list = []
        prev_dim = input_dim
        
        for hidden_dim in hidden_dims:
            layers_list.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout_rate)
            ])
            prev_dim = hidden_dim
        
        self.feature_extractor = nn.Sequential(*layers_list)
        self.classifier = nn.Linear(prev_dim, num_classes)
    
    def forward(self, x):
        features = self.feature_extractor(x)
        output = self.classifier(features)
        return output


def simulate_cancer_genomics_data(n_samples=500, n_genes=2000, n_cancer_types=5, random_state=42):
    """
    Simulate multi-omics cancer genomics data for analysis demonstration.
    
    Parameters:
    -----------
    n_samples : int
        Number of samples (patients)
    n_genes : int
        Number of genes in the expression matrix
    n_cancer_types : int
        Number of different cancer types
    random_state : int
        Random seed for reproducibility
    
    Returns:
    --------
    gene_expression : pd.DataFrame
        Simulated gene expression data (RNA-seq like)
    mutations : pd.DataFrame
        Binary mutation matrix
    clinical_data : pd.DataFrame
        Clinical and survival information
    """
    np.random.seed(random_state)
    
    # Simulate gene expression data (log-normal distribution, typical for RNA-seq)
    base_expression = np.random.lognormal(mean=5, sigma=2, size=(n_samples, n_genes))
    
    # Add cancer-type-specific expression patterns
    cancer_types = [f"Cancer_Type_{i}" for i in range(1, n_cancer_types + 1)]
    labels = np.random.choice(cancer_types, size=n_samples)
    
    # Introduce cancer-type-specific differentially expressed genes
    for cancer_type in cancer_types:
        mask = labels == cancer_type
        n_de_genes = n_genes // 10  # 10% of genes are differentially expressed
        de_indices = np.random.choice(n_genes, size=n_de_genes, replace=False)
        base_expression[mask][:, de_indices] *= np.random.uniform(2, 5, size=(mask.sum(), n_de_genes))
    
    gene_expression = pd.DataFrame(
        base_expression,
        columns=[f"Gene_{i}" for i in range(1, n_genes + 1)]
    )
    
    # Simulate mutation data (binary matrix)
    mutation_rate = 0.05  # 5% mutation rate per gene-sample pair
    mutations = pd.DataFrame(
        np.random.binomial(1, mutation_rate, size=(n_samples, n_genes // 10)),
        columns=[f"Gene_{i}_mut" for i in range(1, n_genes // 10 + 1)]
    )
    
    # Simulate clinical data
    ages = np.random.normal(60, 15, n_samples)
    ages = np.clip(ages, 18, 90)  # Realistic age range
    
    # Survival times (exponential distribution with cancer-type-specific rates)
    survival_times = []
    for cancer_type in labels:
        rate = np.random.uniform(0.01, 0.05)  # Different rates per cancer type
        survival_times.append(np.random.exponential(1/rate))
    survival_times = np.array(survival_times)
    
    clinical_data = pd.DataFrame({
        'sample_id': [f"Sample_{i}" for i in range(1, n_samples + 1)],
        'cancer_type': labels,
        'age': ages,
        'survival_time': survival_times,
        'stage': np.random.choice(['I', 'II', 'III', 'IV'], size=n_samples, p=[0.2, 0.3, 0.3, 0.2])
    })
    
    return gene_expression, mutations, clinical_data


def preprocess_data(gene_expr, clinical, n_top_genes=500):
    """
    Preprocess gene expression data for machine learning.
    
    Parameters:
    -----------
    gene_expr : pd.DataFrame
        Raw gene expression data
    clinical : pd.DataFrame
        Clinical data with cancer type labels
    n_top_genes : int
        Number of top variable genes to select
    
    Returns:
    --------
    X_train : np.ndarray
        Training features
    X_test : np.ndarray
        Test features
    y_train : np.ndarray
        Training labels
    y_test : np.ndarray
        Test labels
    scaler : StandardScaler
        Fitted scaler for inverse transformation if needed
    label_encoder : LabelEncoder
        Fitted label encoder
    feature_names : list
        Names of selected features
    """
    # Normalization: Standard scaling
    scaler = StandardScaler()
    gene_expr_scaled = pd.DataFrame(
        scaler.fit_transform(gene_expr),
        columns=gene_expr.columns,
        index=gene_expr.index
    )
    
    # Feature selection: Select top variable genes (variance-based)
    gene_variance = gene_expr_scaled.var(axis=0).sort_values(ascending=False)
    top_genes = gene_variance.head(n_top_genes).index
    gene_expr_selected = gene_expr_scaled[top_genes]
    
    print(f"Selected {n_top_genes} most variable genes from {len(gene_expr.columns)} total genes")
    
    # Encode labels
    label_encoder = LabelEncoder()
    y_encoded = label_encoder.fit_transform(clinical['cancer_type'])
    
    # Train-test split with stratification
    X_train, X_test, y_train, y_test = train_test_split(
        gene_expr_selected.values,
        y_encoded,
        test_size=0.2,
        random_state=42,
        stratify=y_encoded
    )
    
    print(f"Train set size: {X_train.shape[0]}, Test set size: {X_test.shape[0]}")
    print(f"Number of features: {X_train.shape[1]}, Number of classes: {len(label_encoder.classes_)}")
    
    return X_train, X_test, y_train, y_test, scaler, label_encoder, top_genes.tolist()


def train_pytorch_model(X_train, X_test, y_train, y_test, num_classes, epochs=50, batch_size=32):
    """Train a PyTorch neural network model."""
    if not TORCH_AVAILABLE:
        print("PyTorch not available. Skipping PyTorch model training.")
        return None, None, None
    
    # Set random seed
    torch.manual_seed(42)
    np.random.seed(42)
    
    # Convert to tensors
    X_train_tensor = torch.FloatTensor(X_train)
    X_test_tensor = torch.FloatTensor(X_test)
    y_train_tensor = torch.LongTensor(y_train)
    y_test_tensor = torch.LongTensor(y_test)
    
    # Initialize model
    model = CancerGenomicsNet(
        input_dim=X_train.shape[1],
        hidden_dims=[256, 128, 64],
        num_classes=num_classes,
        dropout_rate=0.3
    )
    
    # Define loss and optimizer
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-5)
    
    # Training loop
    train_losses = []
    train_accuracies = []
    
    print(f"Training PyTorch model for {epochs} epochs...")
    
    for epoch in range(epochs):
        model.train()
        epoch_loss = 0
        correct = 0
        total = 0
        
        # Mini-batch training
        for i in range(0, len(X_train_tensor), batch_size):
            batch_X = X_train_tensor[i:i+batch_size]
            batch_y = y_train_tensor[i:i+batch_size]
            
            optimizer.zero_grad()
            outputs = model(batch_X)
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()
            
            epoch_loss += loss.item()
            _, predicted = torch.max(outputs.data, 1)
            total += batch_y.size(0)
            correct += (predicted == batch_y).sum().item()
        
        avg_loss = epoch_loss / (len(X_train_tensor) // batch_size + 1)
        accuracy = 100 * correct / total
        train_losses.append(avg_loss)
        train_accuracies.append(accuracy)
        
        if (epoch + 1) % 10 == 0:
            print(f"Epoch [{epoch+1}/{epochs}], Loss: {avg_loss:.4f}, Accuracy: {accuracy:.2f}%")
    
    # Evaluation
    model.eval()
    with torch.no_grad():
        test_outputs = model(X_test_tensor)
        _, test_predicted = torch.max(test_outputs.data, 1)
        test_accuracy = (test_predicted == y_test_tensor).sum().item() / len(y_test_tensor) * 100
        test_probs = torch.softmax(test_outputs, dim=1).numpy()
    
    print(f"Test Accuracy: {test_accuracy:.2f}%")
    
    return model, test_predicted.numpy(), test_probs


def train_tensorflow_model(X_train, X_test, y_train, y_test, num_classes, epochs=50):
    """Train a TensorFlow/Keras neural network model."""
    if not TF_AVAILABLE:
        print("TensorFlow not available. Skipping TensorFlow model training.")
        return None, None, None
    
    tf.random.set_seed(42)
    np.random.seed(42)
    
    # Build model
    inputs = keras.Input(shape=(X_train.shape[1],))
    
    x = layers.Dense(256, activation='relu')(inputs)
    x = layers.BatchNormalization()(x)
    x = layers.Dropout(0.3)(x)
    
    x = layers.Dense(128, activation='relu')(x)
    x = layers.BatchNormalization()(x)
    x = layers.Dropout(0.3)(x)
    
    x = layers.Dense(64, activation='relu')(x)
    x = layers.BatchNormalization()(x)
    x = layers.Dropout(0.3)(x)
    
    outputs = layers.Dense(num_classes, activation='softmax')(x)
    
    model = keras.Model(inputs=inputs, outputs=outputs)
    
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=0.001),
        loss='sparse_categorical_crossentropy',
        metrics=['accuracy']
    )
    
    print(f"Training TensorFlow/Keras model for {epochs} epochs...")
    
    history = model.fit(
        X_train, y_train,
        batch_size=32,
        epochs=epochs,
        validation_split=0.2,
        verbose=1
    )
    
    # Evaluate on test set
    test_loss, test_accuracy = model.evaluate(X_test, y_test, verbose=0)
    print(f"Test Accuracy: {test_accuracy*100:.2f}%")
    
    # Predictions
    test_predictions = model.predict(X_test, verbose=0)
    test_pred_classes = np.argmax(test_predictions, axis=1)
    
    return model, test_pred_classes, test_predictions


def evaluate_model(y_test, y_pred, y_pred_proba, class_names):
    """Comprehensive model evaluation."""
    accuracy = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred, average='weighted')
    
    # Multi-class ROC AUC
    if len(class_names) > 2:
        roc_auc = roc_auc_score(y_test, y_pred_proba, multi_class='ovr', average='weighted')
    else:
        roc_auc = roc_auc_score(y_test, y_pred_proba[:, 1])
    
    print("\n" + "="*50)
    print("Model Evaluation Metrics:")
    print("="*50)
    print(f"Accuracy: {accuracy:.4f} ({accuracy*100:.2f}%)")
    print(f"F1-Score (weighted): {f1:.4f}")
    print(f"ROC-AUC (weighted, one-vs-rest): {roc_auc:.4f}")
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred, target_names=class_names))
    
    return {
        'accuracy': accuracy,
        'f1_score': f1,
        'roc_auc': roc_auc
    }


def compute_feature_importance(model, X_test, y_test, feature_names, model_type='pytorch'):
    """Compute feature importance using permutation importance."""
    print("Computing permutation importance (this may take a few moments)...")
    
    def model_predict(X):
        if model_type == 'pytorch':
            model.eval()
            with torch.no_grad():
                X_tensor = torch.FloatTensor(X)
                outputs = model(X_tensor)
                _, predicted = torch.max(outputs, 1)
                return predicted.numpy()
        elif model_type == 'tensorflow':
            pred = model.predict(X, verbose=0)
            return np.argmax(pred, axis=1)
    
    perm_importance = permutation_importance(
        model_predict, X_test, y_test, n_repeats=10, random_state=42, n_jobs=-1
    )
    
    feature_importance = perm_importance.importances_mean
    
    # Get top features
    top_n = 20
    top_indices = np.argsort(feature_importance)[-top_n:][::-1]
    top_features = [feature_names[i] for i in top_indices]
    top_importance = feature_importance[top_indices]
    
    print(f"\nTop {top_n} Most Important Features:")
    for i, (feat, imp) in enumerate(zip(top_features, top_importance), 1):
        print(f"{i:2d}. {feat}: {imp:.6f}")
    
    return feature_importance, top_features, top_importance


def visualize_results(X_train, y_train, y_test, y_pred, y_pred_proba, 
                     class_names, feature_importance=None, top_features=None):
    """Create comprehensive visualizations."""
    
    # Set style
    plt.style.use('seaborn-v0_8-darkgrid')
    sns.set_palette("husl")
    
    # 1. Confusion Matrix
    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(10, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                xticklabels=class_names, yticklabels=class_names)
    plt.title('Confusion Matrix')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.tight_layout()
    plt.savefig('confusion_matrix.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. ROC Curves
    from sklearn.preprocessing import label_binarize
    from sklearn.metrics import roc_curve, auc
    
    y_test_binarized = label_binarize(y_test, classes=range(len(class_names)))
    
    plt.figure(figsize=(10, 8))
    for i in range(len(class_names)):
        fpr, tpr, _ = roc_curve(y_test_binarized[:, i], y_pred_proba[:, i])
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, label=f'{class_names[i]} (AUC = {roc_auc:.3f})')
    
    plt.plot([0, 1], [0, 1], 'k--', label='Random')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves for Multi-Class Classification')
    plt.legend(loc='lower right')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('roc_curves.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. PCA Visualization
    pca_2d = PCA(n_components=2)
    X_train_2d = pca_2d.fit_transform(X_train)
    
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(X_train_2d[:, 0], X_train_2d[:, 1],
                         c=y_train, cmap='tab10', alpha=0.6, s=50)
    plt.colorbar(scatter, label='Cancer Type')
    plt.xlabel(f'PC1 ({pca_2d.explained_variance_ratio_[0]*100:.1f}% variance)')
    plt.ylabel(f'PC2 ({pca_2d.explained_variance_ratio_[1]*100:.1f}% variance)')
    plt.title('PCA Visualization: First Two Principal Components')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('pca_visualization.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Feature Importance (if available)
    if feature_importance is not None and top_features is not None:
        top_n = len(top_features)
        top_importance_vals = [feature_importance[np.where([f == feat for f in top_features])[0][0]] 
                              for feat in top_features]
        
        plt.figure(figsize=(12, 8))
        plt.barh(range(len(top_features)), top_importance_vals[::-1])
        plt.yticks(range(len(top_features)), top_features[::-1])
        plt.xlabel('Feature Importance')
        plt.title(f'Top {top_n} Most Important Features (Potential Biomarkers)')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig('feature_importance.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    print("\nVisualizations saved successfully!")


def main():
    """Main analysis pipeline."""
    print("="*60)
    print("Deep Learning in Cancer Genomics: Comprehensive Analysis")
    print("="*60)
    print(f"Analysis started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Set random seeds
    np.random.seed(42)
    if TORCH_AVAILABLE:
        torch.manual_seed(42)
    if TF_AVAILABLE:
        tf.random.set_seed(42)
    
    # 1. Generate/load data
    print("Step 1: Generating simulated cancer genomics data...")
    gene_expr, mut_data, clinical = simulate_cancer_genomics_data(
        n_samples=500, n_genes=2000, n_cancer_types=5
    )
    print(f"Data generated: {gene_expr.shape[0]} samples, {gene_expr.shape[1]} genes\n")
    
    # 2. Preprocess data
    print("Step 2: Preprocessing data...")
    X_train, X_test, y_train, y_test, scaler, label_encoder, feature_names = preprocess_data(
        gene_expr, clinical, n_top_genes=500
    )
    class_names = label_encoder.classes_
    print()
    
    # 3. Train models
    num_classes = len(class_names)
    
    # Try TensorFlow first, then PyTorch, then Random Forest
    model = None
    y_pred = None
    y_pred_proba = None
    model_type = None
    
    if TF_AVAILABLE:
        print("Step 3: Training TensorFlow/Keras model...")
        model, y_pred, y_pred_proba = train_tensorflow_model(
            X_train, X_test, y_train, y_test, num_classes, epochs=50
        )
        model_type = 'tensorflow'
    elif TORCH_AVAILABLE:
        print("Step 3: Training PyTorch model...")
        model, y_pred, y_pred_proba = train_pytorch_model(
            X_train, X_test, y_train, y_test, num_classes, epochs=50
        )
        model_type = 'pytorch'
    else:
        print("Step 3: Training Random Forest model (no deep learning frameworks available)...")
        rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
        rf_model.fit(X_train, y_train)
        y_pred = rf_model.predict(X_test)
        y_pred_proba = rf_model.predict_proba(X_test)
        model = rf_model
        model_type = 'sklearn'
    
    print()
    
    # 4. Evaluate model
    print("Step 4: Evaluating model...")
    metrics = evaluate_model(y_test, y_pred, y_pred_proba, class_names)
    print()
    
    # 5. Feature importance
    feature_importance = None
    top_features = None
    if model is not None and model_type in ['pytorch', 'tensorflow']:
        print("Step 5: Computing feature importance...")
        feature_importance, top_features, _ = compute_feature_importance(
            model, X_test, y_test, feature_names, model_type=model_type
        )
    elif model_type == 'sklearn':
        feature_importance = model.feature_importances_
        top_n = 20
        top_indices = np.argsort(feature_importance)[-top_n:][::-1]
        top_features = [feature_names[i] for i in top_indices]
    print()
    
    # 6. Visualizations
    print("Step 6: Creating visualizations...")
    visualize_results(X_train, y_train, y_test, y_pred, y_pred_proba,
                     class_names, feature_importance, top_features)
    
    # 7. Save results
    print("Step 7: Saving results...")
    results = {
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'metrics': metrics,
        'model_type': model_type,
        'n_samples': len(y_test),
        'n_features': X_test.shape[1],
        'n_classes': num_classes,
        'top_features': top_features[:10] if top_features else None
    }
    
    with open('analysis_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "="*60)
    print("Analysis completed successfully!")
    print("="*60)
    print(f"\nResults saved to:")
    print("  - analysis_results.json (metrics and summary)")
    print("  - confusion_matrix.png")
    print("  - roc_curves.png")
    print("  - pca_visualization.png")
    if top_features:
        print("  - feature_importance.png")


if __name__ == "__main__":
    main()

