import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from sklearn.model_selection import train_test_split
import numpy as np
import os
import matplotlib.pyplot as plt
import tarfile
import io
from Bio import SeqIO
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.preprocessing import StandardScaler
from sklearn.feature_extraction.text import CountVectorizer
from scipy.cluster.hierarchy import dendrogram, linkage
from tqdm import tqdm
import random
import gc
import seaborn as sns
import warnings
import joblib
warnings.filterwarnings('ignore')

# ------------------- SETUP: Device and Seed -------------------
def set_random_seeds(seed=42):
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
set_random_seeds(42)

# ------------------- DATA ANALYSIS -------------------
def analyze_sequences(sequences, labels=None):
    """Analyze sequence characteristics to understand the data"""
    print("Analyzing sequence characteristics...")
    
    # Calculate GC content
    gc_content = []
    for seq in sequences:
        gc_count = seq.count('G') + seq.count('C')
        gc_content.append(gc_count / len(seq))
    
    # Calculate sequence length
    seq_lengths = [len(seq) for seq in sequences]
    
    # Check for motifs
    promoter_motifs = ["TATAAT", "TTGACA", "GGGCGG", "CCGCCC", "ATGCAT", "GCCGCC", "CGCGCG", "TATATA", "GCGCGC", "ATATAT"]
    terminator_motifs = ["AATAAA", "GCGCGC", "TTTTTT", "CCCCCC", "GGGGGG", "ATATAT", "TATATA", "CGCGCG", "GCCGCC", "CATCAT"]
    class_motifs = ["AGCTAGCT", "TCGATCGA", "GCTAGCTA", "CTAGCTAG", "ATCGATCG", "TAGCTAGC", "GATCGATC", "CGATCGAT", "AATTCCGG", "GGCCAATT"]
    
    promoter_counts = {motif: 0 for motif in promoter_motifs}
    terminator_counts = {motif: 0 for motif in terminator_motifs}
    class_counts = {motif: 0 for motif in class_motifs}
    
    for seq in sequences:
        for motif in promoter_motifs:
            if motif in seq[:50]:  # Check first 50 bases for promoter
                promoter_counts[motif] += 1
        for motif in terminator_motifs:
            if motif in seq[-50:]:  # Check last 50 bases for terminator
                terminator_counts[motif] += 1
        for motif in class_motifs:
            if motif in seq:
                class_counts[motif] += 1
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # GC content distribution
    axes[0, 0].hist(gc_content, bins=30, alpha=0.7)
    axes[0, 0].set_title('GC Content Distribution')
    axes[0, 0].set_xlabel('GC Content')
    axes[0, 0].set_ylabel('Frequency')
    
    # Sequence length distribution
    axes[0, 1].hist(seq_lengths, bins=30, alpha=0.7)
    axes[0, 1].set_title('Sequence Length Distribution')
    axes[0, 1].set_xlabel('Sequence Length')
    axes[0, 1].set_ylabel('Frequency')
    
    # Promoter motifs
    axes[1, 0].bar(promoter_counts.keys(), promoter_counts.values())
    axes[1, 0].set_title('Promoter Motif Counts')
    axes[1, 0].set_ylabel('Count')
    plt.setp(axes[1, 0].xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Class motifs
    axes[1, 1].bar(class_counts.keys(), class_counts.values())
    axes[1, 1].set_title('Class-Specific Motif Counts')
    axes[1, 1].set_ylabel('Count')
    plt.setp(axes[1, 1].xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig('sequence_analysis.png')
    plt.close()
    
    print(f"GC Content: {np.mean(gc_content):.3f} ± {np.std(gc_content):.3f}")
    print(f"Sequence Length: {np.mean(seq_lengths):.1f} ± {np.std(seq_lengths):.1f}")
    print("Promoter motifs:", promoter_counts)
    print("Class-specific motifs:", class_counts)
    
    return gc_content, seq_lengths

# ------------------- HYBRID FEATURE EXTRACTION -------------------
def extract_hybrid_features(sequences, vectorizer=None, k=6):
    """
    Extract both k-mer and sequence-based features
    """
    # 1. K-mer features
    if vectorizer is None:
        vectorizer = CountVectorizer(analyzer='char', ngram_range=(k, k))
        kmer_features = vectorizer.fit_transform(sequences).toarray()
    else:
        kmer_features = vectorizer.transform(sequences).toarray()
    
    # 2. Sequence-based features
    seq_features = []
    for seq in sequences:
        gc_content = (seq.count('G') + seq.count('C')) / len(seq)
        at_content = (seq.count('A') + seq.count('T')) / len(seq)
        purine_content = (seq.count('A') + seq.count('G')) / len(seq)
        pyrimidine_content = (seq.count('C') + seq.count('T')) / len(seq)
        
        # Check for specific motifs
        promoter_motifs = ["TATAAT", "TTGACA", "GGGCGG", "CCGCCC", "ATGCAT"]
        terminator_motifs = ["AATAAA", "GCGCGC", "TTTTTT", "CCCCCC", "GGGGGG"]
        class_motifs = ["AGCTAGCT", "TCGATCGA", "GCTAGCTA", "CTAGCTAG", "ATCGATCG"]
        
        motif_features = []
        for motif in promoter_motifs:
            motif_features.append(1 if motif in seq[:50] else 0)
        for motif in terminator_motifs:
            motif_features.append(1 if motif in seq[-50:] else 0)
        for motif in class_motifs:
            motif_features.append(1 if motif in seq else 0)
        
        seq_features.append([
            gc_content, at_content, purine_content, pyrimidine_content,
            *motif_features
        ])
    
    seq_features = np.array(seq_features)
    
    # Combine features
    hybrid_features = np.hstack([kmer_features, seq_features])
    
    print(f"Generated {hybrid_features.shape[1]} hybrid features")
    return hybrid_features, vectorizer

# ------------------- DATA LOADING & ENCODING -------------------
def process_data(tar_path, max_seq_len=500, is_labeled=True, kmer_size=6, vectorizer=None, scaler=None):
    """Process data from tar file using hybrid features"""
    print(f"Loading and preprocessing reads from {tar_path}...")
    all_sequences = []
    all_labels = []
    label_map = {}
    
    try:
        with tarfile.open(tar_path, "r") as tar:
            for member in tar.getmembers():
                if member.name.endswith('.fasta'):
                    f = tar.extractfile(member)
                    if f is None:
                        continue
                    with io.TextIOWrapper(f, encoding='utf-8') as text_stream:
                        for record in SeqIO.parse(text_stream, "fasta"):
                            seq = str(record.seq).upper()
                            # Pad or truncate to max_seq_len
                            if len(seq) > max_seq_len:
                                seq = seq[:max_seq_len]
                            else:
                                seq = seq + 'N' * (max_seq_len - len(seq))
                            
                            all_sequences.append(seq)
                            
                            if is_labeled:
                                # Extract label from header (format: "Class_X_read_Y")
                                label_name = record.id.split('_')[0] + '_' + record.id.split('_')[1]
                                if label_name not in label_map:
                                    label_map[label_name] = len(label_map)
                                all_labels.append(label_map[label_name])
    except Exception as e:
        print(f"Warning: Could not process {tar_path}. Reason: {e}")
        return None, None, None, None, None, None

    if not all_sequences:
        raise ValueError("No valid sequences were extracted. Check your .tar files and contents.")

    # Analyze sequences
    analyze_sequences(all_sequences, all_labels if is_labeled else None)
    
    # Convert to hybrid features
    hybrid_features, vectorizer = extract_hybrid_features(all_sequences, vectorizer=vectorizer, k=kmer_size)
    
    # Standardize features
    if scaler is None:
        scaler = StandardScaler()
        hybrid_features = scaler.fit_transform(hybrid_features)
    else:
        hybrid_features = scaler.transform(hybrid_features)
    
    X_tensor = torch.from_numpy(hybrid_features).float()
    
    if is_labeled:
        y = np.asarray(all_labels, dtype=np.int64)
        y_tensor = torch.from_numpy(y)
        print(f"Data loaded: {X_tensor.shape[0]} sequences, {len(label_map)} classes.")
        return X_tensor, y_tensor, len(label_map), vectorizer, scaler
    else:
        print(f"Data loaded: {X_tensor.shape[0]} novel sequences.")
        return X_tensor, None, None, vectorizer, scaler

# ------------------- IMPROVED MODEL DEFINITION -------------------
class EnhancedClassifier(nn.Module):
    def __init__(self, input_size, num_classes, hidden_sizes=[512, 256, 128]):
        super().__init__()
        layers = []
        
        # Create hidden layers
        prev_size = input_size
        for hidden_size in hidden_sizes:
            layers.extend([
                nn.Linear(prev_size, hidden_size),
                nn.BatchNorm1d(hidden_size),
                nn.ReLU(),
                nn.Dropout(0.3)
            ])
            prev_size = hidden_size
        
        # Output layer
        layers.append(nn.Linear(prev_size, num_classes))
        
        self.classifier = nn.Sequential(*layers)
        
    def forward(self, x):
        return self.classifier(x)

# ------------------- TRAINING -------------------
def train_model(model, train_loader, val_loader, num_epochs=50, patience=10):
    print("Starting supervised training...")
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-5)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'max', factor=0.5, patience=5)
    
    best_accuracy = 0.0
    epochs_no_improve = 0
    best_model_path = ""
    
    train_losses = []
    val_losses = []
    train_accuracies = []
    val_accuracies = []
    
    for epoch in range(num_epochs):
        # Training phase
        model.train()
        running_loss = 0.0
        correct_train = 0
        total_train = 0
        
        for inputs, labels in tqdm(train_loader, desc=f"Epoch {epoch+1}/{num_epochs}"):
            inputs, labels = inputs.to(DEVICE), labels.to(DEVICE)
            
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            
            running_loss += loss.item() * inputs.size(0)
            _, predicted = torch.max(outputs, 1)
            total_train += labels.size(0)
            correct_train += (predicted == labels).sum().item()
        
        epoch_loss = running_loss / len(train_loader.dataset)
        epoch_acc = 100 * correct_train / total_train
        train_losses.append(epoch_loss)
        train_accuracies.append(epoch_acc)
        
        # Validation phase
        model.eval()
        val_loss = 0.0
        correct_val = 0
        total_val = 0
        
        with torch.no_grad():
            for inputs, labels in val_loader:
                inputs, labels = inputs.to(DEVICE), labels.to(DEVICE)
                outputs = model(inputs)
                loss = criterion(outputs, labels)
                
                val_loss += loss.item() * inputs.size(0)
                _, predicted = torch.max(outputs, 1)
                total_val += labels.size(0)
                correct_val += (predicted == labels).sum().item()
        
        val_epoch_loss = val_loss / len(val_loader.dataset)
        val_epoch_acc = 100 * correct_val / total_val
        val_losses.append(val_epoch_loss)
        val_accuracies.append(val_epoch_acc)
        
        print(f"Epoch [{epoch+1}/{num_epochs}] | Train Loss: {epoch_loss:.4f} | Train Acc: {epoch_acc:.2f}% | "
              f"Val Loss: {val_epoch_loss:.4f} | Val Acc: {val_epoch_acc:.2f}%")
        
        # Save best model
        if val_epoch_acc > best_accuracy:
            if os.path.exists(best_model_path):
                try:
                    os.remove(best_model_path)
                except OSError:
                    pass
            
            best_accuracy = val_epoch_acc
            best_model_path = f"best_model_epoch_{epoch+1}_acc_{best_accuracy:.2f}.pth"
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'accuracy': best_accuracy,
            }, best_model_path)
            
            print(f"New best model saved to: {best_model_path}")
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1
            if epochs_no_improve >= patience:
                print(f"Early stopping (no improvement in {patience} epochs).")
                break
        
        scheduler.step(val_epoch_acc)
        torch.cuda.empty_cache() if torch.cuda.is_available() else None
        gc.collect()
    
    # Plot training history
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.plot(train_losses, label='Training Loss')
    plt.plot(val_losses, label='Validation Loss')
    plt.title('Loss Over Time')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    
    plt.subplot(1, 2, 2)
    plt.plot(train_accuracies, label='Training Accuracy')
    plt.plot(val_accuracies, label='Validation Accuracy')
    plt.title('Accuracy Over Time')
    plt.xlabel('Epoch')
    plt.ylabel('Accuracy (%)')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('training_history.png')
    plt.close()
    
    print(f"Training finished. Best Val Acc: {best_accuracy:.2f}%")
    return best_model_path

# ------------------- COMPREHENSIVE CLUSTERING ANALYSIS -------------------
def comprehensive_clustering_analysis(features, max_k=15):
    """
    Perform comprehensive clustering analysis using multiple methods
    """
    print("Performing comprehensive clustering analysis...")
    
    # Standardize features
    scaler = StandardScaler()
    features_std = scaler.fit_transform(features)
    
    # Reduce dimensionality for visualization
    pca = PCA(n_components=min(50, features_std.shape[1]))
    features_pca = pca.fit_transform(features_std)
    
    # 1. Determine optimal number of clusters using multiple methods
    print("Determining optimal number of clusters...")
    
    # Elbow method
    inertias = []
    k_range = range(2, max_k + 1)
    
    for k in tqdm(k_range, desc="Elbow Method"):
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        kmeans.fit(features_std)
        inertias.append(kmeans.inertia_)
    
    # Silhouette score
    silhouette_scores = []
    for k in tqdm(k_range, desc="Silhouette Score"):
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(features_std)
        if len(set(labels)) > 1:  # Need at least 2 clusters for silhouette score
            silhouette_scores.append(silhouette_score(features_std, labels))
        else:
            silhouette_scores.append(0)
    
    # Calinski-Harabasz score
    ch_scores = []
    for k in tqdm(k_range, desc="Calinski-Harabasz Score"):
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(features_std)
        if len(set(labels)) > 1:
            ch_scores.append(calinski_harabasz_score(features_std, labels))
        else:
            ch_scores.append(0)
    
    # Davies-Bouldin score
    db_scores = []
    for k in tqdm(k_range, desc="Davies-Bouldin Score"):
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(features_std)
        if len(set(labels)) > 1:
            db_scores.append(davies_bouldin_score(features_std, labels))
        else:
            db_scores.append(1)
    
    # Plot all metrics
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Elbow plot
    axes[0, 0].plot(k_range, inertias, marker='o')
    axes[0, 0].set_title('Elbow Method')
    axes[0, 0].set_xlabel('Number of clusters')
    axes[0, 0].set_ylabel('Inertia')
    axes[0, 0].grid(True)
    
    # Silhouette score
    axes[0, 1].plot(k_range, silhouette_scores, marker='o')
    axes[0, 1].set_title('Silhouette Score')
    axes[0, 1].set_xlabel('Number of clusters')
    axes[0, 1].set_ylabel('Score')
    axes[0, 1].grid(True)
    
    # Calinski-Harabasz score
    axes[1, 0].plot(k_range, ch_scores, marker='o')
    axes[1, 0].set_title('Calinski-Harabasz Score')
    axes[1, 0].set_xlabel('Number of clusters')
    axes[1, 0].set_ylabel('Score')
    axes[1, 0].grid(True)
    
    # Davies-Bouldin score
    axes[1, 1].plot(k_range, db_scores, marker='o')
    axes[1, 1].set_title('Davies-Bouldin Score (lower is better)')
    axes[1, 1].set_xlabel('Number of clusters')
    axes[1, 1].set_ylabel('Score')
    axes[1, 1].grid(True)
    
    plt.tight_layout()
    plt.savefig('clustering_metrics.png', dpi=300)
    plt.close()
    
    # Find optimal k based on different metrics
    optimal_k_elbow = find_elbow_point(inertias, k_range)
    optimal_k_silhouette = k_range[np.argmax(silhouette_scores)]
    optimal_k_ch = k_range[np.argmax(ch_scores)]
    optimal_k_db = k_range[np.argmin(db_scores)]
    
    print(f"Optimal k from elbow method: {optimal_k_elbow}")
    print(f"Optimal k from silhouette score: {optimal_k_silhouette}")
    print(f"Optimal k from Calinski-Harabasz: {optimal_k_ch}")
    print(f"Optimal k from Davies-Bouldin: {optimal_k_db}")
    
    # 2. Try different clustering algorithms
    clustering_results = {}
    
    # K-means with different k values
    for k in [optimal_k_elbow, optimal_k_silhouette, 3, 10]:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(features_std)
        clustering_results[f'KMeans_k={k}'] = labels
        
        # Evaluate
        if len(set(labels)) > 1:
            sil_score = silhouette_score(features_std, labels)
            print(f"KMeans (k={k}) Silhouette Score: {sil_score:.3f}")
    
    # DBSCAN with different parameters
    dbscan_params = [
        {'eps': 0.5, 'min_samples': 5},
        {'eps': 1.0, 'min_samples': 10},
        {'eps': 2.0, 'min_samples': 15}
    ]
    
    for i, params in enumerate(dbscan_params):
        dbscan = DBSCAN(eps=params['eps'], min_samples=params['min_samples'])
        labels = dbscan.fit_predict(features_std)
        clustering_results[f'DBSCAN_{i}'] = labels
        
        # Count clusters (excluding noise)
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        print(f"DBSCAN (eps={params['eps']}, min_samples={params['min_samples']}) found {n_clusters} clusters")
    
    # Hierarchical clustering
    for k in [3, 5, 10]:
        agg = AgglomerativeClustering(n_clusters=k)
        labels = agg.fit_predict(features_std)
        clustering_results[f'Agglomerative_k={k}'] = labels
        
        if len(set(labels)) > 1:
            sil_score = silhouette_score(features_std, labels)
            print(f"Agglomerative (k={k}) Silhouette Score: {sil_score:.3f}")
    
    # Gaussian Mixture Model
    for k in [3, 5, 10]:
        gmm = GaussianMixture(n_components=k, random_state=42)
        labels = gmm.fit_predict(features_std)
        clustering_results[f'GMM_k={k}'] = labels
        
        if len(set(labels)) > 1:
            sil_score = silhouette_score(features_std, labels)
            print(f"GMM (k={k}) Silhouette Score: {sil_score:.3f}")
    
    # 3. Visualize different clustering results
    plot_clustering_results(features_pca, clustering_results)
    
    # 4. Create dendrogram for hierarchical clustering
    print("Creating dendrogram...")
    linked = linkage(features_std, 'ward')
    
    plt.figure(figsize=(12, 8))
    dendrogram(linked, orientation='top', distance_sort='descending', show_leaf_counts=True)
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('Sample index')
    plt.ylabel('Distance')
    plt.tight_layout()
    plt.savefig('dendrogram.png', dpi=300)
    plt.close()
    
    # 5. Compare cluster characteristics
    compare_cluster_characteristics(features, clustering_results)
    
    return clustering_results

def find_elbow_point(inertias, k_range):
    """Find the elbow point in the inertia curve"""
    # Calculate the second derivative
    derivatives = np.diff(inertias)
    second_derivatives = np.diff(derivatives)
    
    # Find the point with the maximum curvature (elbow)
    elbow_point = np.argmin(second_derivatives) + 2  # +2 because we start from k=2
    
    return min(elbow_point, len(k_range))

def plot_clustering_results(features_2d, clustering_results):
    """Plot results from different clustering algorithms"""
    n_methods = len(clustering_results)
    n_cols = 3
    n_rows = (n_methods + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
    axes = axes.flatten()
    
    for i, (method, labels) in enumerate(clustering_results.items()):
        scatter = axes[i].scatter(features_2d[:, 0], features_2d[:, 1], c=labels, 
                                 cmap='viridis', s=10, alpha=0.7)
        axes[i].set_title(f'{method} (Clusters: {len(set(labels))})')
        axes[i].set_xlabel('PC1')
        axes[i].set_ylabel('PC2')
        plt.colorbar(scatter, ax=axes[i])
    
    # Hide any unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('clustering_comparison.png', dpi=300)
    plt.close()

def compare_cluster_characteristics(features, clustering_results):
    """Compare characteristics of clusters from different methods"""
    # Select a few methods to compare
    methods_to_compare = ['KMeans_k=3', 'KMeans_k=10', 'DBSCAN_0', 'Agglomerative_k=3']
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    axes = axes.flatten()
    
    for i, method in enumerate(methods_to_compare):
        if method in clustering_results:
            labels = clustering_results[method]
            
            # Calculate cluster sizes
            unique, counts = np.unique(labels, return_counts=True)
            
            axes[i].bar(unique, counts)
            axes[i].set_title(f'{method} - Cluster Sizes')
            axes[i].set_xlabel('Cluster ID')
            axes[i].set_ylabel('Number of Sequences')
            axes[i].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('cluster_size_comparison.png', dpi=300)
    plt.close()

# ------------------- FEATURE EXTRACTION -------------------
def extract_features(model, data_loader, device):
    """Extract features from the model before the classification layer"""
    model.eval()
    features = []
    
    with torch.no_grad():
        for inputs in tqdm(data_loader, desc="Feature Extraction"):
            if isinstance(inputs, (list, tuple)):
                inputs = inputs[0]
            inputs = inputs.to(device)
            
            # Get features from the layer before the classification head
            # Assuming the model has a classifier attribute with sequential layers
            x = inputs
            for layer in model.classifier[:-1]:  # All layers except the last one
                x = layer(x)
            
            features.append(x.cpu().numpy())
    
    return np.vstack(features)

# ------------------- MAIN EXECUTION -------------------
if __name__ == '__main__':
    SEQUENCE_LENGTH = 500
    KMER_SIZE = 6
    data_dir = "sample_data"
    novel_data_dir = "novel_data"
    
    # Load training data and fit transformers
    try:
        X_tensor, y_tensor, num_classes, vectorizer, scaler = process_data(
            os.path.join(data_dir, "sample_reads.tar"), 
            max_seq_len=SEQUENCE_LENGTH, 
            is_labeled=True,
            kmer_size=KMER_SIZE
        )
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}")
        exit()
    
    # Save the vectorizer and scaler for future use
    joblib.dump(vectorizer, 'dna_vectorizer.joblib')
    joblib.dump(scaler, 'dna_scaler.joblib')
    print("Saved vectorizer and scaler to disk.")
    
    # Split data
    X_train, X_val, y_train, y_val = train_test_split(
        X_tensor, y_tensor, test_size=0.2, random_state=42, stratify=y_tensor
    )
    
    # Create data loaders
    batch_size = 128
    train_loader = DataLoader(
        TensorDataset(X_train, y_train), 
        batch_size=batch_size, 
        shuffle=True, 
        pin_memory=True
    )
    val_loader = DataLoader(
        TensorDataset(X_val, y_val), 
        batch_size=batch_size, 
        shuffle=False, 
        pin_memory=True
    )
    
    # Initialize model
    model = EnhancedClassifier(
        input_size=X_tensor.shape[1], 
        num_classes=num_classes,
        hidden_sizes=[512, 256, 128]
    ).to(DEVICE)
    
    print(f"Model has {sum(p.numel() for p in model.parameters()):,} parameters")
    
    # Train model
    best_model_path = train_model(model, train_loader, val_loader, num_epochs=50, patience=10)
    
    # Clean up memory
    del X_train, X_val, y_train, y_val, train_loader, val_loader
    gc.collect()
    torch.cuda.empty_cache() if torch.cuda.is_available() else None

    # ---- Unsupervised Clustering for Novel Taxa Discovery ----
    if not best_model_path or not os.path.exists(best_model_path):
        print("No valid model was trained or saved. Exiting clustering step.")
        exit()
    
    print("\nLoading the best model for feature extraction...")
    # Load the best model
    checkpoint = torch.load(best_model_path, map_location=DEVICE)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    
    # Load novel data using the same vectorizer and scaler
    try:
        novel_sequences_tensor, _, _, _, _ = process_data(
            os.path.join(novel_data_dir, "sample_reads.tar"), 
            max_seq_len=SEQUENCE_LENGTH, 
            is_labeled=False,
            kmer_size=KMER_SIZE,
            vectorizer=vectorizer,  # Use the same vectorizer
            scaler=scaler           # Use the same scaler
        )
    except (ValueError, FileNotFoundError) as e:
        print(f"Error processing novel data: {e}")
        exit()
    
    # Create data loader for novel sequences
    novel_loader = DataLoader(
        TensorDataset(novel_sequences_tensor), 
        batch_size=batch_size, 
        shuffle=False, 
        pin_memory=True
    )
    
    # Extract features
    print("Extracting features from novel sequences...")
    features = extract_features(model, novel_loader, DEVICE)
    
    # Clean up memory
    del novel_sequences_tensor, novel_loader
    gc.collect()
    torch.cuda.empty_cache() if torch.cuda.is_available() else None
    
    # Perform comprehensive clustering analysis
    if features.shape[0] > 10:
        clustering_results = comprehensive_clustering_analysis(features, max_k=15)
        
        # Also try t-SNE for better visualization
        print("Running t-SNE for better visualization...")
        tsne = TSNE(n_components=2, random_state=42, perplexity=30)
        tsne_result = tsne.fit_transform(features)
        
        # Plot t-SNE results for the best clustering method
        best_method = None
        best_score = -1
        
        for method, labels in clustering_results.items():
            if len(set(labels)) > 1:
                score = silhouette_score(features, labels)
                if score > best_score:
                    best_score = score
                    best_method = method
        
        if best_method:
            best_labels = clustering_results[best_method]
            
            plt.figure(figsize=(12, 10))
            scatter = plt.scatter(tsne_result[:, 0], tsne_result[:, 1], 
                                 c=best_labels, cmap='viridis', s=20, alpha=0.7)
            plt.title(f't-SNE Visualization - {best_method} (Silhouette: {best_score:.3f})')
            plt.xlabel('t-SNE Dimension 1')
            plt.ylabel('t-SNE Dimension 2')
            plt.colorbar(scatter, label='Cluster')
            plt.grid(True, linestyle='--', alpha=0.4)
            plt.tight_layout()
            plt.savefig("tsne_best_clustering.png", dpi=300)
            plt.close()
            
            print(f"Best clustering method: {best_method} with silhouette score {best_score:.3f}")
            
            # Analyze cluster composition
            unique, counts = np.unique(best_labels, return_counts=True)
            print("Cluster sizes for best method:")
            for cluster_id, count in zip(unique, counts):
                print(f"  Cluster {cluster_id}: {count} sequences")
    else:
        print("Not enough novel sequences to perform clustering.")