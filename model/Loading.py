import torch
import torch.nn as nn
import numpy as np
import joblib  # Add this import
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import plotly.express as px
import plotly.io as pio
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import StandardScaler
from Bio import SeqIO
import tarfile
import io

# --- 1. Load the Trained Model and Define the Feature Extractor ---
def load_and_create_feature_extractor(model_path, input_size, num_classes):
    """
    Loads full model (architecture + weights) and returns a feature extractor for clustering.
    This is for the EnhancedClassifier model architecture.
    """
    # Model definition matching the training code
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

    # Load the model checkpoint
    checkpoint = torch.load(model_path, map_location=torch.device("cpu"))
    
    # Create model instance
    model = EnhancedClassifier(input_size=input_size, num_classes=num_classes)
    
    # Load state dict
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()

    # Create feature extractor (all layers except the final classification layer)
    feature_extractor = nn.Sequential(*list(model.classifier.children())[:-1])
    feature_extractor.eval()
    
    return feature_extractor, model

# --- 2. Hybrid Feature Extraction (Matches Training) ---
def extract_hybrid_features(sequences, vectorizer, scaler, k=6):
    """
    Extract both k-mer and sequence-based features using pre-fitted transformers
    """
    # 1. K-mer features
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
    
    # Standardize using the pre-fitted scaler
    hybrid_features = scaler.transform(hybrid_features)
    
    print(f"Generated {hybrid_features.shape[1]} hybrid features")
    return hybrid_features

# --- 3. Load Raw eDNA Reads ---
def load_reads_from_tar(tar_path, max_reads=10000, max_seq_len=500):
    """
    Loads reads from a tar file containing FASTA files
    """
    sequences = []
    try:
        with tarfile.open(tar_path, "r") as tar:
            for member in tar.getmembers():
                if member.name.endswith('.fasta'):
                    f = tar.extractfile(member)
                    if f is None:
                        continue
                    with io.TextIOWrapper(f, encoding='utf-8') as text_stream:
                        for i, record in enumerate(SeqIO.parse(text_stream, "fasta")):
                            if i >= max_reads:
                                break
                            seq = str(record.seq).upper()
                            # Pad or truncate to max_seq_len
                            if len(seq) > max_seq_len:
                                seq = seq[:max_seq_len]
                            else:
                                seq = seq + 'N' * (max_seq_len - len(seq))
                            sequences.append(seq)
    except Exception as e:
        print(f"Error loading from tar: {e}")
    
    print(f"Loaded {len(sequences)} sequences")
    return sequences

# --- 4. Clustering, Biodiversity Estimation & Visualization ---
def cluster_and_visualize_features(features, n_clusters=10, sequences=None):
    """
    Clusters features, estimates biodiversity (number of clusters), and creates interactive visualization.
    """
    print("Clustering features with KMeans...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10).fit(features)
    cluster_labels = kmeans.labels_

    print("Running PCA for visualization...")
    pca = PCA(n_components=3)
    reduced = pca.fit_transform(features)

    # Interactive 3D plot
    fig = px.scatter_3d(
        x=reduced[:, 0], y=reduced[:, 1], z=reduced[:, 2],
        color=cluster_labels.astype(str),
        title=f"eDNA Sequence Clusters (K={n_clusters})",
        labels={"color": "Cluster"},
        hover_data={"Sequence": sequences} if sequences is not None else None
    )
    pio.write_html(fig, "interactive_edna_clusters.html")
    print("Interactive cluster plot saved to interactive_edna_clusters.html")

    # 2D plot for backup
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(reduced[:, 0], reduced[:, 1], c=cluster_labels, cmap='viridis', alpha=0.7)
    plt.colorbar(scatter, label='Cluster')
    plt.title(f'eDNA Sequence Clusters (K={n_clusters})')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.savefig('edna_clusters_2d.png', dpi=300)
    plt.close()
    print("2D cluster plot saved to edna_clusters_2d.png")

    # Biodiversity estimation (number of clusters found)
    unique_clusters = np.unique(cluster_labels)
    print(f"Estimated biodiversity (number of clusters): {len(unique_clusters)}")
    
    # Cluster sizes
    cluster_sizes = [np.sum(cluster_labels == i) for i in range(n_clusters)]
    print("Cluster sizes:", cluster_sizes)
    
    return cluster_labels, reduced

# --- 5. Pipeline Driver ---
if __name__ == "__main__":
    # Update these paths to match your files
    MODEL_PATH = "best_model_epoch_2_acc_98.50.pth"  # Trained model from your supervised pipeline
    TAR_PATH = "novel_data/sample_reads.tar"         # Raw eDNA reads file (tar format)
    VECTORIZER_PATH = "dna_vectorizer.joblib"        # Saved vectorizer
    SCALER_PATH = "dna_scaler.joblib"                # Saved scaler
    
    # Parameters that match your training setup
    INPUT_SIZE = 4115        # Number of hybrid features (from your output: "Generated 4115 hybrid features")
    NUM_CLASSES = 10         # Number of classes in your trained model
    KMER_SIZE = 6            # k-mer size used in training
    N_CLUSTERS = 10          # Number of clusters to find (based on your analysis)
    MAX_READS = 5000         # Maximum number of reads to process

    # Load the pre-fitted vectorizer and scaler
    print("Loading pre-fitted vectorizer and scaler...")
    vectorizer = joblib.load(VECTORIZER_PATH)
    scaler = joblib.load(SCALER_PATH)

    # Load feature extractor
    print("Loading model and creating feature extractor...")
    feature_extractor, model = load_and_create_feature_extractor(
        MODEL_PATH, INPUT_SIZE, NUM_CLASSES
    )

    # Load raw, unlabeled eDNA reads
    print("Loading sequences from tar file...")
    sequences = load_reads_from_tar(TAR_PATH, max_reads=MAX_READS, max_seq_len=500)
    
    # Extract hybrid features using the pre-fitted transformers
    print("Extracting hybrid features with pre-fitted transformers...")
    hybrid_features = extract_hybrid_features(sequences, vectorizer, scaler, k=KMER_SIZE)
    
    # Convert to tensor
    features_tensor = torch.tensor(hybrid_features, dtype=torch.float32)
    
    # Feature extraction using the model
    print("Extracting deep features from hybrid features...")
    with torch.no_grad():
        # Get features from the model (output of the penultimate layer)
        deep_features = feature_extractor(features_tensor).numpy()
    
    # Clustering, biodiversity, and interactive visualization
    cluster_labels, reduced = cluster_and_visualize_features(
        deep_features, n_clusters=N_CLUSTERS, sequences=sequences
    )

    print("Pipeline complete. Explore your results in 'interactive_edna_clusters.html'.")