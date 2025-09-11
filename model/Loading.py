import os
import torch
import torch.nn as nn
import numpy as np
import joblib
import csv  # Import CSV module
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
import sys


# --- 1. Load the Trained Model and Define the Feature Extractor ---
def load_and_create_feature_extractor(model_path, input_size, num_classes):
    """
    Loads full model (architecture + weights) and returns both the full model and feature extractor.
    This is for the EnhancedClassifier model architecture.
    """
    # Check if model file exists
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found at: {model_path}")
    
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
def cluster_and_visualize_features(features, n_clusters=10, sequences=None, title_suffix=""):
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
        title=f"eDNA Sequence Clusters {title_suffix} (K={n_clusters})",
        labels={"color": "Cluster"},
        hover_data={"Sequence": sequences} if sequences is not None else None
    )
    # filename = f"interactive_edna_clusters{title_suffix.replace(' ', '_')}.html"
    # pio.write_html(fig, filename)
    # print(f"Interactive cluster plot saved to {filename}")
    #
    # # 2D plot for backup
    # plt.figure(figsize=(10, 8))
    # scatter = plt.scatter(reduced[:, 0], reduced[:, 1], c=cluster_labels, cmap='viridis', alpha=0.7)
    # plt.colorbar(scatter, label='Cluster')
    # plt.title(f'eDNA Sequence Clusters {title_suffix} (K={n_clusters})')
    # plt.xlabel('Principal Component 1')
    # plt.ylabel('Principal Component 2')
    # filename = f"edna_clusters_2d{title_suffix.replace(' ', '_')}.png"
    # plt.savefig(filename, dpi=300)
    # plt.close()
    # print(f"2D cluster plot saved to {filename}")
# create safe suffix (remove parentheses and replace spaces/dashes)
    safe_suffix = title_suffix.strip()
    safe_suffix = safe_suffix.replace(' ', '').replace('(', '').replace(')', '').replace('-', '')

    # File names
    html_filename = f"interactive_edna_clusters_{safe_suffix}.html"
    png_filename = f"edna_clusters_2d_{safe_suffix}.png"

    # Save files
    pio.write_html(fig, html_filename)
    print(f"Interactive cluster plot saved to {html_filename}")

    plt.savefig(png_filename, dpi=300)
    plt.close()
    print(f"2D cluster plot saved to {png_filename}")
    # Biodiversity estimation (number of clusters found)
    unique_clusters = np.unique(cluster_labels)
    print(f"Estimated biodiversity (number of clusters): {len(unique_clusters)}")
    
    # Cluster sizes
    cluster_sizes = [np.sum(cluster_labels == i) for i in range(n_clusters)]
    print("Cluster sizes:", cluster_sizes)
    
    return cluster_labels, reduced

# --- 5. Save results to CSV ---
def save_results_to_csv(filename, data, headers):
    """
    Save results to a CSV file
    """
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(data)
    print(f"Results saved to {filename}")

# --- 6. Pipeline Driver ---
if __name__ == "__main__":
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Update these paths to match your files
    MODEL_PATH = os.path.join(script_dir, "best_model_epoch_2_acc_98.50.pth")  # Trained model from your supervised pipeline
    # TAR_PATH = os.path.join(script_dir, "novel_data", "sample_reads.tar")      # Raw eDNA reads file (tar format)
    # TAR_PATH = "mixed_test_data/mixed_sample_reads.tar"
    # Default fallback (keeps backward compatibility)
    DEFAULT_TAR = os.path.join(script_dir, "mixed_test_data", "mixed_sample_reads.tar")
    # If server passes uploaded path as arg, use it; otherwise fall back to default.
    TAR_PATH = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_TAR

    VECTORIZER_PATH = os.path.join(script_dir, "dna_vectorizer.joblib")        # Saved vectorizer
    SCALER_PATH = os.path.join(script_dir, "dna_scaler.joblib")                # Saved scaler
    
    # Parameters that match your training setup
    INPUT_SIZE = 4115        # Number of hybrid features (from your output: "Generated 4115 hybrid features")
    NUM_CLASSES = 10         # Number of classes in your trained model
    KMER_SIZE = 6            # k-mer size used in training
    N_CLUSTERS = 10          # Number of clusters to find (based on your analysis)
    MAX_READS = 5000         # Maximum number of reads to process
    CONFIDENCE_THRESHOLD = 0.8  # Threshold for considering a prediction as "known"

    # Check if required files exist
    for file_path, description in [
        (MODEL_PATH, "Model file"),
        (TAR_PATH, "Tar file with sequences"),
        (VECTORIZER_PATH, "Vectorizer file"),
        (SCALER_PATH, "Scaler file")
    ]:
        if not os.path.exists(file_path):
            print(f"Warning: {description} not found at: {file_path}")

    # Load the pre-fitted vectorizer and scaler
    print("Loading pre-fitted vectorizer and scaler...")
    try:
        vectorizer = joblib.load(VECTORIZER_PATH)
        scaler = joblib.load(SCALER_PATH)
    except Exception as e:
        print(f"Error loading vectorizer or scaler: {e}")
        exit(1)

    # Load feature extractor and full model
    print("Loading model and creating feature extractor...")
    try:
        feature_extractor, model = load_and_create_feature_extractor(
            MODEL_PATH, INPUT_SIZE, NUM_CLASSES
        )
    except FileNotFoundError as e:
        print(e)
        exit(1)
    except Exception as e:
        print(f"Error loading model: {e}")
        exit(1)

    # Load raw, unlabeled eDNA reads
    print("Loading sequences from tar file...")
    sequences = load_reads_from_tar(TAR_PATH, max_reads=MAX_READS, max_seq_len=500)
    
    # Extract hybrid features using the pre-fitted transformers
    print("Extracting hybrid features with pre-fitted transformers...")
    hybrid_features = extract_hybrid_features(sequences, vectorizer, scaler, k=KMER_SIZE)
    
    # Convert to tensor
    features_tensor = torch.tensor(hybrid_features, dtype=torch.float32)
    
    # Get predictions from the full model
    print("Getting predictions from the supervised model...")
    with torch.no_grad():
        # Get class predictions and confidence scores
        outputs = model(features_tensor)
        probabilities = torch.softmax(outputs, dim=1)
        confidence, predicted_classes = torch.max(probabilities, 1)
    
    # Convert to numpy arrays
    confidence = confidence.numpy()
    predicted_classes = predicted_classes.numpy()
    
    # Separate known and unknown sequences based on confidence threshold
    known_indices = confidence >= CONFIDENCE_THRESHOLD
    unknown_indices = confidence < CONFIDENCE_THRESHOLD
    
    known_sequences = [seq for i, seq in enumerate(sequences) if known_indices[i]]
    unknown_sequences = [seq for i, seq in enumerate(sequences) if unknown_indices[i]]
    
    known_classes = predicted_classes[known_indices]
    known_confidence = confidence[known_indices]
    
    print(f"Found {len(known_sequences)} known sequences (confidence >= {CONFIDENCE_THRESHOLD})")
    print(f"Found {len(unknown_sequences)} unknown sequences (confidence < {CONFIDENCE_THRESHOLD})")
    
    # Save known sequences to CSV
    known_data = []
    for i, (seq, cls, conf) in enumerate(zip(known_sequences, known_classes, known_confidence)):
        known_data.append([i+1, seq, cls, conf])
    
    save_results_to_csv(
        "known_sequences.csv",
        known_data,
        ["Index", "Sequence", "Predicted_Class", "Confidence"]
    )
    
    # Feature extraction using the model for unknown sequences
    print("\nExtracting deep features for unknown sequences...")
    unknown_features_tensor = features_tensor[unknown_indices]
    
    with torch.no_grad():
        # Get features from the model (output of the penultimate layer)
        deep_features = feature_extractor(unknown_features_tensor).numpy()
    
    # Clustering, biodiversity, and interactive visualization for unknown sequences
    if len(unknown_sequences) > 0:
        cluster_labels, reduced = cluster_and_visualize_features(
            deep_features, n_clusters=N_CLUSTERS, sequences=unknown_sequences, title_suffix="(Unknown Sequences)"
        )
        
        # Save unknown sequences with cluster assignments to CSV
        unknown_data = []
        for i, (seq, cluster) in enumerate(zip(unknown_sequences, cluster_labels)):
            unknown_data.append([i+1, seq, cluster])
        
        save_results_to_csv(
            "unknown_sequences.csv",
            unknown_data,
            ["Index", "Sequence", "Cluster"]
        )
    else:
        print("No unknown sequences to cluster.")

    print("Pipeline complete. Results saved to CSV files.")
