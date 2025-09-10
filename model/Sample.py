import os
import tarfile
import random
import shutil
import numpy as np

def set_random_seeds(seed=42):
    random.seed(seed)
    np.random.seed(seed)

def generate_sequence_with_gc_content(length, gc_percent):
    """Generates a random DNA sequence with a specified GC content."""
    gc_count = int(round(length * gc_percent))
    at_count = length - gc_count
    sequence_list = (
        ['A'] * (at_count // 2) + ['T'] * (at_count - at_count // 2) +
        ['G'] * (gc_count // 2) + ['C'] * (gc_count - gc_count // 2)
    )
    random.shuffle(sequence_list)
    return ''.join(sequence_list)

def add_mutations(sequence, mutation_rate=0.01):
    """Adds realistic point mutations (substitutions) to a sequence."""
    dna_bases = "ATCG"
    mutated_sequence = list(sequence)
    for i in range(len(mutated_sequence)):
        if random.random() < mutation_rate:
            original_base = mutated_sequence[i]
            possible_bases = dna_bases.replace(original_base, '')
            mutated_sequence[i] = random.choice(possible_bases)
    return "".join(mutated_sequence)

def generate_unique_sequences(generator_func, n, *args, **kwargs):
    """Ensure each sequence generated is unique for a class."""
    seen = set()
    unique = []
    attempts = 0
    while len(unique) < n and attempts < n * 10:
        seq = generator_func(*args, **kwargs)
        if seq not in seen:
            seen.add(seq)
            unique.append(seq)
        attempts += 1
    if len(unique) < n:
        print(f"Warning: Only generated {len(unique)} unique sequences out of {n} requested.")
    return unique

def generate_complex_data(output_dir, num_classes, reads_per_file, seq_len, mutation_rate=0.02, random_seed=42):
    set_random_seeds(random_seed)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    print(f"Generating complex biological data with {num_classes} classes in '{output_dir}'...")

    # Define regions with different lengths
    promoter_len, terminator_len = 50, 50  # Increased from 15 to 50
    gene_body_len = seq_len - promoter_len - terminator_len
    
    # Expanded motif sets with more distinctive patterns
    promoter_motifs = [
        "TATAAT", "TTGACA", "GGGCGG", "CCGCCC", "ATGCAT", 
        "GCCGCC", "CGCGCG", "TATATA", "GCGCGC", "ATATAT"
    ]
    
    terminator_motifs = [
        "AATAAA", "GCGCGC", "TTTTTT", "CCCCCC", "GGGGGG",
        "ATATAT", "TATATA", "CGCGCG", "GCCGCC", "CATCAT"
    ]
    
    # Class-specific motifs for additional distinctiveness
    class_specific_motifs = [
        "AGCTAGCT", "TCGATCGA", "GCTAGCTA", "CTAGCTAG",
        "ATCGATCG", "TAGCTAGC", "GATCGATC", "CGATCGAT",
        "AATTCCGG", "GGCCAATT", "TTAAAGCT", "AGCTTTAAG"
    ]

    # Assign a wider range of GC content percentages to each class
    gc_percentages = np.linspace(0.2, 0.8, num_classes)  # Wider range: 20% to 80%

    for i in range(num_classes):
        label = f"Class_{i}"
        gc_content = gc_percentages[i]
        file_path = os.path.join(output_dir, f"{label}_reads.fasta")
        
        # Assign a class-specific motif
        class_motif = class_specific_motifs[i % len(class_specific_motifs)]
        
        with open(file_path, "w") as f:
            # Generate unique gene bodies for this class
            gene_bodies = generate_unique_sequences(
                generate_sequence_with_gc_content, reads_per_file, gene_body_len, gc_content
            )
            
            for j in range(reads_per_file):
                # Promoter with class-specific motif
                promoter_base = generate_sequence_with_gc_content(promoter_len, 0.5)
                promoter_motif = promoter_motifs[i % len(promoter_motifs)]
                promoter_start = random.randint(0, promoter_len - len(promoter_motif))
                promoter_final = (
                    promoter_base[:promoter_start] + promoter_motif +
                    promoter_base[promoter_start + len(promoter_motif):]
                )
                
                # Add class-specific motif to promoter
                class_motif_pos = random.randint(0, promoter_len - len(class_motif))
                if class_motif_pos != promoter_start:  # Avoid overlapping with the main motif
                    promoter_final = (
                        promoter_final[:class_motif_pos] + class_motif +
                        promoter_final[class_motif_pos + len(class_motif):]
                    )
                
                # Gene body (unique per read)
                gene_body = gene_bodies[j] if j < len(gene_bodies) else generate_sequence_with_gc_content(gene_body_len, gc_content)
                
                # Add some class-specific patterns to the gene body
                if j % 5 == 0:  # Add pattern to some sequences
                    pattern_pos = random.randint(0, gene_body_len - len(class_motif))
                    gene_body = (
                        gene_body[:pattern_pos] + class_motif +
                        gene_body[pattern_pos + len(class_motif):]
                    )
                
                # Terminator
                terminator_base = generate_sequence_with_gc_content(terminator_len, 0.5)
                terminator_motif = terminator_motifs[i % len(terminator_motifs)]
                terminator_start = random.randint(0, terminator_len - len(terminator_motif))
                terminator_final = (
                    terminator_base[:terminator_start] + terminator_motif +
                    terminator_base[terminator_start + len(terminator_motif):]
                )
                
                # Assemble and mutate
                sequence = promoter_final + gene_body + terminator_final
                final_sequence = add_mutations(sequence, mutation_rate=mutation_rate)
                f.write(f">{label}_read_{j}\n{final_sequence}\n")

    tar_path = os.path.join(output_dir, "sample_reads.tar")
    with tarfile.open(tar_path, "w") as tar:
        for file in os.listdir(output_dir):
            if file.endswith(".fasta"):
                tar.add(os.path.join(output_dir, file), arcname=file)
    print(f"Successfully created '{tar_path}' containing sequences for {num_classes} classes.")

# Run the data generation script for supervised and unsupervised data
if __name__ == "__main__":
    # Generate data with fewer classes but more distinctive features
    generate_complex_data(
        output_dir="sample_data", num_classes=10, reads_per_file=1000,  # Fewer classes, more samples
        seq_len=500, mutation_rate=0.02, random_seed=42  # Longer sequences
    )
    generate_complex_data(
        output_dir="novel_data", num_classes=10, reads_per_file=500,    # Fewer classes
        seq_len=500, mutation_rate=0.02, random_seed=1337  # Longer sequences
    )