# import random

# def random_dna_seq(length, alphabet="ATCGN"):
#     return ''.join(random.choice(alphabet) for _ in range(length))

# def write_fasta(filename, num_reads=50, seq_len=150):
#     with open(filename, "w") as f:
#         for i in range(num_reads):
#             header = f">test_seq_{i}"
#             seq = random_dna_seq(seq_len)
#             f.write(f"{header}\n{seq}\n")

# if __name__ == "__main__":
#     write_fasta("novel_sequences.fasta", num_reads=50, seq_len=150)
#     print("novel_sequences.fasta created with 50 random test sequences.")