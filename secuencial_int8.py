import numpy as np
import time
from Bio import SeqIO
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from tqdm import tqdm

def merge_sequences_from_fasta(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return "".join(sequences)

def sequence_to_int8(sequence):
    base_to_int = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    return np.array([base_to_int[base] for base in sequence], dtype=np.int8)

def dotplot(seq1, seq2, window_size):
    matrix = lil_matrix((len(seq1), len(seq2)), dtype=np.int8)
    
    for i in tqdm(range(len(seq1) - window_size + 1), desc="Dotplot Progress"):
        for j in range(len(seq2) - window_size + 1):
            if np.array_equal(seq1[i:i+window_size], seq2[j:j+window_size]):
                matrix[i, j] = 1
    
    return matrix.tocsc()

file_path_1 = "./dotplot_files/E_coli.fna"
file_path_2 = "./dotplot_files/Salmonella.fna"

start_time = time.time()

merged_sequence_1 = merge_sequences_from_fasta(file_path_1)
merged_sequence_2 = merge_sequences_from_fasta(file_path_2)

seq1_int8 = sequence_to_int8(merged_sequence_1)
seq2_int8 = sequence_to_int8(merged_sequence_2)

data_load_time = time.time() - start_time

print(len(seq1_int8))
print(len(seq2_int8))

window_size = 10

start_time = time.time()

dotplot_matrix = dotplot(seq1_int8, seq2_int8, window_size)

computation_time = time.time() - start_time

print(f"Data load time: {data_load_time} seconds")
print(f"Computation time: {computation_time} seconds")

start_time = time.time()

plt.figure(figsize=(5, 5))
plt.imshow(dotplot_matrix.toarray(), cmap='Greys', aspect='auto')
plt.ylabel("Sequence 1")
plt.xlabel("Sequence 2")
plt.savefig('dotplot.png')

image_generation_time = time.time() - start_time

print(f"Image generation time: {image_generation_time} seconds")
