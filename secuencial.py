import numpy as np
from tqdm import tqdm
from Bio import SeqIO
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix

def merge_sequences_from_fasta(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return "".join(sequences)

def sequence_to_int8(sequence):
    return np.array([ord(char) for char in sequence], dtype=np.int8)

def generate_dotplot(seq1, seq2):
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    
    dotplot_matrix = lil_matrix((len_seq1, len_seq2), dtype=np.int8)
    
    for i in tqdm(range(len_seq1)):
        for j in range(len_seq2):
            if seq1[i] == seq2[j]:
                dotplot_matrix[i, j] = 1
    
    return dotplot_matrix

def plot_dotplot(dotplot_matrix, seq1, seq2):
    plt.figure(figsize=(10, 10))
    plt.spy(dotplot_matrix, markersize=1)
    plt.title('Dotplot')
    plt.xlabel('Sequence 1')
    plt.ylabel('Sequence 2')
    plt.show()

file_path_1 = "./dotplot_files/E_coli.fna"
file_path_2 = "./dotplot_files/Salmonella.fna"

merged_sequence_1 = merge_sequences_from_fasta(file_path_1)
merged_sequence_2 = merge_sequences_from_fasta(file_path_2)

seq1_int8 = sequence_to_int8(merged_sequence_1)
seq2_int8 = sequence_to_int8(merged_sequence_2)

print("Length of Sequence 1:", len(seq1_int8))
print("Length of Sequence 2:", len(seq2_int8))

dotplot_matrix = generate_dotplot(seq1_int8, seq2_int8)

plot_dotplot(dotplot_matrix, seq1_int8, seq2_int8)
