import numpy as np
import time
from tqdm import tqdm
from Bio import SeqIO
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix

def merge_sequences_from_fasta(file_path):
    sequences = []  # List to store all sequences
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return "".join(sequences)

def sequence_to_int8(sequence):
    # Convert the sequence to np.int8
    return np.array([ord(char) for char in sequence], dtype=np.int8)

file_path_1 = "./dotplot_files/E_coli.fna"
file_path_2 = "./dotplot_files/Salmonella.fna"

merged_sequence_1 = merge_sequences_from_fasta(file_path_1) # estas son las secuencias que se van a utilizar para el dotplot
merged_sequence_2 = merge_sequences_from_fasta(file_path_2)

seq1_int8 = sequence_to_int8(merged_sequence_1)
seq2_int8 = sequence_to_int8(merged_sequence_2)

print(len(seq1_int8))
print(len(seq2_int8))

begin = time.time()
dotplot = lil_matrix((len(seq1_int8), len(seq2_int8)), dtype=np.int8)
print("La matriz de resultado tiene tamaño: ", dotplot.shape)

for i in tqdm(range(dotplot.shape[0])):
    for j in range(dotplot.shape[1]):
        if seq1_int8[i] == seq2_int8[j]:
            dotplot[i, j] = 1

print(f"\n El código se ejecutó en: {time.time() - begin} segundos")

def draw_dotplot(matrix, fig_name='dotplot.svg'):
    plt.figure(figsize=(5, 5))
    plt.imshow(matrix.toarray(), cmap='Greys', aspect='auto')
    plt.ylabel("Secuencia 1")
    plt.xlabel("Secuencia 2")
    plt.savefig(fig_name)

draw_dotplot(dotplot)
