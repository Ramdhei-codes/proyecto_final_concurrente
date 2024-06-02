import numpy as np
import time
from tqdm import tqdm 
from Bio import SeqIO
import matplotlib.pyplot as plt

def merge_sequences_from_fasta(file_path):
    sequences = []  # List to store all sequences
    for record in SeqIO.parse(file_path, "fasta"):
        # `record.seq` gives the sequence
        sequences.append(str(record.seq))
    return "".join(sequences)

file_path_1 = "./dotplot_files/E_coli.fna"
file_path_2 = "./dotplot_files/Salmonella.fna"

merged_sequence_1 = merge_sequences_from_fasta(file_path_1) # estas son las secuencias que se van a utilizar para el dotplot
merged_sequence_2 = merge_sequences_from_fasta(file_path_2)

print(len(merged_sequence_1))
print(len(merged_sequence_2))

begin = time.time()
dotplot = np.empty([len(merged_sequence_1),len(merged_sequence_2)])
print("La matriz de resultado tiene tamaño: ", dotplot.shape)

for i in tqdm(range(dotplot.shape[0])):
  for j in range(dotplot.shape[1]):
    if merged_sequence_1[i] == merged_sequence_2[j]:
      dotplot[i,j] = 1
    else:
      dotplot[i,j] = 0

print(f"\n El código se ejecutó en: {time.time() - begin} segundos")


def draw_dotplot(matrix, fig_name='dotplot.svg'):
  plt.figure(figsize=(5,5))
  plt.imshow(matrix, cmap='Greys',aspect='auto')

  plt.ylabel("Secuencia 1")
  plt.xlabel("Secuencia 2")
  plt.savefig(fig_name)