import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.sparse import coo_matrix, save_npz

def read_fasta(file_path, max_length=None):
    """Lee una secuencia de un archivo FASTA y devuelve la secuencia."""
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        if max_length:
            sequence = sequence[:max_length]
        return sequence

def generate_dotplot(seq1, seq2, window_size=1):
    """Genera una matriz dispersa de dotplot para dos secuencias."""
    len1, len2 = len(seq1), len(seq2)
    rows, cols = [], []

    for i in tqdm(range(len1 - window_size + 1), desc="Generando dotplot"):
        for j in range(len2 - window_size + 1):
            if seq1[i:i+window_size] == seq2[j:j+window_size]:
                rows.append(i)
                cols.append(j)
    
    dotplot = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len1, len2), dtype=int)
    return dotplot.tocsr()

def plot_dotplot(dotplot, output_file):
    """Dibuja y guarda la imagen del dotplot."""
    plt.imshow(dotplot.toarray(), cmap='Greys', interpolation='none')
    plt.savefig(output_file, format='png')
    plt.close()

def main(file1, file2, output_file, max_length):
    seq1 = read_fasta(file1, max_length)
    seq2 = read_fasta(file2, max_length)
    
    print(f"Longitud de la secuencia 1: {len(seq1)}")
    print(f"Longitud de la secuencia 2: {len(seq2)}")
    
    dotplot = generate_dotplot(seq1, seq2)
    plot_dotplot(dotplot, output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generar dotplot de dos secuencias FASTA.")
    parser.add_argument("--file1", required=True, help="Archivo FASTA de la primera secuencia.")
    parser.add_argument("--file2", required=True, help="Archivo FASTA de la segunda secuencia.")
    parser.add_argument("--output", required=True, help="Archivo de salida para la imagen del dotplot.")
    parser.add_argument("--max_length", type=int, default=None, help="Número máximo de caracteres a procesar de cada secuencia.")
    
    args = parser.parse_args()
    
    main(args.file1, args.file2, args.output, args.max_length)
