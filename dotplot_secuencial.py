import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.sparse import coo_matrix

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

    seq1_array = np.array(list(seq1))
    seq2_array = np.array(list(seq2))

    for i in tqdm(range(len1 - window_size + 1), desc="Generando dotplot"):
        sub_seq1 = seq1_array[i:i + window_size]
        matches = np.where((seq2_array[:len2 - window_size + 1] == sub_seq1[:, None]).all(axis=0))[0]
        rows.extend([i] * len(matches))
        cols.extend(matches)
    
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