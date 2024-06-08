import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.sparse import coo_matrix
import sys
import multiprocessing as mp

def read_fasta(file_path, max_length=None):
    """Lee una secuencia de un archivo FASTA y devuelve la secuencia."""
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        if max_length:
            sequence = sequence[:max_length]
        return sequence

def compare_subsequences(seq1_array, seq2_array, window_size, start, end, output_queue):
    """Compara sub-secuencias y guarda las coincidencias en una cola."""
    rows, cols = [], []
    len2 = len(seq2_array)
    print(f"Proceso iniciado: rango {start}-{end}")  # Depuración

    for i in range(start, end):
        sub_seq1 = seq1_array[i:i + window_size]
        matches = np.where((seq2_array[:len2 - window_size + 1] == sub_seq1[:, None]).all(axis=0))[0]
        rows.extend([i] * len(matches))
        cols.extend(matches)

    output_queue.put((rows, cols))
    print(f"Proceso terminado: rango {start}-{end}, filas: {len(rows)}, columnas: {len(cols)}")  # Depuración

def generate_dotplot(seq1, seq2, window_size=1, num_processes=4):
    """Genera una matriz dispersa de dotplot para dos secuencias utilizando multiprocessing."""
    len1, len2 = len(seq1), len(seq2)
    seq1_array = np.array(list(seq1))
    seq2_array = np.array(list(seq2))

    chunk_size = len1 // num_processes
    processes = []
    output_queue = mp.Queue()

    try:
        for i in range(num_processes):
            start = i * chunk_size
            end = len1 if i == num_processes - 1 else (i + 1) * chunk_size
            p = mp.Process(target=compare_subsequences, args=(seq1_array, seq2_array, window_size, start, end, output_queue))
            processes.append(p)
            p.start()

        rows, cols = [], []
        for _ in range(num_processes):
            r, c = output_queue.get()
            rows.extend(r)
            cols.extend(c)

        for p in processes:
            p.join()

    except MemoryError:
        print("Error de memoria: No es posible generar el dotplot con las secuencias dadas debido a limitaciones de memoria.")
        sys.exit(1)

    print(f"Total filas: {len(rows)}, Total columnas: {len(cols)}")  # Depuración final
    dotplot = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len1, len2), dtype=int)
    return dotplot.tocsr()

def plot_dotplot(dotplot, output_file):
    """Dibuja y guarda la imagen del dotplot."""
    plt.imshow(dotplot.toarray(), cmap='Greys', interpolation='none')
    plt.savefig(output_file, format='png')
    plt.close()

def main(file1, file2, output_file, max_length, num_processes):
    try:
        seq1 = read_fasta(file1, max_length)
        seq2 = read_fasta(file2, max_length)
        
        print(f"Longitud de la secuencia 1: {len(seq1)}")
        print(f"Longitud de la secuencia 2: {len(seq2)}")
        
        dotplot = generate_dotplot(seq1, seq2, num_processes=num_processes)
        if dotplot is not None:
            plot_dotplot(dotplot, output_file)
    except MemoryError:
        print("Error de memoria durante la ejecución principal.")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generar dotplot de dos secuencias FASTA.")
    parser.add_argument("--file1", required=True, help="Archivo FASTA de la primera secuencia.")
    parser.add_argument("--file2", required=True, help="Archivo FASTA de la segunda secuencia.")
    parser.add_argument("--output", required=True, help="Archivo de salida para la imagen del dotplot.")
    parser.add_argument("--max_length", type=int, default=None, help="Número máximo de caracteres a procesar de cada secuencia.")
    parser.add_argument("--num_processes", type=int, default=4, help="Número de procesos a utilizar.")
    
    args = parser.parse_args()
    
    main(args.file1, args.file2, args.output, args.max_length, args.num_processes)

# Ejecutar el script con los argumentos necesarios
# python multiprocessing-code.py --file1=./dotplot_files/E_coli.fna --file2=./dotplot_files/Salmonella.fna --output=dotplot_multiptocessing.png --max_length=1000