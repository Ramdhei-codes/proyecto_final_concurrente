import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import time
import argparse
from Bio import SeqIO

def read_fasta(file_path, max_length=None):
    """Read a sequence from a FASTA file and return the sequence."""
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        if max_length:
            sequence = sequence[:max_length]
        return sequence

def worker(args):
    i, Secuencia1, Secuencia2, len_secuencia2 = args
    result = np.zeros(len_secuencia2, dtype=bool)
    
    for j in range(len_secuencia2):
        result[j] = Secuencia1[i] == Secuencia2[j]
        
    return result

def parallel_dotplot(Secuencia1, Secuencia2, threads=mp.cpu_count()):
    len_secuencia1 = len(Secuencia1)
    len_secuencia2 = len(Secuencia2)
    
    # Using shared memory for sequences
    with mp.Pool(processes=threads) as pool:
        results = pool.map(worker, [(i, Secuencia1, Secuencia2, len_secuencia2) for i in range(len_secuencia1)])
    
    return np.array(results)

def draw_dotplot(matrix, fig_name='dotplot_multiprocessing.png'):
    plt.figure(figsize=(5, 5))
    plt.imshow(matrix, cmap='Greys', interpolation='none')
    plt.ylabel("Secuencia 1")
    plt.xlabel("Secuencia 2")
    plt.savefig(fig_name)

def main(file1, file2, output_file, max_length, num_processes=mp.cpu_count()):
    merged_sequence_1 = read_fasta(file1, max_length)
    merged_sequence_2 = read_fasta(file2, max_length)

    print(len(merged_sequence_1))
    print(len(merged_sequence_2))

    begin = time.time()
    dotplot = parallel_dotplot(merged_sequence_1, merged_sequence_2, num_processes)
    print("La matriz de resultado tiene tamaño: ", dotplot.shape)
    print(f"\n El código se ejecutó en: {time.time() - begin} segundos")

    draw_dotplot(dotplot, fig_name=output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generar dotplot de dos secuencias utilizando multiprocessing")
    parser.add_argument("--file1", required=True, help="Archivo FASTA de la primera secuencia.")
    parser.add_argument("--file2", required=True, help="Archivo FASTA de la segunda secuencia.")
    parser.add_argument("--output", required=True, help="Archivo de salida para la imagen del dotplot.")
    parser.add_argument("--max_length", type=int, default=1000, help="Número máximo de caracteres a procesar de cada secuencia.")
    parser.add_argument("--num_processes", type=int, default=4, help="Número de procesos para procesar la secuencia (4 por defecto).")

    args = parser.parse_args()

    main(args.file1, args.file2, args.output, args.max_length, args.num_processes)
#python prueba_multiprocessing.py --file1=./dotplot_files/E_coli.fna --file2=./dotplot_files/Salmonella.fna --output=dotplot_multiprocessing.png --max_length=50000 --num_processes 16