import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.sparse import coo_matrix, vstack
import sys
import multiprocessing as mp
import time

def read_fasta(file_path, max_length=None):
    """Lee una secuencia de un archivo FASTA y devuelve la secuencia."""
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        if max_length:
            sequence = sequence[:max_length]
        return sequence

def compare_subsequences(seq1_array, seq2_array, window_size, start, end):
    """Compara sub-secuencias y devuelve las coincidencias."""
    rows, cols = [], []
    len2 = len(seq2_array)

    for i in range(start, end, window_size):
        sub_seq1 = seq1_array[i:i + window_size]
        matches = np.where((seq2_array[:len2 - window_size + 1] == sub_seq1[:, None]).all(axis=0))[0]
        rows.extend([i] * len(matches))
        cols.extend(matches)

    return rows, cols

def generate_dotplot(seq1, seq2, window_size=1, num_processes=4):
    """Genera una matriz dispersa de dotplot para dos secuencias utilizando multiprocessing."""
    len1, len2 = len(seq1), len(seq2)
    seq1_array = np.array(list(seq1))
    seq2_array = np.array(list(seq2))

    chunk_size = max(1, len1 // (num_processes * 10))  # Tamaño de chunk más pequeño para mejor balance de carga
    with mp.Pool(processes=num_processes) as pool:
        tasks = [(seq1_array, seq2_array, window_size, start, min(len1, start + chunk_size))
                 for start in range(0, len1, chunk_size)]

        rows, cols = [], []
        for result in tqdm(pool.imap_unordered(compare_subsequences, tasks), total=len(tasks)):
            r, c = result
            rows.extend(r)
            cols.extend(c)

    print(f"Total filas: {len(rows)}, Total columnas: {len(cols)}")  # Depuración final
    dotplot = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len1, len2), dtype=int)
    return dotplot.tocsr()

def plot_dotplot(dotplot, output_file):
    """Dibuja y guarda la imagen del dotplot."""
    start_time = time.time()  # Tiempo inicial para la generación de la imagen
    plt.imshow(dotplot.toarray(), cmap='Greys', interpolation='none')
    plt.savefig(output_file, format='png')
    plt.close()
    end_time = time.time()  # Tiempo final para la generación de la imagen
    print(f"Tiempo para generar y guardar la imagen: {end_time - start_time:.2f} segundos")

def main(file1, file2, output_file, max_length, num_processes):
    start_time = time.time()  # Tiempo inicial para la ejecución del programa
    try:
        seq1 = read_fasta(file1, max_length)
        seq2 = read_fasta(file2, max_length)
        
        print(f"Longitud de la secuencia 1: {len(seq1)}")
        print(f"Longitud de la secuencia 2: {len(seq2)}")
        
        calc_start_time = time.time()  # Tiempo inicial para los cálculos
        dotplot = generate_dotplot(seq1, seq2, num_processes=num_processes)
        calc_end_time = time.time()  # Tiempo final para los cálculos
        
        if dotplot is not None:
            print(f"Tiempo de cálculo para generar el dotplot: {calc_end_time - calc_start_time:.2f} segundos")
            plot_dotplot(dotplot, output_file)
        else:
            print("No se pudo generar el dotplot debido a un error de memoria.")
    except MemoryError:
        print("Error de memoria durante la ejecución principal.")
        sys.exit(1)

    end_time = time.time()  # Tiempo final para la ejecución del programa
    print(f"Tiempo total de ejecución del programa: {end_time - start_time:.2f} segundos")

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
# python multiprocessing-code.py --file1=./dotplot_files/E_coli.fna --file2=./dotplot_files/Salmonella.fna --output=dotplot_multiprocessing.png --max_length=1000
