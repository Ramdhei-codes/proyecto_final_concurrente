import numpy as np
import time
from tqdm import tqdm
from Bio import SeqIO
import matplotlib.pyplot as plt
import threading
from queue import Queue
import argparse

def merge_sequences_from_fasta(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return "".join(sequences)

def compute_dotplot_section(start_i, end_i, start_j, end_j, merged_sequence_1, merged_sequence_2, result_queue):
    local_dotplot = np.zeros((end_i - start_i, end_j - start_j), dtype=np.uint8)
    for i in tqdm(range(end_i - start_i), f'Índices {start_i} - {end_i}, {start_j} - {end_j} {threading.current_thread().name}'):
        for j in range(end_j - start_j):
            if merged_sequence_1[start_i + i] == merged_sequence_2[start_j + j]:
                local_dotplot[i, j] = 1
    result_queue.put((start_i, start_j, local_dotplot))

def worker(merged_sequence_1, merged_sequence_2, task_queue, result_queue):
    while not task_queue.empty():
        start_i, end_i, start_j, end_j = task_queue.get()
        compute_dotplot_section(start_i, end_i, start_j, end_j, merged_sequence_1, merged_sequence_2, result_queue)
        task_queue.task_done()

def draw_dotplot(matrix, fig_name):
    start_time = time.time()  # Tiempo inicial para la generación de la imagen
    plt.imshow(matrix, cmap='Greys', interpolation='none')
    plt.ylabel("Secuencia 1")
    plt.xlabel("Secuencia 2")
    plt.savefig(fig_name)
    plt.close()
    end_time = time.time()  # Tiempo final para la generación de la imagen
    print(f"Tiempo para generar y guardar la imagen: {end_time - start_time:.2f} segundos")

def main(file1, file2, output_file, max_length, num_threads=4):
    start_time = time.time()  # Tiempo inicial para la ejecución del programa

    merged_sequence_1 = merge_sequences_from_fasta(file1)[0:max_length]
    merged_sequence_2 = merge_sequences_from_fasta(file2)[0:max_length]

    print(len(merged_sequence_1))
    print(len(merged_sequence_2))

    block_size = max_length // num_threads

    task_queue = Queue()
    result_queue = Queue()

    for i in range(0, len(merged_sequence_1), block_size):
        for j in range(0, len(merged_sequence_2), block_size):
            end_i = min(i + block_size, len(merged_sequence_1))
            end_j = min(j + block_size, len(merged_sequence_2))
            task_queue.put((i, end_i, j, end_j))

    threads = []
    for _ in range(num_threads):
        thread = threading.Thread(target=worker, args=(merged_sequence_1, merged_sequence_2, task_queue, result_queue))
        threads.append(thread)
        thread.start()

    task_queue.join()

    for thread in threads:
        thread.join()

    calc_end_time = time.time()  # Tiempo final para los cálculos

    dotplot = np.zeros([len(merged_sequence_1), len(merged_sequence_2)], dtype=np.uint8)

    while not result_queue.empty():
        start_i, start_j, local_dotplot = result_queue.get()
        end_i = start_i + local_dotplot.shape[0]
        end_j = start_j + local_dotplot.shape[1]
        dotplot[start_i:end_i, start_j:end_j] = local_dotplot

    print(f"Tiempo de cálculo para generar el dotplot: {calc_end_time - start_time:.2f} segundos")

    draw_dotplot(dotplot, output_file)

    end_time = time.time()  # Tiempo final para la ejecución del programa
    print(f"Tiempo total de ejecución del programa: {end_time - start_time:.2f} segundos")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generar dotplot de dos secuencias utilizando hilos")
    parser.add_argument("--file1", required=True, help="Archivo FASTA de la primera secuencia.")
    parser.add_argument("--file2", required=True, help="Archivo FASTA de la segunda secuencia.")
    parser.add_argument("--output", required=True, help="Archivo de salida para la imagen del dotplot.")
    parser.add_argument("--max_length", type=int, default=1000, help="Número máximo de caracteres a procesar de cada secuencia.")
    parser.add_argument("--num_threads", type=int, default=4, help="Número de hilos para procesar la secuencia (4 por defecto).")

    args = parser.parse_args()

    main(args.file1, args.file2, args.output, args.max_length, args.num_threads)

# ejecutar con:  
# python .\dotplot_hilos.py --file1 .\dotplot_files\Salmonella.fna --file2 .\dotplot_files\E_coli.fna --max_length 10000 --output .\dotplot_hilos.png
