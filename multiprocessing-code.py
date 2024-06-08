import numpy as np
import time
from Bio import SeqIO
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, csc_matrix, vstack
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import argparse

def merge_sequences_from_fasta(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return "".join(sequences)

def sequence_to_int8(sequence):
    base_to_int = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    return np.array([base_to_int[base] for base in sequence], dtype=np.int8)

def dotplot_chunk_optimized(args):
    seq1, seq2, i_start, i_end, window_size = args
    chunk_matrix = lil_matrix((i_end - i_start, len(seq2) - window_size + 1), dtype=np.int8)

    for i in tqdm(range(i_start, i_end), desc=f"Process {i_start}-{i_end} progress"):
        windows_seq1 = np.lib.stride_tricks.sliding_window_view(seq1[i:i + window_size], window_shape=window_size)
        windows_seq2 = np.lib.stride_tricks.sliding_window_view(seq2, window_shape=window_size)
        match_matrix = (windows_seq1[:, None, :] == windows_seq2).all(axis=-1)
        chunk_matrix[i - i_start, :] = match_matrix.any(axis=1).astype(np.int8)

    return chunk_matrix.tocsr()

def parallel_dotplot_optimized(seq1, seq2, window_size, num_processes):
    pool = Pool(processes=num_processes)
    chunk_size = len(seq1) // num_processes
    chunks = [(seq1, seq2, i, min(i + chunk_size, len(seq1)), window_size)
              for i in range(0, len(seq1), chunk_size)]

    chunk_results = list(tqdm(pool.imap(dotplot_chunk_optimized, chunks), total=len(chunks)))
    pool.close()
    pool.join()

    matrix = vstack(chunk_results)
    return matrix

def main(file1, file2, output_file, max_length, num_processes=cpu_count()):
    start_time = time.time()

    merged_sequence_1 = merge_sequences_from_fasta(file1)[0:max_length]
    merged_sequence_2 = merge_sequences_from_fasta(file2)[0:max_length]

    seq1_int8 = sequence_to_int8(merged_sequence_1)
    seq2_int8 = sequence_to_int8(merged_sequence_2)

    data_load_time = time.time() - start_time

    window_size = 10 

    start_time = time.time()

    dotplot_matrix = parallel_dotplot_optimized(seq1_int8, seq2_int8, window_size, num_processes)

    computation_time = time.time() - start_time

    print(f"Data load time: {data_load_time} seconds")
    print(f"Computation time: {computation_time} seconds")

    start_time = time.time()

    plt.figure(figsize=(10, 10))
    plt.imshow(dotplot_matrix.toarray(), cmap='Greys', aspect='auto')
    plt.ylabel("Sequence 1")
    plt.xlabel("Sequence 2")
    plt.savefig(output_file)

    image_generation_time = time.time() - start_time

    print(f"Image generation time: {image_generation_time} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generar dotplot de dos secuencias utilizando hilos")
    parser.add_argument("--file1", required=True, help="Archivo FASTA de la primera secuencia.")
    parser.add_argument("--file2", required=True, help="Archivo FASTA de la segunda secuencia.")
    parser.add_argument("--output", required=True, help="Archivo de salida para la imagen del dotplot.")
    parser.add_argument("--max_length", type=int, default=1000, help="Número máximo de caracteres a procesar de cada secuencia.")
    parser.add_argument("--num_processes", type=int, default=cpu_count(), help="Número de procesos para procesar la secuencia (Los hilos de la CPU por defecto).")

    args = parser.parse_args()

    main(args.file1, args.file2, args.output, args.max_length, args.num_processes)
