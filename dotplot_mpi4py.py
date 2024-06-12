import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.sparse import coo_matrix
from mpi4py import MPI
import time

def read_fasta(file_path, max_length=None):
    """Lee una secuencia de un archivo FASTA y devuelve la secuencia."""
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        if max_length:
            sequence = sequence[:max_length]
        return sequence

def generate_dotplot_parallel(seq1, seq2, window_size=1, comm=None):
    """Genera una matriz dispersa de dotplot para dos secuencias utilizando MPI."""
    len1, len2 = len(seq1), len(seq2)
    rank = comm.Get_rank()
    size = comm.Get_size()

    rows, cols = [], []

    seq1_array = np.array(list(seq1))
    seq2_array = np.array(list(seq2))

    for i in tqdm(range(rank, len1 - window_size + 1, size), desc=f"Proceso {rank}"):
        sub_seq1 = seq1_array[i:i + window_size]
        matches = np.where((seq2_array[:len2 - window_size + 1] == sub_seq1[:, None]).all(axis=0))[0]
        rows.extend([i] * len(matches))
        cols.extend(matches)

    local_dotplot = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len1, len2), dtype=int)

    # Gather all local dotplots into one at root process
    gathered_dotplots = comm.gather(local_dotplot, root=0)

    if rank == 0:
        # Combine all gathered dotplots into one
        total_rows, total_cols = [], []
        for dp in gathered_dotplots:
            total_rows.extend(dp.row)
            total_cols.extend(dp.col)
        dotplot = coo_matrix((np.ones(len(total_rows)), (total_rows, total_cols)), shape=(len1, len2), dtype=int)
        return dotplot.tocsr()
    else:
        return None

def plot_dotplot(dotplot, output_file):
    """Dibuja y guarda la imagen del dotplot."""
    start_time = time.time()
    plt.imshow(dotplot.toarray(), cmap='Greys', interpolation='none')
    plt.savefig(output_file, format='png')
    plt.close()
    end_time = time.time()
    print(f"Tiempo para generar y guardar la imagen: {end_time - start_time:.2f} segundos")

def main(file1, file2, output_file, max_length, num_processes):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    start_time = time.time()
    
    if rank == 0:
        seq1 = read_fasta(file1, max_length)
        seq2 = read_fasta(file2, max_length)
        print(f"Longitud de la secuencia 1: {len(seq1)}")
        print(f"Longitud de la secuencia 2: {len(seq2)}")
    else:
        seq1 = None
        seq2 = None

    seq1 = comm.bcast(seq1, root=0)
    seq2 = comm.bcast(seq2, root=0)

    calc_start_time = time.time()
    dotplot = generate_dotplot_parallel(seq1, seq2, comm=comm)
    calc_end_time = time.time()

    if rank == 0:
        if dotplot is not None:
            print(f"Tiempo de cálculo para generar el dotplot: {calc_end_time - calc_start_time:.2f} segundos")
            plot_dotplot(dotplot, output_file)
        else:
            print("No se pudo generar el dotplot debido a un error de memoria.")

    end_time = time.time()
    if rank == 0:
        print(f"Tiempo total de ejecución del programa: {end_time - start_time:.2f} segundos")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generar dotplot de dos secuencias FASTA utilizando MPI.")
    parser.add_argument("--file1", required=True, help="Archivo FASTA de la primera secuencia.")
    parser.add_argument("--file2", required=True, help="Archivo FASTA de la segunda secuencia.")
    parser.add_argument("--output", required=True, help="Archivo de salida para la imagen del dotplot.")
    parser.add_argument("--max_length", type=int, default=None, help="Número máximo de caracteres a procesar de cada secuencia.")
    parser.add_argument("--num_processes", type=int, default=4, help="Número de procesos a utilizar.")
    
    args = parser.parse_args()

    main(args.file1, args.file2, args.output, args.max_length, args.num_processes)



# Ejecutar con varios procesos
# mpiexec -n 4 python dotplot_mpi4py.py --file1=./dotplot_files/E_coli.fna --file2=./dotplot_files/Salmonella.fna --output=dotplot_mpi.png --max_length=10000
