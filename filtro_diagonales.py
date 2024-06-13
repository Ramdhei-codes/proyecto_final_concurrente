import cv2
import numpy as np
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt

def load_image(file_path):
    """Carga la imagen en escala de grises."""
    image = cv2.imread(file_path, cv2.IMREAD_GRAYSCALE)
    return image

def filter_diagonals(image):
    """Aplica un filtro para resaltar las diagonales."""
    kernel = np.array([[1, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1]], dtype=np.uint8)
    filtered_image = cv2.filter2D(image, -1, kernel)
    return filtered_image

def process_chunk(chunk):
    """Procesa un fragmento de la imagen para identificar las diagonales."""
    return filter_diagonals(chunk)

def divide_chunks(image, num_chunks):
    """Divide la imagen en fragmentos para el procesamiento en paralelo."""
    chunk_size = image.shape[0] // num_chunks
    return [image[i*chunk_size:(i+1)*chunk_size, :] for i in range(num_chunks)]

def merge_chunks(chunks):
    """Combina los fragmentos procesados en una sola imagen."""
    return np.vstack(chunks)

def main(file_path, top_crop, bottom_crop, left_crop, right_crop):
    image = load_image(file_path)
    num_chunks = cpu_count()
    chunks = divide_chunks(image, num_chunks)
    
    with Pool(processes=num_chunks) as pool:
        processed_chunks = pool.map(process_chunk, chunks)
    
    filtered_image = merge_chunks(processed_chunks)
    
    # Recortar la imagen
    filtered_image = filtered_image[top_crop:filtered_image.shape[0]-bottom_crop, left_crop:filtered_image.shape[1]-right_crop]
    
    plt.imshow(filtered_image, cmap='gray')
    plt.title("Diagonales Identificadas en el Dotplot")
    plt.xlabel("Secuencia 2")
    plt.ylabel("Secuencia 1")
    plt.show()

if __name__ == "__main__":
    file_path = "./dotplot_secuencial.png"
    top_crop = 57      # Ajusta este valor según sea necesario
    bottom_crop = 50   # Ajusta este valor según sea necesario
    left_crop = 142     # Ajusta este valor según sea necesario
    right_crop = 130    # Ajusta este valor según sea necesario
    main(file_path, top_crop, bottom_crop, left_crop, right_crop)
