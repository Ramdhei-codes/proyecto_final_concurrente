import cv2
import numpy as np
from matplotlib import pyplot as plt

# Load the image
image_path = './dotplot_hilos.png'
image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)

# Check if the image was loaded correctly
if image is None:
    raise ValueError("Image not loaded correctly")

# Define diagonal Sobel filters
# Kernel for 45 degree diagonal (top-left to bottom-right)
sobel_45 = np.array([[2, -1, -1],
                     [-1, 2, -1],
                     [-1, -1, 2]])

# Kernel for 135 degree diagonal (top-right to bottom-left)
sobel_135 = np.array([[-1, -1, 2],
                      [-1, 2, -1],
                      [2, -1, -1]])

# Apply the filters using convolution
filtered_45 = cv2.filter2D(image, -1, sobel_45)
filtered_135 = cv2.filter2D(image, -1, sobel_135)

# Combine the results
combined_filtered = cv2.addWeighted(filtered_45, 0.5, filtered_135, 0.5, 0)

# Display the original and filtered images
plt.figure(figsize=(15, 5))
plt.subplot(1, 3, 1)
plt.title("Original Image")
plt.imshow(image, cmap='gray')
plt.subplot(1, 3, 2)
plt.title("Diagonal Edge Detection (45 degrees)")
plt.imshow(filtered_45, cmap='gray')
plt.subplot(1, 3, 3)
plt.title("Diagonal Edge Detection (135 degrees)")
plt.imshow(filtered_135, cmap='gray')
plt.figure(figsize=(5, 5))
plt.title("Combined Diagonal Edge Detection")
plt.imshow(combined_filtered, cmap='gray')
plt.show()
