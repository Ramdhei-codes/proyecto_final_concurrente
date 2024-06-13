import matplotlib.pyplot as plt
import numpy as np

# Datos para gráficas de Multiprocessing
processes = [2, 4, 8, 16]

# Tiempos de ejecución para 10,000 caracteres con distintos números de procesos
times_total = [5.9, 3.9, 3.19, 3.71]
times_calculation = [5.22, 3.15, 2.56, 3.06]
times_image = [0.61, 0.69, 0.57, 0.59]


x = np.arange(len(processes))  # las etiquetas de los ejes x para los procesos
width = 0.25  # ancho de las barras

fig, ax = plt.subplots()
rects1 = ax.bar(x - width, times_total, width, label='Tiempo Total')
rects2 = ax.bar(x, times_calculation, width, label='Tiempo de Cálculo')
rects3 = ax.bar(x + width, times_image, width, label='Tiempo de Guardar Imagen')

# Añadir algunas etiquetas de texto para títulos y etiquetas de ejes
ax.set_xlabel('Número de Procesos')
ax.set_ylabel('Tiempo (s)')
ax.set_title('Tiempos de Ejecución para Multiprocessing con 10,000 Caracteres')
ax.set_xticks(x)
ax.set_xticklabels([f'{p} procesos' for p in processes])
ax.legend()

fig.tight_layout()

# Guardar gráfico
plt.savefig('./assets/imgs/tiempo_ejecucion_multiprocessing_10000.png')
plt.close()

# Datos para gráficas de Multiprocessing para 25,000 caracteres
processes = [2, 4, 8, 16]

# Tiempos de ejecución con 25,000 caracteres
times_total_25000 = [32.59, 18.63, 12.88, 11.89]
times_calculation_25000 = [29.20, 15.95, 10.27, 9.25]
times_image_25000 = [3.33, 2.62, 2.55, 2.58]

x = np.arange(len(processes))  # las etiquetas de los ejes x para los procesos
width = 0.25  # ancho de las barras

fig, ax = plt.subplots()
rects1 = ax.bar(x - width, times_total_25000, width, label='Tiempo Total')
rects2 = ax.bar(x, times_calculation_25000, width, label='Tiempo de Cálculo')
rects3 = ax.bar(x + width, times_image_25000, width, label='Tiempo de Guardar Imagen')

# Añadir algunas etiquetas de texto para títulos y etiquetas de ejes
ax.set_xlabel('Número de Procesos')
ax.set_ylabel('Tiempo (s)')
ax.set_title('Tiempos de Ejecución para Multiprocessing con 25,000 Caracteres')
ax.set_xticks(x)
ax.set_xticklabels([f'{p} procesos' for p in processes])
ax.legend()

fig.tight_layout()

# Guardar gráfico
plt.savefig('./assets/imgs/tiempo_ejecucion_multiprocessing_25000.png')
plt.close()

# Datos para gráficas de Multiprocessing para 50,000 caracteres
processes = [2, 4, 8, 16]

# Tiempos de ejecución con 50,000 caracteres
times_total_50000 = [159.79, 105.11, 62.06, 72.16]
times_calculation_50000 = [116.75, 62.35, 38.94, 33.86]
times_image_50000 = [42.98, 42.67, 23.05, 38.23]

x = np.arange(len(processes))  # las etiquetas de los ejes x para los procesos
width = 0.25  # ancho de las barras

fig, ax = plt.subplots()
rects1 = ax.bar(x - width, times_total_50000, width, label='Tiempo Total')
rects2 = ax.bar(x, times_calculation_50000, width, label='Tiempo de Cálculo')
rects3 = ax.bar(x + width, times_image_50000, width, label='Tiempo de Guardar Imagen')

# Añadir algunas etiquetas de texto para títulos y etiquetas de ejes
ax.set_xlabel('Número de Procesos')
ax.set_ylabel('Tiempo (s)')
ax.set_title('Tiempos de Ejecución para Multiprocessing con 50,000 Caracteres')
ax.set_xticks(x)
ax.set_xticklabels([f'{p} procesos' for p in processes])
ax.legend()

fig.tight_layout()

# Guardar gráfico
plt.savefig('./assets/imgs/tiempo_ejecucion_multiprocessing_50000.png')
plt.close()



# Datos de tiempo total secuencial obtenidos anteriormente
sizes = ['10,000', '25,000', '50,000']
sequential_times = [5.96, 36.84, 666.69]  # Tiempos secuenciales obtenidos para 10,000 y 25,000 caracteres

# Mejores tiempos usando multiprocessing
best_multiprocessing_times_10000 = min([5.9, 3.9, 3.19, 3.71])
best_multiprocessing_times_25000 = min([32.59, 18.63, 12.88, 11.89])
best_multiprocessing_times_50000 = min([159.79, 105.11, 62.06, 72.16])

best_multiprocessing_times = [best_multiprocessing_times_10000, best_multiprocessing_times_25000, best_multiprocessing_times_50000]

x = np.arange(len(sizes))

# Configuración de la gráfica
plt.figure(figsize=(10, 6))
plt.bar(x - 0.2, sequential_times, width=0.4, label='Secuencial', color='b')
plt.bar(x + 0.2, best_multiprocessing_times, width=0.4, label='Mejor tiempo con multiprocessing', color='r')

plt.xlabel('Longitud de la secuencia')
plt.ylabel('Tiempo (s)')
plt.title('Comparación de Tiempos Secuencial vs Multiprocessing')
plt.xticks(x, sizes)
plt.legend()

# Guardar el gráfico
plt.savefig("./assets/imgs/comparacion_secuencial_multiprocessing.png")
plt.close()
