import matplotlib.pyplot as plt
import numpy as np

# # Datos
# sizes = ['10,000', '25,000', '50,000']
# threads = [2, 4, 8, 16]

# times_10000 = {
#     2: [1.37, 0.36, 0.97],
#     4: [1.16, 0.25, 0.85],
#     8: [1.19, 0.28, 0.85],
#     16: [1.29, 0.35, 0.86]
# }

# times_25000 = {
#     2: [5.75, 1.58, 3.42],
#     4: [4.14, 1.44, 2.31],
#     8: [4.32, 1.06, 2.58],
#     16: [4.01, 1.21, 2.38]
# }

# times_50000 = {
#     2: [39.01, 6.31, 31.76],
#     4: [43.42, 5.80, 36.63],
#     8: [33.70, 5.15, 27.53],
#     16: [41.63, 3.94, 36.65]
# }

# # Función para generar y guardar diagramas de barras
# def plot_data(threads, data, title, ylabel, size_label, file_name):
#     fig, ax = plt.subplots()
#     bar_width = 0.2
#     index = np.arange(len(data[threads[0]]))  # Asegurar el tamaño correcto del índice basado en los datos
    
#     for i, thread in enumerate(threads):
#         times = [data[thread][j] for j in range(len(data[thread]))]
#         ax.bar(index + bar_width*i, times, bar_width, label=f'{thread} hilos')
    
#     ax.set_xlabel('Metricas')
#     ax.set_ylabel(ylabel)
#     ax.set_title(f'Métricas de tiempo para el tamaño de la secuencia {size_label}')
#     ax.set_xticks(index + bar_width*(len(threads)-1)/2)
#     ax.set_xticklabels(['Tiempo total', 'Tiempo de cálculo', 'Tiempo de imágen'])
#     ax.legend()

#     plt.savefig(f"./assets/imgs/{file_name}.png")  # Guardar el gráfico como un archivo PNG
#     plt.close()

# # Ejemplo de llamada a la función para graficar y guardar los datos
# plot_data(threads, times_10000, '10,000 caracteres', 'Time (s)', '10,000', 'tiempo_ejecucion_hilos_10000')
# plot_data(threads, times_25000, '25,000 caracteres', 'Time (s)', '25,000', 'tiempo_ejecucion_hilos_25000')
# plot_data(threads, times_50000, '50,000 caracteres', 'Time (s)', '50,000', 'tiempo_ejecucion_hilos_50000')


import matplotlib.pyplot as plt
import numpy as np

sizes = ['10,000', '25,000', '50,000']
sequential_times = [5.96, 36.84, 39.01]  # Tiempos secuenciales 
best_thread_times = [1.16, 4.01, 33.70]  # Mejores tiempos multihilo

x = np.arange(len(sizes))

plt.figure(figsize=(10, 6))
plt.bar(x - 0.2, sequential_times, width=0.4, label='Secuencial', color='b')
plt.bar(x + 0.2, best_thread_times, width=0.4, label='Mejor tiempo con multihilos', color='r')

plt.xlabel('Longitud de la secuencia')
plt.ylabel('Tiempo (s)')
plt.title('Secuencial vs Multihilos')
plt.xticks(x, sizes)
plt.legend()

plt.savefig('/home/erika/20182.pcomp/proyecto_final_concurrente/assets/imgs/secuencial_vs_multihilos.png')  # Guarda la imagen en lugar de mostrarla
plt.close()

