import matplotlib.pyplot as plt
import numpy as np

# Datos para gráficas de MPI
processes = [2, 4, 8, 16]

# Datos para 10,000 caracteres
times_total_10000 = [4.67, 3.80, 4.97, 10.01]
times_calculation_10000 = [8.1, 2.04, 1.62, 1.73]
times_image_10000 = [3.21, 1.70, 3.30, 8.22]

# Datos para 25,000 caracteres
times_total_25000 = [35.73, 29.26, 81.97, 84.66]
times_calculation_25000 = [30.32, 16.48, 11.01, 81.53]  
times_image_25000 = [5.35, 12.72, 70.9, 3.06]  

x = np.arange(len(processes))  # las etiquetas de los ejes x para los procesos
width = 0.3  # ancho de las barras

def create_chart(times_total, times_calculation, times_image, title, file_name):
    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width, times_total, width, label='Tiempo total')
    rects2 = ax.bar(x, times_calculation, width, label='Tiempo de cálculo')
    rects3 = ax.bar(x + width, times_image, width, label='Tiempo para guardar imagen')

    # Añadir algunas etiquetas de texto para títulos y etiquetas de ejes
    ax.set_xlabel('Número de procesos')
    ax.set_ylabel('Tiempo (s)')
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{p} procesos' for p in processes])
    ax.legend()

    fig.tight_layout()
    plt.savefig(f"./assets/imgs/{file_name}.png")
    plt.close()

# Generar gráficas para 10,000 caracteres
create_chart(times_total_10000, times_calculation_10000, times_image_10000,
             'Métricas de tiempo para MPI con tamaño de secuencia de 10,000 caracteres', 'tiempo_ejecucion_mpi_10000')

# Generar gráficas para 25,000 caracteres
create_chart(times_total_25000, times_calculation_25000, times_image_25000,
             'Métricas de tiempo para MPI con tamaño de secuencia de 25,000 caracteres', 'tiempo_ejecucion_mpi_25000')


# Sequential times for 10,000 and 25,000 characters
sequential_times = [5.96, 36.84]  # Adjust these if you have a specific time for 50,000 characters

# Minimum MPI times extracted from your provided data
min_mpi_times = [3.80, 29.26]  # Minimum times for 10,000 and 25,000 characters from MPI

# Labels for the plot
sizes = ['10,000', '25,000']  # You can add '50,000' if you have that data

# Creating the plot
fig, ax = plt.subplots()
x = np.arange(len(sizes))  # label locations
width = 0.35  # the width of the bars

rects1 = ax.bar(x - width/2, sequential_times, width, label='Secuencial')
rects2 = ax.bar(x + width/2, min_mpi_times, width, label='MPI Mejor tiempo')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_xlabel('Tamaño de la Secuencia (caracteres)')
ax.set_ylabel('Tiempo (s)')
ax.set_title('Comparación de tiempos de ejecución: Secuencial vs MPI')
ax.set_xticks(x)
ax.set_xticklabels(sizes)
ax.legend()

# Function to add a label above each bar, showing the time
def autolabel(rects):
    """Attach a text label above each bar displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)

fig.tight_layout()

# Save the plot to the specified directory instead of showing it
plt.savefig("./assets/imgs/mpi_vs_secuencial.png")
plt.close()
