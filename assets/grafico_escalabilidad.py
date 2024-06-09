import matplotlib.pyplot as plt

# Datos de ejemplo
# Tiempos de ejecución para diferentes números de procesadores (strong scalability)
times_strong = [100, 60, 35, 20, 15]  # tiempos medidos en segundos
processors_strong = [1, 2, 4, 8, 16]  # número de procesadores correspondientes

# Tiempos de ejecución para diferentes tamaños de problema y números de procesadores (weak scalability)
times_weak = [100, 110, 115, 120, 125]  # tiempos medidos en segundos
processors_weak = [1, 2, 4, 8, 16]  # número de procesadores correspondientes

# Función para calcular speedup (strong scalability)
def calculate_speedup(times, processors):
    T1 = times[0]
    speedup = [T1 / t for t in times]
    return speedup

# Función para calcular eficiencia (weak scalability)
def calculate_efficiency(times, processors):
    T1 = times[0]
    efficiency = [T1 / (t * p) for t, p in zip(times, processors)]
    return efficiency

# Calcular speedup (strong scalability)
speedup_strong = calculate_speedup(times_strong, processors_strong)

# Calcular eficiencia (weak scalability)
efficiency_weak = calculate_efficiency(times_weak, processors_weak)

# Graficar resultados
plt.figure(figsize=(12, 5))

# Strong Scalability
plt.subplot(1, 2, 1)
plt.plot(processors_strong, speedup_strong, marker='o')
plt.xlabel('Number of Processors')
plt.ylabel('Speedup')
plt.title('Strong Scalability')
plt.grid(True)

# Weak Scalability
plt.subplot(1, 2, 2)
plt.plot(processors_weak, efficiency_weak, marker='o')
plt.xlabel('Number of Processors')
plt.ylabel('Efficiency')
plt.title('Weak Scalability')
plt.grid(True)

plt.tight_layout()
plt.show()
