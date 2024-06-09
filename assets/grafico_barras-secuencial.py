import matplotlib.pyplot as plt

# Datos
tamanos_secuencia = ['10,000', '25,000', '50,000']
tiempos_ejecucion = [5.96, 36.84, 666.69]  # Tiempos totales de ejecución para cada tamaño

# Crear el gráfico de barras
plt.figure(figsize=(10, 6))
plt.bar(tamanos_secuencia, tiempos_ejecucion, color='blue')

# Añadir título y etiquetas
plt.title('Tiempo Total de Ejecución vs. Tamaño de Secuencia')
plt.xlabel('Tamaño de la Secuencia (caracteres)')
plt.ylabel('Tiempo Total de Ejecución (segundos)')

# Mostrar valores en las barras
for i, tiempo in enumerate(tiempos_ejecucion):
    plt.text(i, tiempo + 10, f'{tiempo:.2f}s', ha = 'center', color = 'black')

# Guardar el gráfico como imagen PNG
plt.tight_layout()
plt.savefig('/home/erika/20182.pcomp/proyecto_final_concurrente/assets/imgs/tiempo_ejecucion_dotplot-secuencial.png')

# Opcional: Cierra el gráfico para liberar memoria
plt.close()
