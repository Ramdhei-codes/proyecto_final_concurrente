# Proyecto Dotplot

Este proyecto genera gráficos de dotplot para comparar secuencias de ADN utilizando diferentes enfoques de programación: secuencial, multiprocessing, hilos, MPI, y PyCUDA.

## Descripción

El objetivo de este proyecto es comparar el rendimiento de diferentes enfoques de programación para generar gráficos de dotplot, que son herramientas visuales utilizadas para comparar secuencias de ADN.

Las imagenes dependiendo de cada una de las ejecuciones se crean sobre la carpte raíz según lo que especifiquemos en el --output= del comando de ejecución de cada archivo.

## Instalación

1. Clona este repositorio:
    ```bash
    git clone https://github.com/Ramdhei-codes/proyecto_final_concurrente.git
    ```

2. Instala las dependencias:
    ```bash
    pip install -r requirements.txt
    ```

## Uso

### Ejecución Secuencial

Para generar un dotplot utilizando el enfoque secuencial:
```bash
python dotplot_secuencial.py --file1=./dotplot_files/E_coli.fna --file2=./dotplot_files/Salmonella.fna --output=dotplot_secuencial.png --max_length=10000
```

### Ejecución Multiprocessing
Para generar un dotplot utilizando multiprocessing:
```bash
python prueba_multiprocessing.py --file1=./dotplot_files/E_coli.fna --file2=./dotplot_files/Salmonella.fna --output=dotplot_multiprocessing.png --max_length=10000 --num_processes=4
```

### Ejecución con Hilos
Para generar un dotplot utilizando hilos:
```bash
python dotplot_hilos.py --file1=./dotplot_files/E_coli.fna --file2=./dotplot_files/Salmonella.fna --output=dotplot_hilos.png --max_length=10000 --num_threads=4

```

### Ejecución con MPI
Para generar un dotplot utilizando MPI:
```bash
mpirun -np 4 python dotplot_mpi4py.py --file1=./dotplot_files/E_coli.fna --file2=./dotplot_files/Salmonella.fna --output=dotplot_mpi.png --max_length=10000
```
- Cabe recalcar que el número de procesos en este caso se especifica en el # que se encuentra despues del -np

### Ejecución con PyCUDA
Para generar un dotplot utilizando PyCUDA:

```python
# Definir argumentos para el entorno
class Args:
    file1 = './dotplot_files/E_coli.fna'
    file2 = './dotplot_files/Salmonella.fna'
    output = 'dotplot_cuda.png'
    max_length = 25000
    block_size = 25000

args = Args()
main(args.file1, args.file2, args.output, args.max_length, args.block_size)
```
- Este dotplot se recomienda ejecutarlo en google colab debido a que pycuda utiliza la GPU y la que nos presta este entorno tiene una buna capacidad de procesamiento,
el enlace al colab en caso de querer probarlo allí es https://colab.research.google.com/drive/1nWXGivroH174YhhpeT9NmKySMeng0j3H?usp=sharing, en otro caso se puede abrir en el navegador el archivo dotplot_pycda.ipynb que se encuentra en la carpeta raíz del proyecto, donde lo unico que debe hacerse para ambos casos de prueba es subir los archivos fasta que se encuentran en la carpeta dotplot_files y corregir las rutas según lo necesario.

### Contribución
- Haz un fork del proyecto.
- Crea una nueva rama (git checkout -b feature/nueva-caracteristica).
- Realiza los cambios y haz commit (git commit -am 'Agrega nueva característica').
- Envía los cambios a tu rama (git push origin feature/nueva-caracteristica).
- Abre un Pull Request.

### Autores
- Erika Fernanda Orrego Giraldo
- Jojhan Perez Arroyave
- Ramdhei Lopez Arcila