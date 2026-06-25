# Ejercicio 4 - Parser de la salida de BLAST
# Lee el reporte blast.out del Ejercicio 2 y busca los hits que en su
# descripcion contienen un Pattern que pasamos como parametro (ej. "Mus musculus").
# Despues, para cada hit encontrado, baja la secuencia completa en formato FASTA
# desde NCBI usando Bio.Entrez (el equivalente al modulo Bio::DB::GenBank).
# Uso: python src/ex4_blast_parser.py "Mus musculus"

import sys
from Bio import SeqIO, Entrez

# NCBI pide un email para saber quien hace las consultas
Entrez.email = "hernanhrib@gmail.com"

archivo_blast = "output/blast.out"        # reporte del Ejercicio 2 (entrada)
archivo_salida = "output/ex4_hits.fasta"  # secuencias de los hits (salida)

# El Pattern es el primer parametro. Si no se pasa ninguno, usamos "Mus musculus".
if len(sys.argv) > 1:
    pattern = sys.argv[1]
else:
    pattern = "Mus musculus"

print('Buscando hits que contengan: "' + pattern + '"')
print()

# Recorremos el archivo linea por linea. Guardamos en que frame estamos y
# nos quedamos con las lineas que contienen el Pattern.
hits = []      # cada elemento es [frame, accession, linea]
frame = "?"

archivo = open(archivo_blast, encoding="utf-8")
for linea in archivo:
    linea = linea.strip()
    if linea.startswith("Frame"):
        # ejemplo de linea: "Frame +2 - 5 hits"  -> nos quedamos con el "+2"
        frame = linea.split()[1]
    elif pattern.lower() in linea.lower():
        # si la linea contiene el Pattern, es un hit que nos sirve.
        # el accession es lo primero de la linea (ej. "P42859.2")
        accession = linea.split()[0]
        hits.append([frame, accession, linea])
archivo.close()

# Si no encontramos nada, avisamos y terminamos
if len(hits) == 0:
    print("No se encontraron hits que coincidan con el pattern.")
    sys.exit()

print("Hits encontrados: " + str(len(hits)))
print()
for frame, accession, linea in hits:
    print("Frame " + frame + " -> " + linea)

# Para cada hit bajamos la secuencia completa en FASTA desde NCBI
print()
print("Bajando las secuencias completas desde NCBI...")
secuencias = []
for frame, accession, linea in hits:
    print("  bajando " + accession + " ...")
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
    registro = SeqIO.read(handle, "fasta")
    handle.close()
    secuencias.append(registro)

# Guardamos todas las secuencias en un solo archivo FASTA
SeqIO.write(secuencias, archivo_salida, "fasta")
print()
print(str(len(secuencias)) + " secuencia(s) guardada(s) en " + archivo_salida)
