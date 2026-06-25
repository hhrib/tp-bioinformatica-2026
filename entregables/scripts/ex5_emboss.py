# Ejercicio 5 - EMBOSS
# Este script usa programas de EMBOSS para analizar la secuencia del mRNA.
# Pasos:
#  1) calcula los ORFs del mRNA y los traduce a secuencias de proteinas (getorf)
#  2) se queda con el ORF mas largo (la proteina que estudiamos)
#  3) baja la base de datos PROSITE (prosite.dat) de dominios/motivos funcionales
#  4) busca esos dominios en la secuencia de aminoacidos (patmatmotifs)
#
# Input:  data/NM_002111.gb
# Output: output/ex5_orfs.fasta, output/ex5_longest_orf.fasta, output/ex5_domains.txt

import os
import subprocess
import urllib.request
from Bio import SeqIO

mrna_genbank = "data/NM_002111.gb"
mrna_fasta   = "output/ex5_mrna.fasta"
orfs_fasta   = "output/ex5_orfs.fasta"
orf_largo    = "output/ex5_longest_orf.fasta"
salida_dominios = "output/ex5_domains.txt"

carpeta_prosite = "data/prosite"
urls_prosite = {
    "prosite.dat": "https://ftp.expasy.org/databases/prosite/prosite.dat",
    "prosite.doc": "https://ftp.expasy.org/databases/prosite/prosite.doc",
}


# funcion para correr un programa de EMBOSS
def correr(comando):
    print("  $ " + " ".join(comando))
    subprocess.run(comando, check=True)


# Paso 1: pasar el mRNA de formato GenBank a FASTA
print("== Paso 1: pasar el mRNA a FASTA ==")
correr(["seqret", "-sequence", mrna_genbank, "-outseq", mrna_fasta, "-auto"])

# Paso 2: calcular los ORFs y traducirlos a proteina
# -find 1 traduce desde el codon de inicio hasta el de stop.
# -minsize 300 descarta los ORFs muy cortos (300 nucleotidos = 100 aminoacidos).
print()
print("== Paso 2: calcular ORFs (getorf) ==")
correr(["getorf", "-sequence", mrna_fasta, "-outseq", orfs_fasta,
        "-find", "1", "-minsize", "300", "-auto"])

orfs = list(SeqIO.parse(orfs_fasta, "fasta"))
print("  getorf encontro " + str(len(orfs)) + " ORFs")

# Paso 3: quedarnos con el ORF mas largo
mas_largo = orfs[0]
for orf in orfs:
    if len(orf.seq) > len(mas_largo.seq):
        mas_largo = orf
SeqIO.write([mas_largo], orf_largo, "fasta")
print()
print("== Paso 3: ORF mas largo ==")
print("  " + mas_largo.id + " (" + str(len(mas_largo.seq)) + " aa)")

# Paso 4: bajar la base de datos PROSITE (solo si no la tenemos ya)
print()
print("== Paso 4: preparar PROSITE ==")
if not os.path.exists(carpeta_prosite):
    os.makedirs(carpeta_prosite)
for nombre in urls_prosite:
    destino = os.path.join(carpeta_prosite, nombre)
    if os.path.exists(destino):
        print("  " + nombre + " ya esta descargado")
    else:
        print("  bajando " + nombre + " ...")
        urllib.request.urlretrieve(urls_prosite[nombre], destino)
# prosextract deja la base PROSITE lista para usarla
correr(["prosextract", "-prositedir", carpeta_prosite, "-auto"])

# Paso 5: buscar los dominios/motivos en la secuencia de aminoacidos
print()
print("== Paso 5: buscar dominios (patmatmotifs) ==")
correr(["patmatmotifs", "-sequence", orf_largo, "-outfile", salida_dominios,
        "-full", "Y", "-auto"])

# contar los motivos que encontro, leyendo el archivo de salida
motivos = []
archivo = open(salida_dominios, encoding="utf-8")
for linea in archivo:
    if linea.strip().startswith("Motif ="):
        motivos.append(linea.split("=")[1].strip())
archivo.close()

print()
print("Dominios/motivos encontrados: " + str(len(motivos)))
for m in motivos:
    print("  - " + m)
print("Resultado completo en " + salida_dominios)
