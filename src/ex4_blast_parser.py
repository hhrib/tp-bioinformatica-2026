# Ejercicio 4 - Parser de salida BLAST
# Parsea el reporte blast.out del Ejercicio 2 e identifica los hits cuya
# descripcion contiene un Pattern dado como parametro (ej. "Mus musculus").
# Punto extra: para cada hit que matchea extrae su ACCESSION y baja la
# secuencia completa del hit en formato FASTA usando Bio.Entrez (equivalente
# al modulo Bio::DB::GenBank de BioPerl), escribiendola a un archivo.
# Uso: python src/ex4_blast_parser.py "Mus musculus"
# Referencia: Cock et al. (2009) Bioinformatics 25(11):1422-1423

import re
import sys
from Bio import SeqIO, Entrez

# email requerido por NCBI para identificar al usuario en consultas remotas
Entrez.email = "hernanhrib@gmail.com"

BLAST_FILE = "output/blast.out"        # reporte legible del Ejercicio 2 (input)
HITS_FASTA = "output/ex4_hits.fasta"   # secuencias completas de los hits (output)

# Pattern de busqueda: primer argumento de la linea de comandos.
# Por defecto "Mus musculus" como en el ejemplo de la consigna.
pattern = sys.argv[1] if len(sys.argv) > 1 else "Mus musculus"

# Cada fila de datos de blast.out arranca con un ACCESSION tipo P42859.2 o
# Q76P24.1: una o mas letras/numeros, un punto y la version. Esta regex lo
# captura y nos sirve para distinguir filas de hit de las lineas decorativas.
ACCESSION_RE = re.compile(r"^([A-Z0-9]+\.\d+)\s")

print(f'Buscando hits cuya descripcion contenga el pattern: "{pattern}"\n')

# Recorrer el reporte. Guardamos el "Frame N" actual para reportar en que
# frame aparecio cada hit, y nos quedamos con las filas que matchean.
matches = []  # lista de (frame, accession, linea_completa)
current_frame = "?"

with open(BLAST_FILE, encoding="utf-8") as f:
    for line in f:
        line = line.rstrip("\n")
        frame_header = re.search(r"Frame\s+(\S+)", line)
        if frame_header:
            current_frame = frame_header.group(1)
            continue
        m = ACCESSION_RE.match(line)
        if not m:
            continue  # no es una fila de hit (separadores, encabezados, etc.)
        # match del pattern contra la fila completa (accession + organismo),
        # case-insensitive para que "mus musculus" y "Mus musculus" coincidan
        if pattern.lower() in line.lower():
            matches.append((current_frame, m.group(1), line.strip()))

# Reporte de los hits encontrados
if not matches:
    print("No se encontraron hits que coincidan con el pattern.")
    sys.exit(0)

print(f"Hits que coinciden con \"{pattern}\": {len(matches)}\n")
print(f"{'Frame':<8} {'Accession':<14} Descripcion")
print("-" * 70)
for frame, accession, line in matches:
    print(f"{frame:<8} {accession:<14} {line}")

# --- Punto extra: bajar la secuencia completa de cada hit en FASTA ---
print(f"\nBajando secuencias completas de NCBI (db=protein)...")
sequences = []
for frame, accession, _ in matches:
    print(f"  Bajando {accession}...", flush=True)
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    sequences.append(record)
    print(f"  {record.description[:70]} ({len(record.seq)} aa)")

with open(HITS_FASTA, "w", encoding="utf-8") as out:
    SeqIO.write(sequences, out, "fasta")

print(f"\n{len(sequences)} secuencia(s) guardada(s) en {HITS_FASTA}")
