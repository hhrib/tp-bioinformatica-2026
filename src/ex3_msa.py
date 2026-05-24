# Ejercicio 3 - Multiple Sequence Alignment (MSA)
# Descarga las secuencias proteicas de los hits encontrados en el BLAST
# y las guarda en un archivo FASTA listo para alinear con Clustal Omega.
# El MSA se realiza online en https://www.ebi.ac.uk/jdispatcher/
# Referencia: Cock et al. (2009) Bioinformatics 25(11):1422-1423

from Bio import SeqIO, Entrez
from Bio.Blast import NCBIXML

# email requerido por NCBI para identificar al usuario en consultas remotas
Entrez.email = "hernanhrib@gmail.com"

BLAST_FILE = "output/blast.out"        # reporte del Ejercicio 2
SEQUENCES_FILE = "output/msa_input.fasta"  # secuencias a alinear

# extraer los IDs de las proteínas encontradas en el BLAST
# el formato del hit_id es: sp|P42858.2|HD_HUMAN → tomamos P42858.2
print("Leyendo resultados del BLAST...")
with open(BLAST_FILE) as f:
    blast_records = list(NCBIXML.parse(f))

hit_ids = []
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        parts = alignment.hit_id.split("|")
        if len(parts) >= 2:
            hit_ids.append(parts[1])

print(f"Hits encontrados: {hit_ids}")

# descargar cada secuencia proteica de NCBI en formato FASTA
print("\nBajando secuencias de NCBI...")
sequences = []
for hit_id in hit_ids:
    print(f"  Bajando {hit_id}...")
    handle = Entrez.efetch(db="protein", id=hit_id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    sequences.append(record)
    print(f"  {record.description[:70]} ({len(record.seq)} aa)")

# guardar todas las secuencias juntas en un solo FASTA para el MSA
with open(SEQUENCES_FILE, "w") as out:
    SeqIO.write(sequences, out, "fasta")

print(f"\n{len(sequences)} secuencias guardadas en {SEQUENCES_FILE}")
print("\nPara hacer el MSA:")
print("  1. Abri https://www.ebi.ac.uk/jdispatcher/")
print("  2. Multiple Sequence Alignment → Clustal Omega")
print("  3. Subi el archivo output/msa_input.fasta")
print("  4. Descarga el resultado y guardalo en output/msa.aln")
