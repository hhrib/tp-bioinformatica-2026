# Ejercicio 2a - BLAST remoto
# Toma la proteína del frame correcto (output del Ej1) y ejecuta blastp
# contra SwissProt en el servidor del NCBI. Guarda el reporte y muestra
# los mejores hits con su E-value, identidad y score.
# Referencia: Cock et al. (2009) Bioinformatics 25(11):1422-1423

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

ORFS_FILE = "output/orfs.fasta"  # output del Ejercicio 1
OUTPUT_FILE = "output/blast.out" # reporte BLAST en formato XML
CORRECT_FRAME = "+2"             # frame identificado en el Ejercicio 1
TOP_HITS = 10                    # cantidad de hits a mostrar

# extraer la secuencia proteica del frame correcto
# se corta en el primer codón de stop (*) para obtener solo la proteína real
protein = None
for record in SeqIO.parse(ORFS_FILE, "fasta"):
    if f"frame_{CORRECT_FRAME}" in record.id:
        protein = str(record.seq).split("*")[0]
        break

if not protein:
    print(f"No se encontro el frame {CORRECT_FRAME} en {ORFS_FILE}")
    exit(1)

print(f"Secuencia del frame {CORRECT_FRAME}: {len(protein)} aa")
print("Ejecutando BLAST remoto contra SwissProt (puede tardar 2-5 minutos)...")

# ejecutar blastp de forma remota usando el servidor del NCBI
# blastp: compara proteína contra base de datos de proteínas
# swissprot: base de datos curada con proteínas de alta calidad
result = NCBIWWW.qblast(
    program="blastp",
    database="swissprot",
    sequence=protein,
)

print("BLAST finalizado. Guardando resultados...")

# guardar el resultado XML completo en disco
raw = result.read()
with open(OUTPUT_FILE, "w") as out:
    out.write(raw)

print(f"Output guardado en {OUTPUT_FILE}")
print(f"\nTop {TOP_HITS} hits:\n")

# parsear el XML y mostrar los mejores hits con sus métricas estadísticas
with open(OUTPUT_FILE) as f:
    blast_records = list(NCBIXML.parse(f))

for blast_record in blast_records:
    for i, alignment in enumerate(blast_record.alignments[:TOP_HITS]):
        hsp = alignment.hsps[0]  # hsp: high-scoring segment pair (mejor alineamiento)
        identity_pct = (hsp.identities / hsp.align_length) * 100
        print(f"{i+1:2}. {alignment.title[:70]}")
        print(f"    E-value: {hsp.expect:.2e}  Identity: {identity_pct:.1f}%  Score: {hsp.score}")
        print()
