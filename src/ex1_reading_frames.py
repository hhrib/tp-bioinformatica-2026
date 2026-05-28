# Ejercicio 1 - Procesamiento de Secuencias
# Lee un archivo GenBank de mRNA, calcula los 6 marcos de lectura posibles,
# traduce cada uno a aminoacidos y escribe los resultados en formato FASTA.
# El frame correcto se determina en el Ejercicio 2 via BLAST.
# Referencia: Cock et al. (2009) Bioinformatics 25(11):1422-1423

from Bio import SeqIO
from Bio.SeqUtils import six_frame_translations

INPUT_FILE = "data/NM_002111.gb"   # mRNA de HTT en formato GenBank
OUTPUT_FILE = "output/orfs.fasta"  # proteinas traducidas en formato FASTA


def get_six_frames(seq):
    frames = {}
    rc = seq.reverse_complement()
    for i in range(3):
        frames[f"+{i+1}"] = seq[i:].translate(to_stop=False)
        frames[f"-{i+1}"] = rc[i:].translate(to_stop=False)
    return frames


records = list(SeqIO.parse(INPUT_FILE, "genbank"))
print(f"Secuencias encontradas en el archivo: {len(records)}")

with open(OUTPUT_FILE, "w") as out:
    for record in records:
        seq = record.seq
        print(f"\nSecuencia: {record.id} ({len(seq)} bp)")

        frames = get_six_frames(seq)

        for label, protein in frames.items():
            out.write(f">{record.id}_frame_{label}\n{protein}\n")
            before_stop = str(protein).split("*")[0]
            print(f"  Frame {label}: {len(before_stop)} aa antes del primer stop")

print(f"\nOutput guardado en {OUTPUT_FILE}")

# --- COMPARACION: six_frame_translations() de Bio.SeqUtils ---
print("\n" + "="*60)
print("COMPARACION: Bio.SeqUtils.six_frame_translations()")
print("="*60)
for record in SeqIO.parse(INPUT_FILE, "genbank"):
    print(six_frame_translations(record.seq))
    break
