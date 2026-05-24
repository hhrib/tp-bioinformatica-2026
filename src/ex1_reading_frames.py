# Ejercicio 1 - Procesamiento de Secuencias
# Lee un archivo GenBank de mRNA, calcula los 6 marcos de lectura posibles,
# traduce cada uno a aminoácidos y escribe los resultados en formato FASTA.
# Referencia: Cock et al. (2009) Bioinformatics 25(11):1422-1423

from Bio import SeqIO
from Bio.Seq import Seq

INPUT_FILE = "data/NM_002111.gb"   # mRNA de HTT en formato GenBank
OUTPUT_FILE = "output/orfs.fasta"  # proteínas traducidas en formato FASTA


def get_six_frames(seq):
    """
    Calcula los 6 marcos de lectura posibles de una secuencia de nucleótidos.
    +1, +2, +3: cadena directa empezando en base 1, 2 o 3
    -1, -2, -3: complemento reverso empezando en base 1, 2 o 3
    """
    frames = {}
    rc = seq.reverse_complement()  # complemento reverso para los frames negativos
    for i in range(3):
        frames[f"+{i+1}"] = seq[i:].translate(to_stop=False)   # frame directo
        frames[f"-{i+1}"] = rc[i:].translate(to_stop=False)    # frame reverso
    return frames


def find_best_frame(frames):
    """
    Identifica el marco de lectura correcto: el que produce la proteína
    más larga antes del primer codón de stop (*).
    """
    best_frame = None
    best_length = 0
    for label, protein in frames.items():
        seq_str = str(protein)
        before_stop = seq_str.split("*")[0]  # tomar solo hasta el primer stop
        if len(before_stop) > best_length:
            best_length = len(before_stop)
            best_frame = label
    return best_frame


# leer el archivo GenBank — contiene la secuencia y metadatos del gen HTT
records = list(SeqIO.parse(INPUT_FILE, "genbank"))
print(f"Secuencias encontradas en el archivo: {len(records)}")

with open(OUTPUT_FILE, "w") as out:
    for record in records:
        seq = record.seq
        print(f"\nSecuencia: {record.id} ({len(seq)} bp)")

        # calcular los 6 marcos de lectura y encontrar el correcto
        frames = get_six_frames(seq)
        best = find_best_frame(frames)

        # escribir cada frame en el archivo FASTA de salida
        for label, protein in frames.items():
            marker = " <-- FRAME CORRECTO" if label == best else ""
            header = f">{record.id}_frame_{label}{marker}"
            out.write(f"{header}\n{protein}\n")
            before_stop = str(protein).split("*")[0]
            print(f"  Frame {label}: {len(before_stop)} aa antes del primer stop{marker}")

        print(f"\n  Frame correcto: {best}")

print(f"\nOutput guardado en {OUTPUT_FILE}")
