# Ejercicio 2a - BLAST remoto
# Corre blastp remoto contra SwissProt para los 6 frames del Ejercicio 1.
# El frame correcto se identifica como el único con hits significativos.
# Genera blast.xml (para el Ejercicio 3) y blast.out (reporte legible).
# Referencia: Cock et al. (2009) Bioinformatics 25(11):1422-1423

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from io import StringIO

ORFS_FILE = "output/orfs.fasta"  # output del Ejercicio 1
OUTPUT_XML  = "output/blast.xml" # XML para consumir en Ejercicio 3
OUTPUT_TEXT = "output/blast.out" # reporte legible para humanos
TOP_HITS = 5

# leer los 6 frames y cortar cada uno en el primer codón de stop
frames = {}
for record in SeqIO.parse(ORFS_FILE, "fasta"):
    label = record.id.split("frame_")[1].replace(" <--", "").strip()
    protein = str(record.seq).split("*")[0]
    frames[label] = protein

print(f"Frames a procesar: {list(frames.keys())}\n")

all_results = []  # lista de (label, xml_string, blast_record)

for label, protein in frames.items():
    print(f"Frame {label}: {len(protein)} aa — ejecutando BLAST remoto...", flush=True)
    result = NCBIWWW.qblast(program="blastp", database="swissprot", sequence=protein)
    xml = result.read()
    blast_record = list(NCBIXML.parse(StringIO(xml)))[0]
    all_results.append((label, xml, blast_record))
    hits = len(blast_record.alignments)
    best = blast_record.alignments[0].hsps[0].expect if hits > 0 else "-"
    print(f"  -> {hits} hits, mejor E-value: {best}\n")

# guardar XML completo (para Ejercicio 3)
with open(OUTPUT_XML, "w") as f:
    f.write("\n".join(xml for _, xml, _ in all_results))

# guardar reporte legible en texto plano
with open(OUTPUT_TEXT, "w") as f:
    for label, _, blast_record in all_results:
        f.write(f"\n{'='*60}\nFrame {label} — {len(blast_record.alignments)} hits\n{'='*60}\n\n")
        if not blast_record.alignments:
            f.write("Sin hits significativos.\n")
            continue
        f.write(f"{'Accesión':<15} {'Organismo':<30} {'Identity':>9} {'E-value':>12} {'Score':>8}\n")
        f.write("-" * 78 + "\n")
        for alignment in blast_record.alignments[:TOP_HITS]:
            hsp = alignment.hsps[0]
            identity_pct = (hsp.identities / hsp.align_length) * 100
            # extraer organismo del título: "... [Homo sapiens]"
            org_match = alignment.title.split("[")[-1].rstrip("]") if "[" in alignment.title else "?"
            # extraer accesión: sp|P42858.2|HD_HUMAN → P42858.2
            accession = alignment.hit_id.split("|")[1] if "|" in alignment.hit_id else alignment.hit_id
            f.write(f"{accession:<15} {org_match:<30} {identity_pct:>8.1f}% {hsp.expect:>12.2e} {hsp.score:>8.0f}\n")

print(f"XML guardado en {OUTPUT_XML}")
print(f"Reporte de texto guardado en {OUTPUT_TEXT}")

# resumen por pantalla
print("\n" + "="*60)
print("RESUMEN por frame")
print("="*60)

for label, _, blast_record in all_results:
    sig = [a for a in blast_record.alignments if a.hsps[0].expect < 0.001]
    marker = " <-- FRAME CORRECTO" if sig else ""
    print(f"\nFrame {label}: {len(sig)} hits significativos{marker}")
    for i, alignment in enumerate(blast_record.alignments[:TOP_HITS]):
        hsp = alignment.hsps[0]
        identity_pct = (hsp.identities / hsp.align_length) * 100
        print(f"  {i+1}. {alignment.title[:65]}")
        print(f"     E-value: {hsp.expect:.2e}  Identity: {identity_pct:.1f}%  Score: {hsp.score}")
