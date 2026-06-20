# Ejercicio 5 - EMBOSS
# 1) Calcula los ORFs del mRNA de HTT y traduce a proteina con `getorf`.
# 2) Selecciona el ORF mas largo (la huntingtina, ~3144 aa).
# 3) Baja la base de datos PROSITE (prosite.dat + prosite.doc) y la prepara
#    con `prosextract`.
# 4) Analiza los dominios/motivos funcionales del ORF con `patmatmotifs`.
#
# Input:  data/NM_002111.gb   (mRNA de HTT en GenBank)
# Output: output/ex5_orfs.fasta     (todos los ORFs traducidos por getorf)
#         output/ex5_longest_orf.fasta  (el ORF mas largo = proteina analizada)
#         output/ex5_domains.txt    (dominios PROSITE encontrados)
#
# Requiere EMBOSS instalado (paquete `emboss` en el Dockerfile).
# Referencia: Rice et al. (2000) Trends in Genetics 16(6):276-277 (EMBOSS)

import os
import subprocess
import urllib.request
from Bio import SeqIO

GENBANK_IN   = "data/NM_002111.gb"
MRNA_FASTA   = "output/ex5_mrna.fasta"
ORFS_FASTA   = "output/ex5_orfs.fasta"
LONGEST_ORF  = "output/ex5_longest_orf.fasta"
DOMAINS_OUT  = "output/ex5_domains.txt"

PROSITE_DIR  = "data/prosite"   # cache en volumen montado (no se re-baja)
PROSITE_URLS = {
    "prosite.dat": "https://ftp.expasy.org/databases/prosite/prosite.dat",
    "prosite.doc": "https://ftp.expasy.org/databases/prosite/prosite.doc",
}


def run(cmd):
    """Ejecuta un comando EMBOSS y aborta si falla."""
    print(f"  $ {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)


# --- Paso 1: mRNA GenBank -> FASTA de nucleotidos ---
print("== Paso 1: convertir mRNA a FASTA (seqret) ==")
run(["seqret", "-sequence", GENBANK_IN, "-outseq", MRNA_FASTA, "-auto"])

# --- Paso 2: calcular ORFs y traducir a proteina (getorf) ---
# -find 1 = traducir las regiones entre un codon START (ATG) y el STOP.
# -minsize en nucleotidos: 300 nt = 100 aa, descarta ORFs cortos espurios.
print("\n== Paso 2: calcular ORFs y traducir (getorf) ==")
run(["getorf", "-sequence", MRNA_FASTA, "-outseq", ORFS_FASTA,
     "-find", "1", "-minsize", "300", "-auto"])

orfs = list(SeqIO.parse(ORFS_FASTA, "fasta"))
print(f"  getorf encontro {len(orfs)} ORF(s) >= 100 aa")

# --- Paso 3: elegir el ORF mas largo (la huntingtina) ---
longest = max(orfs, key=lambda r: len(r.seq))
SeqIO.write([longest], LONGEST_ORF, "fasta")
print(f"\n== Paso 3: ORF mas largo ==")
print(f"  {longest.id}  ({len(longest.seq)} aa)  -> {LONGEST_ORF}")

# --- Paso 4: bajar y preparar PROSITE ---
print(f"\n== Paso 4: preparar PROSITE ==")
os.makedirs(PROSITE_DIR, exist_ok=True)
for name, url in PROSITE_URLS.items():
    dest = os.path.join(PROSITE_DIR, name)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"  {name} ya esta en cache ({os.path.getsize(dest)//1024} KB)")
        continue
    print(f"  bajando {name} ...", flush=True)
    urllib.request.urlretrieve(url, dest)
    print(f"  -> {dest} ({os.path.getsize(dest)//1024} KB)")

# prosextract indexa prosite.dat/.doc en el directorio de datos de EMBOSS,
# dejando la base lista para que patmatmotifs la consulte.
print("  indexando con prosextract ...")
run(["prosextract", "-prositedir", PROSITE_DIR, "-auto"])

# --- Paso 5: analisis de dominios/motivos (patmatmotifs) ---
print(f"\n== Paso 5: analisis de dominios (patmatmotifs) ==")
run(["patmatmotifs", "-sequence", LONGEST_ORF, "-outfile", DOMAINS_OUT,
     "-full", "Y", "-auto"])

# resumen: contar y listar los motivos encontrados
motifs = []
with open(DOMAINS_OUT, encoding="utf-8") as f:
    for line in f:
        if line.strip().startswith("Motif ="):
            motifs.append(line.split("=", 1)[1].strip())

print(f"\nDominios/motivos PROSITE encontrados: {len(motifs)}")
for m in motifs:
    print(f"  - {m}")
print(f"\nResultado completo en {DOMAINS_OUT}")
