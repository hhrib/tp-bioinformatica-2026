---
marp: true
theme: default
paginate: true
size: 16:9
header: 'TP Bioinformática 2026 — Enfermedad de Huntington / gen HTT'
style: |
  section { font-size: 24px; }
  h1 { color: #1a5276; }
  h2 { color: #1a5276; }
  table { font-size: 20px; }
  section.lead h1 { font-size: 48px; }
  .small { font-size: 18px; }
---

<!-- _class: lead -->
<!-- _paginate: false -->
<!-- _header: '' -->

# Enfermedad de Huntington
## Pipeline bioinformático del gen *HTT*

Introducción a la Bioinformática (K5061) — UTN FRBA · 2026
Profesor: Patricio Yankilevich

**Integrantes:** Tomás Guillermo León · Gustavo Di Peppe · Alejo Gurfein · Hernán Hrib

---

## Agenda

1. **Introducción** — la enfermedad y el gen
2. **Métodos** — stack y reproducibilidad
3. **Resultados** por ejercicio
   - Ej1 marcos de lectura · Ej2 BLAST · Ej3 MSA
   - Ej4 parser BLAST · Ej5 EMBOSS/PROSITE
   - Ej6 bases de datos biológicas
4. **Conclusiones**

> 10 min · ~15 slides

---

## La enfermedad de Huntington

- Enfermedad **neurodegenerativa hereditaria**, autosómica **dominante** (OMIM #143100).
- Deterioro progresivo: **movimiento, cognición y conducta**. Sin cura.
- **No hay portadores sanos:** quien hereda el alelo mutado desarrolla la enfermedad.

**Por qué la elegimos:** causa molecular única y muy clara → caso ideal para recorrer
todo el pipeline, de la secuencia a la variante patogénica.

---

## Causa molecular: expansión CAG

- Expansión inestable del triplete **CAG** en el exón 1 de *HTT*.
- Se traduce en un tracto de **poliglutaminas (polyQ)** tóxico para las neuronas.

| Nº repeticiones CAG | Fenotipo |
|---------------------|----------|
| ≤ 35 | Normal |
| 36–39 | Penetrancia incompleta |
| ≥ 40 | **Penetrancia completa — HD** |
| ≥ 60 | Forma juvenil |

> A más repeticiones → inicio más temprano y curso más grave.

---

## El gen *HTT*

- Cromosoma **4p16.3**, ~180 kb, **67 exones**.
- Codifica la **huntingtina** (~3144 aa): proteína de **andamiaje** del transporte axonal,
  transporte de BDNF, autofagia, desarrollo neuronal.
- Secuencia de referencia del TP: **NM_002111.8** (NCBI Nucleotide).

<!-- Captura sugerida: ficha NCBI Gene 3064 -->

---

## Métodos — ambiente de trabajo

- **Lenguaje:** Python 3.11 + **BioPython**
- **BLAST+** y **EMBOSS** (getorf, patmatmotifs)
- **Docker** sobre Windows → reproducibilidad
- **Git/GitHub** → control de versiones

```bash
docker compose run tp python src/ex1_reading_frames.py
docker compose run tp python src/ex2_blast.py
docker compose run tp python src/ex3_msa.py
docker compose run tp python src/ex4_blast_parser.py "Mus musculus"
docker compose run tp python src/ex5_emboss.py
```

---

## Ej1 — Marcos de lectura

Leer GenBank → 6 marcos (+1/+2/+3 y -1/-2/-3) → traducir → FASTA.

| Frame | aa antes del 1er stop | Observación |
|-------|----------------------|-------------|
| +1 | 131 | stop prematuro |
| **+2** | **3192** | **frame correcto** |
| +3 | 27 | stop prematuro |
| -1 / -2 / -3 | 23 / 13 / 5 | stop prematuro |

> El frame **+2** (3192 aa) coincide con el tamaño conocido (~3144 aa). La diferencia
> con el resto es enorme — pero el frame se **confirma** en el Ej2 con BLAST.

---

## Ej2 — BLAST (conservación evolutiva)

`blastp` remoto contra SwissProt. Solo el frame **+2** dio hits significativos.

| # | Proteína | Organismo | Identity | E-value |
|---|----------|-----------|----------|---------|
| 1 | Huntingtin | *Homo sapiens* | 99.9% | 0.0 |
| 2 | Huntingtin | *Mus musculus* | 91.2% | 0.0 |
| 3 | Huntingtin | *Rattus norvegicus* | 90.8% | 0.0 |
| 4 | Huntingtin | *Takifugu rubripes* | 69.7% | 0.0 |
| 5 | HD homolog | *Dictyostelium discoideum* | 28.8% | 8e-19 |

> Homólogo hasta en una **ameba** → gen **muy antiguo** (eucariotas en general).
> También se corrió BLAST **local** contra SwissProt (482.697 secuencias).

---

## Ej3 — Alineamiento múltiple (MSA)

- MSA de la huntingtina en **5 especies** (humano, ratón, rata, pez globo, *Dictyostelium*).
- Las posiciones **conservadas en todas** las especies son funcionalmente críticas.
- Salida: `output/msa.aln`.

<!-- Captura sugerida: fragmento de msa.aln mostrando columnas conservadas -->

---

## Ej4 — Parser de salida BLAST

- Script que **parsea** `blast.out` y filtra hits cuya descripción contiene un **Pattern**
  (ej. `"Mus musculus"`, case-insensitive).
- Punto extra: extrae el **ACCESSION** de cada hit y baja la secuencia completa en FASTA
  con `Bio.Entrez` (equivalente a `Bio::DB::GenBank`).

```bash
docker compose run tp python src/ex4_blast_parser.py "Mus musculus"
# → output/ex4_hits.fasta
```

---

## Ej5 — EMBOSS: ORFs + dominios PROSITE

- `getorf` → ORFs del mRNA; se elige el más largo = **huntingtina (3144 aa)**.
- `prosextract` + `patmatmotifs` → análisis de dominios contra **PROSITE**.

**Motivos encontrados (HitCount: 4):**

| Motivo | Posición |
|--------|----------|
| AMIDATION | 2545–2548 |
| AMIDATION | 1532–1535 |
| LEUCINE_ZIPPER | 1446–1467 |
| TYR_PHOSPHO_SITE_2 | 2716–2723 |

---

## Ej6 — Bases de datos biológicas (1/2)

| Pregunta | Resultado (verificado jun 2026) |
|----------|----------------------------------|
| Homólogos | HomoloGene (legacy, ~9–11 spp) vs **Ensembl: 412 ortólogos** |
| Transcriptos | RefSeq **2 reviewed** vs **Ensembl 24** (10 coding) |
| Interacciones | HAP1, HIP1, HIP14, DCTN1… → patrón: **transporte sobre microtúbulos** |
| GO | citoplasma/vesículas · transporte axonal · *protein binding* |
| Pathways | KEGG **hsa05016** "Huntington disease" · Reactome |

> Diferencia metodológica clave: RefSeq curado/conservador vs Ensembl automático/exhaustivo.

---

## Ej6 — La variante patogénica (2/2)

- No es un SNP clásico: es una **expansión de microsatélite** (CAG).
- **ClinVar Variation ID 409** — *Pathogenic*, review *practice guideline*.

**Frecuencia / etnia:**
- Mayor en **ascendencia europea** (~5–15/100.000), haplogrupo A — efecto fundador.
- Mucho menor en Asia oriental y África subsahariana.
- Caso emblemático: **Lago de Maracaibo** (Venezuela) → clonado del gen en 1993.

<!-- Captura sugerida: entrada ClinVar Variation ID 409 -->

---

## Conclusiones

- *HTT* es un gen **antiguo y muy conservado** (de amebas a humanos).
- El pipeline completo —de la secuencia a la variante— confirmó la **causa molecular**:
  expansión **CAG / polyQ**.
- Las herramientas coinciden y se complementan: BLAST, MSA, EMBOSS/PROSITE y bases de datos
  cuentan la misma historia desde ángulos distintos.
- Diferencias entre bases (RefSeq vs Ensembl, HomoloGene legacy) → importa **saber qué
  base usar** según el objetivo.

---

<!-- _class: lead -->
<!-- _paginate: false -->

# ¡Gracias!

¿Preguntas?

<span class="small">Repo: github.com/hhrib/tp-bioinformatica-2026</span>
