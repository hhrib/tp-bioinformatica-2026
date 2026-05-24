# TP Bioinformática 2026

**Materia:** Introducción a la Bioinformática — UTN FRBA Ingeniería en Sistemas  
**Entrega informe de avance:** lunes 1 de junio 2026

---

## Enfermedad y gen

**Enfermedad:** Enfermedad de Huntington  
**Gen:** HTT (Huntingtin) — `NM_002111`  
**OMIM:** [#143100](https://omim.org/entry/143100)

La enfermedad de Huntington es una enfermedad neurodegenerativa hereditaria autosómica dominante. Una repetición anormal del triplete CAG en el gen HTT (más de 36 repeticiones) causa inevitablemente la enfermedad. Quien hereda el gen mutado la desarrolla — no hay portadores sanos.

---

## Stack

- Python 3.11 + BioPython
- BLAST+
- Docker

---

## Cómo ejecutar

```bash
docker compose up --build
docker compose run tp python src/ex1_reading_frames.py
docker compose run tp python src/ex2_blast.py
docker compose run tp python src/ex3_msa.py
```

---

## Estructura

```
tp-bioinformatica-2026/
├── Dockerfile
├── docker-compose.yml
├── data/
│   └── NM_002111.gbk        # mRNA de HTT en formato GenBank (input Ej1)
├── output/
│   ├── orfs.fasta            # 6 reading frames traducidos (output Ej1)
│   ├── blast.out             # reporte BLAST (output Ej2)
│   └── msa.aln               # alineamiento múltiple (output Ej3)
├── src/
│   ├── ex1_reading_frames.py
│   ├── ex2_blast.py
│   └── ex3_msa.py
└── informe/
    └── informe_avance.pdf
```

---

## Progreso

- [x] Fase 0 — Setup: repo, Dockerfile, estructura de carpetas
- [x] Fase 1 — Datos: descargar `NM_002111.gb` de NCBI → `data/`
- [x] Fase 2 — Ejercicio 1: 6 reading frames → `output/orfs.fasta` (frame correcto: +2, 3192 aa)
- [x] Fase 3 — Ejercicio 2a: BLAST remoto → `output/blast.out` (5 hits, top hit 99.9% identity humano, hit 5 Dictyostelium discoideum 28.8% — gen de ~1000M años de antigüedad)
- [x] Fase 4 — Ejercicio 2b: interpretación de resultados BLAST
- [x] Fase 5 — Ejercicio 3: MSA de HTT en 5 especies (humano, ratón, rata, pez globo, Dictyostelium) → `output/msa.aln`
- [ ] Fase 6 — Informe de avance PDF

---

## Ejercicios

### Ejercicio 1 — Reading Frames

Lee el archivo GenBank de HTT, calcula los 6 marcos de lectura posibles (+1/+2/+3 en cadena directa y -1/-2/-3 en complemento reverso), traduce cada uno a aminoácidos y escribe los resultados en formato FASTA. Identifica cuál es el marco de lectura correcto (el que produce la proteína más larga sin codones de stop prematuros).

```bash
docker compose run tp python src/ex1_reading_frames.py
# output: output/orfs.fasta
```

### Ejercicio 2a — BLAST

Toma la proteína del frame correcto y ejecuta blastp contra SwissProt. Busca proteínas similares en otras especies para entender qué tan conservada está la huntingtina evolutivamente.

```bash
docker compose run tp python src/ex2_blast.py
# output: output/blast.out
```

### Ejercicio 2b — Interpretación BLAST

| Valor | Significado |
|-------|------------|
| E-value | Probabilidad de que el hit sea por azar. Menor = más confiable |
| Identity % | Porcentaje de aminoácidos idénticos entre las secuencias |
| Bit score | Calidad del alineamiento. Mayor = mejor |
| Query coverage | Qué porción de nuestra secuencia cubre el hit |

### Ejercicio 3 — MSA

Descarga las secuencias de HTT de 3+ especies encontradas en el BLAST y realiza un alineamiento múltiple. Las posiciones conservadas en todas las especies son funcionalmente críticas.

```bash
docker compose run tp python src/ex3_msa.py
# output: output/msa.aln
```
