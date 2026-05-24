# Trabajo Práctico Parte 1 — Informe de Avance
**Introducción a la Bioinformática — UTN FRBA**  
**Fecha de entrega: 1 de junio de 2026**

---

## 1. Introducción

### Enfermedad elegida: Enfermedad de Huntington

La Enfermedad de Huntington (EH) es una enfermedad neurodegenerativa hereditaria autosómica dominante, catalogada en OMIM bajo el número #143100. Se caracteriza por el deterioro progresivo de las funciones motoras, cognitivas y psiquiátricas, llevando inevitablemente a la muerte entre 10 y 30 años después de la aparición de los síntomas.

A diferencia de la mayoría de las enfermedades, la EH no tiene portadores sanos: cualquier persona que herede el alelo mutado del gen HTT desarrollará la enfermedad a lo largo de su vida. Esto la convierte en un caso de estudio central en genética y bioinformática.

**Causa molecular:** la enfermedad es causada por una expansión anormal de repeticiones del triplete CAG (citosina-adenina-guanina) en el exón 1 del gen HTT. En individuos sanos, el número de repeticiones es menor a 36. Cuando supera ese umbral, la proteína huntingtina resultante se vuelve tóxica para las neuronas. A mayor número de repeticiones, mayor severidad y menor edad de inicio de la enfermedad.

### Gen elegido: HTT (Huntingtin)

El gen HTT, ubicado en el cromosoma 4 (4p16.3), codifica para la proteína huntingtina, una proteína de gran tamaño (~3144 aminoácidos) cuya función normal incluye el transporte intracelular, la supervivencia neuronal y la regulación de la transcripción génica. La secuencia de referencia utilizada en este trabajo es el transcripto **NM_002111.8** obtenido de la base de datos NCBI Nucleotide.

---

## 2. Ambiente de trabajo

El trabajo fue desarrollado utilizando:

- **Lenguaje:** Python 3.11 con la librería BioPython
- **Ambiente:** contenedor Docker sobre Windows, garantizando reproducibilidad
- **Control de versiones:** repositorio Git en GitHub
- **BLAST:** ejecutado de forma remota contra la base de datos SwissProt del NCBI

Para reproducir el ambiente y ejecutar los scripts:

```bash
docker compose run tp bash
python src/ex1_reading_frames.py
python src/ex2_blast.py
python src/ex3_msa.py
```

---

## 3. Ejercicio 1 — Procesamiento de Secuencias

### Descripción

Se desarrolló el script `src/ex1_reading_frames.py` que realiza las siguientes operaciones:

1. Lee el archivo `data/NM_002111.gb` en formato GenBank usando `Bio.SeqIO`
2. Extrae la secuencia de nucleótidos del mRNA maduro de HTT (13.498 bp)
3. Calcula los **6 marcos de lectura posibles**: +1, +2, +3 sobre la cadena directa y -1, -2, -3 sobre el complemento reverso
4. Traduce cada marco a su secuencia de aminoácidos usando la tabla de código genético estándar
5. Escribe los 6 resultados en el archivo `output/orfs.fasta` en formato FASTA
6. Identifica el marco de lectura correcto como aquel que produce la proteína más larga antes del primer codón de stop

### Resultados

| Frame | Aminoácidos antes del primer stop | Observación |
|-------|----------------------------------|-------------|
| +1 | 131 | Stop prematuro |
| **+2** | **3192** | **Frame correcto** |
| +3 | 27 | Stop prematuro |
| -1 | 23 | Stop prematuro |
| -2 | 13 | Stop prematuro |
| -3 | 5 | Stop prematuro |

**El frame correcto es el +2**, con 3192 aminoácidos antes del primer codón de stop, coherente con el tamaño conocido de la huntingtina humana (~3144 aa). Los otros 5 marcos producen proteínas truncadas de menos de 132 aminoácidos, lo que indica que no corresponden al marco de lectura real.

---

## 4. Ejercicio 2a — BLAST

### Descripción

Se desarrolló el script `src/ex2_blast.py` que:

1. Extrae la secuencia proteica del frame +2 del archivo `output/orfs.fasta`
2. Ejecuta una búsqueda **blastp remota** contra la base de datos SwissProt del NCBI usando `Bio.Blast.NCBIWWW`
3. Guarda el resultado completo en `output/blast.out` en formato XML
4. Parsea el resultado y muestra los 10 mejores hits con E-value, porcentaje de identidad y score

### Resultados

| # | Proteína | Organismo | Identity | E-value | Score |
|---|----------|-----------|----------|---------|-------|
| 1 | Huntingtin (HD_HUMAN) | *Homo sapiens* | 99.9% | 0.0 | 16720 |
| 2 | Huntingtin (HD_MOUSE) | *Mus musculus* | 91.2% | 0.0 | 14691 |
| 3 | Huntingtin (HD_RAT) | *Rattus norvegicus* | 90.8% | 0.0 | 14620 |
| 4 | Huntingtin (HD_TAKRU) | *Takifugu rubripes* | 69.7% | 0.0 | 11538 |
| 5 | HD protein homolog | *Dictyostelium discoideum* | 28.8% | 8.22e-19 | 244 |

---

## 5. Ejercicio 2b — Interpretación del resultado BLAST

### Significado de los valores estadísticos

**E-value (Expect value):** representa la probabilidad de encontrar un hit con esa puntuación por azar en una base de datos del tamaño de SwissProt. Un E-value de 0.0 significa que la probabilidad es tan pequeña que Python la redondea a cero — el hit es real con absoluta certeza. El hit 5 con E-value de 8.22e-19 también es estadísticamente significativo: implica que hay 1 en 10^19 posibilidades de que esa similitud sea producto del azar.

**Identity %:** porcentaje de aminoácidos idénticos en el alineamiento. Los hits 2, 3 y 4 (ratón, rata y pez globo) presentan alta identidad con la huntingtina humana, lo que indica que el gen HTT está fuertemente conservado en vertebrados.

**Bit score:** medida normalizada de la calidad del alineamiento, independiente del tamaño de la base de datos. A mayor bit score, mejor es el alineamiento.

### Interpretación biológica

Los primeros 4 hits corresponden a la huntingtina de otros vertebrados, con identidades superiores al 69%. Esto demuestra que HTT es un gen altamente conservado en vertebrados, lo que implica que cumple funciones celulares fundamentales más allá de la patología asociada en humanos.

El hit más revelador es el número 5: *Dictyostelium discoideum*, un organismo unicelular (ameba) que divergió de los animales hace aproximadamente 1.000 millones de años. La presencia de un homólogo de huntingtina con 28.8% de identidad en este organismo primitivo indica que el gen HTT tiene una antigüedad evolutiva extraordinaria y cumple funciones celulares básicas conservadas desde los inicios de la vida eucariótica.

---

## 6. Ejercicio 3 — Multiple Sequence Alignment (MSA)

### Descripción

Se desarrolló el script `src/ex3_msa.py` que descarga automáticamente las secuencias proteicas de los 5 hits encontrados en el BLAST desde NCBI usando `Bio.Entrez`. Las secuencias fueron guardadas en `output/msa_input.fasta` y el alineamiento múltiple fue realizado con **Clustal Omega** (versión online, EMBL-EBI). El resultado se guardó en `output/msa.aln`.

Las especies incluidas en el MSA son:

| Código | Especie | Identidad con HTT humana |
|--------|---------|--------------------------|
| HD_HUMAN | *Homo sapiens* | referencia |
| HD_MOUSE | *Mus musculus* | 91.2% |
| HD_RAT | *Rattus norvegicus* | 90.8% |
| HD_TAKRU | *Takifugu rubripes* (pez globo) | 69.7% |
| HD_DICDI | *Dictyostelium discoideum* (ameba) | 28.8% |

### Interpretación del alineamiento

**Conservación general:** los símbolos `*` en el alineamiento indican posiciones idénticas en todas las especies. Las regiones con alta densidad de `*` corresponden a dominios funcionalmente críticos de la proteína huntingtina.

**Observación clave — región poliglutamina (polyQ):**

```
HD_HUMAN   MATLEKLMKAFESLKSFQQQQQQQQQQQQQQQQQQQPPPPP...
HD_MOUSE   MATLEKLMKAFESLKSFQQQQQQQQPP...
HD_RAT     -------MKAFESLKSFQQQQQQQQP...
HD_TAKRU   MATMEKLMKAFESLKSFQQQQG...
HD_DICDI   ----------------------------------------------------------MD
```

La región de repeticiones de glutamina (Q) en el extremo N-terminal de la proteína es directamente observable en el alineamiento. Los humanos presentan una región polyQ significativamente más extensa que el resto de los vertebrados, y es en esta misma región donde la expansión patológica de repeticiones CAG causa la Enfermedad de Huntington. Esta observación conecta directamente el análisis bioinformático con el mecanismo molecular de la enfermedad.

**Conservación entre vertebrados:** las secuencias de humano, ratón, rata y pez globo muestran amplias regiones de alta conservación (`*` y `:`), especialmente en la mitad C-terminal de la proteína, donde se encuentran los dominios HEAT repeats responsables de las interacciones proteína-proteína.

**Divergencia con Dictyostelium:** el alineamiento con la ameba muestra menor conservación general, aunque persisten bloques conservados distribuidos a lo largo de toda la secuencia, corroborando el origen evolutivo antiguo del gen.

---

## 7. Conclusión

El análisis bioinformático del gen HTT permitió:

1. Identificar el marco de lectura correcto (+2) de entre 6 posibles, produciendo una proteína de 3192 aminoácidos coherente con la huntingtina humana conocida.
2. Confirmar mediante BLAST que la huntingtina está conservada en vertebrados (69-91% de identidad) y que tiene homólogos incluso en organismos unicelulares como *Dictyostelium discoideum*, lo que sugiere una función celular fundamental de más de 1.000 millones de años de antigüedad.
3. Visualizar directamente en el MSA la región polyQ que causa la Enfermedad de Huntington, observando cómo esta región es más extensa en humanos que en otras especies.

Estos resultados ilustran cómo las herramientas bioinformáticas permiten conectar la información genómica con el mecanismo molecular de una enfermedad hereditaria.

---

## Anexo — Cómo ejecutar los scripts

**Requisitos:** Docker Desktop instalado y corriendo.

```bash
# Clonar el repositorio
git clone https://github.com/[usuario]/tp-bioinformatica-2026

# Entrar al contenedor
docker compose run tp bash

# Ejecutar ejercicios
python src/ex1_reading_frames.py   # Genera output/orfs.fasta
python src/ex2_blast.py            # Genera output/blast.out
python src/ex3_msa.py              # Descarga secuencias para MSA
```

Los archivos de input (`data/NM_002111.gb`) y output (`output/`) se encuentran en el repositorio.
