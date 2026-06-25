---
css: |
  body { font-family: Arial, Helvetica, sans-serif; }
  h1, h2, h3, h4, h5, h6, p, li, td, th, blockquote { font-family: Arial, Helvetica, sans-serif; }
  .page-break { page-break-after: always; }
---

# Trabajo Práctico Parte 1 — Informe de Avance
**Introducción a la Bioinformática — UTN FRBA**  
**Fecha de entrega: 1 de junio de 2026**

---

## 1. Introducción

### Enfermedad elegida: Enfermedad de Huntington

La Enfermedad de Huntington (EH) es una enfermedad hereditaria catalogada en OMIM bajo el número #143100. La elegimos porque tiene una causa genética muy clara y bien documentada, lo que la hace ideal para trabajar con herramientas bioinformáticas. Se caracteriza por un deterioro progresivo que afecta el movimiento, la cognición y el comportamiento, y no tiene cura.

Lo que más nos llamó la atención al investigarla es que no existen portadores sanos: cualquier persona que herede el alelo mutado va a desarrollar la enfermedad en algún momento de su vida. Eso la hace diferente a la mayoría de las enfermedades hereditarias.

**Causa molecular:**la elegimos por ser una enfermedad conocida y hereditaria

### Gen elegido: HTT (Huntingtin)

El gen HTT está en el cromosoma 4 y codifica para la proteína huntingtina. Según la base de datos SwissProt, la proteína tiene alrededor de 3144 aminoácidos. La secuencia de referencia que usamos en este trabajo es el transcripto **NM_002111.8** obtenido de NCBI Nucleotide.

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

**Nota de implementación:** Al investigar las herramientas disponibles se encontró que BioPerl ofrece `Bio::SeqUtils->translate_6frames()`, función que calcula los 6 marcos de lectura incluyendo los frames negativos (-1, -2, -3) sobre el complemento reverso de la secuencia, y devuelve 6 objetos de secuencia directamente utilizables. BioPython cuenta con `Bio.SeqUtils.six_frame_translations()` que realiza el mismo cálculo pero devuelve un string formateado para visualización, no objetos procesables. Por este motivo se implementó la función `get_six_frames()` utilizando las primitivas `seq.reverse_complement()` y `seq.translate()`, reproduciendo el comportamiento de la función de BioPerl.

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

### Conclusión del ejercicio

Lo que más me llamó la atención al ver los resultados fue la diferencia de longitudes: el frame +2 produce 3192 aminoácidos antes del primer stop, y el segundo mejor llega solo a 131. No hacía falta saber de antemano cuál era el correcto — la diferencia es tan grande que se ve sola. Igual, en este ejercicio no determinamos el frame correcto: eso lo confirmamos recién en el ejercicio 2 con BLAST, que es el criterio más sólido.

---

## 4. Ejercicio 2a — BLAST

### Descripción

Se desarrolló el script `src/ex2_blast.py` que:

1. Lee los 6 frames del archivo `output/orfs.fasta` (output del Ejercicio 1)
2. Ejecuta una búsqueda **blastp remota** contra SwissProt del NCBI para cada frame usando `Bio.Blast.NCBIWWW`
3. Guarda todos los resultados en `output/blast.xml` (para el Ejercicio 3) y en `output/blast.out` (reporte legible en texto plano)
4. Muestra un resumen por frame indicando cuántos hits significativos obtuvo cada uno

### Resultados

El BLAST remoto se corrió sobre los 6 frames para identificar el correcto sin conocimiento previo. Solo el frame +2 produjo hits significativos:

| Frame | aa antes del stop | Hits significativos | Mejor E-value |
|-------|-------------------|---------------------|---------------|
| +1 | 131 | 0 | — |
| -1 | 23 | 0 (E-value: 1.40) | — |
| **+2** | **3192** | **5** | **0.0** |
| -2 | 13 | 0 | — |
| +3 | 27 | 0 | — |
| -3 | 5 | 0 | — |

Los 5 hits del frame +2:

| # | Proteína | Organismo | Identity | E-value | Score |
|---|----------|-----------|----------|---------|-------|
| 1 | Huntingtin (HD_HUMAN) | *Homo sapiens* | 99.9% | 0.0 | 16720 |
| 2 | Huntingtin (HD_MOUSE) | *Mus musculus* | 91.2% | 0.0 | 14691 |
| 3 | Huntingtin (HD_RAT) | *Rattus norvegicus* | 90.8% | 0.0 | 14620 |
| 4 | Huntingtin (HD_TAKRU) | *Takifugu rubripes* | 69.7% | 0.0 | 11538 |
| 5 | HD protein homolog | *Dictyostelium discoideum* | 28.8% | 8.22e-19 | 244 |

### BLAST local

Además del BLAST remoto, se ejecutó un BLAST local contra una copia descargada de SwissProt (482.697 secuencias). Se descargó la base de datos en formato FASTA desde el FTP del NCBI y se la formateó con `makeblastdb`. El query utilizado fue el mismo que en el BLAST remoto: la secuencia proteica del frame +2 extraída de `output/orfs.fasta` (3192 aa). El comando ejecutado fue:

```
blastp.exe -db C:\Users\herna\facu\bioinformatica-2026\tp-bioinformatica-2026\data\swissprotdb\swissprot -query C:\Users\herna\facu\bioinformatica-2026\tp-bioinformatica-2026\output\orfs.fasta -out C:\Users\herna\facu\bioinformatica-2026\tp-bioinformatica-2026\output\local-ncbi-blast-repos.txt
```

El resumen de hits significativos producido por BLAST para el frame +2:

```
Query= NM_002111.8_frame_+2

Sequences producing significant alignments:                          (Bits)  Value

P42858.2 RecName: Full=Huntingtin [Homo sapiens]                     6397    0.0
P42859.2 RecName: Full=Huntingtin [Mus musculus]                     5616    0.0
P51111.1 RecName: Full=Huntingtin [Rattus norvegicus]                5589    0.0
P51112.1 RecName: Full=Huntingtin [Takifugu rubripes]                4424    0.0
Q76P24.1 RecName: Full=HD protein homolog [Dictyostelium discoideum]   97.4    3e-18
```

Los resultados son consistentes con el BLAST remoto, confirmando ambos análisis:

| # | Accesión | Proteína | Organismo | Identidad | E-value | Bit score |
|---|----------|----------|-----------|-----------|---------|-----------|
| 1 | P42858.2 | Huntingtin (HD_HUMAN) | *Homo sapiens* | 99% (3142/3144) | 0.0 | 6397 |
| 2 | P42859.2 | Huntingtin (HD_MOUSE) | *Mus musculus* | 91% (2792/3063) | 0.0 | 5616 |
| 3 | P51111.1 | Huntingtin (HD_RAT) | *Rattus norvegicus* | 91% (2781/3063) | 0.0 | 5589 |
| 4 | P51112.1 | Huntingtin (HD_TAKRU) | *Takifugu rubripes* | 70% (2245/3223) | 0.0 | 4424 |
| 5 | Q76P24.1 | HD protein homolog | *Dictyostelium discoideum* | 29% | 3e-18 | 97.4 |

El hit 1 muestra 99% (3142/3144) y no 100% porque la traducción del frame +2 incluye 2 aminoácidos extra en el N-terminal respecto a la secuencia curada de SwissProt (P42858.2, 3142 aa), lo que genera 2 gaps en el alineamiento. Esto es esperable: el frame +2 traduce desde la posición 2 del mRNA crudo, mientras que la secuencia SwissProt representa la proteína procesada y anotada manualmente.


---

## 5. Ejercicio 2b — Interpretación del resultado BLAST

### Significado de los valores estadísticos

**E-value (Expect value):** es la cantidad de hits que esperarías encontrar por azar en una base de datos del tamaño de SwissProt con una puntuación igual o mejor. Un E-value de 0.0 significa que esa probabilidad es tan pequeña que Python la redondea a cero — el hit es real. El hit 5 (Dictyostelium) con E-value de 8.22e-19 también es significativo: hay 1 en 10^19 chances de que esa similitud sea producto del azar. En general, se considera significativo cualquier E-value menor a 0.001.

**Identity %:** porcentaje de aminoácidos idénticos en el alineamiento. Los hits 2, 3 y 4 (ratón, rata y pez globo) presentan alta identidad con la huntingtina humana, lo que indica que el gen HTT está fuertemente conservado en vertebrados.

**Score:** mide la calidad del alineamiento sumando puntos por cada posición — aminoácidos idénticos suman más, similares suman menos, gaps restan. Existen dos variantes: el **raw score** depende de los parámetros de configuración usados y no es comparable entre distintas corridas; el **bit score** es una versión normalizada del raw score que sí es comparable entre cualquier corrida de BLAST. Por eso el BLAST remoto y el local reportan scores distintos para los mismos hits: el remoto usa raw score y el local bit score.

### Interpretación biológica

Los primeros 4 hits corresponden a la huntingtina de otros vertebrados, con identidades superiores al 69%. Esto demuestra que HTT es un gen altamente conservado en vertebrados, lo que implica que cumple funciones celulares fundamentales más allá de la patología asociada en humanos.

El hit más sorprendente es el número 5: *Dictyostelium discoideum*, que es una ameba unicelular. No esperábamos encontrar un homólogo de huntingtina en un organismo tan distinto a los vertebrados. Con 28.8% de identidad todavía tiene regiones similares, lo que sugiere que este gen existe desde hace muchísimo tiempo y probablemente cumple alguna función básica en la célula.

### Conclusión del ejercicio

Correr BLAST sobre los 6 frames sin saber cuál era el correcto resultó ser la forma más directa de identificarlo: solo el frame +2 dio hits reales, todos los demás no encontraron nada. Eso ya responde la pregunta del ejercicio 1.

Lo que no esperaba era la ameba. Que exista un homólogo en *Dictyostelium discoideum* me hizo entender que la huntingtina no es solo una proteína asociada a una enfermedad humana — debe tener alguna función más básica que se conservó a lo largo de la evolución. También fue útil comparar el BLAST remoto con el local: los dos dieron los mismos 5 hits en el mismo orden, lo que me da más confianza en los resultados.

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

Al abrir el archivo de alineamiento lo primero que se ve son los símbolos `*`: indican posiciones donde todas las especies tienen el mismo aminoácido. Hay muchos en la parte central y final de la secuencia, lo que significa que esas regiones se mantuvieron iguales a lo largo de la evolución.

Lo más llamativo está al principio de la secuencia:

```
HD_HUMAN   MATLEKLMKAFESLKSFQQQQQQQQQQQQQQQQQQQPPPPP...
HD_MOUSE   MATLEKLMKAFESLKSFQQQQQQQQPP...
HD_RAT     -------MKAFESLKSFQQQQQQQQP...
HD_TAKRU   MATMEKLMKAFESLKSFQQQQG...
HD_DICDI   ----------------------------------------------------------MD
```

La secuencia humana tiene una cadena mucho más larga de Q repetidas que el resto de las especies. Esa región es exactamente donde ocurre la expansión de repeticiones CAG que causa la enfermedad — cada Q en la proteína corresponde a un CAG en el ADN. Verlo directamente en el alineamiento fue bastante impactante porque conecta todo el análisis con la causa molecular de la enfermedad.

En el resto de la secuencia, ratón y rata son muy similares a la humana. El pez globo tiene más diferencias pero todavía comparte bastantes posiciones. La ameba es la que más gaps tiene y menos `*` muestra, aunque aparecen algunos bloques conservados en distintos lugares de la secuencia.

### Conclusión del ejercicio

Al ver el alineamiento lo primero que se nota es que ratón y rata son casi iguales a la secuencia humana, lo cual tiene sentido porque son mamíferos. El pez globo ya diverge más pero todavía comparte regiones. La ameba tiene la mayor cantidad de gaps y es la más distinta, pero que aparezca en el alineamiento con bloques conservados es lo que más me sorprendió de todo el trabajo. También fue interesante poder ver directamente en el alineamiento la región de glutaminas (Q) que causa la enfermedad — los humanos tienen una cadena mucho más larga que el resto.

![Alineamiento múltiple de huntingtina](../img/alineamiento.png)

*Resultado completo del MSA disponible en: https://www.ebi.ac.uk/jdispatcher/msa/clustalo/summary?jobId=clustalo-I20260526-034619-0482-90639946-p1m*

---

## 7. Conclusión

Este trabajo me resultó más interesante de lo que esperaba. Empezamos con un archivo de texto en formato GenBank y terminamos viendo en el alineamiento la región exacta que causa la enfermedad.

Lo que más me quedó fue la progresión: el ejercicio 1 traduce los 6 frames sin saber cuál es el correcto, el ejercicio 2 lo identifica a través de BLAST, y el ejercicio 3 permite visualizar las diferencias entre especies. Cada paso tiene sentido en función del siguiente.

La sorpresa fue la ameba. Que un organismo tan distinto a nosotros tenga una versión del mismo gen sugiere que la proteína cumple alguna función muy básica — algo que no habría notado sin hacer el análisis. Y ver en el alineamiento que la cadena de Q repetidas es más larga en humanos que en cualquier otra especie cierra el círculo con la causa molecular de la enfermedad.

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

<div class="page-break"></div>

# Trabajo Práctico Parte 2 — Informe

**Introducción a la Bioinformática — UTN FRBA**
**Fecha de entrega: 26 de junio de 2026**

Esta segunda parte continúa el trabajo sobre el gen **HTT** (Huntingtin) y la enfermedad de
Huntington iniciado en la Parte 1. Acá se desarrollan los ejercicios 4 (parser de la salida de
BLAST), 5 (análisis con EMBOSS) y 6 (trabajo con bases de datos biológicas). El ambiente de trabajo
es el mismo de la Parte 1: Python 3.11 con BioPython, EMBOSS y BLAST+ dentro de un contenedor Docker,
con el código versionado en Git.

---

## 1. Ejercicio 4 — Parser de la salida de BLAST

### Descripción

Se desarrolló el script `src/ex4_blast_parser.py`, que toma como entrada el reporte `output/blast.out`
generado en el Ejercicio 2 y un **Pattern** pasado como parámetro por línea de comandos (por defecto
`"Mus musculus"`, como en el ejemplo de la consigna). El script recorre el reporte e identifica los
hits cuya descripción contiene ese Pattern, de forma **case-insensitive** (para que `"mus musculus"`
y `"Mus musculus"` coincidan).

Para distinguir las filas de hit de las líneas decorativas del reporte se usa una expresión regular
que detecta el ACCESSION al inicio de cada fila (formato tipo `P42859.2`). Además, mientras recorre
el archivo, el script va guardando el frame actual (`Frame +2`, etc.) para poder indicar en qué marco
de lectura apareció cada hit.

**Punto extra:** para cada hit que coincide con el Pattern, el script extrae su ACCESSION y descarga
la secuencia completa de la proteína en formato FASTA desde NCBI usando `Bio.Entrez` (el equivalente
en BioPython al módulo `Bio::DB::GenBank` de BioPerl que menciona la consigna), y la escribe a
`output/ex4_hits.fasta`.

```bash
docker compose run tp python src/ex4_blast_parser.py "Mus musculus"
# input:  output/blast.out + Pattern
# output: output/ex4_hits.fasta
```

### Resultados

Ejecutando con el Pattern `"Mus musculus"` sobre el `blast.out` del Ejercicio 2, el script encontró
**1 hit** que coincide:

| Frame | Accession | Descripción |
|-------|-----------|-------------|
| +2 | P42859.2 | Huntingtin — *Mus musculus* (ratón) |

A continuación descargó la secuencia completa de ese hit desde NCBI y la guardó en
`output/ex4_hits.fasta`:

```
>sp|P42859.2|HD_MOUSE RecName: Full=Huntingtin; AltName: Full=Huntington
disease protein homolog; Short=HD protein homolog; ...
```

Si se ejecuta con otro Pattern (por ejemplo `"Rattus"` o `"Takifugu"`) el script devuelve los hits
correspondientes a esa otra especie, lo que confirma que el filtro funciona de forma genérica.

### Conclusión del ejercicio

El ejercicio muestra cómo procesar de forma programática la salida de una herramienta como BLAST en
lugar de leerla a mano. Lo más interesante fue el punto extra: con solo el ACCESSION del hit se puede
recuperar automáticamente la secuencia original completa desde NCBI, encadenando el resultado del
BLAST con una nueva consulta a la base de datos. Esto deja todo listo para análisis posteriores (por
ejemplo, un alineamiento) sin pasos manuales.

---

## 2. Ejercicio 5 — EMBOSS (ORFs y dominios PROSITE)

### Descripción

Se desarrolló el script `src/ex5_emboss.py`, que arma un pequeño pipeline llamando a varios programas
del paquete **EMBOSS** sobre el mRNA de HTT (`data/NM_002111.gb`):

1. **`seqret`** convierte el mRNA del formato GenBank a FASTA de nucleótidos.
2. **`getorf`** calcula los ORFs (marcos de lectura abiertos) del mRNA y los traduce a proteína. Se
   usó `-find 1` (traducir las regiones entre un codón de inicio ATG y el de stop) y `-minsize 300`
   (300 nucleótidos = 100 aa, para descartar ORFs cortos espurios).
3. El script selecciona el **ORF más largo**, que corresponde a la huntingtina.
4. Descarga la base de datos **PROSITE** (`prosite.dat` + `prosite.doc`) desde el FTP de ExPASy y la
   indexa con **`prosextract`**.
5. Finalmente, **`patmatmotifs`** analiza los dominios/motivos funcionales conocidos de la proteína
   contra PROSITE y escribe el resultado en `output/ex5_domains.txt`.

```bash
docker compose run tp python src/ex5_emboss.py
# input:  data/NM_002111.gb
# output: output/ex5_orfs.fasta, output/ex5_longest_orf.fasta, output/ex5_domains.txt
```

### Resultados

`getorf` encontró **8 ORFs** de al menos 100 aa. El más largo es la huntingtina:

```
>NM_002111_1 [146 - 9577] Homo sapiens huntingtin (HTT), transcript variant 2, mRNA
Longitud: 3144 aa
```

Este tamaño (3144 aa) coincide con el de la huntingtina humana conocida y con el frame +2 que ya
habíamos identificado en la Parte 1.

El análisis de dominios con `patmatmotifs` contra PROSITE encontró **4 motivos funcionales**:

| Motivo | Posición (aa) | Qué es |
|--------|---------------|--------|
| AMIDATION | 1532–1535 | Sitio de amidación |
| AMIDATION | 2545–2548 | Sitio de amidación |
| LEUCINE_ZIPPER | 1446–1467 | Patrón "cremallera de leucina" |
| TYR_PHOSPHO_SITE_2 | 2716–2723 | Sitio de fosforilación de tirosina |

### Conclusión del ejercicio

Lo más útil de este ejercicio fue ver cómo se pueden encadenar varias herramientas de EMBOSS en un
mismo pipeline: a partir del mRNA en GenBank se obtienen los ORFs, se elige automáticamente el más
largo y sobre esa proteína se buscan dominios funcionales contra una base de datos externa (PROSITE).
Que `getorf` recupere la huntingtina completa de 3144 aa de forma independiente confirma el resultado
de la Parte 1 por otro camino. Los motivos PROSITE encontrados (sitios de fosforilación, amidación y
una cremallera de leucina) dan una primera pista de las regiones funcionales de la proteína.

---

## 3. Ejercicio 6 — Trabajo con Bases de Datos Biológicas

**Gen:** HTT (Huntingtin) · **Enfermedad:** Enfermedad de Huntington (OMIM #143100)
**NCBI Gene ID:** 3064 · **UniProt:** P42858 · **Ensembl:** ENSG00000197386

*Datos verificados contra las bases en vivo en junio de 2026. Como las bases se actualizan
periódicamente, los conteos podrían variar en el futuro.*

### a) Gen / proteína de interés — NCBI Gene

**Link Entrez Gene:** https://www.ncbi.nlm.nih.gov/gene/3064

Entramos a la página del gen en NCBI Gene y de ahí sacamos los datos principales. El gen se llama
**HTT** (Huntingtin) y tiene el ID **3064**. Está en el **cromosoma 4** (posición 4p16.3) y es un
gen bastante grande: ocupa unas 180 kb y está formado por **67 exones**.

Este gen tiene las instrucciones para fabricar una proteína llamada **huntingtina**. Por lo que
leímos, la huntingtina aparece sobre todo en el cerebro y es importante para que las neuronas
funcionen bien. Cuando el gen tiene un error (una parte de su secuencia, el "CAG", se repite de más),
la proteína sale defectuosa y termina dañando a las neuronas, lo que produce la enfermedad de
Huntington.

![Ficha del gen HTT en NCBI Gene](capturas/ex6_a_ncbi_gene.png)

![Genomic context del gen HTT (ubicación, assembly y 67 exones)](capturas/ex6_a_genomic_context.png)

**Por qué la elegimos:** es una enfermedad **conocida y reconocible** —aparece, por ejemplo, como
diagnóstico en la serie *Dr. House*—, lo que nos motivó a investigarla. Además es **hereditaria**,
con una causa genética clara y bien documentada, lo que la hace ideal para trabajar con herramientas
bioinformáticas.

### b) Genes / proteínas homólogas en otros organismos

Buscamos en dos bases distintas para ver en cuántos otros organismos existe un gen parecido a HTT.

Primero, una aclaración nuestra para entender qué estábamos buscando. Un **ortólogo** es el mismo gen,
pero en otra especie. Lo vimos clarísimo en el Ejercicio 2: cuando comparamos la proteína humana de
HTT contra una base de datos, encontramos proteínas casi iguales en otros animales —en el ratón
coincidía un **91%**, en la rata un **91%** y hasta en un pez un **70%**—. No son idénticas a la
humana, pero son tan parecidas y hacen lo mismo que claramente son "la misma" proteína en cada animal.
Eso es un ortólogo: la versión de un gen en otra especie. Por eso, contar cuántos ortólogos tiene HTT
es ver en cuántos animales aparece este mismo gen; y como aparece en muchísimos, sabemos que es un gen
importante que se mantuvo a lo largo de la evolución.

**HomoloGene (NCBI):** quisimos usarla, pero descubrimos que **NCBI la dio de baja**. Hoy, al entrar,
te redirige automáticamente a **NCBI Datasets / Gene**. Según el anuncio oficial de NCBI esto pasó el
**30 de enero de 2024**, porque la reemplazaron por una herramienta nueva (NCBI Orthologs) que cubre
más genes y más organismos ([anuncio](https://ncbiinsights.ncbi.nlm.nih.gov/2024/01/30/homologene-redirects-ncbi-datasets-gene/)).
O sea que HomoloGene quedó vieja y ya no se puede usar directamente — un buen ejemplo de cómo las
bases biológicas se van actualizando y reemplazando.

![HomoloGene redirige a NCBI Datasets / Gene (discontinuada)](capturas/ex6_b_homologene_retired.png)

**NCBI Orthologs (el reemplazo de HomoloGene):** https://www.ncbi.nlm.nih.gov/datasets/gene/3064/#orthologs
Acá sí encontramos los datos. NCBI lista **849 genes ortólogos** de HTT, es decir, versiones de este
mismo gen repartidas en muchísimas especies distintas.

![Ortólogos de HTT en NCBI Datasets (849 genes)](capturas/ex6_b_ncbi_orthologs.png)

**Ensembl (Comparative Genomics / orthologues):**
https://www.ensembl.org/Homo_sapiens/Gene/Compara_Ortholog?g=ENSG00000197386
Esta base sí está activa. Comparó el gen HTT contra **199 especies** y encontró que tiene un gen
equivalente (ortólogo) en **190 de ellas**:

| Tipo de relación | Qué significa | Especies |
|------------------|---------------|----------|
| 1:1 | un gen humano = un gen en la otra especie | **180** |
| 1:varios | un gen humano = varios en la otra (duplicaciones) | 10 |
| sin ortólogo | no se encontró equivalente | 9 |

![Ortólogos de HTT en Ensembl Compara (199 especies)](capturas/ex6_b_ensembl.png)

**Diferencia entre las bases:** HomoloGene era una base vieja y estática, y NCBI ya la discontinuó.
Sus dos reemplazos —**NCBI Orthologs** y **Ensembl**— están actualizados y comparan automáticamente
contra muchísimas especies. Una aclaración importante: los números no se comparan directamente porque
cada base cuenta cosas distintas. NCBI cuenta **genes** (849, contando varias especies y a veces más
de un gen por especie), mientras que Ensembl cuenta **especies** (199 comparadas, con ortólogo en 190).
Más allá del número exacto, las dos coinciden en lo mismo: HTT aparece en una enorme cantidad de
organismos.

**Qué tan común es el gen y grupos taxonómicos:** el resultado de Ensembl (ortólogo en 190 de 199
especies, casi todos 1:1) muestra que HTT está **muy conservado en todos los vertebrados**. Además,
en el Ejercicio 2 (BLAST) habíamos encontrado un homólogo hasta en *Dictyostelium discoideum* (una
ameba, 28.8% de identidad), lo que indica que el gen es **muy antiguo**: existe desde mucho antes de
los vertebrados y aparece en eucariotas en general, no sólo en animales.

### c) Transcriptos y splicing alternativo

Antes que nada, qué es esto en simple: un **gen** es como una receta. A veces, de la misma receta, la
célula puede armar **versiones distintas** del plato (cambiando o salteando algunos pasos). Cada
versión es un **transcripto**, y esa "edición" de la receta se llama **splicing alternativo**. Acá lo
que hicimos fue ver cuántas versiones distintas de HTT figuran en cada base.

**NCBI (RefSeq):** https://www.ncbi.nlm.nih.gov/gene/3064 (sección *RefSeq transcripts*)
En NCBI encontramos solo **2 transcriptos**: `NM_001388492.1` y `NM_002111.8` (este último es el que
usamos en el TP). Son pocos porque NCBI los **revisa a mano**, así que solo deja los que están bien
confirmados.

![Transcriptos RefSeq de HTT en NCBI (2)](capturas/ex6_c_ncbi.png)

**Ensembl:** https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000197386
Ensembl, en cambio, lista **24 transcriptos** para HTT. Muchos más que NCBI, porque esta base los
detecta de forma automática e incluye también versiones más raras o poco usadas.

![Transcriptos de HTT en Ensembl (24)](capturas/ex6_c_ensembl.png)

| Base | Nº de transcriptos | Cómo los arma |
|------|--------------------|---------------|
| NCBI RefSeq | **2** | Revisados a mano, alta confianza |
| Ensembl | **24** | Automático, capta más variantes |

**¿Por qué dan números tan distintos?** Como nos llamó la atención la diferencia, lo buscamos y
encontramos la explicación en un foro de bioinformática ([Biostars](https://www.biostars.org/p/72845/)):
RefSeq (NCBI) es una colección **curada** de transcriptos —menos, pero bien confirmados—, mientras que
Ensembl es más **inclusiva** e incorpora muchas variantes, incluso algunas con poco respaldo. O sea:
NCBI prioriza estar seguro de lo que muestra, y Ensembl prioriza mostrar todo lo posible.

**¿Cuáles se expresan y cuál base es más precisa?** La versión principal es la que fabrica la
huntingtina completa, que es la forma importante. Muchas de las 24 de Ensembl son versiones cortas o
poco frecuentes. Por eso **las dos bases sirven para cosas distintas**: NCBI es más confiable (pocas
pero seguras, ideal para uso médico), y Ensembl es más completa (muestra toda la variedad posible,
aunque algunas sean dudosas). Es la misma diferencia que vimos en el punto b).

### d) Interacciones proteína–proteína

Qué es esto en simple: las proteínas no trabajan solas, se "enganchan" con otras para hacer su
tarea. Acá miramos con cuántas otras proteínas se relaciona la huntingtina, comparando dos bases:
**NCBI Gene** y **UniProt**.

- **NCBI Gene** (https://www.ncbi.nlm.nih.gov/gene/3064, sección *Interactions*): nos dio una tabla
  con **560 interacciones**.
- **UniProt** (https://www.uniprot.org/uniprotkb/P42858/entry#interaction): nos dio una tabla con
  **860 interacciones**.

![Interacciones de HTT en NCBI Gene (560)](capturas/ex6_d_ncbi.png)

![Interacciones de HTT en UniProt (860)](capturas/ex6_d_uniprot.png)

O sea, la huntingtina se relaciona con **muchísimas** otras proteínas. Algunas de las más conocidas:

| Proteína | Para qué |
|----------|----------|
| **HAP1** | Transporte de cargas dentro de la neurona |
| **HIP1** | Ayuda a meter cosas dentro de la célula (endocitosis) |
| **HIP14** | Modifica químicamente a la huntingtina |
| **DCTN1** | Parte del "motor" que mueve cargas dentro de la célula |

**¿Hay un patrón?** Sí: la mayoría de estas proteínas tienen que ver con **mover cosas de un lado a
otro dentro de la neurona** (transporte interno). Eso encaja con la función principal de la
huntingtina, que es justamente ayudar en ese transporte.

**¿Por qué los números no coinciden (560 vs 860)?** Otra vez, como en b) y c): cada base junta los
datos de fuentes distintas y con criterios distintos, así que una lista más interacciones que la otra.
UniProt llegó a un número más alto (860) que NCBI (560), pero las dos coinciden en lo importante: HTT
es una proteína **muy conectada**, lo que muestra que cumple un rol central en la neurona.

### e) Gene Ontology (GO): componente, proceso y función

Qué es esto en simple: **Gene Ontology (GO)** es un "etiquetado estándar" que describe cada proteína
respondiendo tres preguntas: **dónde está**, **en qué procesos participa** y **qué hace**. Sacamos
estos datos de la página de UniProt (P42858). El listado completo está en las capturas; acá lo
resumimos agrupado por tema.

**Dónde está (Componente celular):** sobre todo en el **citoplasma** (el interior de la célula) y en
sus **vesículas** (las "cajitas" de transporte), también en el **endosoma** y en el **núcleo**.

![Ubicación de la huntingtina en la célula (UniProt)](capturas/ex6_e_subcellullar_location.png)

**En qué procesos participa (Proceso biológico):** los ~20 términos se agrupan en pocos temas:
- **Desarrollo del sistema nervioso** (neurogénesis, desarrollo del cerebro).
- **Transporte dentro de la célula** (mover vesículas a lo largo de los microtúbulos, transporte en
  las sinapsis, organización del Golgi).
- **Limpieza/reciclado celular** (autofagia: aggrephagy, lipophagy, mitophagy).
- **Muerte celular programada** (apoptosis) y **señalización** (calcio, cascada CAMKK-AMPK).

![Procesos biológicos (GO) de la huntingtina en UniProt](capturas/ex6_e_uniprot.png)

**Qué hace (Función molecular):** principalmente **se une a otras proteínas**. Y casi todas esas
proteínas son del **sistema de transporte** de la célula: se une a *beta-tubulina*, *dinactina* y
*dineína* (las piezas de los "rieles y motores" que mueven cargas dentro de la neurona), además de
unirse a otras proteínas reguladoras (p53, kinasas, proteínas de choque térmico).

![Función molecular (GO) de la huntingtina en UniProt](capturas/ex6_e_mollecular_function.png)

**En resumen:** la huntingtina vive en el **citoplasma y las vesículas**, trabaja sobre todo en el
**transporte interno de la neurona y el desarrollo del cerebro**, y lo hace **uniéndose a las
proteínas del sistema de transporte**. Las tres categorías de GO cuentan la misma historia que ya
veníamos viendo en los puntos anteriores.

### f) Vías metabólicas / pathways

Qué es esto en simple: un **pathway** es una "cadena de pasos" donde varias proteínas trabajan juntas
para lograr algo en la célula, como una línea de producción. Acá vimos en qué cadenas participa HTT.

**KEGG:** https://www.genome.jp/dbget-bin/www_bget?hsa:3064
KEGG ubica a HTT en dos pathways, los dos relacionados con la enfermedad:
- **hsa05016 — Huntington disease** (la vía propia de la enfermedad)
- **hsa05022 — Pathways of neurodegeneration** (vía general de enfermedades neurodegenerativas)

Lo interesante es que KEGG además muestra **cómo la huntingtina dañada (mutada) rompe muchas cosas
distintas** dentro de la neurona. Algunos ejemplos que lista:
- frena el **transporte de cargas** dentro de la neurona (axonal),
- altera la **lectura de genes** (transcripción: p53, CREB, REST),
- daña las **mitocondrias** (la "central de energía" de la célula),
- desregula la **autofagia** y la **eliminación de basura** celular (proteasoma),
- dispara la **muerte de la neurona** (apoptosis).

Es decir, un solo gen roto termina afectando un montón de procesos a la vez — por eso la enfermedad
es tan grave. KEGG también muestra que ya hay **medicamentos** apuntando a este gen (Tominersen,
Votoplam), un dato interesante.

![Pathways de HTT en KEGG (hsa05016)](capturas/ex6_f_kegg.png)

**Patrón:** todos los pathways apuntan a lo mismo que veníamos viendo en GO e interacciones —
**transporte dentro de la neurona, autofagia y neurodegeneración**.

### g) Variante genética (dbSNP / ClinVar)

Qué es esto en simple: una **variante** es un "cambio" en el ADN respecto de lo normal. Acá buscamos
la variante que causa la enfermedad de Huntington y qué se sabe de ella.

**La variante:** entramos a **ClinVar** (https://www.ncbi.nlm.nih.gov/clinvar/?term=HTT%5Bgene%5D) y,
filtrando por *Pathogenic*, encontramos **132 variantes** que causan enfermedad en este gen. La
principal —la que causa Huntington— es la **expansión del "CAG"**: en vez de un cambio de una sola
letra, lo que pasa es que un pedacito del ADN (las letras CAG) **se repite de más**. En ClinVar es la
entrada **Variation ID 409**, clasificada como **Pathogenic** para la enfermedad de Huntington, con el
nivel de confianza más alto que da la base (*practice guideline*).

![Variante de la expansión CAG en ClinVar (Variation ID 409)](capturas/ex6_g_clinvar.png)

Lo interesante de esta variante es que **cuántas veces se repite el CAG decide si te enfermás y qué
tan grave**:

| Veces que se repite el CAG | Qué pasa |
|----------------------------|----------|
| 35 o menos | Normal (no te enfermás) |
| 36 a 39 | Zona de riesgo (puede o no aparecer) |
| 40 o más | **Aparece la enfermedad de Huntington** |
| 60 o más | Forma juvenil (empieza en la infancia/adolescencia) |

**¿A quiénes afecta más?** Para esto usamos **MedlinePlus** (la base de divulgación de NCBI que
sugiere la consigna). En su sección de frecuencia dice textualmente:

> La enfermedad de Huntington afecta a un estimado de 3 a 7 de cada 100.000 personas de ascendencia
> europea. El trastorno parece ser menos común en algunas otras poblaciones, incluidas las personas
> de ascendencia japonesa, china y africana.

Es decir, la enfermedad es claramente **más frecuente en personas de ascendencia europea** y menos
común en personas de ascendencia japonesa, china y africana.

![Frecuencia de la enfermedad de Huntington en MedlinePlus](capturas/ex6_g_medlineplus.png)

### Fuentes del Ejercicio 6

- [NCBI Gene 3064 (HTT)](https://www.ncbi.nlm.nih.gov/gene/3064) · [OMIM #143100](https://omim.org/entry/143100) · [OMIM *613004 (HTT)](https://omim.org/entry/613004)
- [Ensembl ENSG00000197386](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000197386)
- [UniProt P42858](https://www.uniprot.org/uniprotkb/P42858/entry) · [GeneCards HTT](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTT)
- [KEGG hsa05016 — Huntington disease](https://www.genome.jp/dbget-bin/www_bget?hsa:3064) · [ClinVar Variation ID 409](https://www.ncbi.nlm.nih.gov/clinvar/variation/409/)
- [MedlinePlus — Huntington disease](https://medlineplus.gov/genetics/condition/huntington-disease/) · [Biostars — RefSeq vs Ensembl](https://www.biostars.org/p/72845/)
- [NCBI Insights — HomoloGene redirige a NCBI Datasets/Gene (2024)](https://ncbiinsights.ncbi.nlm.nih.gov/2024/01/30/homologene-redirects-ncbi-datasets-gene/)
- [Wikipedia — Huntingtin](https://en.wikipedia.org/wiki/Huntingtin)

---

## Anexo — Cómo ejecutar los scripts

**Requisitos:** Docker Desktop instalado y corriendo.

```bash
# Entrar al contenedor
docker compose run tp bash

# Ejercicios de la Parte 2
python src/ex4_blast_parser.py "Mus musculus"   # Genera output/ex4_hits.fasta
python src/ex5_emboss.py                         # Genera output/ex5_*.fasta y ex5_domains.txt
```
