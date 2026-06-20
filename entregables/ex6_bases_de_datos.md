# Ejercicio 6 — Trabajo con Bases de Datos Biológicas
**Gen:** HTT (Huntingtin) · **Enfermedad:** Enfermedad de Huntington (OMIM #143100)
**NCBI Gene ID:** 3064 · **UniProt:** P42858 · **Ensembl:** ENSG00000197386

> [!warning] Verificar en vivo antes de entregar
> Este documento tiene los datos correctos a junio 2026, pero **las bases se actualizan**.
> Antes de cerrar el TP, entrá a cada link, confirmá los números (conteos de homólogos,
> transcriptos, interacciones) y **sacá las capturas de pantalla** para la presentación (Ej7).
> Los puntos marcados con 🔎 son los que conviene capturar.

---

## a) Gen / proteína de interés — NCBI Gene

🔎 **Link Entrez Gene:** https://www.ncbi.nlm.nih.gov/gene/3064

El gen **HTT** está en el cromosoma **4p16.3** (GRCh38: chr4:3,074,510–3,243,960), abarca ~180 kb
y tiene **67 exones**. Codifica la **huntingtina**, una proteína grande (~3144 aa) de expresión
ubicua, con mayor expresión en cerebro adulto y fetal.

**Qué hace la proteína:** la huntingtina es una proteína de andamiaje (*scaffold*) involucrada en:
- **Transporte axonal** de vesículas y organelas a lo largo de microtúbulos (vía su asociación
  con HAP1, dineína/dinactina y kinesina).
- **Transporte de BDNF** (factor neurotrófico) hacia las neuronas estriatales — clave para su
  supervivencia.
- Regulación de **autofagia**, **transcripción** y **apoptosis**.
- Desarrollo del sistema nervioso (su knockout es letal embrionario en ratón).

**Por qué la elegimos:** la enfermedad de Huntington tiene una causa molecular única y muy clara
—la expansión del triplete **CAG** en el exón 1, que produce un tracto de poliglutaminas (polyQ)
tóxico— y no hay portadores sanos: quien hereda el alelo expandido desarrolla la enfermedad. Eso
la vuelve un caso ideal para recorrer todo el pipeline bioinformático, desde la secuencia hasta
la variante patogénica.

---

## b) Genes / proteínas homólogas en otros organismos

**HomoloGene (NCBI):** 🔎 buscar "HTT" en https://www.ncbi.nlm.nih.gov/homologene
La huntingtina aparece en un grupo de homólogos que abarca el rango de **vertebrados modelo**
(humano, chimpancé, perro, vaca, ratón, rata, gallina, *Xenopus*, pez cebra). Son ~9–11 organismos.

> **Importante para el ejercicio:** HomoloGene es una base **legacy** — NCBI dejó de actualizarla
> (último build 68, 2014). Por eso su lista es corta y fija. Conviene mencionarlo como diferencia
> metodológica clave.

**Ensembl (Comparative Genomics / orthologues):** 🔎
https://www.ensembl.org/Homo_sapiens/Gene/Compara_Ortholog?g=ENSG00000197386
Ensembl Compara reporta ortólogos de HTT en **>200 especies** (todo el linaje de vertebrados y
más allá), calculados con árboles filogenéticos de genes.

**Diferencia entre ambas:** HomoloGene da un puñado de organismos modelo curados de forma estática;
Ensembl da cientos de ortólogos calculados automáticamente y actualizados. Ensembl es mucho más
exhaustivo; HomoloGene quedó congelado.

**Qué tan común es el gen y grupos taxonómicos:** HTT está **fuertemente conservado en todos los
vertebrados**. De hecho, en el Ejercicio 2 (BLAST) encontramos un homólogo hasta en *Dictyostelium
discoideum* (una ameba, 28.8% identidad), lo que indica que el gen es muy antiguo y existe desde
mucho antes de los vertebrados — pertenece a eucariotas en general, no sólo a animales.

---

## c) Transcriptos y splicing alternativo

**NCBI (RefSeq):** 🔎 https://www.ncbi.nlm.nih.gov/gene/3064 (sección *mRNA and Protein(s)* / RefSeq)
RefSeq mantiene pocas variantes **curadas manualmente**. La que usamos en el TP es
**NM_002111.8** (*transcript variant 2*).

**Ensembl:** 🔎 https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000197386
Ensembl lista **23 transcriptos** para HTT (entre codificantes, con retención de intrón y
procesados), bastantes más que RefSeq.

| Base | Nº de transcriptos | Criterio |
|------|--------------------|----------|
| NCBI RefSeq | pocos (curados) | Revisión manual, alta confianza |
| Ensembl | 23 | Automático + evidencia, más exhaustivo |

**¿Cuáles se expresan / funciones alternativas?** El transcripto principal codifica la huntingtina
completa de ~3144 aa, que es la forma funcional dominante. Varias de las isoformas de Ensembl son
fragmentos cortos o transcriptos con retención de intrón de expresión baja o incierta.

**¿Por qué difieren y cuál es más precisa?** Difieren por la metodología: **RefSeq es manual y
conservador** (menos transcriptos pero de alta confianza), mientras que **Ensembl es automático y
amplio** (capta más variantes, incluyendo algunas dudosas). Para una referencia clínica, RefSeq es
más preciso/confiable; para explorar la diversidad de splicing, Ensembl es más completo.

---

## d) Interacciones proteína–proteína

**NCBI Gene** (sección *Interactions*, datos de BioGRID) 🔎 y **UniProt P42858**
(sección *Interaction*) 🔎 https://www.uniprot.org/uniprotkb/P42858/entry#interaction

Interactores destacados de la huntingtina:

| Interactor | Rol | Aparece en |
|------------|-----|------------|
| **HAP1** | Transporte axonal; puente con dineína/dinactina y kinesina | NCBI + UniProt |
| **HIP1** (huntingtin-interacting protein 1) | Endocitosis mediada por clatrina | NCBI + UniProt |
| **HIP14 / ZDHHC17** | Palmitoilación de HTT (vía dominio ankyrin) | NCBI + UniProt |
| **HIP14L / ZDHHC13** | Palmitoilación de HTT | NCBI |
| **HAP40 / F8A1** | Estabiliza la conformación de HTT | UniProt |
| **DCTN1** (dinactina p150Glued) | Motor de transporte retrógrado | NCBI |

**¿Hay un patrón?** Sí: la mayoría de las interacciones giran en torno al **transporte intracelular
sobre microtúbulos** (HAP1, dinactina, kinesina) y al **tráfico de vesículas/endocitosis** (HIP1,
HIP14). Esto refuerza el rol de la huntingtina como proteína de andamiaje del transporte neuronal.

**Comparación de tablas:** UniProt tiende a listar interacciones **binarias curadas** (con
referencia a IntAct), mientras que NCBI Gene agrega datos de **BioGRID** (incluye ensayos de alto
rendimiento). 🔎 Compará ambas listas y anotá: probablemente haya interactores presentes sólo en
una de las dos (p. ej. HIP14L/ZDHHC13 suele verse en NCBI/BioGRID; HAP40/F8A1 destaca en UniProt).

---

## e) Gene Ontology (GO): componente, proceso y función

Fuentes: NCBI Gene, UniProt P42858 (sección *Function / GO*) y AmiGO 🔎
http://amigo.geneontology.org/amigo/gene_product/UniProtKB:P42858

| Aspecto GO | Términos principales |
|------------|----------------------|
| **Componente celular** | citoplasma, núcleo, citosol, retículo endoplásmico, aparato de Golgi, mitocondria, axón/dendrita, membrana de vesícula, autofagosoma |
| **Proceso biológico** | desarrollo del sistema nervioso, **transporte axonal** (anterógrado y retrógrado), transporte de vesículas sobre microtúbulos, **transporte de BDNF**, regulación de la autofagia, regulación de la apoptosis |
| **Función molecular** | *protein binding*, unión a **dinactina**, unión a **kinesina**, unión a microtúbulos, *identical protein binding* (dimerización) |

En resumen: la huntingtina es **citoplasmática/asociada a organelas y vesículas**, participa en el
**transporte intracelular y el desarrollo neuronal**, y su función molecular es de **proteína de
andamiaje** que une motores y cargos del transporte sobre microtúbulos.

---

## f) Vías metabólicas / pathways

**KEGG:** 🔎 https://www.genome.jp/dbget-bin/www_bget?hsa:3064
HTT es el gen central de la vía **hsa05016 "Huntington disease"**, y también figura en
**hsa05022 "Pathways of neurodegeneration - multiple diseases"**. La vía KEGG de Huntington
muestra cómo la huntingtina mutada afecta el transporte axonal, la transcripción (CBP, SP1, p53),
la función mitocondrial y la maquinaria de degradación de proteínas.

**Reactome:** 🔎 https://reactome.org/content/detail/R-HSA-... (buscar "HTT" / "P42858")
Reactome ubica a HTT en eventos de **transporte de vesículas mediado por microtúbulos**, **autofagia**
y vías de respuesta al estrés. 🔎 Capturá la lista de *pathways* exacta que devuelve Reactome para
P42858.

Patrón: todos los pathways apuntan a lo mismo que vimos en GO e interacciones — **transporte
neuronal, autofagia y neurodegeneración**.

---

## g) Variante genética (dbSNP / ClinVar)

🔎 ClinVar (filtrar por HTT + Pathogenic): https://www.ncbi.nlm.nih.gov/clinvar/?term=HTT%5Bgene%5D
🔎 dbSNP: https://www.ncbi.nlm.nih.gov/snp/?term=HTT

**La variante:** la causa de la enfermedad de Huntington es la **expansión del repeat CAG en el
exón 1 de HTT**. No es un SNP bialélico clásico sino una **expansión de microsatélite**
(repetición de trinucleótido), catalogada en ClinVar como variante **patogénica**. Umbrales:

| Nº repeticiones CAG | Fenotipo |
|---------------------|----------|
| ≤ 35 | Normal (no afectado) |
| 36–39 | Penetrancia incompleta |
| ≥ 40 | **Penetrancia completa — HD** |
| ≥ 60 | Forma juvenil (inicio en infancia/adolescencia) |

**Frecuencia en la población y grupo étnico más afectado:** como es una expansión, no tiene una
"frecuencia alélica" tipo SNP; se mide por **prevalencia** y por el **haplogrupo** sobre el que
ocurre la expansión:

- Mayor prevalencia en poblaciones de **ascendencia europea**: ~5–15 / 100.000, sobre el
  **haplogrupo A (variantes A1/A2)**, predisponente a la expansión.
- Mucho menor en **Asia oriental** (China, Japón: ~1–2 / 100.000 o menos), donde las expansiones
  ocurren sobre el **haplogrupo C** y faltan las variantes A1/A2 de alto riesgo.
- Baja también en **África subsahariana** (~1–2 / 100.000), sobre el **haplogrupo B**.

→ **El grupo más afectado es el de ascendencia europea**, por un efecto fundador ligado al
haplogrupo A. (Caso emblemático: el conglomerado de HD en la región del Lago de Maracaibo,
Venezuela, que permitió clonar el gen en 1993.)

> 🔎 Para la entrega: capturá la entrada concreta de ClinVar de la expansión CAG patogénica de HTT
> (anotá su *Variation ID*) y, si querés mostrar un rs concreto, podés citar un **modificador**
> documentado como **rs79727797** (repeat en *TCERG1*, retrasa la edad de inicio ~2.3 años por alelo).

---

## Fuentes

- [NCBI Gene 3064 (HTT)](https://www.ncbi.nlm.nih.gov/gene/3064) · [OMIM #143100](https://omim.org/entry/143100) · [OMIM *613004 (HTT)](https://omim.org/entry/613004)
- [Ensembl ENSG00000197386](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000197386)
- [UniProt P42858](https://www.uniprot.org/uniprotkb/P42858/entry) · [GeneCards HTT](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTT)
- [KEGG hsa05016 — Huntington disease](https://www.genome.jp/pathway/hsa05016) · [AmiGO P42858](http://amigo.geneontology.org/amigo/gene_product/UniProtKB:P42858)
- Warby et al. *HTT haplotypes contribute to differences in HD prevalence between Europe and East Asia*, Eur J Hum Genet (2011) — [nature.com/articles/ejhg2010229](https://www.nature.com/articles/ejhg2010229)
- [Wikipedia — Huntingtin](https://en.wikipedia.org/wiki/Huntingtin)
