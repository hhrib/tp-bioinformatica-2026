**Fecha de entrega del TP completo (Partes 1 y 2 en PDF): Viernes 26 de Junio 2026** 

**Introducción a la Bioinformática** 

**Trabajo Práctico – Parte 2** 

**Ejercicio 4 – BLAST OUTPUT.** Escribir un script para analizar (*parsear*) un reporte de salida de blast que  identifique los hits que en su descripción aparezca un Pattern determinado que le damos como parámetro  de entrada. El pattern puede ser una palabra. Punto extra: pueden a su vez parsear cuál es el ACCESSION  del hit identificado (donde hay una coincidencia del Pattern) y con el módulo Bio::DB::GenBank obtener la  secuencia completa del hit en formato FASTA y escribirla a un archivo, es decir, levantar las secuencias originales completas de los hits seleccionados. 

− **Input: Reporte Blast (blast.out del ej. 2\) y un Pattern (por ej. “Mus Musculus”).** − **Output: Lista de los hits que coincidan con el pattern (por ej. solo los hits de Ratones).** Deben entregar el script Ex4.pm y su input file con una breve descripción. 

**Ejercicio 5 \- EMBOSS.** Instalar EMBOSS. Escribir un script que llame a uno o más programas EMBOSS para  hacer algún análisis sobre la secuencia de nucleótidos del mRNA y/o sobre la secuencia de amino ácidos de  la proteína investigada. 

Por ejemplo, pueden correr un programa que calcule los ORFs y obtenga las secuencias de proteínas  posibles o algún otro programa EMBOSS que les resulte de interés.  

Luego deben bájense la base de datos PROSITE de dominios/motivos funcionales conocidos (archivo  prosite.dat), y por medio del llamado a otro programa EMBOSS realizar un análisis de dominios de las secuencias de aminoácidos obtenidas y escribir los resultados en un archivo de salida. 

− **Input: Archivo de secuencias Fasta (por ej. Xxxxx.fas con una o más secuencias de aa).** − **Output: Archivo de resultados de dominios funcionales encontrados en las secuencias de aa.** 

**Ejercicio 6\. Trabajo con Bases de Datos Biológicas.** 

a) A partir del gen o proteína de interés para ustedes dar su link a NCBI-Gene como una entrada de Entrez,  por ej.: http://www.ncbi.nlm.nih.gov/gene/3630 

Expliquen brevemente lo que hace la proteína y por qué la eligieron. 

b) ¿Cuántos genes / proteínas homólogas se conocen en otros organismos? Utilicen la información que está 

1   
en la base de datos de HomoloGene y en las bases de datos Ensembl. Describan los resultados en ambas  bases de datos, y en qué se diferencian. Mencionen sobre qué tan común creen son estos genes o  proteínas y a qué grupos taxonómicos pertenecen (sólo en las bacterias, en los vertebrados, etc.) 

c) ¿Cuántos transcriptos y cuántas formas alternativas de *splicing* son conocidos para este gen / proteína?  ¿Cuáles de estos *splicing* alternativos se expresan? ¿Tienen funciones alternativas? Buscar evidencia de  esto en las bases de datos de NCBI y en los transcriptos de Ensembl ¿Cómo el número de *splicings* alternativos diferente entre las dos bases de datos y cuál piensan que es más precisa y por qué? 

d) ¿Con cuántas otras proteínas interactúa el producto génico de su gen? ¿Existe un patrón o relación entre las interacciones? Mencione las interacciones interesantes o inusuales. Usted encontrará las interacciones  de su gene/proteína tanto en la base de datos NCBI Gene como en la base de datos UniProt. Compare las  dos tablas entre sí. ¿Hay proteínas que interactúan únicas para cada tabla? 

e) Expliquen brevemente de qué componente celular forma parte su proteína (pista: se puede estudiar la  información de Gene Ontology \- GO), ¿A qué procesos biológicos pertenece (pista idem)? y ¿En qué función  molecular trabaja esta proteína? Los términos ontológicos de genes los pueden encontrar tanto en NCBI  Gene y en la base de datos UniProt como haciendo una búsqueda en AmiGO. 

f) Discutan brevemente en qué estructura o vías metabólicas específicas (*pathways*) estaría participando su  gen / proteína? (Reactome, KEGG son algunas bases de datos de pathways). 

g) Entrar en la base de datos de variantes genéticas dbSNP e intentar interpretar o encontrar info sobre  alguna variante (reference SNP \- rsXXXX) asociada con la patología investigada en su gen de interés. ¿Qué  variante es? ¿Hay información sobre la frecuencia que tiene esta variante en la población? ¿Qué grupo  étnico parece ser el más afectado?  

NOTA: Para hacer este ejercicio les pueden servir algunas otras bases de datos como estas (entre otras): http://www.genecards.org 

https://www.ncbi.nlm.nih.gov/snp/ (para obtener información de la variante en la población) http://www.ncbi.nlm.nih.gov/clinvar/ (para obtener información clínica del gen y sus variantes) https://ghr.nlm.nih.gov 

**Ejercicio 7\. Armar una presentación donde expliquen la enfermedad que investigaron, lo que hicieron y  los resultados que fueron obteniendo en los ejercicios del TP.**  

Los integrantes de cada grupo tendrán un máximo de 10 minutos (\~15 diapositivas) para exponer el trabajo  práctico. Comentar sobre sus investigaciones: Intro \- Métodos (código implementado, muy brevemente) \- Resultados. La correcta presentación del trabajo realizado es también parte importante de la evaluación (solo enviar el pdf de la presentación). 

2 