# SOURCES.md

Tenemos tres carpetas principales:
    - analysis: scripts de bash para ejecutar los distintos pasos del análisis
    - simulations: scripts de R para simular datos de una célula única y calcular la métricas
    - src: scripts de R que implementan los cálculos de las métricas

En esencial bruto del proyecto es:

    - SingleCheck.sh que depende de src/*.R
    - analysis/*.sh que son wrappers de herramientas externas (picard, samtools, etc) para analizar los datos

## Propuestas

### Algunas mejoras en la calidad de vida:
    - Mejorar la cohesión del pipeline (formato homogéneo, nombres de archivos consistentes, etc)
    - Pasar a un modelo jerárquico más claro (leemos una sola vez los datos y los pasamos a las etapas de forma controlada, evitando relecturas innecesarias)
    - Las diferentes etapas de un pipeline deben poder ejecutarse en paralelo

### Testers sólidos:
    - Vamos a mantener, en la medida de lo posible la modularidad, para cada módulo vamos a generar unos tests automatizados para
    asegurarnos de que nuestra versión es correcta. Aquí haremos el esfuerzo extra y usaremos TDD.

### Eficiencia:

    - La eficiencia finalmente va a depender de la capacidad de paralelización del pipeline como un todo (partir el dataset)
    - Actualmente están intentando cargar todo un dataset en un único nodo de 60gb.
    - Las herramientas se pueden mejorar marginalmente, son código muy simples ya de por sí.

## Descripción para trabajo interno de cada fichero fuente

.
├── CreateInputApp
├── README.md
├── ShinyApp
│   ├── server.R
│   └── ui.R
├── SingleCheck
├── Workflow-SingleCheck.png
├── analysis
│   ├── CheckSeqDepth.sh
│   ├── CollectInsertSizeMetrics.sh
│   ├── CreatePlotsHUANG.sh
│   ├── CreatePlotsWANG.sh
│   ├── DownsampleSam.sh
│   ├── GetSeqDepth.sh
│   ├── GetSeqDepthFixingReadLength.sh
│   ├── RunExample.sh
│   ├── RunHUANG.sh
│   ├── RunWANG.sh
│   ├── SimulateSingleCellReads.sh
│   ├── SingleCheckArray
│   └── SortSam.sh
├── simulations
│   ├── Autocorrelation.R
│   ├── CoefficientOfVariation.R
│   ├── GiniIndex.R
│   ├── MAD.R
│   └── SmallSimulations.R
└── src
    ├── Autocorrelation.R
    ├── CoefficientOfVariation.R
    ├── GiniIndex.R
    └── MAD.R

### scripts en SRC:
    Son todos muy similares:

    Todos leen un fichero de entrada en argv1 con un sufijo concreto (.shiftedcov.txt o .contiguous.txt, etc)
    lo pasan a una tabla con 2 ó 3 columnas (depth, count) ó (depth, depth_fwd, count) y en algunos casos reciben
    un parametro extra (argv2 es el lag en Autocorrelation.R)

    Procesan un calculo sencillo con dplyr y generan un fichero de salida con los datos. (nombre script).arg1.txt

### scripts en simulations

    Los que llevan el mismo nombre que en SRC son versiones equivalentes a los scripts de src, pero preparados para
    ejecutarse desde un entorno probablemente interactivo (u otro script), no importan la librería dplyr (esperan que
    la importemos nosotros) y no leen los datos de un fichero, sino que reciben un data.frame como argumento.

    El fichero SmallSimulations.R parece que importa las librerías necesarias y llama a las simulaciones tomando como
    datos de entrada unos ficheros hardcodeados. (que no existen, claro está). 

    En general no tienen mucho interés

### Scripts en ShinyApp

    Es un visualizador web simple, fuera del alcance de este proyecto salvo que lo solicite el equipo de biología.

### Scripts en analysis

    En esencia estos scripts solo transforman con herramientas sencillas (wc, awk, cat, etc) los ficheros de entrada
    y en varios casos llaman a las librerías picard tools y samtools para hacer las operaciones más complejas.

    - CheckSeqDepth.sh: Wrapper de samtools y picard tools (¿Reimplementamos esto?)
    - CollectInsertSizeMetrics.sh: Wrapper de picard tools
    - CreatePlotsHUANG.sh y CreatePlotsWANG.sh: Manipulación simple de datos
    - DownsampleSam.sh: De nuevo un wrapper de picard
    - GetSeqDepth.sh: wrapper de picard + samtools
    - GetSeqDepthFixingReadLength.sh: similar al anterior (solo picard)
    - SortSam.sh: wrapper de picard

    También hay varios scripts que llaman a estos anteriores:

    - RunExample.sh, RuHuang.sh, RunWang.sh son los más sencillos
    - SimulateSingleCellReads.sh y SingleCheckArray son más complejos

### Ficheros en la raíz

    Solo hay dos ficheros relevantes: CreateInputApp y SingleCheck

    - CreateInputApp: Script bash que une varios ficheros de entrada
    - SingleCheck: Script principal (slurm) que llama a los ficheros SRC y realiza todo el trabajo