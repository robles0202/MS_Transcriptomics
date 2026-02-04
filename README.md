# Análisis Transcriptómico Comparativo de Subtipos de Esclerosis Múltiple

**Autor:** Rubén Ballester Robles\
**Máster:** Bioinformática - Universidad de Valencia\
**Asignatura:** Estudios in silico en Biomedicina\
**Curso:** 2025-2026

------------------------------------------------------------------------

## Descripción del Proyecto

Análisis transcriptómico del dataset GSE17048 para identificar genes diferencialmente expresados entre los tres subtipos clínicos de Esclerosis Múltiple (RRMS, SPMS, PPMS) y controles sanos, con el objetivo de caracterizar perfiles moleculares específicos e identificar biomarcadores que permitan distinguir entre formas recurrentes y progresivas.

## Estructura del Proyecto

```         
MS_Transcriptomics/
├── data/
│   ├── raw/              # Datos crudos de GEO (GSE17048)
│   ├── processed/        # Datos procesados y normalizados
│   └── metadata/         # Información de muestras y fenotipos
├── scripts/
│   ├── 01_setup/         # Configuración inicial y descarga de datos
│   ├── 02_qc/            # Control de calidad y exploración
│   ├── 03_differential_expression/  # Análisis de expresión diferencial
│   ├── 04_functional_analysis/      # Enriquecimiento funcional y redes
│   └── 05_visualization/ # Generación de figuras
├── results/
│   ├── figures/          # Figuras generadas (PCA, volcano, heatmaps...)
│   └── tables/           # Tablas de genes diferenciales
└── README.md            # Este archivo
```

------------------------------------------------------------------------

## Dataset

**Accession:** GSE17048\
**Repositorio:** Gene Expression Omnibus (GEO)\
**URL:** <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17048>

### Características del Dataset:

-   **N total:** 144 muestras
    -   36 RRMS (Remitente-Recurrente)
    -   20 SPMS (Secundaria Progresiva)
    -   43 PPMS (Primaria Progresiva)
    -   45 Controles sanos
-   **Plataforma:** Illumina HumanHT-12 V3.0 expression beadchip (GPL6947)
-   **Tipo de muestra:** Sangre total (whole blood)
-   **Población:** Pacientes australianos sin tratamiento previo

------------------------------------------------------------------------

## Pipeline de Análisis

### 1. Preparación de Datos

-   Descarga de matriz de expresión normalizada de GEO
-   Extracción y organización de metadatos
-   Anotación de sondas

### 2. Control de Calidad y Exploración

-   Análisis de Componentes Principales (PCA)
-   Clustering jerárquico
-   Detección de outliers
-   Evaluación de distribuciones (boxplots, densidad)

### 3. Análisis de Expresión Diferencial (limma)

Seis comparaciones principales: 1. RRMS vs Control 2. SPMS vs Control 3. PPMS vs Control 4. RRMS vs SPMS 5. RRMS vs PPMS 6. SPMS vs PPMS

**Criterios de significancia:** - \|log2FC\| \> 0.3 - FDR \< 0.1 (corrección Benjamini-Hochberg)

### 4. Identificación de Genes Específicos

-   Diagramas de Venn
-   Genes únicos por subtipo
-   Genes de transición RRMS→SPMS

### 5. Análisis Funcional

-   Enriquecimiento Gene Ontology (BP, MF, CC)
-   Análisis de rutas KEGG/Reactome

### 6. Redes de Interacción

-   Análisis PPI con STRING
-   Identificación de genes hub
-   Análisis de módulos funcionales

### 7. Validación

-   Comparación con literatura

## Requisitos del Sistema

### Software Requerido:

-   **R:** ≥ 4.5.0

### Paquetes de R Necesarios:

``` r
# Bioconductor
BiocManager::install(c(
  "GEOquery",      # Descarga de datos de GEO
  "limma",         # Expresión diferencial
  "biomaRt",       # Anotación de genes
  "clusterProfiler", # Enriquecimiento funcional
  "enrichplot",    # Visualización de enriquecimiento
  "org.Hs.eg.db",  # Base de datos de genes humanos
  "ComplexHeatmap" # Heatmaps avanzados
))

# CRAN
install.packages(c(
  "tidyverse",     # Manipulación de datos
  "pheatmap",      # Heatmaps
  "ggplot2",       # Visualización
  "ggrepel",       # Etiquetas en gráficos
  "VennDiagram",   # Diagramas de Venn
  "UpSetR",        # Gráficos UpSet
  "RColorBrewer",  # Paletas de colores
  "gridExtra"      # Organización de gráficos
))
```

## Instrucciones de Uso

### 1. Clonar el repositorio

``` bash
git clone https://github.com/robles0202/MS_Transcriptomics.git
cd MS_Transcriptomics_TFM
```

### 2. Instalar dependencias

``` r
source("scripts/01_setup/install_packages.R")
```

### 3. Descargar datos

``` r
source("scripts/01_setup/download_data.R")
```

### 4. Ejecutar análisis

Los scripts deben ejecutarse en orden numérico:

``` r
# Control de calidad
source("scripts/02_qc/01_quality_control.R")
source("scripts/02_qc/02_exploratory_analysis.R")

# Expresión diferencial
source("scripts/03_differential_expression/01_limma_analysis.R")

# Y así sucesivamente...
```

## Contacto

**Rubén Ballester Robles**\
Máster en Bioinformática\
Universidad de Valencia\
rubaro\@alumni.uv.es

**Última actualización:** Diciembre 2025
