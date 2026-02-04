# Este archivo contiene todos los parámetros y configuraciones del proyecto

# INFORMACIÓN DEL DATASET
# =============================================================================
GEO_ACCESSION <- "GSE17048"
PLATFORM <- "GPL6947"
ORGANISM <- "Homo sapiens"

# Información de muestras
SAMPLE_INFO <- list(RRMS = 36,SPMS = 20,PPMS = 43,Control = 45)
TOTAL_SAMPLES <- sum(unlist(SAMPLE_INFO))

# DIRECTORIOS DEL PROYECTO
# =============================================================================
# Directorio raíz del proyecto
PROJECT_ROOT <- here::here()

# Directorios de datos
DATA_DIR <- file.path(PROJECT_ROOT, "data")
DATA_RAW_DIR <- file.path(DATA_DIR, "raw")
DATA_PROCESSED_DIR <- file.path(DATA_DIR, "processed")
DATA_METADATA_DIR <- file.path(DATA_DIR, "metadata")

# Directorio de scripts
SCRIPTS_DIR <- file.path(PROJECT_ROOT, "scripts")

# Directorios de resultados
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
TABLES_DIR <- file.path(RESULTS_DIR, "tables")

# Directorio de documentación
DOCS_DIR <- file.path(PROJECT_ROOT, "docs")

# PARÁMETROS DE ANÁLISIS DE EXPRESIÓN DIFERENCIAL
# =============================================================================

# Criterios de significancia
LOG2FC_THRESHOLD <- 1
FDR_THRESHOLD <- 0.05
P_VALUE_THRESHOLD <- 0.05

# Método de ajuste p-valor
P_ADJUST_METHOD <- "BH"

# Comparaciones a realizar
COMPARISONS <- list(
  "RRMS_vs_Control" = c("RRMS", "Control"),
  "SPMS_vs_Control" = c("SPMS", "Control"),
  "PPMS_vs_Control" = c("PPMS", "Control"),
  "RRMS_vs_SPMS" = c("RRMS", "SPMS"),
  "RRMS_vs_PPMS" = c("RRMS", "PPMS"),
  "SPMS_vs_PPMS" = c("SPMS", "PPMS")
)

# PARÁMETROS DE CONTROL DE CALIDAD
# =============================================================================

# PCA
PCA_N_COMPONENTS <- 5  # Número de componentes principales a calcular

# Detección de outliers
OUTLIER_METHOD <- "IQR"
OUTLIER_THRESHOLD <- 3   # Umbral para detección

# Correlación entre muestras
CORRELATION_METHOD <- "pearson"


# PARÁMETROS DE FILTRADO DE GENES
# =============================================================================

# Filtrado por expresión
MIN_EXPRESSION <- 6.1     # Expresión mínima
MIN_SAMPLES <- 0.2      # Proporción mínima de muestras con expresión

# Variabilidad
MIN_VARIANCE <- 0.1     # Varianza mínima

# PARÁMETROS DE ANÁLISIS FUNCIONAL
# =============================================================================

# Enriquecimiento
ENRICHMENT_FDR <- 0.05
ENRICHMENT_P_VALUE <- 0.05
MIN_GENE_SET_SIZE <- 10
MAX_GENE_SET_SIZE <- 500

# Bases de datos de enriquecimiento
ENRICHMENT_DBS <- c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome")

# Gene Ontology
GO_ONTOLOGIES <- c("BP", "MF", "CC")
GO_LEVEL <- 3

# KEGG
KEGG_ORGANISM <- "hsa"

# PARÁMETROS DE REDES PPI
# =============================================================================

# STRING database
STRING_VERSION <- "12.0"
STRING_SPECIES <- 9606  # Homo sapiens
STRING_SCORE_THRESHOLD <- 400  # Confianza mínima

# Network analysis
NETWORK_ALGORITHM <- "louvain"
MIN_MODULE_SIZE <- 5

# Hub genes
HUB_DEGREE_PERCENTILE <- 90  # Percentil para considerar hub (top 10%)

# PARÁMETROS DE VISUALIZACIÓN
# =============================================================================

# Configuración general de figuras
FIG_WIDTH <- 10
FIG_HEIGHT <- 8
FIG_DPI <- 300           # Resolución
FIG_FORMAT <- c("png", "pdf")

# Paleta de colores para subtipos
COLOR_PALETTE <- list(RRMS = "#E41A1C",SPMS = "#377EB8",PPMS = "#4DAF4A",
                      Control = "#999999")

# Paleta para expresión diferencial
EXPRESSION_COLORS <- list(upregulated = "#D62728",downregulated = "#1F77B4",
                          not_significant = "#CCCCCC")

# Heatmaps
HEATMAP_COLORS <- circlize::colorRamp2(c(-3, 0, 3),c("blue", "white", "red"))

# Número de genes a mostrar en figuras
TOP_GENES_HEATMAP <- 50
TOP_GENES_VOLCANO <- 20 
TOP_PATHWAYS <- 20

# CONFIGURACIÓN DE REPRODUCIBILIDAD
# =============================================================================

# Semilla para análisis aleatorios
RANDOM_SEED <- 42

# Número de cores para paralelización
N_CORES <- parallel::detectCores() - 1  # Dejar 1 core libre

# OPCIONES DE SESIÓN R
# =============================================================================
options(stringsAsFactors = FALSE)
options(scipen = 999)  # Evitar notación científica

# FUNCIONES AUXILIARES
# =============================================================================

create_dir_if_not_exists <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Directorio creado: ", dir_path)
  }
}

check_directories <- function() {
  dirs <- c(DATA_DIR, DATA_RAW_DIR, DATA_PROCESSED_DIR, DATA_METADATA_DIR,
            RESULTS_DIR, FIGURES_DIR, TABLES_DIR)
  
  for (dir in dirs) {
    create_dir_if_not_exists(dir)
  }
  
  message("Estructura de directorios verificada")
}

# Función para guardar figura
save_figure <- function(plot_obj, filename, width = FIG_WIDTH, 
                        height = FIG_HEIGHT, dpi = FIG_DPI) {
  base_path <- file.path(FIGURES_DIR, tools::file_path_sans_ext(filename))
  
  for (format in FIG_FORMAT) {
    filepath <- paste0(base_path, ".", format)
    
    if (format == "png") {
      ggplot2::ggsave(filepath, plot_obj, width = width, height = height, dpi = dpi)
    } else if (format == "pdf") {
      ggplot2::ggsave(filepath, plot_obj, width = width, height = height)
    }
  }
  
  message("Figura guardada: ", filename)
}

# Función para guardar tabla
save_table <- function(data, filename) {
  filepath <- file.path(TABLES_DIR, filename)
  write.csv(data, filepath, row.names = FALSE)
  message("Tabla guardada: ", filename)
}

message("Archivo de configuración cargado correctamente")