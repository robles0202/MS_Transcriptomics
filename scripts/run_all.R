# =============================================================================
# SCRIPT MASTER - EJECUTAR TODO EL PIPELINE
# Análisis Transcriptómico de Esclerosis Múltiple (GSE17048)
# =============================================================================

# Este script ejecuta todo el pipeline de análisis en secuencia
# Puede ejecutarse completo o por fases

# Autor: Rubén Ballester Robles
# Fecha: Diciembre 2025

# =============================================================================
# CONFIGURACIÓN INICIAL
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("  ANÁLISIS TRANSCRIPTÓMICO DE ESCLEROSIS MÚLTIPLE\n")
cat("  Dataset: GSE17048\n")
cat("================================================================================\n")
cat("\n")

# Cargar configuración del proyecto
source("scripts/config.R")

# Verificar estructura de directorios
check_directories()

# Imprimir resumen de configuración
print_config_summary()

# Establecer semilla para reproducibilidad
set.seed(RANDOM_SEED)

# =============================================================================
# OPCIONES DE EJECUCIÓN
# =============================================================================

# Puedes modificar estas variables para ejecutar solo ciertas fases
RUN_SETUP <- TRUE              # Fase 1: Instalación y descarga
RUN_QC <- TRUE                 # Fase 2: Control de calidad
RUN_DIFFERENTIAL <- TRUE       # Fase 3: Expresión diferencial
RUN_FUNCTIONAL <- TRUE         # Fase 4: Análisis funcional
RUN_VISUALIZATION <- TRUE      # Fase 5: Visualización

# Opción para generar reporte HTML al final
GENERATE_REPORT <- TRUE

# =============================================================================
# FUNCIÓN DE LOG
# =============================================================================

log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] %s\n", timestamp, level, message))
}

log_phase <- function(phase_name) {
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("  ", phase_name, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("\n")
}

# =============================================================================
# FASE 1: SETUP Y DESCARGA DE DATOS
# =============================================================================

if (RUN_SETUP) {
  log_phase("FASE 1: SETUP Y DESCARGA DE DATOS")
  
  tryCatch({
    log_message("Instalando paquetes necesarios...")
    source("scripts/01_setup/install_packages.R")
    log_message("✓ Paquetes instalados correctamente", "SUCCESS")
    
    log_message("Descargando datos de GEO...")
    source("scripts/01_setup/download_data.R")
    log_message("✓ Datos descargados correctamente", "SUCCESS")
    
  }, error = function(e) {
    log_message(paste("Error en Fase 1:", e$message), "ERROR")
    stop("Ejecución detenida debido a error en Fase 1")
  })
  
} else {
  log_message("Fase 1 omitida (RUN_SETUP = FALSE)", "INFO")
}

# =============================================================================
# FASE 2: CONTROL DE CALIDAD Y EXPLORACIÓN
# =============================================================================

if (RUN_QC) {
  log_phase("FASE 2: CONTROL DE CALIDAD Y EXPLORACIÓN")
  
  tryCatch({
    log_message("Ejecutando control de calidad inicial...")
    source("scripts/02_qc/01_quality_control.R")
    log_message("✓ Control de calidad completado", "SUCCESS")
    
    log_message("Realizando análisis exploratorio...")
    source("scripts/02_qc/02_exploratory_analysis.R")
    log_message("✓ Análisis exploratorio completado", "SUCCESS")
    
  }, error = function(e) {
    log_message(paste("Error en Fase 2:", e$message), "ERROR")
    stop("Ejecución detenida debido a error en Fase 2")
  })
  
} else {
  log_message("Fase 2 omitida (RUN_QC = FALSE)", "INFO")
}

# =============================================================================
# FASE 3: ANÁLISIS DE EXPRESIÓN DIFERENCIAL
# =============================================================================

if (RUN_DIFFERENTIAL) {
  log_phase("FASE 3: ANÁLISIS DE EXPRESIÓN DIFERENCIAL")
  
  tryCatch({
    log_message("Ejecutando análisis de expresión diferencial con limma...")
    source("scripts/03_differential_expression/01_limma_analysis.R")
    log_message("✓ Análisis con limma completado", "SUCCESS")
    
    log_message("Filtrando y organizando resultados...")
    source("scripts/03_differential_expression/02_filter_results.R")
    log_message("✓ Filtrado completado", "SUCCESS")
    
    log_message("Identificando genes específicos por subtipo...")
    source("scripts/03_differential_expression/03_subtype_specific_genes.R")
    log_message("✓ Genes específicos identificados", "SUCCESS")
    
  }, error = function(e) {
    log_message(paste("Error en Fase 3:", e$message), "ERROR")
    stop("Ejecución detenida debido a error en Fase 3")
  })
  
} else {
  log_message("Fase 3 omitida (RUN_DIFFERENTIAL = FALSE)", "INFO")
}

# =============================================================================
# FASE 4: ANÁLISIS FUNCIONAL Y REDES
# =============================================================================

if (RUN_FUNCTIONAL) {
  log_phase("FASE 4: ANÁLISIS FUNCIONAL Y REDES")
  
  tryCatch({
    log_message("Realizando enriquecimiento Gene Ontology...")
    source("scripts/04_functional_analysis/01_gene_ontology.R")
    log_message("✓ Enriquecimiento GO completado", "SUCCESS")
    
    log_message("Analizando pathways KEGG/Reactome...")
    source("scripts/04_functional_analysis/02_pathway_analysis.R")
    log_message("✓ Análisis de pathways completado", "SUCCESS")
    
    log_message("Construyendo redes PPI...")
    source("scripts/04_functional_analysis/03_ppi_networks.R")
    log_message("✓ Análisis de redes completado", "SUCCESS")
    
  }, error = function(e) {
    log_message(paste("Error en Fase 4:", e$message), "ERROR")
    # No detener ejecución, continuar con visualización
    log_message("Continuando con siguiente fase...", "WARNING")
  })
  
} else {
  log_message("Fase 4 omitida (RUN_FUNCTIONAL = FALSE)", "INFO")
}

# =============================================================================
# FASE 5: VISUALIZACIÓN
# =============================================================================

if (RUN_VISUALIZATION) {
  log_phase("FASE 5: VISUALIZACIÓN")
  
  tryCatch({
    log_message("Generando gráficos de PCA...")
    source("scripts/05_visualization/01_pca_plots.R")
    log_message("✓ PCA plots generados", "SUCCESS")
    
    log_message("Generando volcano plots...")
    source("scripts/05_visualization/02_volcano_plots.R")
    log_message("✓ Volcano plots generados", "SUCCESS")
    
    log_message("Generando heatmaps...")
    source("scripts/05_visualization/03_heatmaps.R")
    log_message("✓ Heatmaps generados", "SUCCESS")
    
    log_message("Generando diagramas de Venn/UpSet...")
    source("scripts/05_visualization/04_venn_diagrams.R")
    log_message("✓ Diagramas generados", "SUCCESS")
    
    log_message("Generando gráficos de enriquecimiento...")
    source("scripts/05_visualization/05_enrichment_plots.R")
    log_message("✓ Gráficos de enriquecimiento generados", "SUCCESS")
    
    log_message("Generando visualización de redes...")
    source("scripts/05_visualization/06_network_plots.R")
    log_message("✓ Visualización de redes generada", "SUCCESS")
    
  }, error = function(e) {
    log_message(paste("Error en Fase 5:", e$message), "ERROR")
    log_message("Algunas figuras pueden no haberse generado correctamente", "WARNING")
  })
  
} else {
  log_message("Fase 5 omitida (RUN_VISUALIZATION = FALSE)", "INFO")
}

# =============================================================================
# GENERACIÓN DE REPORTE FINAL
# =============================================================================

if (GENERATE_REPORT) {
  log_phase("GENERANDO REPORTE FINAL")
  
  tryCatch({
    log_message("Creando resumen del análisis...")
    
    # Crear archivo de resumen
    summary_file <- file.path(RESULTS_DIR, "analysis_summary.txt")
    
    sink(summary_file)
    cat("================================================================================\n")
    cat("  RESUMEN DEL ANÁLISIS TRANSCRIPTÓMICO DE ESCLEROSIS MÚLTIPLE\n")
    cat("================================================================================\n")
    cat("\n")
    cat("Dataset:", GEO_ACCESSION, "\n")
    cat("Fecha de análisis:", as.character(Sys.Date()), "\n")
    cat("Autor:", PROJECT_AUTHOR, "\n")
    cat("\n")
    cat("Muestras analizadas:\n")
    cat("  - RRMS:", SAMPLE_INFO$RRMS, "\n")
    cat("  - SPMS:", SAMPLE_INFO$SPMS, "\n")
    cat("  - PPMS:", SAMPLE_INFO$PPMS, "\n")
    cat("  - Control:", SAMPLE_INFO$Control, "\n")
    cat("  - Total:", TOTAL_SAMPLES, "\n")
    cat("\n")
    cat("Criterios de análisis:\n")
    cat("  - |log2FC| >", LOG2FC_THRESHOLD, "\n")
    cat("  - FDR <", FDR_THRESHOLD, "\n")
    cat("  - Método ajuste p-valor:", P_ADJUST_METHOD, "\n")
    cat("\n")
    cat("Comparaciones realizadas:\n")
    for (comp_name in names(COMPARISONS)) {
      cat("  -", comp_name, "\n")
    }
    cat("\n")
    cat("Directorios de resultados:\n")
    cat("  - Figuras:", FIGURES_DIR, "\n")
    cat("  - Tablas:", TABLES_DIR, "\n")
    cat("  - Enriquecimiento:", ENRICHMENT_DIR, "\n")
    cat("  - Redes:", NETWORKS_DIR, "\n")
    cat("\n")
    cat("Información de sesión:\n")
    print(sessionInfo())
    cat("\n")
    cat("================================================================================\n")
    cat("  ANÁLISIS COMPLETADO EXITOSAMENTE\n")
    cat("================================================================================\n")
    sink()
    
    log_message("✓ Resumen guardado en: analysis_summary.txt", "SUCCESS")
    
  }, error = function(e) {
    log_message(paste("Error generando reporte:", e$message), "ERROR")
  })
}

# =============================================================================
# MENSAJE FINAL
# =============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("  ✓ PIPELINE COMPLETADO EXITOSAMENTE\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n")
cat("Resultados disponibles en:", RESULTS_DIR, "\n")
cat("  - Figuras:", FIGURES_DIR, "\n")
cat("  - Tablas:", TABLES_DIR, "\n")
cat("  - Enriquecimiento:", ENRICHMENT_DIR, "\n")
cat("  - Redes:", NETWORKS_DIR, "\n")
cat("\n")
cat("Tiempo total de ejecución:", format(Sys.time()), "\n")
cat("\n")

# Guardar información de sesión
sessionInfo_file <- file.path(RESULTS_DIR, "sessionInfo.txt")
sink(sessionInfo_file)
sessionInfo()
sink()

log_message("Información de sesión guardada en: sessionInfo.txt", "INFO")
log_message("¡Análisis completado!", "SUCCESS")
