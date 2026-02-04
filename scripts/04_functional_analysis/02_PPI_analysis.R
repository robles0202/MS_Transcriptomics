# Este script analiza las redes de interacción entre los genes
# diferencialmente expresados usando STRING database

# Cargar configuración
source("scripts/config.R")

# Cargar librerías
library(tidyverse)
library(STRINGdb)
library(igraph)

set.seed(RANDOM_SEED)

# Cargar datos
cat("Paso 1/5: Cargando datos...\n")

all_results <- readRDS(file.path(DATA_PROCESSED_DIR, "limma_results.rds"))
gene_annotation <- read.csv(file.path(DATA_RAW_DIR, "annotation.csv"))

cat("Columnas en anotación:\n")
print(names(gene_annotation))

cat("\n")
cat("Datos cargados\n\n")

# Preparar genes para STRING
cat("Paso 2/5: Preparando datos para STRING...\n")

# Función para preparar genes
prepare_genes_for_string <- function(deg_results, subtype_name) {
  degs <- deg_results %>%
    filter(regulation != "NS")
  annotation_subset <- gene_annotation %>%
    dplyr::select(probe_id, gene_symbol) %>%
    dplyr::rename(
      illumina_id = probe_id,
      symbol = gene_symbol
    )
  degs <- degs %>%
    left_join(annotation_subset, 
              by = c("gene_id" = "illumina_id"))
  degs <- degs %>%
    filter(!is.na(symbol) & symbol != "" & symbol != "NA")
  
  cat("  ", subtype_name, ":", nrow(degs), "genes con símbolo\n")
  
  return(degs)
}

rrms_genes <- prepare_genes_for_string(all_results[["RRMS_vs_Control"]], "RRMS")
ppms_genes <- prepare_genes_for_string(all_results[["PPMS_vs_Control"]], "PPMS")

cat("\n Datos preparados\n\n")


# Análisis STRING
cat("Paso 3/5: Analizando red PPI con STRING...\n")

# Función para análisis STRING
analyze_string_network <- function(genes_df, subtype_name, score_threshold = 400) {
  
  cat("  Conectando a STRING database para", subtype_name, "...\n")
  
  string_db <- STRINGdb$new(
    version = "12.0",
    species = 9606,  # Homo sapiens
    score_threshold = score_threshold,
    network_type = "full",
    input_directory = ""
  )
  
  genes_mapped <- string_db$map(
    genes_df, 
    "symbol", 
    removeUnmappedRows = TRUE
  )
  
  cat("    Genes mapeados:", nrow(genes_mapped), "\n")
  
  # Obtener interacciones
  interactions <- string_db$get_interactions(genes_mapped$STRING_id)
  
  cat("    Interacciones encontradas:", nrow(interactions), "\n")

  # Crear grafo
  if(nrow(interactions) > 0) {
    g <- graph_from_data_frame(
      interactions[, c("from", "to")], 
      directed = FALSE,
      vertices = unique(c(interactions$from, interactions$to))
    )
    metrics <- data.frame(
      node = V(g)$name,
      degree = degree(g),
      betweenness = betweenness(g),
      closeness = closeness(g),
      eigenvector = eigen_centrality(g)$vector
    )
    metrics <- metrics %>%
      left_join(
        genes_mapped %>% select(STRING_id, symbol),
        by = c("node" = "STRING_id")
      )
    
    return(list(
      string_db = string_db,
      genes_mapped = genes_mapped,
      interactions = interactions,
      graph = g,
      metrics = metrics
    ))
    
  }
  return(NULL)
}

# Análisis para RRMS
rrms_network <- NULL
if(nrow(rrms_genes) >= 5) {
  rrms_network <- analyze_string_network(rrms_genes, "RRMS", score_threshold = 400)
  
  if(!is.null(rrms_network)) {
    # Guardar métricas
    write.csv(rrms_network$metrics, 
              file.path(TABLES_DIR, "ppi_metrics_rrms.csv"), 
              row.names = FALSE)
    
    # Visualizar red
    png(file.path(FIGURES_DIR, "ppi_network_rrms.png"), 
        width = 12, height = 10, units = "in", res = FIG_DPI)
    
    rrms_network$string_db$plot_network(rrms_network$genes_mapped$STRING_id)
    
    dev.off()
    
    # Top hubs
    top_hubs <- rrms_network$metrics %>%
      arrange(desc(degree)) %>%
      head(10)
    
    cat("\n  Top 10 hubs en RRMS:\n")
    print(top_hubs[, c("symbol", "degree", "betweenness")])
    cat("\n")
  }
}

# Análisis para PPMS
ppms_network <- NULL
if(nrow(ppms_genes) >= 5) {
  ppms_network <- analyze_string_network(ppms_genes, "PPMS", score_threshold = 400)
  
  if(!is.null(ppms_network)) {
    # Guardar métricas
    write.csv(ppms_network$metrics, 
              file.path(TABLES_DIR, "ppi_metrics_ppms.csv"), 
              row.names = FALSE)
    
    # Visualizar red
    png(file.path(FIGURES_DIR, "ppi_network_ppms.png"), 
        width = 12, height = 10, units = "in", res = FIG_DPI)
    
    ppms_network$string_db$plot_network(ppms_network$genes_mapped$STRING_id)
    
    dev.off()
    
    # Top hubs
    top_hubs <- ppms_network$metrics %>%
      arrange(desc(degree)) %>%
      head(10)
    
    cat("\n  Top 10 hubs en PPMS:\n")
    print(top_hubs[, c("symbol", "degree", "betweenness")])
    cat("\n")
  }
}

cat("Análisis STRING completado\n\n")


# Análisis de clusters
cat("Paso 4/5: Detectando módulos en la red...\n")

detect_modules <- function(network_result, subtype_name) {
  
  if(is.null(network_result)) return(NULL)
  
  g <- network_result$graph
  
  # Usar algoritmo Louvain
  communities <- cluster_louvain(g)
  
  # Añadir comunidad a las métricas
  metrics_with_modules <- network_result$metrics
  metrics_with_modules$module <- membership(communities)[metrics_with_modules$node]
  
  cat("  ", subtype_name, ":", length(unique(membership(communities))), "módulos detectados\n")
  
  write.csv(metrics_with_modules, 
            file.path(TABLES_DIR, paste0("ppi_modules_", tolower(subtype_name), ".csv")), 
            row.names = FALSE)
  
  return(metrics_with_modules)
}

rrms_modules <- detect_modules(rrms_network, "RRMS")
ppms_modules <- detect_modules(ppms_network, "PPMS")


cat("Detección de módulos completada\n\n")

# Resumen
cat("Paso 5/5: Generando resumen...\n")

summary_data <- data.frame(
  Subtipo = character(),
  Genes_Totales = integer(),
  Genes_Mapeados = integer(),
  Interacciones = integer(),
  Nodos_Red = integer(),
  Modulos = integer(),
  stringsAsFactors = FALSE
)

summary_data <- rbind(summary_data, data.frame(
  Subtipo = "RRMS",
  Genes_Totales = nrow(rrms_genes),
  Genes_Mapeados = nrow(rrms_network$genes_mapped),
  Interacciones = nrow(rrms_network$interactions),
  Nodos_Red = vcount(rrms_network$graph),
  Modulos = if(exists("rrms_modules")) length(unique(rrms_modules$module)) else 0
))


summary_data <- rbind(summary_data, data.frame(
  Subtipo = "PPMS",
  Genes_Totales = nrow(ppms_genes),
  Genes_Mapeados = nrow(ppms_network$genes_mapped),
  Interacciones = nrow(ppms_network$interactions),
  Nodos_Red = vcount(ppms_network$graph),
  Modulos = if(exists("ppms_modules")) length(unique(ppms_modules$module)) else 0
))

write.csv(summary_data, 
          file.path(TABLES_DIR, "ppi_network_summary.csv"), 
          row.names = FALSE)

cat("\n")
print(summary_data)
cat("\n")

cat("Resumen generado\n\n")

cat("Archivos generados:\n")
cat("  - Tablas:\n")
cat("    * ppi_metrics_rrms.csv\n")
cat("    * ppi_modules_rrms.csv\n")
cat("    * ppi_metrics_ppms.csv\n")
cat("    * ppi_modules_ppms.csv\n")
cat("    * ppi_network_summary.csv\n")
cat("  - Figuras:\n")
cat("    * ppi_network_rrms.png\n")
cat("    * ppi_network_ppms.png\n")
cat("\n")

cat("Para continuar ejecuta: source('scripts/05_visualization/01_venn_diagrams.R')\n")