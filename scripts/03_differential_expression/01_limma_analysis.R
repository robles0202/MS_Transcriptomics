# Este script realiza un análisis de expresión diferencial usando el paquete 
# limma

# Cargar configuración
source("scripts/config.R")

# Cargar librerías
library(limma)
library(tidyverse)

# Establecer semilla para reproducibilidad
set.seed(RANDOM_SEED)

# Cargar datos
cat("Paso 1/7: Cargando datos filtrados...\n")

expression_filtered <- readRDS(file.path(DATA_PROCESSED_DIR, 
                                         "expression_matrix_filtered.rds"))
sample_metadata <- read.csv(file.path(DATA_METADATA_DIR, "sample_metadata.csv"))

cat("Datos cargados\n")

# Preparar datos para limma
cat("Paso 2/7: Preparando diseño experimental...\n")

# Asegurar que el orden de muestras coincida
sample_metadata <- sample_metadata %>%
  filter(geo_accession %in% colnames(expression_filtered)) %>%
  arrange(match(geo_accession, colnames(expression_filtered)))

# Crear factor de grupo
group_factor <- factor(sample_metadata$group, 
                       levels = c("Control", "RRMS", "SPMS", "PPMS"))

cat("Diseño preparado\n")
cat("Distribución de grupos:\n")
print(table(group_factor))
cat("\n")

# Crear matriz de diseño
cat("Paso 3/7: Creando matriz de diseño...\n")

# Sin intercepto para facilitar los contrastes
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)

cat("Matriz de diseño creada\n")
cat("Dimensiones:", nrow(design), "x", ncol(design), "\n\n")

# Ajustar modelo lineal
cat("Paso 4/7: Ajustando modelo lineal...\n")

# Ajustar modelo lineal
fit <- lmFit(expression_filtered, design)

cat("Modelo lineal ajustado\n\n")

# Definir contrastes para las 6 comparaciones
cat("Paso 5/7: Definiendo contrastes...\n")

contrast_matrix <- makeContrasts(
  RRMS_vs_Control = RRMS - Control,
  SPMS_vs_Control = SPMS - Control,
  PPMS_vs_Control = PPMS - Control,
  RRMS_vs_SPMS = RRMS - SPMS,
  RRMS_vs_PPMS = RRMS - PPMS,
  SPMS_vs_PPMS = SPMS - PPMS,
  levels = design
)

cat("Contrastes definidos:\n")
print(colnames(contrast_matrix))
cat("\n")

# Ajustar contrastes y aplicar Bayes
cat("Paso 6/7: Ajustando contrastes y moderación bayesiana...\n")

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

cat("Análisis estadístico completado\n\n")

# Extraer resultados
cat("Paso 7/7: Extrayendo resultados...\n\n")

all_results <- list()

for (contrast_name in colnames(contrast_matrix)) {
  cat("  Procesando:", contrast_name, "\n")
  results <- topTable(fit2, 
                      coef = contrast_name, 
                      number = Inf,  # Todos los genes
                      adjust.method = "BH",
                      sort.by = "P")
  
  # Añadir columna de contraste y de gene_id
  results$contrast <- contrast_name
  results$gene_id <- rownames(results)
  
  # Clasificar genes como UP, DOWN o NS (not significant)
  results$regulation <- ifelse(
    results$adj.P.Val < 0.1 & results$logFC > 0.3, "UP",
    ifelse(results$adj.P.Val < 0.1 & results$logFC < -0.3, "DOWN", "NS")
  )
  
  # Guardar en lista
  all_results[[contrast_name]] <- results
  
  # Resumen
  n_up <- sum(results$regulation == "UP")
  n_down <- sum(results$regulation == "DOWN")
  n_total <- n_up + n_down
  
  cat("    Genes UP:", n_up, "\n")
  cat("    Genes DOWN:", n_down, "\n")
  cat("    Total DEGs:", n_total, "\n\n")
  
  # Guardar tabla individual
  output_file <- file.path(TABLES_DIR, 
                           paste0("degs_", gsub(" ", "_", contrast_name), ".csv"))
  write.csv(results, output_file, row.names = FALSE)
}

# Crear tabla resumen
cat("Creando tabla resumen de todas las comparaciones...\n")

# Resumen de DEGs por contraste
summary_df <- data.frame(
  contrast = names(all_results),
  total_genes = sapply(all_results, nrow),
  deg_up = sapply(all_results, function(x) sum(x$regulation == "UP")),
  deg_down = sapply(all_results, function(x) sum(x$regulation == "DOWN")),
  deg_total = sapply(all_results, function(x) sum(x$regulation != "NS"))
)

cat("\n")
cat("GENES DIFERENCIALMENTE EXPRESADOS:\n")
print(summary_df)
cat("\n")

# Guardar resumen
write.csv(summary_df, 
          file.path(TABLES_DIR, "differential_expression_summary.csv"), 
          row.names = FALSE)

# Combinar todos los resultados en una tabla
all_results_combined <- bind_rows(all_results)

write.csv(all_results_combined, 
          file.path(TABLES_DIR, "all_comparisons_combined.csv"), 
          row.names = FALSE)

# Guardar objetos importantes
saveRDS(fit2, file.path(DATA_PROCESSED_DIR, "limma_fit.rds"))
saveRDS(all_results, file.path(DATA_PROCESSED_DIR, "limma_results.rds"))

cat("Ejecuta para seguir: source('scripts/03_differential_expression/
    02_subtype_specific_genes.R')\n")


# Guardar workspace
save.image(file.path(DATA_PROCESSED_DIR, "limma_workspace.RData"))