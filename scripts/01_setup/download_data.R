# Este script descarga los datos de la serie GSE17048 del respositorio GEO con 
# los que vamos a trabajar 

# Cargar configuración del proyecto
source("scripts/config.R")


# Cargar librerías
library(GEOquery)
library(Biobase)
library(tidyverse)

# Descargar datos de GEO
cat("Paso 1/5: Descargando datos de GEO...\n")

# Descargar GSE17048
gse <- tryCatch({
  getGEO(GEO_ACCESSION, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = DATA_RAW_DIR)
}, error = function(e) {
  cat("Error descargando datos:", e$message, "\n")
  cat("Intentando nuevamente...\n")
  Sys.sleep(5)  # Esperar 5 segundos
  getGEO(GEO_ACCESSION, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = DATA_RAW_DIR)
})

gse <- gse[[1]]

cat("Datos descargados correctamente\n")
cat("Muestras:", ncol(gse), "\n")
cat("Genes:", nrow(gse), "\n\n")

# Extraer matriz de expresión
cat("Paso 2/5: Extrayendo matriz de expresión...\n")
expression_matrix <- exprs(gse)

cat("Rango original:", round(range(expression_matrix, na.rm = TRUE), 2), "\n")
cat("Valores negativos:", sum(expression_matrix < 0, na.rm = TRUE), "\n")
cat("Procesando datos...\n")

# Reemplazar valores negativos por valor pequeño positivo
expression_matrix[expression_matrix < 0] <- 0.1

# Transformar todo a log2 
expression_matrix <- log2(expression_matrix)

# Verificar resultado
cat("Datos procesados\n")
cat("Nuevo rango:", round(range(expression_matrix, na.rm = TRUE), 2), "\n")
cat("NAs generados:", sum(is.na(expression_matrix)), "\n\n")

cat("Matriz extraída:", nrow(expression_matrix), "genes x", 
    ncol(expression_matrix), "muestras\n\n")

cat("Matriz de expresión extraída\n")

# Extraer metadatos 
cat("Paso 3/5: Procesando metadatos de muestras...\n")

# Obtener información fenotípica
phenoData <- pData(gse)
print(names(phenoData))
print(phenoData$characteristics_ch1)

# Extraer columnas importantes
sample_info <- data.frame(
  sample_id = rownames(phenoData),
  title = phenoData$title,
  geo_accession = phenoData$geo_accession,
  source = phenoData$source_name_ch1,
  characteristics = phenoData$characteristics_ch1
)

# Para extraer el subtipo de MS se hace uso de la variable "characteristics"
sample_info$group <- ifelse(grepl("tissue:", sample_info$characteristics), "Control",
                            ifelse(grepl("RRMS", sample_info$characteristics), "RRMS",
                                   ifelse(grepl("SPMS", sample_info$characteristics), "SPMS", 
                                          ifelse(grepl("PPMS", sample_info$characteristics), "PPMS", "Unknown"))))

# Verificar
cat("Metadatos procesados\n")
print(table(sample_info$group))
cat("\n")

# Obtener features
cat("Paso 4/5: Procesando anotación de genes...\n")

feature_data <- fData(gse)
print(names(feature_data))

annotation_table <- data.frame(
  probe_id = rownames(feature_data),
  gene_symbol = feature_data$`Gene symbol`,
  gene_name = feature_data$`Gene title`,
  gene_id = feature_data$`Gene ID`,
  stringsAsFactors = FALSE
)

cat("Anotación procesada\n")
cat("Total de sondas:", nrow(annotation_table), "\n")
cat("Sondas con símbolo:", sum(!is.na(annotation_table$gene_symbol) & 
                                 annotation_table$gene_symbol != ""), "\n\n")

# Guardar datos
cat("Paso 5/5: Guardando datos procesados...\n")

# Guardar matriz de expresión
saveRDS(expression_matrix, file = file.path(DATA_RAW_DIR, "expression_matrix.rds"))
# Guardar metadatos
write.csv(sample_info, file = file.path(DATA_METADATA_DIR, "sample_metadata.csv"), 
          row.names = FALSE)
# Guardar anotación
write.csv(annotation_table, file = file.path(DATA_RAW_DIR, "annotation.csv"), 
          row.names = FALSE)
# Guardar ExpressionSet
saveRDS(gse, file = file.path(DATA_RAW_DIR, "gse_object.rds"))


cat("Para continuar ejecuta: source('scripts/02_qc/01_quality_control.R')\n")