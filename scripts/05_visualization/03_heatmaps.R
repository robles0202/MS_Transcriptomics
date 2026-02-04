# Este script se encarga de la generación de  heatmaps

# Cargar configuración
source("scripts/config.R")

# Cargar librerías
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

# Cargar datos
cat("Paso 1/5: Cargando datos...\n")

expression_filtered <- readRDS(file.path(DATA_PROCESSED_DIR, "expression_matrix_filtered.rds"))
sample_metadata <- read.csv(file.path(DATA_METADATA_DIR, "sample_metadata.csv"))
all_results <- readRDS(file.path(DATA_PROCESSED_DIR, "limma_results.rds"))

cat("Datos cargados\n\n")

# Preparar anotaciones de muestras
cat("Paso 2/5: Preparando anotaciones...\n\n")

# Ordenar por grupo
sample_metadata <- sample_metadata %>%
  arrange(group) %>%
  filter(geo_accession %in% colnames(expression_filtered))

# Reordenar matriz de expresión
expression_ordered <- expression_filtered[, sample_metadata$geo_accession]

# Crear anotaciones
sample_annotation <- sample_metadata %>%
  select(geo_accession, group) %>%
  column_to_rownames("geo_accession")
annotation_colors <- list(
  group = unlist(COLOR_PALETTE)
)

# Generación de heatmaps
cat("Paso 3/4: Generando heatmap de Top genes diferenciales...\n")

# Obtener top genes más significativos de todas las comparaciones
all_degs <- bind_rows(all_results) %>%
  filter(regulation != "NS") %>%
  arrange(adj.P.Val) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  head(TOP_GENES_HEATMAP)

# Extraer expresión de estos genes
heatmap_data <- expression_ordered[all_degs$gene_id, ]

# Escalar datos por filas (z-score)
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# Crear anotación columnas
ha_column <- HeatmapAnnotation(
  Group = sample_annotation$group,
  col = annotation_colors,
  show_legend = TRUE,
  show_annotation_name = TRUE
)

ht1 <- Heatmap(
  heatmap_data_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  top_annotation = ha_column,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title = paste("Top", TOP_GENES_HEATMAP, "Genes Diferencialmente Expresados"),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2")
  )
)

png(file.path(FIGURES_DIR, "heatmap_top_degs.png"),
    width = 12, height = 10, units = "in", res = FIG_DPI)
draw(ht1)
dev.off()

cat("Heatmap guardado: heatmap_top_degs.png\n\n")

cat("Paso 4/4: Generando heatmap de genes específicos por subtipo...\n")

# Cargar listas de genes específicos
gene_lists <- readRDS(file.path(DATA_PROCESSED_DIR, "subtype_specific_gene_lists.rds"))

# Seleccionar top genes de cada categoría
n_per_category <- 15

genes_to_plot <- c(
  head(gene_lists$rrms_specific, n_per_category),
  head(gene_lists$ppms_specific, n_per_category)
)

# Filtrar genes que existen en la matriz
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(expression_ordered)]

  
# Extraer expresión
heatmap_data2 <- expression_ordered[genes_to_plot, ]
heatmap_data2_scaled <- t(scale(t(heatmap_data2)))

# Crear anotación de filas
gene_category <- data.frame(
  gene = genes_to_plot,
  category = c(
    rep("RRMS", min(n_per_category, length(gene_lists$rrms_specific))),
    rep("PPMS", min(n_per_category, length(gene_lists$ppms_specific)))
  )[1:length(genes_to_plot)]
)

gene_annotation <- gene_category %>%
  column_to_rownames("gene")

ha_row <- rowAnnotation(
  Category = gene_annotation$category,
  col = list(Category = c(
    "RRMS" = COLOR_PALETTE$RRMS,
    "PPMS" = COLOR_PALETTE$PPMS
  )),
  show_legend = TRUE
)

ht2 <- Heatmap(
  heatmap_data2_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  top_annotation = ha_column,
  left_annotation = ha_row,
  cluster_rows = FALSE,  # Mantener agrupación por categoría
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title = "Genes Específicos por Subtipo de MS",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)
png(file.path(FIGURES_DIR, "heatmap_subtype_specific.png"),
    width = 12, height = 12, units = "in", res = FIG_DPI)
draw(ht2)
dev.off()

cat("Heatmap guardado: heatmap_subtype_specific.png\n\n")


cat("Figuras generadas:\n")
cat("  - heatmap_top_degs.png\n")
cat("  - heatmap_subtype_specific.png\n")
cat("\n")