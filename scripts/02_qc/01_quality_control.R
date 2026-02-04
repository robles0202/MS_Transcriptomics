# Este script realiza un control de calidad inicial a los datos que acabamos de
# cargar

# Cargar configuración
source("scripts/config.R")

# Cargar librerías
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

# Establecer semilla para reproducibilidad
set.seed(RANDOM_SEED)

# Cargar datos
cat("Paso 1/6: Cargando datos...\n")

# Cargar matriz de expresión
expression_matrix <- readRDS(file.path(DATA_RAW_DIR, "expression_matrix.rds"))

# Cargar metadatos
sample_metadata <- read.csv(file.path(DATA_METADATA_DIR, "sample_metadata.csv"))

cat("Datos cargados\n")
cat("Genes/sondas:", nrow(expression_matrix), "\n")
cat("Muestras:", ncol(expression_matrix), "\n\n")

# Distribución de expresión por muestra
cat("Paso 2/6: Generando boxplot...\n")

# Preparar datos para boxplot
expr_df <- as.data.frame(expression_matrix)
expr_df$probe <- rownames(expr_df)

expr_long <- expr_df %>%
  pivot_longer(cols = -probe, names_to = "sample", values_to = "expression")

# Añadir información de grupo
expr_long <- expr_long %>%
  left_join(sample_metadata %>% select(geo_accession, group), 
            by = c("sample" = "geo_accession"))

cat("Grupos asignados:\n")
print(table(expr_long$group))
cat("\nPrimeros registros de expr_long:\n")
print(head(expr_long %>% select(sample, probe, expression, group)))


# Boxplot de control de calidad
p_boxplot <- ggplot(expr_long, aes(x = sample, y = expression)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, 
               fill = "lightblue", color = "gray30", linewidth = 0.2) +
  labs(
    title = "Distribución de Expresión por Muestra",
    subtitle = paste0("Dataset: ", GEO_ACCESSION, " (n=144 muestras)"),
    x = "Muestras",
    y = "Expresión (log2)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(FIGURES_DIR, "qc_boxplot.png"), p_boxplot, 
       width = 14, height = 6, dpi = FIG_DPI)

# Guardar figura
ggsave(file.path(FIGURES_DIR, "qc_boxplot.png"), p_boxplot, 
       width = 16, height = 6, dpi = FIG_DPI)

cat("Boxplot guardado: results/figures/qc_boxplot.png\n\n")


# Gráfico de densidad
cat("Paso 3/6: Generando gráfico de densidad...\n")

expr_density_df <- expr_df %>%
  pivot_longer(cols = -probe, names_to = "sample", values_to = "expression") %>%
  left_join(sample_metadata %>% select(geo_accession, group), 
            by = c("sample" = "geo_accession"))

p_density <- ggplot(expr_density_df, aes(x = expression, group = sample)) +
  geom_density(alpha = 0.6, linewidth = 0.5, color = "black") +  # linewidth más fino
  facet_wrap(~group, nrow = 2) +
  labs(
    title = "Densidad de Expresión por Grupo",
    subtitle = paste0("Todas las muestras (n=", ncol(expression_matrix), ")"),
    x = "Expresión (log2)",
    y = "Densidad"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(file.path(FIGURES_DIR, "qc_density.png"), p_density, 
       width = 10, height = 6, dpi = FIG_DPI)
cat("Gráfico de densidad guardado: results/figures/qc_density.png\n\n")

# Heatmap
cat("Paso 4/6: Generando heatmap...\n")

correlation_matrix <- cor(expression_matrix, method = "pearson")

# Preparar anotaciones
annotation_df <- sample_metadata %>%
  select(geo_accession, group) %>%
  column_to_rownames("geo_accession")

annotation_colors <- list(
  group = unlist(COLOR_PALETTE)
)

# Crear heatmap
png(file.path(FIGURES_DIR, "qc_heatmap.png"), 
    width = 12, height = 12, units = "in", res = FIG_DPI)
pheatmap(correlation_matrix,
         annotation_col = annotation_df,
         annotation_colors = annotation_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Heatmap (n=144)",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(min(correlation_matrix), max(correlation_matrix), 
                      length.out = 101))
dev.off()
cat("Heatmap guardado: results/figures/qc_heatmap.png\n\n")

# Estadísticas
cat("Paso 5/6: Calculando estadísticas\n")

# Calcular estadísticas por muestra
qc_stats <- data.frame(
  sample = colnames(expression_matrix),
  median_expr = apply(expression_matrix, 2, median),
  mean_expr = apply(expression_matrix, 2, mean),
  sd_expr = apply(expression_matrix, 2, sd),
  q25 = apply(expression_matrix, 2, function(x) quantile(x, 0.25)),
  q75 = apply(expression_matrix, 2, function(x) quantile(x, 0.75)),
  iqr = apply(expression_matrix, 2, IQR)
)

# Añadir información de grupo
qc_stats <- qc_stats %>%
  left_join(sample_metadata %>% select(geo_accession, group), 
            by = c("sample" = "geo_accession"))

# Detectar outliers basados en la mediana
median_mean <- mean(qc_stats$median_expr)
median_sd <- sd(qc_stats$median_expr)
qc_stats$outlier <- abs(qc_stats$median_expr - median_mean) > 3 * median_sd

# Resumen
cat("  Estadísticas generales:\n")
cat("    Mediana de expresión promedio:", round(mean(qc_stats$median_expr), 2), "\n")
cat("    Desviación estándar promedio:", round(mean(qc_stats$sd_expr), 2), "\n")
cat("    Outliers detectados:", sum(qc_stats$outlier), "\n")
cat("\n")

# Guardar estadísticas
write.csv(qc_stats, file.path(TABLES_DIR, "qc_statistics.csv"), row.names = FALSE)
cat("Estadísticas guardadas: results/tables/qc_statistics.csv\n\n")

# Graficar estadísticas 
cat("Paso 6/6: Generando gráficos por grupo...\n")

# Boxplot
p_median_by_group <- ggplot(qc_stats, aes(x = group, y = median_expr, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = unlist(COLOR_PALETTE)) +
  labs(
    title = "Mediana de Expresión por Grupo",
    x = "Grupo",
    y = "Mediana(log2)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(FIGURES_DIR, "qc_median.png"), p_median_by_group, 
       width = 8, height = 6, dpi = FIG_DPI)

cat("Gráfico por grupo guardado: results/figures/qc_median.png\n\n")

cat("Ejecuta para seguir: source('scripts/02_qc/02_exploratory_analysis.R')\n")
cat("\n")

# Guardar workspace para siguiente análisis
save.image(file.path(DATA_PROCESSED_DIR, "qc_workspace.RData"))
cat("Workspace guardado: data/processed/qc_workspace.RData\n\n")