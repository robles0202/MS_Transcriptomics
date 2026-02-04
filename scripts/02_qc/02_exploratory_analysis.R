# Este script realiza un análisis exploratorio de nuestro estudio mediante la 
# realización de un PCA y clustering


# Cargar configuración
source("scripts/config.R")

# Cargar librerías
library(tidyverse)
library(ggplot2)
library(reshape)
library(dendextend)

# Establecer semilla para reproducibilidad
set.seed(RANDOM_SEED)

# Cargar datos
cat("Paso 1/5: Cargando datos...\n")

expression_matrix <- readRDS(file.path(DATA_RAW_DIR, "expression_matrix.rds"))
sample_metadata <- read.csv(file.path(DATA_METADATA_DIR, "sample_metadata.csv"))

cat("Datos cargados\n\n")

# Filtrado inicial
cat("Paso 2/5: Filtrando genes de baja expresión...\n")

# Filtrar genes con expresión muy baja manteniendo genes que tienen una 
# expresión mayor a 6.1 en escala logaritmica en al menos el 20% de las muestras

n_samples_threshold <- ceiling(ncol(expression_matrix) * 0.2)

genes_keep <- apply(expression_matrix, 1, function(x) {
  sum(x > 6.1) >= n_samples_threshold
})

expression_filtered <- expression_matrix[genes_keep, ]

cat("  Filtrado completado\n")
cat("  Genes originales:", nrow(expression_matrix), "\n")
cat("  Genes después del filtro:", nrow(expression_filtered), "\n")

# Guardar matriz filtrada
saveRDS(expression_filtered, 
        file.path(DATA_PROCESSED_DIR, "expression_matrix_filtered.rds"))

# PCA
cat("Paso 3/5: Realizando PCA...\n")

# Transponer: muestras en filas, genes en columnas
expr_t <- t(expression_filtered)

# PCA
pcaexpr <- prcomp(expr_t)

# Explorar el objeto pca
names(pcaexpr)
summary(pcaexpr)

# Crear dataframe
pcadf <- data.frame(pcaexpr$x, 
                    group = sample_metadata$group)
print(dim(pcadf))
print(head(pcadf))

# Calcular varianza explicada
variance_explained <- (pcaexpr$sdev^2) / sum(pcaexpr$sdev^2) * 100

# Graficar PCA 
p1 <- ggplot(data = pcadf, aes(x = PC1, y = PC2)) +
  geom_point(size = 2) +
  theme_minimal() +
  ggtitle("PCA de Muestras MS")

ggsave(file.path(FIGURES_DIR, "pca.png"), p1, 
       width = 8, height = 6, dpi = FIG_DPI)

cat("PCA guardado: results/figures/pca.png\n")

# Graficar PCA coloreado por "group"
p2 <- ggplot(data = pcadf, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = group), size = 3) +
  theme_minimal() +
  labs(
    title = "PCA de Muestras MS coloreado por grupo",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)")
  ) +
  scale_color_manual(values = unlist(COLOR_PALETTE))

ggsave(file.path(FIGURES_DIR, "pca_group.png"), p2, 
       width = 10, height = 7, dpi = FIG_DPI)

cat("PCA coloreado por grupo guardado: results/figures/pca_group.png\n")

# PC1 vs PC3
p3 <- ggplot(data = pcadf, aes(x = PC1, y = PC3)) +
  geom_point(aes(color = group), size = 3) +
  theme_minimal() +
  labs(
    title = "PCA de Muestras MS: PC1 vs PC3",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC3 (", round(variance_explained[3], 1), "%)")
  ) +
  scale_color_manual(values = unlist(COLOR_PALETTE))

ggsave(file.path(FIGURES_DIR, "pc1_pc3.png"), p3, 
       width = 10, height = 7, dpi = FIG_DPI)

cat("PC1 vs PC3 guardado: results/figures/pc1_pc3.png\n")

# PC2 vs PC3
p4 <- ggplot(data = pcadf, aes(x = PC2, y = PC3)) +
  geom_point(aes(color = group), size = 3) +
  theme_minimal() +
  labs(
    title = "PCA de Muestras MS: PC2 vs PC3",
    x = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
    y = paste0("PC3 (", round(variance_explained[3], 1), "%)")
  ) +
  scale_color_manual(values = unlist(COLOR_PALETTE))

ggsave(file.path(FIGURES_DIR, "pc2_pc3.png"), p4, 
       width = 10, height = 7, dpi = FIG_DPI)

cat("PC2 vs PC3 guardado: results/figures/pc2_pc3.png\n")


cat("PCA completado\n\n")

# Guardar resultados
saveRDS(pcaexpr, file.path(DATA_PROCESSED_DIR, "pca_result.rds"))
write.csv(pcadf, file.path(TABLES_DIR, "pca_coordinates.csv"), row.names = FALSE)

# Clustering
cat("Paso 4/5: Realizando clustering...\n")
cat("  Clustering por distancia euclidiana...\n")

# Calculamos distancias
distexpr <- dist(x = expr_t, method = "euclidean")

# Hacemos el clustering
hcexpr <- hclust(d = distexpr, method = "complete")
order <- hcexpr$order
hcexpr <- as.dendrogram(hcexpr)

# Graficar el clustering
png(file.path(FIGURES_DIR, "euclidean_simple.png"), 
    width = 14, height = 6, units = "in", res = FIG_DPI)
plot(hcexpr, main = "Clustering Euclidiano")
dev.off()

cat("  Dendograma guardado: results/figures/euclidean_simple.png\n")

# Colorear por "group"
colgroup <- unlist(COLOR_PALETTE)
names(colgroup) <- names(COLOR_PALETTE)
labels_colors(hcexpr) <- colgroup[sample_metadata$group][order]

png(file.path(FIGURES_DIR, "euclidean_colored.png"), 
    width = 14, height = 6, units = "in", res = FIG_DPI)
plot(hcexpr, main = "Clustering Euclidiano - Coloreado por Grupo")
dev.off()

cat("  Dendograma coloreado guardado: results/figures/euclidean_colored.png\n")

cat("  Clustering euclidiano completado\n")

cat("  Clustering por correlación...\n")

# Calcular las distancias por correlación
correlacion <- cor(expression_filtered)
distancia <- as.dist((1 - correlacion) / 2)

# Hacemos el clustering
hccor <- hclust(d = distancia, method = "complete")
ordercor <- hccor$order
hccor <- as.dendrogram(hccor)

# Graficar el clustering
png(file.path(FIGURES_DIR, "correlation_simple.png"), 
    width = 14, height = 6, units = "in", res = FIG_DPI)
plot(hccor, main = "Clustering por Correlación")
dev.off()

cat("  Dendograma guardado: results/figures/correlation_simple.png\n")

# Colorear por "group"
labels_colors(hccor) <- colgroup[sample_metadata$group][ordercor]

png(file.path(FIGURES_DIR, "correlation_colored.png"), 
    width = 14, height = 6, units = "in", res = FIG_DPI)
plot(hccor, main = "Clustering por Correlación - Coloreado por Grupo")
dev.off()

cat("  Dendograma coloreado guardado: results/figures/correlation_colored.png\n")

cat("  Clustering por correlación completado\n\n")

# Guardar resultados de clustering
saveRDS(hclust(distexpr, method = "complete"), 
        file.path(DATA_PROCESSED_DIR, "hclust_euclidean.rds"))
saveRDS(hclust(distancia, method = "complete"), 
        file.path(DATA_PROCESSED_DIR, "hclust_correlation.rds"))

# Resumen por grupos
cat("Paso 5/5: Generando resumen por grupos...\n")

# Calcular estadísticas por grupo
group_summary <- pcadf %>%
  group_by(group) %>%
  summarise(
    n_samples = n(),
    PC1_mean = mean(PC1),
    PC1_sd = sd(PC1),
    PC2_mean = mean(PC2),
    PC2_sd = sd(PC2),
    PC3_mean = mean(PC3),
    PC3_sd = sd(PC3)
  )

cat("\nResumen PCA por grupo:\n")
print(group_summary)
cat("\n")

write.csv(group_summary, file.path(TABLES_DIR, "pca_summary.csv"), 
          row.names = FALSE)

cat("Resumen por grupo guardado\n\n")


cat("Ejecuta para seguir: source('scripts/03_differential_expression/
    01_limma_analysis.R')\n")

# Guardar workspace
save.image(file.path(DATA_PROCESSED_DIR, "exploratory_workspace.RData"))