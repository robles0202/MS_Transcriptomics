# Este script se encarga de la generación de  volcano-pots

# Cargar configuración
source("scripts/config.R")

# Cargar librerías
library(tidyverse)
library(ggplot2)
library(ggrepel)

# Cargar datos
cat("Paso 1/4: Cargando resultados de expresión diferencial...\n")

all_results <- readRDS(file.path(DATA_PROCESSED_DIR, "limma_results.rds"))

cat("Resultados cargados\n\n")

# Crear volcano plots
cat("Paso 2/4: Creando volcano-plots...\n")

volcano_plot <- function(results_df, contrast_name, top_n = 20) {
  
  # Preparar datos
  plot_data <- results_df %>%
    mutate(
      neg_log10_pval = -log10(adj.P.Val),
      significant = regulation != "NS",
      label = ifelse(
        rank(adj.P.Val) <= top_n & regulation != "NS",
        gene_id,
        ""
      )
    )
  
  # Contar genes
  n_up <- sum(results_df$regulation == "UP")
  n_down <- sum(results_df$regulation == "DOWN")
  
  # Crear plot
  p <- ggplot(plot_data, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(
      aes(color = regulation),
      alpha = 0.6,
      size = 1.5
    ) +
    scale_color_manual(
      values = c(
        "UP" = "#D62728",
        "DOWN" = "#1F77B4",
        "NS" = "#CCCCCC"
      ),
      labels = c(
        "UP" = paste0("UP (", n_up, ")"),
        "DOWN" = paste0("DOWN (", n_down, ")"),
        "NS" = "NS"
      )
    ) +
    geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), 
               linetype = "dashed", color = "gray30") +
    geom_hline(yintercept = -log10(FDR_THRESHOLD), 
               linetype = "dashed", color = "gray30") +
    geom_text_repel(
      aes(label = label),
      size = 2.5,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.3
    ) +
    labs(
      title = paste("Volcano Plot:", contrast_name),
      subtitle = paste0(
        "Criterios: |log2FC| > ", LOG2FC_THRESHOLD, 
        ", FDR < ", FDR_THRESHOLD
      ),
      x = "log2 Fold Change",
      y = "-log10 (FDR)",
      color = "Regulación"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Generar volcano-plots
cat("Paso 3/4: Generando volcano plots...\n\n")

control_contrasts <- c("RRMS_vs_Control", "SPMS_vs_Control", "PPMS_vs_Control")

for (contrast_name in control_contrasts) {
  
  cat("Generando:", contrast_name, "\n")
  
  # Crear plot
  p <- volcano_plot(all_results[[contrast_name]], contrast_name)
  
  # Guardar
  filename <- paste0("volcano_", gsub(" ", "_", contrast_name), ".png")
  ggsave(
    file.path(FIGURES_DIR, filename),
    p,
    width = 10,
    height = 8,
    dpi = FIG_DPI
  )
  
  cat("Guardado:", filename, "\n")
}

cat("\n")

# Crear panel de volcano-plots
cat("Paso 4/4: Creando panel de comparaciones vs Control...\n")

library(cowplot)

p1 <- volcano_plot(all_results[["RRMS_vs_Control"]], "RRMS vs Control")
p2 <- volcano_plot(all_results[["SPMS_vs_Control"]], "SPMS vs Control")
p3 <- volcano_plot(all_results[["PPMS_vs_Control"]], "PPMS vs Control")

combined_panel <- plot_grid(
  p1, p2, p3,
  ncol = 1,
  labels = c("A", "B", "C"),
  label_size = 14
)

ggsave(
  file.path(FIGURES_DIR, "volcano_combined_vs_control.png"),
  combined_panel,
  width = 12,
  height = 18,
  dpi = FIG_DPI
)

cat("Panel combinado guardado: volcano_combined_vs_control.png\n\n")


cat("Figuras generadas:\n")
cat("  - volcano_RRMS_vs_Control.png\n")
cat("  - volcano_SPMS_vs_Control.png\n")
cat("  - volcano_PPMS_vs_Control.png\n")
cat("  - volcano_combined_vs_control.png (panel combinado)\n")
cat("\n")

cat("Para continuar ejecuta: source('scripts/05_visualization/03_heatmaps.R')\n")