# Este script genera diagramas de Venn para visualizar el overlap
# entre genes diferencialmente expresados de diferentes subtipos

# Cargar configuración
source("scripts/config.R")

# Cargar librerías
library(tidyverse)
library(VennDiagram)
library(RColorBrewer)
library(grid)

set.seed(RANDOM_SEED)

# Cargar datos
cat("Paso 1/5: Cargando resultados de expresión diferencial...\n")

all_results <- readRDS(file.path(DATA_PROCESSED_DIR, "limma_results.rds"))

cat("Resultados cargados\n\n")

# Extraer DEGs por subtipo
cat("Paso 2/5: Extrayendo genes diferencialmente expresados...\n")

# Extraer genes DEG de cada subtipo vs Control
rrms_genes <- all_results[["RRMS_vs_Control"]] %>% 
  filter(regulation != "NS") %>% 
  pull(gene_id)

ppms_genes <- all_results[["PPMS_vs_Control"]] %>% 
  filter(regulation != "NS") %>% 
  pull(gene_id)

cat("Genes por subtipo:\n")
cat("  RRMS:", length(rrms_genes), "\n")
cat("  PPMS:", length(ppms_genes), "\n")

cat("Genes extraídos\n\n")

# Crear diagramas de Venn
cat("Paso 3/5: Creando diagramas de Venn...\n")

# Crear lista de genes para Venn
gene_lists <- list()
gene_lists[["RRMS"]] <- rrms_genes
gene_lists[["PPMS"]] <- ppms_genes

if(length(gene_lists) >= 2) {
  
  cat("  Generando diagrama con", length(gene_lists), "grupos...\n")
  
  # Configurar colores
  if(length(gene_lists) == 2) {
    venn_colors <- c("red2", "green3")
    cat_dist_val <- c(0.05, 0.05)
  } else if(length(gene_lists) == 3) {
    venn_colors <- c("red2", "green3", "blue2")
    cat_dist_val <- c(0.05, 0.05, 0.05)
  }
  
  # Crear Venn diagram principal
  venn.plot <- venn.diagram(
    x = gene_lists,
    filename = NULL,
    category.names = names(gene_lists),
    fill = venn_colors,
    alpha = 0.5,
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.dist = cat_dist_val,
    cat.fontfamily = "sans",
    main = "Genes Diferencialmente Expresados por Subtipo",
    main.cex = 2,
    main.fontface = "bold",
    height = 3000,
    width = 3000,
    resolution = 300,
    compression = "lzw",
    lwd = 2,
    lty = 'blank'
  )
  
  png(file.path(FIGURES_DIR, "venn_diagram_all_subtypes.png"),
      width = 10, height = 8, units = "in", res = FIG_DPI)
  grid.draw(venn.plot)
  dev.off()
  
  cat("Diagrama guardado: venn_diagram_all_subtypes.png\n")
  
}

# Venn para up/down regulated
cat("Paso 4/5: Creando diagramas para up/down-regulated...\n")

# Extraer genes UP
rrms_up <- all_results[["RRMS_vs_Control"]] %>% 
  filter(regulation == "UP") %>% 
  pull(gene_id)

ppms_up <- all_results[["PPMS_vs_Control"]] %>% 
  filter(regulation == "UP") %>% 
  pull(gene_id)

# Extraer genes DOWN
rrms_down <- all_results[["RRMS_vs_Control"]] %>% 
  filter(regulation == "DOWN") %>% 
  pull(gene_id)

ppms_down <- all_results[["PPMS_vs_Control"]] %>% 
  filter(regulation == "DOWN") %>% 
  pull(gene_id)

# Venn para genes UP
if(length(rrms_up) > 0 && length(ppms_up) > 0) {
  venn_up <- venn.diagram(
    x = list(RRMS = rrms_up, PPMS = ppms_up),
    filename = NULL,
    category.names = c("RRMS", "PPMS"),
    fill = c("red2", "green3"),
    alpha = 0.5,
    cex = 2,
    fontface = "bold",
    cat.cex = 2,
    cat.fontface = "bold",
    main = "Genes Sobreexpresados (UP)",
    main.cex = 2,
    lwd = 2,
    lty = 'blank'
  )
  
  png(file.path(FIGURES_DIR, "venn_diagram_upregulated.png"),
      width = 8, height = 6, units = "in", res = FIG_DPI)
  grid.draw(venn_up)
  dev.off()
  
  cat("Diagrama guardado: venn_diagram_upregulated.png\n")
}

# Venn para genes DOWN
if(length(rrms_down) > 0 && length(ppms_down) > 0) {
  venn_down <- venn.diagram(
    x = list(RRMS = rrms_down, PPMS = ppms_down),
    filename = NULL,
    category.names = c("RRMS", "PPMS"),
    fill = c("red2", "green3"),
    alpha = 0.5,
    cex = 2,
    fontface = "bold",
    cat.cex = 2,
    cat.fontface = "bold",
    main = "Genes Infraexpresados (DOWN)",
    main.cex = 2,
    lwd = 2,
    lty = 'blank'
  )
  
  png(file.path(FIGURES_DIR, "venn_diagram_downregulated.png"),
      width = 8, height = 6, units = "in", res = FIG_DPI)
  grid.draw(venn_down)
  dev.off()
  
  cat("Diagrama guardado: venn_diagram_downregulated.png\n")
}


cat("Archivos generados:\n")
cat("Figuras:\n")
cat("  - venn_diagram_all_subtypes.png\n")
cat("  - venn_diagram_upregulated.png\n")
cat("  - venn_diagram_downregulated.png\n")
cat("\n")

cat("Para continuar ejecuta: source('scripts/05_visualization/02_volcano_plots.R')\n")