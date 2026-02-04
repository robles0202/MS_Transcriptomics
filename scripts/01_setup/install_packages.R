# Este script instala todos los paquetes necesarios para el análisis

cat("Paso 1/4: Instalando BiocManager...\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  cat("BiocManager instalado\n")
} else {
  cat("BiocManager ya está instalado\n")
}

# Usar BiocManager e la versión 3.22
BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)

cat("\nPaso 2/4: Instalando paquetes de Bioconductor...\n")

bioc_packages <- c(
  "GEOquery",
  "Biobase",
  "limma",
  "biomaRt",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "clusterProfiler",
  "enrichplot",
  "ReactomePA",
  "DOSE",
  "pathview",
  "ComplexHeatmap",
  "EnhancedVolcano"
)

for (x in bioc_packages) {
  if (!requireNamespace(x, quietly = TRUE)) {
    BiocManager::install(x, ask = FALSE, update = FALSE)
  } else {
    cat(x, "ya está instalado\n")
  }
}

cat("\nPaso 3/4: Instalando paquetes de CRAN...\n")

cran_packages <- c(
  "tidyverse",        
  "dplyr",            
  "tidyr",            
  "readr",            
  "data.table",       
  "ggplot2",          
  "pheatmap",         
  "ggrepel",          
  "VennDiagram",      
  "UpSetR",           
  "RColorBrewer",     
  "viridis",          
  "gridExtra",        
  "cowplot",          
  "ggpubr",           
  "igraph",           
  "STRINGdb",         
  "here",             
  "knitr",            
  "rmarkdown",        
  "openxlsx",         
  "writexl"
)

for (x in cran_packages) {
  if (!requireNamespace(x, quietly = TRUE)) {
    install.packages(x, repos = "https://cloud.r-project.org", 
                     dependencies = TRUE, quiet = TRUE)
  } else {
    cat(x, "ya está instalado\n")
  }
}

cat("\nPaso 4/4: Verificando instalación...\n\n")

all_packages <- c(bioc_packages, cran_packages)
installed <- sapply(all_packages, requireNamespace, quietly = TRUE)

if (all(installed)) {
  cat("Todos los paquetes se han instalado correctamente.\n")
  cat("Para continuar con la descarga de datos ejecuta: source('scripts/01_setup/download_data.R')\n")
} else {
  cat("Atención: Algunos paquetes no se instalaron correctamente:\n")
  failed <- names(installed[!installed])
  for (x in failed) {
    cat(x, "\n")
  }
  cat("\nProbar instalación manual:\n")
  cat("BiocManager::install(c('", paste(failed, collapse = "', '"), "'))\n")
}