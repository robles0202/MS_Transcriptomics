# Este script permite identificar genes específicos de cada subtipo de MS 

# Cargar configuración
source("scripts/config.R")

# Cargar librerías
library(tidyverse)

# Cargar resultados de expresión diferencial
cat("Paso 1/4: Cargando resultados de limma...\n")

all_results <- readRDS(file.path(DATA_PROCESSED_DIR, "limma_results.rds"))

cat("Resultados cargados\n\n")

# Extraer DEGs
cat("Paso 2/4: Extrayendo genes diferencialmente expresados...\n")

# Crear listas de genes DEGs para cada comparación relevante
degs_rrms <- all_results[["RRMS_vs_Control"]] %>%
  filter(regulation != "NS") %>%
  pull(gene_id)

degs_spms <- all_results[["SPMS_vs_Control"]] %>%
  filter(regulation != "NS") %>%
  pull(gene_id)

degs_ppms <- all_results[["PPMS_vs_Control"]] %>%
  filter(regulation != "NS") %>%
  pull(gene_id)

cat("DEGs en RRMS vs Control:", length(degs_rrms), "\n")
cat("DEGs en SPMS vs Control:", length(degs_spms), "\n")
cat("DEGs en PPMS vs Control:", length(degs_ppms), "\n\n")

# Identificar genes
cat("Paso 3/4: Identificando genes específicos y compartidos...\n")

# Genes específicos de cada subtipo
rrms_specific <- setdiff(degs_rrms, union(degs_spms, degs_ppms))
spms_specific <- setdiff(degs_spms, union(degs_rrms, degs_ppms))
ppms_specific <- setdiff(degs_ppms, union(degs_rrms, degs_spms))

# Genes compartidos
all_three <- Reduce(intersect, list(degs_rrms, degs_spms, degs_ppms))
rrms_spms <- setdiff(intersect(degs_rrms, degs_spms), degs_ppms)
rrms_ppms <- setdiff(intersect(degs_rrms, degs_ppms), degs_spms)
spms_ppms <- setdiff(intersect(degs_spms, degs_ppms), degs_rrms)

# Genes de formas progresivas (SPMS + PPMS)
progressive_specific <- setdiff(union(degs_spms, degs_ppms), degs_rrms)

cat("Genes identificados:\n")
cat("RRMS específico:", length(rrms_specific), "\n")
cat("SPMS específico:", length(spms_specific), "\n")
cat("PPMS específico:", length(ppms_specific), "\n")
cat("Compartidos entre los tres subtipos:", length(all_three), "\n")
cat("RRMS + SPMS:", length(rrms_spms), "\n")
cat("RRMS + PPMS:", length(rrms_ppms), "\n")
cat("SPMS + PPMS:", length(spms_ppms), "\n")
cat("Progresivos específicos:", length(progressive_specific), "\n\n")


# Crear tablas
cat("Paso 4/4: Generando tablas con información detallada...\n")

# Función para obtener información de genes
get_gene_info <- function(gene_list, contrast_name) {
  if (length(gene_list) == 0) {
    return(data.frame())
  }
  
  all_results[[contrast_name]] %>%
    filter(gene_id %in% gene_list) %>%
    select(gene_id, logFC, AveExpr, t, P.Value, adj.P.Val, regulation)
}

# RRMS específico
if (length(rrms_specific) > 0) {
  rrms_specific_df <- get_gene_info(rrms_specific, "RRMS_vs_Control")
  write.csv(rrms_specific_df, 
            file.path(TABLES_DIR, "genes_rrms_specific.csv"), 
            row.names = FALSE)
  cat("Tabla RRMS específico guardada\n")
}

# SPMS específico (no generada dado que no hay)
# if (length(spms_specific) > 0) {
#   spms_specific_df <- get_gene_info(spms_specific, "SPMS_vs_Control")
#   write.csv(spms_specific_df, 
#             file.path(TABLES_DIR, "genes_spms_specific.csv"), 
#             row.names = FALSE)
#   cat("Tabla SPMS específico guardada\n")
# }

# PPMS específico
if (length(ppms_specific) > 0) {
  ppms_specific_df <- get_gene_info(ppms_specific, "PPMS_vs_Control")
  write.csv(ppms_specific_df, 
            file.path(TABLES_DIR, "genes_ppms_specific.csv"), 
            row.names = FALSE)
  cat("Tabla PPMS específico guardada\n")
}

# Genes compartidos por todos(tampoco se genera)
# if (length(all_three) > 0) {
#   all_three_df <- get_gene_info(all_three, "RRMS_vs_Control")
#   write.csv(all_three_df, 
#             file.path(TABLES_DIR, "genes_shared_all_subtypes.csv"), 
#             row.names = FALSE)
#   cat("Tabla genes compartidos guardada\n")
# }

# Genes de formas progresivas (no se genera dado que es la misma que PPMS)
# if (length(progressive_specific) > 0) {
#   progressive_df <- get_gene_info(progressive_specific, "SPMS_vs_Control")
#   write.csv(progressive_df, 
#             file.path(TABLES_DIR, "genes_progressive_specific.csv"), 
#             row.names = FALSE)
#   cat("Tabla genes progresivos guardada\n")
# }

# Crear tabla resumen
summary_subtype <- data.frame(
  category = c("RRMS específico", "SPMS específico", "PPMS específico",
               "Compartido (todos)", "RRMS + SPMS", "RRMS + PPMS", 
               "SPMS + PPMS", "Progresivos específicos"),
  n_genes = c(length(rrms_specific), length(spms_specific), length(ppms_specific),
              length(all_three), length(rrms_spms), length(rrms_ppms),
              length(spms_ppms), length(progressive_specific))
)

write.csv(summary_subtype, 
          file.path(TABLES_DIR, "subtype_specific_genes_summary.csv"), 
          row.names = FALSE)

# Guardar listas de genes
gene_lists <- list(
  rrms_specific = rrms_specific,
  spms_specific = spms_specific,
  ppms_specific = ppms_specific,
  all_three = all_three,
  rrms_spms = rrms_spms,
  rrms_ppms = rrms_ppms,
  spms_ppms = spms_ppms,
  progressive_specific = progressive_specific
)

saveRDS(gene_lists, file.path(DATA_PROCESSED_DIR, "subtype_specific_gene_lists.rds"))

# RESUMEN FINAL
cat("Resumen:\n")
print(summary_subtype)
cat("\n")
cat("Archivos generados:\n")
cat("  - genes_rrms_specific.csv\n")
cat("  - genes_ppms_specific.csv\n")
cat("  - subtype_specific_gene_lists.rds\n")
cat("\n")
cat("Ejecuta pra seguir: source('scripts/04_functional_analysis/01_GO_KEGG.R')\n")