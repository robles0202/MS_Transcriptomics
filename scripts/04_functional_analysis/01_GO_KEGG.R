# Este script lleva a cabo un enrriquecimiento funcional mediante las bases de 
# datos GO y KEGG para los DEGs 

# Cargar configuración
source("scripts/config.R")

# Cargar librerías
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

set.seed(RANDOM_SEED)

# Cargar datos y anotación
cat("Paso 1/6: Cargando datos...\n")

all_results <- readRDS(file.path(DATA_PROCESSED_DIR, "limma_results.rds"))
# Cargar anotación de genes
gene_annotation <- read.csv(file.path(DATA_RAW_DIR, "annotation.csv"), 
                            stringsAsFactors = FALSE)

# Verificar columnas disponibles
cat("Columnas en anotación:\n")
print(names(gene_annotation))

cat("\n")
cat("Datos cargados\n\n")

# Mapear ids de Illumina a Entrez
cat("Paso 2/6: Mapeando IDs de genes...\n")

# Función para mapear Illumina IDs a Entrez IDs
map_to_entrez <- function(deg_results) {
  degs <- deg_results %>%
    filter(regulation != "NS")
  annotation_subset <- gene_annotation[, c("probe_id", "gene_symbol", "gene_id")]
  degs <- degs %>%
    left_join(annotation_subset, 
              by = c("gene_id" = "probe_id"))
  degs <- degs %>%
    mutate(entrez = as.character(gene_id.y)) %>%
    filter(!is.na(entrez) & entrez != "" & entrez != "NA")
  
  return(degs)
}

# Mapear para cada subtipo
rrms_mapped <- map_to_entrez(all_results[["RRMS_vs_Control"]])
ppms_mapped <- map_to_entrez(all_results[["PPMS_vs_Control"]])
spms_mapped <- map_to_entrez(all_results[["SPMS_vs_Control"]])

cat("Genes mapeados a Entrez ID:\n")
cat("  RRMS:", nrow(rrms_mapped), "genes\n")
cat("  PPMS:", nrow(ppms_mapped), "genes\n")
cat("  SPMS:", nrow(spms_mapped), "genes\n\n")

cat("Mapeo completado\n\n")

# Analisis GO - RRMS
cat("Paso 3/6: Realizando análisis GO...\n")

run_go_analysis <- function(genes_df, subtype_name) {

  cat("  Analizando", subtype_name, "(", nrow(genes_df), "genes)...\n")
  
  entrez_ids <- unique(genes_df$entrez)
  
  # GO Biological Process
  go_bp <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  
  # GO Molecular Function
  go_mf <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  
  # GO Cellular Component
  go_cc <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  
  return(list(BP = go_bp, MF = go_mf, CC = go_cc))
}

# Análisis GO para RRMS
go_rrms <- run_go_analysis(rrms_mapped, "RRMS")

write.csv(as.data.frame(go_rrms$BP), 
          file.path(TABLES_DIR, "go_bp_rrms.csv"), 
          row.names = FALSE)
    
p <- dotplot(go_rrms$BP, showCategory = 15) + 
  ggtitle("GO Biological Process - RRMS") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(FIGURES_DIR, "go_bp_rrms.png"), 
       p, width = 10, height = 8, dpi = FIG_DPI)
  
write.csv(as.data.frame(go_rrms$MF), 
          file.path(TABLES_DIR, "go_mf_rrms.csv"), 
          row.names = FALSE)
  
write.csv(as.data.frame(go_rrms$CC), 
          file.path(TABLES_DIR, "go_cc_rrms.csv"), 
          row.names = FALSE)

# Análisis GO para PPMS
go_ppms <- run_go_analysis(ppms_mapped, "PPMS")
  
write.csv(as.data.frame(go_ppms$BP), 
          file.path(TABLES_DIR, "go_bp_ppms.csv"), 
          row.names = FALSE)

p <- dotplot(go_ppms$BP, showCategory = 15) + 
  ggtitle("GO Biological Process - PPMS") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(FIGURES_DIR, "go_bp_ppms.png"), 
       p, width = 10, height = 8, dpi = FIG_DPI)

write.csv(as.data.frame(go_ppms$MF), 
          file.path(TABLES_DIR, "go_mf_ppms.csv"), 
          row.names = FALSE)

write.csv(as.data.frame(go_ppms$CC), 
          file.path(TABLES_DIR, "go_cc_ppms.csv"), 
          row.names = FALSE)

cat(" Análisis GO completado\n\n")


# Análisis KEGG
cat("Paso 4/6: Realizando análisis KEGG...\n")

run_kegg_analysis <- function(genes_df, subtype_name) {
  
  cat("  Analizando", subtype_name, "(", nrow(genes_df), "genes)...\n")
  
  entrez_ids <- unique(genes_df$entrez)
  
  # KEGG pathway enrichment
  kegg <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )
  
  return(kegg)
}

# Análisis KEGG para RRMS
kegg_rrms <- run_kegg_analysis(rrms_mapped, "RRMS")

if(!is.null(kegg_rrms) && nrow(as.data.frame(kegg_rrms)) > 0) {
  
  write.csv(as.data.frame(kegg_rrms), 
            file.path(TABLES_DIR, "kegg_rrms.csv"), 
            row.names = FALSE)
  # Visualizar
  p <- dotplot(kegg_rrms, showCategory = 15) + 
    ggtitle("KEGG Pathways - RRMS") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(FIGURES_DIR, "kegg_rrms.png"), 
         p, width = 10, height = 8, dpi = FIG_DPI)
}

# Análisis KEGG para PPMS
kegg_ppms <- run_kegg_analysis(ppms_mapped, "PPMS")

if(!is.null(kegg_ppms) && nrow(as.data.frame(kegg_ppms)) > 0) {
  
  write.csv(as.data.frame(kegg_ppms), 
            file.path(TABLES_DIR, "kegg_ppms.csv"), 
            row.names = FALSE)
      
  # Visualizar
  p <- dotplot(kegg_ppms, showCategory = 15) + 
    ggtitle("KEGG Pathways - PPMS") +
    theme(plot.title = element_text(hjust = 0.5))
      
  ggsave(file.path(FIGURES_DIR, "kegg_ppms.png"), 
         p, width = 10, height = 8, dpi = FIG_DPI)
}

cat("Análisis KEGG completado\n\n")

# Comparación GO entre RRMS y PPMS
cat("Paso 5/6: Comparando enriquecimientos entre subtipos...\n")

go_comparison <- list(RRMS = go_rrms$BP, PPMS = go_ppms$BP)

comparison_plot <- compareCluster(
  geneClusters = list(
    RRMS = unique(rrms_mapped$entrez),
    PPMS = unique(ppms_mapped$entrez)
    ),
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
  )
  
p <- dotplot(comparison_plot, showCategory = 10) + 
  ggtitle("Comparación GO BP: RRMS vs PPMS") +
  theme(plot.title = element_text(hjust = 0.5))
      
ggsave(file.path(FIGURES_DIR, "go_comparison_rrms_ppms.png"), 
       p, width = 12, height = 8, dpi = FIG_DPI)

cat("Comparación completada\n\n")

# Resumen
cat("Paso 6/6: Generando resumen...\n")

summary_list <- list()

summary_list$RRMS_GO_BP <- nrow(as.data.frame(go_rrms$BP))
summary_list$PPMS_GO_BP <- nrow(as.data.frame(go_ppms$BP))
summary_list$PPMS_KEGG <- nrow(as.data.frame(kegg_ppms))

summary_df <- data.frame(
  Categoria = names(summary_list),
  Terminos_Enriquecidos = unlist(summary_list)
)

cat("\n")
print(summary_df)
cat("\n")
cat("Archivos generados:\n")
cat("  - Tablas:\n")
cat("    go_bp_rrms.csv\n")
cat("    go_mf_rrms.csv\n")
cat("    go_cc_rrms.csv\n")
cat("    go_bp_ppms.csv\n")
cat("   go_mf_ppms.csv\n")
cat("    go_cc_ppms.csv\n")
cat("    kegg_ppms.csv\n")
cat("  - Figuras:\n")
cat("    go_bp_rrms.png\n")
cat("    go_bp_ppms.png\n")
cat("    kegg_ppms.png\n")
cat("    go_comparison_rrms_ppms.png\n")
cat("\n")

cat("Para continuar ejecuta: source('scripts/04_functional_analysis/02_network_analysis.R')\n")