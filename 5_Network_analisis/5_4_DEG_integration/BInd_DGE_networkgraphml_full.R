############################################################
## SCRIPT: 5_2_5_Subnetwork_to_GraphML_with_DE.R
##
## Objetivo:
##   - Cargar subnetworks A7/A9
##   - Anotar nodos con:
##       * logFC
##       * adj.P.Val
##       * direction ("up"/"down")
##       * is_signif (adj.P.Val < 0.05)
##   - Exportar a GraphML para Cytoscape.
############################################################

setwd("/home/joelr/CESC_Network/4_Network_analisis/5_1_Subnetwork")

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(readr)
})

## =========================
## 1) Rutas de archivos
## =========================
# Subnetworks (ya los tienes en este directorio)
file_sub_A7 <- "5_2_4_Subnetwork_A7_max.rds"
file_sub_A9 <- "5_2_4_Subnetwork_A9_max.rds"

# RESULTADOS LIMMA
# OJO: AJUSTA ESTAS RUTAS/NOMBRES SI NO COINCIDEN CON TU PIPELINE
file_limma_A7 <- "/home/joelr/CESC_Network/4_Network_analisis/4_3_limma_A7_vs_SNT_allGenes_annot_PR_DEG.tsv"
file_limma_A9 <- "/home/joelr/CESC_Network/4_Network_analisis/4_3_limma_A9_vs_SNT_allGenes_annot_PR_DEG.tsv"

## =========================
## 2) Leer subnetworks
## =========================
g_A7 <- readRDS(file_sub_A7)
g_A9 <- readRDS(file_sub_A9)

## =========================
## 3) Leer resultados limma
## =========================
deg_A7 <- read_tsv(file_limma_A7, show_col_types = FALSE)
deg_A9 <- read_tsv(file_limma_A9, show_col_types = FALSE)

## =========================
## 4) Helper para detectar columnas
## =========================
detect_col <- function(df, candidates, label) {
  found <- intersect(candidates, colnames(df))
  if (length(found) == 0) {
    stop("No encontré ninguna columna de tipo ", label,
         ". Revisa los nombres en la tabla limma.")
  }
  found[1]
}

## =========================
## 5) Función para anotar un grafo
## =========================
annotate_graph_with_DE <- function(g, deg_tbl, label = "") {
  message("=== Anotando grafo ", label, " ===")
  
  # Detectar columnas clave
  gene_col  <- detect_col(deg_tbl,
                          c("SYMBOL", "Gene.symbol", "Symbol", "gene", "Gene"),
                          "GENE/SYMBOL")
  logfc_col <- detect_col(deg_tbl,
                          c("logFC", "log2FoldChange", "log_fc"),
                          "logFC")
  padj_col  <- detect_col(deg_tbl,
                          c("adj.P.Val", "padj", "FDR", "qvalue"),
                          "adj.P.Val / padj / FDR")
  
  message("Usando columnas: gene = ", gene_col,
          ", logFC = ", logfc_col,
          ", adj = ", padj_col)
  
  # Tabla reducida de DE
  deg_small <- deg_tbl %>%
    select(
      gene  = all_of(gene_col),
      logFC = all_of(logfc_col),
      adjP  = all_of(padj_col)
    ) %>%
    distinct(gene, .keep_all = TRUE)
  
  # Nombre de los nodos del grafo
  v_names <- igraph::vertex_attr(g, "name")
  if (is.null(v_names)) {
    stop("El grafo no tiene atributo de vértice 'name'.")
  }
  
  # Mapear DE a nodos
  df_map <- tibble(name = v_names) %>%
    left_join(deg_small, by = c("name" = "gene")) %>%
    mutate(
      direction = dplyr::case_when(
        is.na(logFC)        ~ NA_character_,
        logFC > 0           ~ "up",
        logFC < 0           ~ "down",
        TRUE                ~ "zero"
      ),
      # criterio de significancia: adjP < 0.05
      is_signif = !is.na(adjP) & adjP < 0.05
    )
  
  # Asignar atributos a los nodos
  V(g)$logFC     <- df_map$logFC
  V(g)$adjP      <- df_map$adjP
  V(g)$direction <- df_map$direction
  V(g)$is_signif <- df_map$is_signif
  
  g
}

## =========================
## 6) Anotar A7 y A9
## =========================
g_A7_annot <- annotate_graph_with_DE(g_A7, deg_A7, label = "A7")
g_A9_annot <- annotate_graph_with_DE(g_A9, deg_A9, label = "A9")

## =========================
## 7) Guardar como RDS y GraphML
## =========================
saveRDS(g_A7_annot, "5_2_4_Subnetwork_A7_max_annot_DGE.rds")
saveRDS(g_A9_annot, "5_2_4_Subnetwork_A9_max_annot_DGE.rds")

write_graph(g_A7_annot,
            file = "5_2_4_Subnetwork_A7_max_annot_DGE.graphml",
            format = "graphml")

write_graph(g_A9_annot,
            file = "5_2_4_Subnetwork_A9_max_annot_DGE.graphml",
            format = "graphml")

message("Listo. Se generaron los GraphML anotados para A7 y A9.")
