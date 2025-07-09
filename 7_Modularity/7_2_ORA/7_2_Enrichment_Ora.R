### Cargar librerías ###
start_time <- Sys.time()

library(igraph)
library(dplyr)
library(clusterProfiler)
library(vroom)
library(org.Hs.eg.db)
library(ggplot2)
library(purrr)
library(future.apply)
library(readr)

plan(multicore, workers = 3)  # Configura el uso de 3 cores
set.seed(123)

### Cargar redes ###
graph_A7 <- readRDS(file = "~/CESC_Network/7_Modularity/7_1_1_Louvain_A7_graph.rds")
graph_A9 <- readRDS(file = "~/CESC_Network/7_Modularity/7_1_1_Louvain_A9_graph.rds")

### Cargar universo de genes ###
universo <- read_rds("~/CESC_Network/2_Prepro_TCGA_GTEx/2_10_Universe_annotation_table.rds") %>% pull(hgnc_symbol)

### Función para realizar ORA ###
perform_ora <- function(genes, universe, community_id, community_name, ontology) {
  if (length(genes) == 0) return(data.frame())
  
  enrich_results <- tryCatch({
    enrichGO(
      gene = genes, 
      universe = universe, 
      OrgDb = org.Hs.eg.db, 
      keyType = "SYMBOL", 
      ont = ontology, 
      pAdjustMethod = "fdr", 
      pvalueCutoff = 0.05
    ) %>% as.data.frame()
  }, error = function(e) data.frame())
  
  if (nrow(enrich_results) == 0) return(enrich_results)
  
  enrich_results <- enrich_results %>%
    mutate(Community = community_id, Community_Name = community_name, Ontology = ontology)
  
  return(enrich_results)
}

### NUEVA versión robusta de extract_communities_and_ora ###
extract_communities_and_ora <- function(graph, universe) {
  communities <- unique(V(graph)$Community)
  ontologies <- c("BP", "CC", "MF")
  
  ora_results <- future_lapply(communities, function(community) {
    # Genes del módulo
    genes <- V(graph)[V(graph)$Community == community]$name
    community_name <- unique(V(graph)[V(graph)$Community == community]$Community_name)
    
    # Calcular atributos extra
    community_size <- length(genes)
    avg_pagerank <- mean(V(graph)[V(graph)$Community == community]$PageRank, na.rm = TRUE)
    
    # ORA para las tres ontologías
    enrich_df <- bind_rows(lapply(ontologies, function(ontology) {
      perform_ora(genes, universe, community, community_name, ontology)
    }))
    
    # Si no hubo enriquecimiento, devolver vacío
    if (nrow(enrich_df) == 0) return(enrich_df)
    
    # Añadir columnas extra
    enrich_df <- enrich_df %>%
      mutate(
        Community_size = community_size,
        Avg_PageRank = avg_pagerank
      )
    
    return(enrich_df)
  }, future.seed = TRUE) %>% bind_rows()
  
  return(ora_results)
}


### Ejecutar ORA para cada red ###
ora_A7 <- extract_communities_and_ora(graph_A7, universo)
ora_A9 <- extract_communities_and_ora(graph_A9, universo)

print(length(unique(ora_A7$Community_Name))) # número de módulos A7
print(length(unique(ora_A9$Community_Name))) # número de módulos A9

### Filtrar por mínimo de genes enriquecidos ###
ora_A7 <- ora_A7 %>% filter(Count >= 3)
ora_A9 <- ora_A9 %>% filter(Count >= 3)

print(length(unique(ora_A7$Community_Name)))
print(length(unique(ora_A9$Community_Name)))

### Guardar resultados ###
write.table(ora_A7, file = "~/CESC_Network/7_Modularity/7_2_1_Ora7_600k.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ora_A9, file = "~/CESC_Network/7_Modularity/7_2_1_Ora9_600k.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

end_time <- Sys.time()
cat("Tiempo total:", round(end_time - start_time, 2), "minutos\n")
