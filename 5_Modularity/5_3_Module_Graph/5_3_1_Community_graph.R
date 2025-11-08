### ========================
### Cargar librerías
### ========================
library(igraph)
library(tidyverse)

### ========================
### Cargar grafos originales
### ========================
graph_A7_louv <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A7_graph_max.rds")
graph_A9_louv <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A9_graph_max.rds")

graph_A7_info <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_1_1_Infomap_A7_graph_max.rds")
graph_A9_info <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_1_1_Infomap_A9_graph_max.rds")

### ========================
### FUNCION para crear metared
### ========================
build_community_network <- function(graph) {
  
  message("✔️ Processing graph with ", vcount(graph), " nodes and ", ecount(graph), " edges.")
  
  # Nodos con su comunidad asignada
  nodes <- data.frame(
    name = V(graph)$name,
    community = V(graph)$Community_name
  )
  
  # Edgelist original con las comunidades
  edges <- as.data.frame(get.edgelist(graph)) %>%
    rename(from = V1, to = V2) %>%
    left_join(nodes, by = c("from" = "name")) %>%
    rename(community_from = community) %>%
    left_join(nodes, by = c("to" = "name")) %>%
    rename(community_to = community)
  
  # Filtrar solo enlaces ENTRE comunidades diferentes y contar
  # Filtrar solo enlaces ENTRE comunidades diferentes y sumar
  community_edges <- edges %>%
    filter(community_from != community_to) %>%
    mutate(
      community_pair = ifelse(
        community_from < community_to,
        paste(community_from, community_to, sep = "|"),
        paste(community_to, community_from, sep = "|")
      )
    ) %>%
    group_by(community_pair) %>%
    summarise(
      community_from = first(community_from),
      community_to = first(community_to),
      weight = n(),
      .groups = "drop"
    ) %>%
    dplyr::select(-community_pair)
  # ======== Atributos de nodos (size y PageRank promedio)
  
  ## Community size
  community_size <- nodes %>% 
    count(community, name = "size")
  
  ## Normalized avg PageRank per community
  pagerank_values <- data.frame(
    name = V(graph)$name,
    pagerank = V(graph)$PageRank
  )
  
  community_pagerank <- nodes %>%
    left_join(pagerank_values, by = "name") %>%
    group_by(community) %>%
    summarise(
      avg_pagerank = mean(pagerank, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      normalized_pagerank = (avg_pagerank - min(avg_pagerank)) / (max(avg_pagerank) - min(avg_pagerank))
    )
  
  # ======== Armar lista de nodos para grafo
  vertices <- community_size %>%
    left_join(community_pagerank, by = c("community")) %>%
    rename(
      name = community,
      size = size,
      avg_pagerank = avg_pagerank,
      normalized_pagerank = normalized_pagerank
    )
  
  # ======== Crear grafo de comunidades
  community_graph <- graph_from_data_frame(
    d = community_edges,
    vertices = vertices,
    directed = FALSE
  )
  
  message("✔️ Built community-level graph with ", vcount(community_graph), " nodes and ", ecount(community_graph), " edges.")
  
  return(community_graph)
}

### ========================
### Construir redes de comunidades
### ========================
community_graph_A7_louv <- build_community_network(graph_A7_louv)
community_graph_A9_louv <- build_community_network(graph_A9_louv)

community_graph_A7_info <- build_community_network(graph_A7_info)
community_graph_A9_info <- build_community_network(graph_A9_info)

### ========================
### Guardar en RDS
### ========================
saveRDS(community_graph_A7_louv, "~/CESC_Network/5_Modularity/5_3_Module_Graph/5_3_1_1_Community_Graph_A7_louv_max.rds")
saveRDS(community_graph_A9_louv, "~/CESC_Network/5_Modularity/5_3_Module_Graph/5_3_1_1_Community_Graph_A9_louv_max.rds")

saveRDS(community_graph_A7_info, "~/CESC_Network/5_Modularity/5_3_Module_Graph/5_3_1_1_Community_Graph_A7_info_max.rds")
saveRDS(community_graph_A9_info, "~/CESC_Network/5_Modularity/5_3_Module_Graph/5_3_1_1_Community_Graph_A9_info_max.rds")

### ========================
### Exportar en GraphML
### ========================
write_graph(community_graph_A7_louv, "/STORAGE/csbig/jruiz/cesc_network/5_3_1_1_Community_Graph_A7_louv_max.graphml", format = "graphml")
write_graph(community_graph_A9_louv, "/STORAGE/csbig/jruiz/cesc_network/5_3_1_1_Community_Graph_A9_louv_max.graphml", format = "graphml")

write_graph(community_graph_A7_info, "/STORAGE/csbig/jruiz/cesc_network/5_3_1_1_Community_Graph_A7_info_max.graphml", format = "graphml")
write_graph(community_graph_A9_info, "/STORAGE/csbig/jruiz/cesc_network/5_3_1_1_Community_Graph_A9_info_max.graphml", format = "graphml")


message("✅ All done!")
