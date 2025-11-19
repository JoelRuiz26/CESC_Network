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
      normalized_pagerank = (avg_pagerank - min(avg_pagerank)) /
        (max(avg_pagerank) - min(avg_pagerank))
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
  
  message("✔️ Built community-level graph with ", vcount(community_graph),
          " nodes and ", ecount(community_graph), " edges.")
  
  return(community_graph)
}

### ========================
### Construir redes de comunidades
### ========================
community_graph_A7_louv <- build_community_network(graph_A7_louv)
#✔️ Processing graph with 10846 nodes and 443707 edges.
#✔️ Built community-level graph with 20 nodes and 190 edges.

community_graph_A9_louv <- build_community_network(graph_A9_louv)
#✔️ Processing graph with 10846 nodes and 443707 edges.
#✔️ Built community-level graph with 34 nodes and 561 edges.

community_graph_A7_info <- build_community_network(graph_A7_info)
community_graph_A9_info <- build_community_network(graph_A9_info)

### ========================
### NUEVO: función para tabla tipo edgelist + atributos
### ========================
build_edge_df <- function(g) {
  # edgelist con peso (weight = # de enlaces entre módulos)
  edges <- igraph::as_data_frame(g, what = "edges")  # from, to, weight
  
  # atributos de cada módulo (nodo)
  verts <- igraph::as_data_frame(g, what = "vertices") %>%
    dplyr::select(name, size, avg_pagerank, normalized_pagerank)
  
  edges %>%
    left_join(
      verts %>%
        rename(
          from = name,
          size_from = size,
          avg_pagerank_from = avg_pagerank,
          norm_pagerank_from = normalized_pagerank
        ),
      by = "from"
    ) %>%
    left_join(
      verts %>%
        rename(
          to = name,
          size_to = size,
          avg_pagerank_to = avg_pagerank,
          norm_pagerank_to = normalized_pagerank
        ),
      by = "to"
    )
}

### ========================
### NUEVO: construir las tablas edgelist + atributos
### ========================
community_A7_louv_df  <- build_edge_df(community_graph_A7_louv)
community_A9_louv_df  <- build_edge_df(community_graph_A9_louv)
community_A7_info_df  <- build_edge_df(community_graph_A7_info)
community_A9_info_df  <- build_edge_df(community_graph_A9_info)


### ========================
### NUEVO: guardar las tablas como RDS (y luego si quieres TSV)
### ========================
saveRDS(community_A7_louv_df,
        "~/CESC_Network/5_Modularity/5_3_Module_Graph/5_3_1_1_Community_Graph_A7_louv_edgelist_df.rds")
saveRDS(community_A9_louv_df,
        "~/CESC_Network/5_Modularity/5_3_Module_Graph/5_3_1_1_Community_Graph_A9_louv_edgelist_df.rds")

saveRDS(community_A7_info_df,
        "~/CESC_Network/5_Modularity/5_3_Module_Graph/5_3_1_1_Community_Graph_A7_info_edgelist_df.rds")
saveRDS(community_A9_info_df,
        "~/CESC_Network/5_Modularity/5_3_Module_Graph/5_3_1_1_Community_Graph_A9_info_edgelist_df.rds")

### (Opcional) si luego quieres TSV:
# readr::write_tsv(community_A7_louv_df, "A7_louv_edgelist_attrib.tsv")

### ========================
### Guardar en RDS (grafos, como ya tenías)
### ========================
saveRDS(community_graph_A7_louv, "~/CESC_Network/5_Modularity/5_3_Module_Graph/5_3_1_1_Community_Graph_A7_louv_max.rds")
saveRDS(community_graph_A9_louv, "~/CESC_Network/5_Modularity/5_3_Module_Graph/5_3_1_1_Community_Graph_A9_louv_max.rds")

saveRDS(community_graph_A7_info, "~/CESC_Network/5_Modularity/5_3_Module_Graph/5_3_1_1_Community_Graph_A7_info_max.rds")
saveRDS(community_graph_A9_info, "~/CESC_Network/5_Modularity/5_3_Module_Graph/5_3_1_1_Community_Graph_A9_info_max.rds")

### ========================
### Exportar en GraphML (igual que antes)
### ========================
write_graph(community_graph_A7_louv, "/STORAGE/csbig/jruiz/cesc_network/5_3_1_1_Community_Graph_A7_louv_max.graphml", format = "graphml")
write_graph(community_graph_A9_louv, "/STORAGE/csbig/jruiz/cesc_network/5_3_1_1_Community_Graph_A9_louv_max.graphml", format = "graphml")

write_graph(community_graph_A7_info, "/STORAGE/csbig/jruiz/cesc_network/5_3_1_1_Community_Graph_A7_info_max.graphml", format = "graphml")
write_graph(community_graph_A9_info, "/STORAGE/csbig/jruiz/cesc_network/5_3_1_1_Community_Graph_A9_info_max.graphml", format = "graphml")

message("✅ All done!")
