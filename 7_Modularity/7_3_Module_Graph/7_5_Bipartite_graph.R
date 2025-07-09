### ========================
### Cargar librerías
### ========================
library(igraph)
library(tidyverse)
library(vroom)

### ========================
### 1️⃣ Cargar metaredes desde RDS
### ========================
community_graph_A7 <- readRDS("~/CESC_Network/7_Modularity/7_4_1_Community_Graph_A7.rds")
community_graph_A9 <- readRDS("~/CESC_Network/7_Modularity/7_4_1_Community_Graph_A9.rds")

### ========================
### 2️⃣ Cargar tablas ORA
### ========================
ora_A7 <- vroom("~/CESC_Network/7_Modularity/7_2_1_Ora7_600k.tsv")
ora_A9 <- vroom("~/CESC_Network/7_Modularity/7_2_1_Ora9_600k.tsv")

### ========================
### 3️⃣ Función para construir red combinada
### ========================
build_combined_enrichment_network <- function(ora_table, community_graph) {
  
  message("✔️ Processing ORA table with ", n_distinct(ora_table$Community_Name), " communities")
  
  # ------------------------------------
  # 1. Edges ENTRE comunidades (del metagrafo)
  # ------------------------------------
  community_edges <- igraph::as_data_frame(community_graph, what = "edges")
  
  message("✔️ Community edges: ", nrow(community_edges))
  
  # ------------------------------------
  # 2. Filtrar top 3 por Ontología por Comunidad para ENRICHMENT
  # ------------------------------------
  top_ora <- ora_table %>%
    group_by(Community_Name) %>%
    arrange(p.adjust) %>%
    slice_head(n = 3) %>%
    ungroup()
  
  # Edges bipartitos (comunidad ↔ GO)
  enrichment_edges <- top_ora %>%
    select(Community_Name, Description) %>%
    distinct() %>%
    rename(from = Community_Name, to = Description) %>%
    mutate(weight = 20,
           normalized_weight = 20)
  
  message("✔️ Enrichment edges: ", nrow(enrichment_edges))
  
  # ------------------------------------
  # 3. Combinar ambos sets de edges
  # ------------------------------------
  all_edges <- bind_rows(
    community_edges,
    enrichment_edges
  )
  
  message("✔️ Total combined edges: ", nrow(all_edges))
  
  # ------------------------------------
  # 4. Nodos de comunidades (con atributos)
  # ------------------------------------
  community_nodes <- data.frame(
    name = V(community_graph)$name,
    size = V(community_graph)$size,
    normalized_pagerank = V(community_graph)$normalized_pagerank
  ) %>%
    mutate(type = "community")
  
  # Filtrar solo comunidades que aparecen en edges finales
  community_nodes <- community_nodes %>%
    filter(name %in% all_edges$from | name %in% all_edges$to)
  
  # ------------------------------------
  # 5. Nodos GO terms (Ontology, size=1)
  # ------------------------------------
  go_nodes <- top_ora %>%
    distinct(Description, Ontology) %>%
    rename(name = Description) %>%
    mutate(
      type = "GO",
      size = 1,
      normalized_pagerank = NA
    )
  
  # ------------------------------------
  # 6. Combinar todos los nodos
  # ------------------------------------
  nodes <- bind_rows(community_nodes, go_nodes)
  
  message("✔️ Total nodes: ", nrow(nodes))
  
  # ------------------------------------
  # 7. Crear grafo combinado
  # ------------------------------------
  graph_combined <- graph_from_data_frame(
    d = all_edges,
    vertices = nodes,
    directed = FALSE
  )
  
  message("✔️ Graph built with ", vcount(graph_combined), " nodes and ", ecount(graph_combined), " edges")
  
  return(graph_combined)
}

### ========================
### 4️⃣ Construir redes combinadas
### ========================
combined_graph_A7 <- build_combined_enrichment_network(ora_A7, community_graph_A7)
combined_graph_A9 <- build_combined_enrichment_network(ora_A9, community_graph_A9)


### ========================
### 5️⃣ Guardar en GraphML
### ========================
write_graph(combined_graph_A7, "~/CESC_Network/7_Modularity/7_5_1_Bipartite_Enrichment_A7.graphml", format = "graphml")
write_graph(combined_graph_A9, "~/CESC_Network/7_Modularity/7_5_1_Bipartite_Enrichment_A9.graphml", format = "graphml")
message("✅ All done! Combined GraphML files created.")
