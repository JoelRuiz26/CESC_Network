### ========================
### Load libraries
### ========================
library(igraph)
library(tidyverse)

### ========================
### Load original graphs
### ========================
graph_A7_louv <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A7_graph_max.rds")
graph_A9_louv <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A9_graph_max.rds")
#graph_A7_info <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_1_1_Infomap_A7_graph_max.rds")
#graph_A9_info <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_1_1_Infomap_A9_graph_max.rds")

out_rds_dir <- "~/CESC_Network/5_Modularity/5_4_Module_Graph/"
out_gml_dir <- "/STORAGE/csbig/jruiz/cesc_network/"

### ========================
### Single function: build community graph + edge table
### ========================
build_comm_graph_and_df <- function(graph) {
  message("Processing graph: ", vcount(graph), " nodes, ", ecount(graph), " edges")
  
  # Chequeo de atributos necesarios
  if (!all(c("Community", "Community_name", "PageRank") %in% vertex_attr_names(graph))) {
    stop("Faltan atributos de vértice: se requieren 'Community', 'Community_name' y 'PageRank'.")
  }
  
  # ----------------------------
  # 1) Info a nivel de gen
  # ----------------------------
  nodes <- tibble(
    gene          = V(graph)$name,           # gen
    comm_id       = V(graph)$Community,      # ID de comunidad (membership)
    comm_name     = V(graph)$Community_name, # gen representante del módulo
    pagerank      = V(graph)$PageRank
  )
  
  # Cada gen debe pertenecer a un único módulo
  if (any(duplicated(nodes$gene))) {
    warning("Hay genes repetidos en V(graph)$name; la red original no debería tener nodos duplicados.")
  }
  
  # ----------------------------
  # 2) Atributos de cada módulo
  #    (agrupamos por comm_id, NO por comm_name)
  # ----------------------------
  vertices <- nodes %>%
    group_by(comm_id) %>%
    summarise(
      community        = first(comm_name),         # nombre del módulo = gen representativo
      size             = n(),                      # tamaño del módulo
      avg_pagerank     = mean(pagerank, na.rm = TRUE),
      .groups          = "drop"
    )
  
  # Normalizar PageRank medio (por módulo)
  vertices <- vertices %>%
    mutate(
      normalized_pagerank = ifelse(
        max(avg_pagerank) == min(avg_pagerank),
        0,
        (avg_pagerank - min(avg_pagerank)) /
          (max(avg_pagerank) - min(avg_pagerank))
      )
    )
  
  # Chequeo clave: ¿un mismo gen (community_name) representando >1 módulo?
  dup_comm_names <- vertices$community[duplicated(vertices$community)]
  if (length(dup_comm_names) > 0) {
    stop(
      "El/los gen(es) siguiente(s) aparecen como 'Community_name' de más de un módulo: ",
      paste(unique(dup_comm_names), collapse = ", "),
      ".\nUn gen no puede representar dos módulos diferentes: revisa la asignación de Community_name."
    )
  }
  
  # Para igraph, 'name' será el nombre de comunidad (gen representante)
  vertices <- vertices %>%
    mutate(
      name = community
    ) %>%
    select(name, comm_id, community, size, avg_pagerank, normalized_pagerank)
  
  # ----------------------------
  # 3) Agregar aristas entre módulos
  #    (trabajamos con comm_id y luego traducimos a community_name)
  # ----------------------------
  edges_by_gene <- as.data.frame(get.edgelist(graph), stringsAsFactors = FALSE) %>%
    rename(gene_from = V1, gene_to = V2) %>%
    left_join(nodes[, c("gene", "comm_id")], by = c("gene_from" = "gene")) %>%
    rename(comm_from = comm_id) %>%
    left_join(nodes[, c("gene", "comm_id")], by = c("gene_to" = "gene")) %>%
    rename(comm_to = comm_id) %>%
    filter(!is.na(comm_from), !is.na(comm_to), comm_from != comm_to) %>%
    mutate(
      comm_min = pmin(comm_from, comm_to),
      comm_max = pmax(comm_from, comm_to)
    ) %>%
    count(comm_min, comm_max, name = "weight")
  
  # Mapear IDs de comunidad → nombres de comunidad (genes)
  edges_named <- edges_by_gene %>%
    left_join(vertices %>% select(comm_id, community), by = c("comm_min" = "comm_id")) %>%
    rename(from = community) %>%
    left_join(vertices %>% select(comm_id, community), by = c("comm_max" = "comm_id")) %>%
    rename(to = community) %>%
    transmute(from, to, weight)
  
  # ----------------------------
  # 4) Construir la metared de módulos
  # ----------------------------
  g_comm <- graph_from_data_frame(
    d = edges_named,
    vertices = vertices,  # 'name' = gene representante del módulo
    directed = FALSE
  )
  message("Built community graph: ", vcount(g_comm), " nodes, ", ecount(g_comm), " edges")
  
  # ----------------------------
  # 5) Dataframe final de aristas con atributos de los módulos
  # ----------------------------
  verts <- igraph::as_data_frame(g_comm, what = "vertices") %>%
    select(name, comm_id, community, size, avg_pagerank, normalized_pagerank)
  
  edge_df <- igraph::as_data_frame(g_comm, what = "edges") %>%
    left_join(
      verts %>%
        rename(
          from  = name,
          size_from          = size,
          avg_pagerank_from  = avg_pagerank,
          norm_pagerank_from = normalized_pagerank,
          comm_id_from       = comm_id
        ),
      by = "from"
    ) %>%
    left_join(
      verts %>%
        rename(
          to  = name,
          size_to          = size,
          avg_pagerank_to  = avg_pagerank,
          norm_pagerank_to = normalized_pagerank,
          comm_id_to       = comm_id
        ),
      by = "to"
    )
  
  list(graph = g_comm, edge_df = edge_df)
}

### ========================
### Helper para guardar
### ========================
process_graph <- function(graph, label) {
  res <- build_comm_graph_and_df(graph)
  
  saveRDS(
    res$graph,
    file.path(out_rds_dir, paste0("5_3_1_1_Community_Graph_", label, "_max.rds"))
  )
  
  saveRDS(
    res$edge_df,
    file.path(out_rds_dir, paste0("5_3_1_1_Community_Graph_", label, "_edgelist_df.rds"))
  )
  
  write_graph(
    res$graph,
    file.path(out_gml_dir, paste0("5_3_1_1_Community_Graph_", label, "_max.graphml")),
    format = "graphml"
  )
  
  invisible(res)
}

### ========================
### Run para los cuatro grafos
### ========================
res_A7_louv <- process_graph(graph_A7_louv, "A7_louv")
res_A9_louv <- process_graph(graph_A9_louv, "A9_louv")
#res_A7_info <- process_graph(graph_A7_info, "A7_info")
#res_A9_info <- process_graph(graph_A9_info, "A9_info")

message("✅ All done!")
