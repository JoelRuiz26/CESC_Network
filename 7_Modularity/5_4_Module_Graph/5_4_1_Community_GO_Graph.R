### ========================
### Build module-level meta-network annotated with GO terms (ORA)
###
### For each network (A7, A9, Louvain communities):
###   1) Collapse the gene network into a module-module graph
###      (nodes = modules, edges = inter-module gene-gene links).
###   2) Attach GO term nodes (BP/CC/MF) from ORA results, with
###      module-GO edges weighted by GeneRatio.
###   3) Final touches: GO node size = 500, Ontology NA -> "Community",
###      avg_degree rescaled to [0.1, 1] for visualization.
###
### Output (per network): RDS graph + GraphML for Cytoscape.
### ========================
library(igraph)
library(tidyverse)
library(vroom)

### ========================
### Paths
### ========================
mod_dir     <- "~/CESC_Network/7_Modularity/5_1_Modules/"
ora_dir     <- "~/CESC_Network/7_Modularity/5_2_ORA/"
out_dir     <- "~/CESC_Network/7_Modularity/5_4_Module_Graph/"
storage_dir <- "/STORAGE/csbig/jruiz/cesc_network"

### ========================
### Helper: GeneRatio ("x/y") -> numeric, to pick the "dominant"
### term when the same GO appears in >1 module
### ========================
geneRatio_to_num <- function(x) {
  parts <- strsplit(x, "/")
  vapply(parts, function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))
}

### ========================
### Step 1: module-module graph + edge table
### ========================
build_module_module_graph <- function(graph) {
  message("Processing graph: ", vcount(graph), " nodes, ", ecount(graph), " edges")

  if (!all(c("Community", "Community_name") %in% vertex_attr_names(graph))) {
    stop("Faltan atributos de vértice: se requieren 'Community' y 'Community_name'.")
  }

  # Info a nivel de gen
  nodes <- tibble(
    gene      = V(graph)$name,
    comm_id   = V(graph)$Community,
    comm_name = V(graph)$Community_name,
    degree    = igraph::degree(graph)
  )

  # Atributos de cada módulo (agrupado por comm_id)
  # avg_degree = media geométrica del grado (mu de la log-normal ajustada en
  # la distribución de grado global), localizada a cada módulo.
  # community  = gen de mayor grado dentro del módulo (reemplaza al gen de
  # mayor PageRank usado antes como Community_name).
  vertices <- nodes %>%
    group_by(comm_id) %>%
    summarise(
      old_name   = first(comm_name),
      community  = gene[which.max(degree)],
      size       = n(),
      avg_degree = exp(mean(log(pmax(degree, 1)), na.rm = TRUE)),
      .groups    = "drop"
    ) %>%
    mutate(
      normalized_degree = ifelse(
        max(avg_degree) == min(avg_degree),
        0,
        (avg_degree - min(avg_degree)) / (max(avg_degree) - min(avg_degree))
      )
    )

  dup_comm_names <- vertices$community[duplicated(vertices$community)]
  if (length(dup_comm_names) > 0) {
    stop(
      "El/los gen(es) siguiente(s) son el de mayor grado en más de un módulo: ",
      paste(unique(dup_comm_names), collapse = ", "),
      ".\nUn gen no puede representar dos módulos diferentes."
    )
  }

  # Mapeo de nombres: Community_name (PageRank, en la tabla ORA) -> gen de
  # mayor grado, para renombrar la tabla ORA sin tocar el .tsv en disco.
  name_map <- vertices %>% select(old_name, community)

  vertices <- vertices %>%
    mutate(name = community) %>%
    select(name, comm_id, community, size, avg_degree, normalized_degree)

  # Aristas entre módulos (vía aristas gen-gen que cruzan módulos)
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

  edges_named <- edges_by_gene %>%
    left_join(vertices %>% select(comm_id, community), by = c("comm_min" = "comm_id")) %>%
    rename(from = community) %>%
    left_join(vertices %>% select(comm_id, community), by = c("comm_max" = "comm_id")) %>%
    rename(to = community) %>%
    transmute(from, to, weight)

  g_comm <- graph_from_data_frame(d = edges_named, vertices = vertices, directed = FALSE)
  message("Built module-module graph: ", vcount(g_comm), " nodes, ", ecount(g_comm), " edges")

  # Tabla de aristas con atributos de los módulos en ambos extremos
  verts <- igraph::as_data_frame(g_comm, what = "vertices") %>%
    select(name, comm_id, community, size, avg_degree, normalized_degree)

  edge_df <- igraph::as_data_frame(g_comm, what = "edges") %>%
    left_join(
      verts %>% rename(from = name, size_from = size, avg_degree_from = avg_degree,
                        norm_degree_from = normalized_degree, comm_id_from = comm_id),
      by = "from"
    ) %>%
    left_join(
      verts %>% rename(to = name, size_to = size, avg_degree_to = avg_degree,
                        norm_degree_to = normalized_degree, comm_id_to = comm_id),
      by = "to"
    )

  list(graph = g_comm, edge_df = edge_df, name_map = name_map)
}

### ========================
### Step 2: extend module-module graph with GO term nodes (ORA)
### ========================
build_mod_GO_network <- function(df_modules, ora_tbl) {

  # Atributos por módulo (de los lados "from" y "to" de las aristas módulo-módulo)
  mod_attr <- df_modules %>%
    transmute(Community_Name = from, size = size_from,
              avg_degree = avg_degree_from, norm_degree = norm_degree_from) %>%
    bind_rows(
      df_modules %>%
        transmute(Community_Name = to, size = size_to,
                  avg_degree = avg_degree_to, norm_degree = norm_degree_to)
    ) %>%
    distinct(Community_Name, .keep_all = TRUE)

  # Aristas módulo-GO (peso fijo = 1)
  df_mod_GO <- ora_tbl %>%
    left_join(mod_attr, by = "Community_Name") %>%
    transmute(
      from             = Community_Name,
      to               = Description,
      weight           = 1L,
      GeneRatio        = GeneRatio,
      size_from        = size,
      avg_degree_from  = avg_degree,
      norm_degree_from = norm_degree,
      size_to          = NA_integer_,
      avg_degree_to    = NA_real_,
      norm_degree_to   = NA_real_,
      Ontology         = Ontology
    )

  # Tabla de aristas extendida: módulo-módulo + módulo-GO
  df_ext <- bind_rows(
    df_modules %>% mutate(Ontology = NA_character_),
    df_mod_GO
  )

  # Nodos de módulo
  module_nodes <- df_ext %>%
    transmute(name = from, size = size_from, avg_degree = avg_degree_from, normalized_degree = norm_degree_from) %>%
    bind_rows(
      df_ext %>%
        transmute(name = to, size = size_to, avg_degree = avg_degree_to, normalized_degree = norm_degree_to)
    ) %>%
    filter(!is.na(size) | !is.na(avg_degree) | !is.na(normalized_degree)) %>%
    distinct(name, .keep_all = TRUE) %>%
    mutate(node_type = "module")

  # Nodos de término GO: si un GO aparece en varios módulos, nos quedamos
  # con el de mayor GeneRatio
  go_nodes <- df_ext %>%
    filter(!is.na(Ontology)) %>%
    transmute(name = to, Ontology = Ontology, GeneRatio = GeneRatio, node_type = "GO_term") %>%
    group_by(name) %>%
    slice_max(order_by = geneRatio_to_num(GeneRatio), with_ties = FALSE) %>%
    ungroup() %>%
    distinct(name, .keep_all = TRUE)

  vertices_all <- bind_rows(module_nodes, go_nodes)

  edges_all <- df_ext %>% dplyr::select(from, to, weight, Ontology, GeneRatio)

  graph_from_data_frame(d = edges_all, vertices = vertices_all, directed = FALSE)
}

### ========================
### Step 3: final touches for visualization
### ========================
finalize_graph <- function(g) {
  # GO term nodes (BP/MF/CC) get size = 500
  V(g)$size <- ifelse(V(g)$Ontology %in% c("BP", "MF", "CC"), 500, V(g)$size)

  # Module nodes have Ontology = NA -> "Community"
  V(g)$Ontology <- ifelse(is.na(V(g)$Ontology), "Community", V(g)$Ontology)

  # Rescale avg_degree to [0.1, 1] for visualization
  ad  <- V(g)$avg_degree
  rng <- range(ad, na.rm = TRUE)
  if (rng[1] == rng[2]) {
    ad_rescaled <- ifelse(is.na(ad), NA_real_, 1)
  } else {
    ad_rescaled <- (ad - rng[1]) / (rng[2] - rng[1]) * (1 - 0.1) + 0.1
  }
  V(g)$avg_degree_rescaled <- ad_rescaled

  g
}

### ========================
### Run pipeline for one network and save outputs
### ========================
process_network <- function(graph_louv, ora_tbl, label) {
  message("=== Processing ", label, " ===")

  comm <- build_module_module_graph(graph_louv)

  # Renombrar módulos en la tabla ORA (PageRank -> gen de mayor grado) vía
  # el mapeo por módulo; la tabla en disco no se modifica.
  map_vec <- setNames(comm$name_map$community, comm$name_map$old_name)
  ora_tbl <- ora_tbl %>% mutate(Community_Name = unname(map_vec[Community_Name]))

  g_mod_go <- build_mod_GO_network(comm$edge_df, ora_tbl)
  g_final  <- finalize_graph(g_mod_go)

  saveRDS(g_final, file.path(out_dir, paste0("5_4_1_Community_Graph_", label, "_with_GO.rds")))

  write_graph(
    g_final,
    file.path(storage_dir, paste0("Community_Graph_", label, "_with_GO.graphml")),
    format = "graphml"
  )

  message("✔ Done: ", label, " (", vcount(g_final), " nodes, ", ecount(g_final), " edges)")
  invisible(g_final)
}

### ========================
### Inputs
### ========================
graph_A7_louv <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
graph_A9_louv <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))

Ora_A7_SIMP <- vroom(file.path(ora_dir, "5_2_2_Ora7_louv_max_SIMPLIFIED_0_2.tsv"))
Ora_A9_SIMP <- vroom(file.path(ora_dir, "5_2_2_Ora9_louv_max_SIMPLIFIED_0_2.tsv"))

### ========================
### Run for A7 and A9
### ========================
res_A7 <- process_network(graph_A7_louv, Ora_A7_SIMP, "A7_louv")
res_A9 <- process_network(graph_A9_louv, Ora_A9_SIMP, "A9_louv")

message("✅ Built modular networks with GO processes for A7 and A9 (RDS + GraphML).")
