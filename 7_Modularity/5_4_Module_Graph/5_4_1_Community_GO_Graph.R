### ========================
### Build module-level meta-networks annotated with functional terms (ORA)
###
### For each network (A7, A9; Louvain modules named M1..Mn by size) build TWO
### meta-networks:
###   (a) modules + GO terms      (from the SIMPLIFIED GO ORA)
###   (b) modules + Hallmark sets  (from the MSigDB Hallmark ORA)
###
### In each: nodes = modules (M1..Mn) + term nodes; module-module edges are
### inter-module gene-gene links; module-term edges come from the ORA table.
### Term nodes get size = 500; module Ontology = "Community"; avg_degree is
### rescaled to [0.1, 1] for visualization.
###
### Output per (network x term type): RDS (local) + GraphML (STORAGE, Cytoscape).
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

geneRatio_to_num <- function(x) {
  parts <- strsplit(x, "/")
  vapply(parts, function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))
}

### ========================
### Step 1: module-module graph (nodes = modules, named M1..Mn)
### ========================
build_module_module_graph <- function(graph) {
  message("Processing graph: ", vcount(graph), " nodes, ", ecount(graph), " edges")
  if (!all(c("Community", "Community_name") %in% vertex_attr_names(graph)))
    stop("Faltan atributos de vértice: se requieren 'Community' y 'Community_name'.")

  nodes <- tibble(
    gene      = V(graph)$name,
    comm_id   = V(graph)$Community,
    comm_name = V(graph)$Community_name,
    degree    = igraph::degree(graph)
  )

  # Module attributes. Node name = Louvain module label (M1..Mn).
  # avg_degree = geometric mean of within-module degree.
  vertices <- nodes %>%
    group_by(comm_id) %>%
    summarise(
      name       = first(comm_name),
      size       = n(),
      avg_degree = exp(mean(log(pmax(degree, 1)), na.rm = TRUE)),
      .groups    = "drop"
    ) %>%
    mutate(
      normalized_degree = ifelse(
        max(avg_degree) == min(avg_degree), 0,
        (avg_degree - min(avg_degree)) / (max(avg_degree) - min(avg_degree))
      )
    )

  id2name <- setNames(vertices$name, as.character(vertices$comm_id))

  # Inter-module edges via gene-gene links that cross modules
  edges_by_comm <- as.data.frame(get.edgelist(graph), stringsAsFactors = FALSE) %>%
    rename(gene_from = V1, gene_to = V2) %>%
    left_join(nodes[, c("gene", "comm_id")], by = c("gene_from" = "gene")) %>%
    rename(comm_from = comm_id) %>%
    left_join(nodes[, c("gene", "comm_id")], by = c("gene_to" = "gene")) %>%
    rename(comm_to = comm_id) %>%
    filter(!is.na(comm_from), !is.na(comm_to), comm_from != comm_to) %>%
    mutate(comm_min = pmin(comm_from, comm_to), comm_max = pmax(comm_from, comm_to)) %>%
    count(comm_min, comm_max, name = "weight")

  edges_named <- edges_by_comm %>%
    transmute(from = id2name[as.character(comm_min)],
              to   = id2name[as.character(comm_max)], weight)

  vertices <- vertices %>% select(name, comm_id, size, avg_degree, normalized_degree)
  g_comm <- graph_from_data_frame(d = edges_named, vertices = vertices, directed = FALSE)
  message("Built module-module graph: ", vcount(g_comm), " nodes, ", ecount(g_comm), " edges")

  verts <- igraph::as_data_frame(g_comm, what = "vertices")
  edge_df <- igraph::as_data_frame(g_comm, what = "edges") %>%
    left_join(verts %>% rename(from = name, size_from = size, avg_degree_from = avg_degree,
                               norm_degree_from = normalized_degree, comm_id_from = comm_id),
              by = "from") %>%
    left_join(verts %>% rename(to = name, size_to = size, avg_degree_to = avg_degree,
                               norm_degree_to = normalized_degree, comm_id_to = comm_id),
              by = "to")

  list(graph = g_comm, edge_df = edge_df)
}

### ========================
### Step 2: extend module-module graph with term nodes (GO or Hallmark)
###   ora_tbl must have columns: Community_Name, Description, GeneRatio, Ontology
###   (Ontology = BP/CC/MF for GO, or "Hallmark" for Hallmark sets)
### ========================
build_mod_term_network <- function(df_modules, ora_tbl) {

  mod_attr <- df_modules %>%
    transmute(Community_Name = from, size = size_from,
              avg_degree = avg_degree_from, norm_degree = norm_degree_from) %>%
    bind_rows(df_modules %>%
              transmute(Community_Name = to, size = size_to,
                        avg_degree = avg_degree_to, norm_degree = norm_degree_to)) %>%
    distinct(Community_Name, .keep_all = TRUE)

  df_mod_term <- ora_tbl %>%
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

  df_ext <- bind_rows(df_modules %>% mutate(Ontology = NA_character_), df_mod_term)

  module_nodes <- df_ext %>%
    transmute(name = from, size = size_from, avg_degree = avg_degree_from, normalized_degree = norm_degree_from) %>%
    bind_rows(df_ext %>%
              transmute(name = to, size = size_to, avg_degree = avg_degree_to, normalized_degree = norm_degree_to)) %>%
    filter(!is.na(size) | !is.na(avg_degree) | !is.na(normalized_degree)) %>%
    distinct(name, .keep_all = TRUE) %>%
    mutate(node_type = "module")

  term_nodes <- df_ext %>%
    filter(!is.na(Ontology)) %>%
    transmute(name = to, Ontology = Ontology, GeneRatio = GeneRatio, node_type = "term") %>%
    group_by(name) %>%
    slice_max(order_by = geneRatio_to_num(GeneRatio), with_ties = FALSE) %>%
    ungroup() %>%
    distinct(name, .keep_all = TRUE)

  vertices_all <- bind_rows(module_nodes, term_nodes)
  edges_all <- df_ext %>% dplyr::select(from, to, weight, Ontology, GeneRatio)

  graph_from_data_frame(d = edges_all, vertices = vertices_all, directed = FALSE)
}

### ========================
### Step 3: final touches for visualization
### ========================
finalize_graph <- function(g) {
  # term nodes (non-module) get size = 500
  V(g)$size <- ifelse(!is.na(V(g)$Ontology), 500, V(g)$size)
  # module nodes have Ontology = NA -> "Community"
  V(g)$Ontology <- ifelse(is.na(V(g)$Ontology), "Community", V(g)$Ontology)
  # rescale avg_degree to [0.1, 1]
  ad  <- V(g)$avg_degree
  rng <- range(ad, na.rm = TRUE)
  V(g)$avg_degree_rescaled <- if (rng[1] == rng[2]) ifelse(is.na(ad), NA_real_, 1) else
    (ad - rng[1]) / (rng[2] - rng[1]) * (1 - 0.1) + 0.1
  g
}

### ========================
### Run pipeline for one network x term type
### ========================
process_network <- function(comm, ora_tbl, label, term_suffix) {
  message("=== ", label, " | ", term_suffix, " ===")
  g_final <- finalize_graph(build_mod_term_network(comm$edge_df, ora_tbl))

  saveRDS(g_final, file.path(out_dir, paste0("5_4_1_Community_Graph_", label, "_", term_suffix, ".rds")))
  write_graph(g_final,
              file.path(storage_dir, paste0("Community_Graph_", label, "_", term_suffix, ".graphml")),
              format = "graphml")

  n_terms <- sum(V(g_final)$Ontology != "Community")
  message("✔ ", label, " ", term_suffix, ": ", vcount(g_final), " nodes / ",
          ecount(g_final), " edges (", n_terms, " term nodes)")
  invisible(g_final)
}

### ========================
### Inputs
### ========================
graph_A7_louv <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A7_graph_max.rds"))
graph_A9_louv <- readRDS(file.path(mod_dir, "5_1_1_Louvain_A9_graph_max.rds"))

# GO (simplified) — already has an Ontology column (BP/CC/MF)
Ora_A7_GO <- vroom(file.path(ora_dir, "5_2_2_Ora7_louv_max_SIMPLIFIED_0_2.tsv"), show_col_types = FALSE)
Ora_A9_GO <- vroom(file.path(ora_dir, "5_2_2_Ora9_louv_max_SIMPLIFIED_0_2.tsv"), show_col_types = FALSE)

# Hallmark — no Ontology column; tag all terms as "Hallmark" and drop the
# "HALLMARK_" prefix from the set name so nodes show only the Hallmark name.
Ora_A7_HM <- vroom(file.path(ora_dir, "5_2_3_Hallmark7_louv_max.tsv"), show_col_types = FALSE) %>%
  mutate(Ontology = "Hallmark", Description = sub("^HALLMARK_", "", Description))
Ora_A9_HM <- vroom(file.path(ora_dir, "5_2_3_Hallmark9_louv_max.tsv"), show_col_types = FALSE) %>%
  mutate(Ontology = "Hallmark", Description = sub("^HALLMARK_", "", Description))

### ========================
### Build the module-module backbone once per network, then annotate twice
### ========================
comm_A7 <- build_module_module_graph(graph_A7_louv)
comm_A9 <- build_module_module_graph(graph_A9_louv)

process_network(comm_A7, Ora_A7_GO, "A7_louv", "with_GO")
process_network(comm_A9, Ora_A9_GO, "A9_louv", "with_GO")
process_network(comm_A7, Ora_A7_HM, "A7_louv", "with_Hallmark")
process_network(comm_A9, Ora_A9_HM, "A9_louv", "with_Hallmark")

message("✅ Built module meta-networks (GO + Hallmark) for A7 and A9 (RDS + GraphML).")
