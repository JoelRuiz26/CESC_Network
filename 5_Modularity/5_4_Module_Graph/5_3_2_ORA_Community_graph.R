### ========================
### Load libraries
### ========================
library(igraph)
library(tidyverse)
library(vroom)

### ========================
### Paths
### ========================
mod_graph_dir <- "~/CESC_Network/5_Modularity/5_4_Module_Graph/"
ora_dir       <- "~/CESC_Network/5_Modularity/5_3_GSEA_module/"
storage_dir   <- "/STORAGE/csbig/jruiz/cesc_network"

### ========================
### Helper function:
### from module–module df + ORA/GSEA to:
###   - extended df (modules + GO edges)
###   - igraph object
###   - RDS (df + graph) + GraphML
### ========================
build_mod_GO_network <- function(df_modules,
                                 ora_tbl,
                                 label) {
  message("=== Processing ", label, " ===")
  
  # --- Module-level attributes (one row per module) ---
  mod_attr <- df_modules %>%
    # attributes from "from" side
    transmute(
      Community_Name = from,
      size          = size_from,
      avg_pagerank  = avg_pagerank_from,
      norm_pagerank = norm_pagerank_from
    ) %>%
    bind_rows(
      # attributes from "to" side
      df_modules %>%
        transmute(
          Community_Name = to,
          size          = size_to,
          avg_pagerank  = avg_pagerank_to,
          norm_pagerank = norm_pagerank_to
        )
    ) %>%
    distinct(Community_Name, .keep_all = TRUE)
  
  # --- Module–GO edges (weight = 1) ---
  df_mod_GO <- ora_tbl %>%
    left_join(mod_attr, by = "Community_Name") %>%
    transmute(
      from               = Community_Name,   # module name
      to                 = Description,      # GO term as node name
      weight             = 1L,               # fixed weight
      NES                = NES,              # NES from GSEA
      size_from          = size,
      avg_pagerank_from  = avg_pagerank,
      norm_pagerank_from = norm_pagerank,
      size_to            = NA_integer_,      # GO node: no size
      avg_pagerank_to    = NA_real_,         # GO node: no PageRank
      norm_pagerank_to   = NA_real_,
      Ontology           = Ontology          # BP / MF / CC
    )
  
  # --- Extend original module–module df with module–GO edges ---
  df_ext <- bind_rows(
    df_modules %>% mutate(Ontology = NA_character_),  # module–module edges
    df_mod_GO                                         # module–GO edges
  )
  
  # --- Build vertex table: modules + GO terms ---
  module_nodes <- df_ext %>%
    transmute(
      name          = from,
      size          = size_from,
      avg_pagerank  = avg_pagerank_from,
      normalized_pr = norm_pagerank_from
    ) %>%
    bind_rows(
      df_ext %>%
        transmute(
          name          = to,
          size          = size_to,
          avg_pagerank  = avg_pagerank_to,
          normalized_pr = norm_pagerank_to
        )
    ) %>%
    filter(!is.na(size) | !is.na(avg_pagerank) | !is.na(normalized_pr)) %>%
    distinct(name, .keep_all = TRUE) %>%
    mutate(node_type = "module")
  
  # --- GO term nodes (one row per GO term) ---
  # Aquí hacemos el cambio sutil:
  #   - Tomamos NES de df_ext (aristas módulo–GO).
  #   - Si un GO aparece en varios módulos, nos quedamos con el NES de mayor |NES|.
  go_nodes <- df_ext %>%
    filter(!is.na(Ontology)) %>%   # only module–GO edges
    transmute(
      name     = to,
      Ontology = Ontology,
      NES      = NES,
      node_type = "GO_term"
    ) %>%
    group_by(name) %>%
    slice_max(order_by = abs(NES), with_ties = FALSE) %>%  # NES "dominante"
    ungroup() %>%
    distinct(name, .keep_all = TRUE)
  
  vertices_all <- bind_rows(module_nodes, go_nodes)
  
  # --- Build igraph object from extended edge table ---
  edges_all <- df_ext %>%
    dplyr::select(from, to, weight, Ontology, NES)  # NES sigue como atributo de arista
  
  g_mod_GO <- graph_from_data_frame(
    d        = edges_all,
    vertices = vertices_all,
    directed = FALSE
  )
  
  # --- Save outputs: df + graph (RDS) + GraphML ---
  # Extended dataframe
  saveRDS(
    df_ext,
    file.path(
      mod_graph_dir,
      paste0("5_3_2_1_Community_Graph_", label, "_with_GO_df.rds")
    )
  )
  
  # Graph as RDS
  saveRDS(
    g_mod_GO,
    file.path(
      mod_graph_dir,
      paste0("5_3_2_1_Community_Graph_", label, "_with_GO.rds")
    )
  )
  
  # GraphML in STORAGE
  write_graph(
    g_mod_GO,
    file.path(
      storage_dir,
      paste0("5_3_2_1_Community_Graph_", label, "_with_GO.graphml")
    ),
    format = "graphml"
  )
  
  message("✔ Done: ", label)
  invisible(list(df = df_ext, graph = g_mod_GO))
}

### ========================
### Load inputs for A7 and A9
### ========================

# A7 Louvain module–module edges
df_modules_louv_A7 <- readRDS(
  file.path(
    mod_graph_dir,
    "5_3_1_1_Community_Graph_A7_louv_edgelist_df.rds"
  )
)

# A9 Louvain module–module edges
df_modules_louv_A9 <- readRDS(
  file.path(
    mod_graph_dir,
    "5_3_1_1_Community_Graph_A9_louv_edgelist_df.rds"
  )
)

# GSEA A7 (simplified)
Ora_A7_SIMP <- vroom(
  file.path(
    ora_dir,
    "5_3_1_GSEA_A7_SIMPLIFIED.tsv"
  )
) %>%
  dplyr::rename(
    Community_Name = Module_Name   # mapeo a nombre esperado
  )

# GSEA A9 (simplified)
Ora_A9_SIMP <- vroom(
  file.path(
    ora_dir,
    "5_3_1_GSEA_A9_SIMPLIFIED.tsv"
  )
) %>%
  dplyr::rename(
    Community_Name = Module_Name
  )

### ========================
### Run for A7 and A9
### ========================
res_A7 <- build_mod_GO_network(
  df_modules = df_modules_louv_A7,
  ora_tbl    = Ora_A7_SIMP,
  label      = "A7_louv"
)

res_A9 <- build_mod_GO_network(
  df_modules = df_modules_louv_A9,
  ora_tbl    = Ora_A9_SIMP,
  label      = "A9_louv"
)

message("✅ Built modular networks with GO processes for A7 and A9 (df + RDS + GraphML).")
