### Cargar librerías ###
start_time <- Sys.time()

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(clusterProfiler)
  library(vroom)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(purrr)
  library(future.apply)
  library(readr)
  library(enrichplot)
})

## 0) Define el corte global ANTES de las funciones
SIMPLIFY_CUTOFF <- 0.2    # <-- muy astringente (colapsa casi todos los términos similares)

plan(multicore, workers = 3)
set.seed(123)

### Cargar redes ###
graph_A7 <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A7_graph_max.rds")
graph_A9 <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A9_graph_max.rds")
length(igraph::V(graph_A7)) # [1] 17481

### Universo ###
universo <- readr::read_rds("~/CESC_Network/2_Prepro_TCGA_GTEx/2_1_Universe_annotation_table.rds") %>%
  dplyr::pull(hgnc_symbol)
length(universo) #[1] 46731
#genes_A7 <- igraph::V(graph_A7)$name
#genes_A9 <- igraph::V(graph_A9)$name
#universo <- union(genes_A7, genes_A9)
#length(universo) #[1] 17481


### ORA con opción de simplificación ###
perform_ora <- function(genes, universe, community_id, community_name, ontology,
                        p_cut = 0.05,
                        q_cut = 0.20,
                        key_type = "SYMBOL",
                        simplify_terms = FALSE,
                        simplify_cutoff = 0.7,       # default interno
                        simplify_by = "p.adjust",
                        simplify_measure = "Wang",
                        simplify_select = min) {
  if (length(genes) == 0) return(data.frame())
  
  eg <- tryCatch({
    enrichGO(
      gene          = genes,
      universe      = universe,
      OrgDb         = org.Hs.eg.db,
      keyType       = key_type,
      ont           = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cut,
      qvalueCutoff  = q_cut,
      readable      = FALSE
    )
  }, error = function(e) NULL)
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) return(data.frame())
  
  if (isTRUE(simplify_terms)) {
    eg <- tryCatch({
      enrichplot::pairwise_termsim(eg, method = simplify_measure)
    }, error = function(e) eg)
    
    eg <- tryCatch({
      clusterProfiler::simplify(
        eg,
        cutoff     = simplify_cutoff,
        by         = simplify_by,
        select_fun = simplify_select,
        measure    = simplify_measure
      )
    }, error = function(e) eg)
  }
  
  df <- as.data.frame(eg)
  if (nrow(df) == 0) return(df)
  
  df %>%
    arrange(p.adjust, desc(Count)) %>%
    mutate(Community = community_id,
           Community_Name = community_name,
           Ontology = ontology)
}

### Extrae comunidades y corre ORA (FULL o SIMPLIFIED) ###
extract_communities_and_ora <- function(graph, universe,
                                        p_cut = 0.05,
                                        q_cut = 0.20,
                                        key_type = "SYMBOL",
                                        simplify_terms = FALSE,
                                        simplify_cutoff = 0.7,
                                        simplify_by = "p.adjust",
                                        simplify_measure = "Wang",
                                        simplify_select = min,
                                        min_count = 1) {
  communities <- unique(V(graph)$Community)
  ontologies  <- c("BP", "CC", "MF")
  
  ora_results <- future_lapply(communities, function(community) {
    genes <- V(graph)[V(graph)$Community == community]$name
    community_name <- unique(V(graph)[V(graph)$Community == community]$Community_name)
    community_size <- length(genes)
    avg_pagerank   <- mean(V(graph)[V(graph)$Community == community]$PageRank, na.rm = TRUE)
    
    enrich_df <- dplyr::bind_rows(lapply(ontologies, function(ontology) {
      perform_ora(
        genes, universe,
        community_id   = community,
        community_name = community_name,
        ontology       = ontology,
        p_cut          = p_cut,
        q_cut          = q_cut,
        key_type       = key_type,
        simplify_terms = simplify_terms,
        simplify_cutoff = simplify_cutoff,
        simplify_by     = simplify_by,
        simplify_measure = simplify_measure,
        simplify_select  = simplify_select
      )
    }))
    
    if (nrow(enrich_df) == 0) return(enrich_df)
    
    enrich_df %>%
      mutate(Community_size = community_size,
             Avg_PageRank   = avg_pagerank)
  }, future.seed = TRUE) %>% dplyr::bind_rows()
  
  if (nrow(ora_results) == 0) return(ora_results)
  ora_results %>% dplyr::filter(Count >= min_count)
}

### Salidas por red: FULL y SIMPLIFIED ###
out_dir <- "~/CESC_Network/5_Modularity/5_2_ORA/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Test rápido en un módulo (opcional)
mod <- unique(V(graph_A7)$Community)[1]
ont <- "BP"
g <- V(graph_A7)[V(graph_A7)$Community == mod]$name
full1 <- perform_ora(g, universo, mod, "demo", ont, simplify_terms = FALSE)
simp1 <- perform_ora(g, universo, mod, "demo", ont, simplify_terms = TRUE, simplify_cutoff = SIMPLIFY_CUTOFF)
cat("Full terms:", nrow(full1), " | Simplified terms:", nrow(simp1), "\n")

# A7 FULL
ora_A7_full <- extract_communities_and_ora(graph_A7, universo,
                                           simplify_terms = FALSE)

# A7 SIMPLIFIED
ora_A7_simpl <- extract_communities_and_ora(graph_A7, universo,
                                            simplify_terms = TRUE,
                                            simplify_cutoff = SIMPLIFY_CUTOFF)

# A9 FULL
ora_A9_full <- extract_communities_and_ora(graph_A9, universo,
                                           simplify_terms = FALSE)

# A9 SIMPLIFIED
ora_A9_simpl <- extract_communities_and_ora(graph_A9, universo,
                                            simplify_terms = TRUE,
                                            simplify_cutoff = SIMPLIFY_CUTOFF)

# Guardar (usa el cutoff real en el nombre)
cutoff_tag <- gsub("\\.", "_", as.character(SIMPLIFY_CUTOFF))

write.table(ora_A7_full,
            file = file.path(out_dir, "5_2_1_Ora7_louv_max_FULL.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ora_A7_simpl,
            file = file.path(out_dir, paste0("5_2_2_Ora7_louv_max_SIMPLIFIED_", cutoff_tag, ".tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(ora_A9_full,
            file = file.path(out_dir, "5_2_1_Ora9_louv_max_FULL.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ora_A9_simpl,
            file = file.path(out_dir, paste0("5_2_2_Ora9_louv_max_SIMPLIFIED_", cutoff_tag, ".tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Resumen rápido
cat("[A7] módulos:", length(unique(ora_A7_full$Community_Name)),
    " | FULL términos únicos:", length(unique(ora_A7_full$Description)),
    " | SIMPL términos únicos:", length(unique(ora_A7_simpl$Description)), "\n")
#[A7] módulos: 12  | FULL términos únicos: 904  | SIMPL términos únicos: 58 

cat("[A9] módulos:", length(unique(ora_A9_full$Community_Name)),
    " | FULL términos únicos:", length(unique(ora_A9_full$Description)),
    " | SIMPL términos únicos:", length(unique(ora_A9_simpl$Description)), "\n")
#[A9] módulos: 23  | FULL términos únicos: 477  | SIMPL términos únicos: 58 

end_time <- Sys.time()
cat("Tiempo total:", round(end_time - start_time, 2), "minutos\n")
