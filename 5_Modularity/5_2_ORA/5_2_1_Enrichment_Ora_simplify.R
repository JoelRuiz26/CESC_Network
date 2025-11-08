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
library(enrichplot) # [ADDED] para pairwise_termsim() usado al mapear términos simplificados

plan(multicore, workers = 3)  # Configura el uso de 3 cores
set.seed(123)

### Cargar redes ###
#graph_A7 <- readRDS(file = "~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A7_graph_max.rds")
#graph_A9 <- readRDS(file = "~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A9_graph_max.rds")

graph_A7 <- readRDS(file = "~/CESC_Network/5_Modularity/5_1_Modules/5_1_1_Infomap_A7_graph_max.rds")
graph_A9 <- readRDS(file = "~/CESC_Network/5_Modularity/5_1_Modules/5_1_1_Infomap_A9_graph_max.rds")


### Cargar universo de genes ###
universo <- read_rds("~/CESC_Network/2_Prepro_TCGA_GTEx/2_1_Universe_annotation_table.rds") %>% pull(hgnc_symbol)

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
    )
  }, error = function(e) NULL) # [MODIFIED] devolver objeto enrichResult (no data.frame) para poder usar simplify()
  
  if (is.null(enrich_results)) return(data.frame()) # [ADDED]
  if (nrow(as.data.frame(enrich_results)) == 0) return(as.data.frame(enrich_results)) # [ADDED]
  
  # --- Simplificación de términos y mapeo por redundancia (cambios mínimos) ---
  # Intentar simplificar con parámetros conservadores (como en tu otro script)
  ego_simpl <- tryCatch(
    simplify(enrich_results,
             cutoff = 0.15,          # [ADDED] umbral de similitud semántica
             by = "p.adjust",        # [ADDED]
             select_fun = min,       # [ADDED]
             measure = "Wang"),      # [ADDED]
    error = function(e) enrich_results # fallback al original si falla
  ) # [ADDED]
  
  # IDs “conservados” por simplify (representantes)
  kept_ids <- tryCatch(ego_simpl@result$ID, error = function(e) character(0)) # [ADDED]
  
  # Similitud par-a-par para mapear cada término al representante conservado
  sim_mat <- tryCatch({
    pairwise_termsim(enrich_results)
    # pairwise_termsim devuelve enrichResult con slot @pairwise_termsim
    enrich_results@pairwise_termsim
  }, error = function(e) NULL) # [ADDED]
  
  df <- as.data.frame(enrich_results) # [MODIFIED] mover as.data.frame aquí
  if (is.null(df) || nrow(df) == 0) return(df) # [ADDED]
  
  # Construir columna "Simplified_Description"
  simplified_desc <- df$Description # por defecto, se deja igual # [ADDED]
  
  if (!is.null(sim_mat) && length(kept_ids) > 0) { # [ADDED]
    # Para cada término original, si no fue conservado por simplify,
    # se mapea al término conservado más similar por encima del cutoff (0.25).
    id_vec <- df$ID
    desc_lookup <- df$Description
    names(desc_lookup) <- id_vec
    
    # Descripciones de los representantes conservados
    kept_desc <- tryCatch({
      kd <- ego_simpl@result[, c("ID", "Description")]
      setNames(kd$Description, kd$ID)
    }, error = function(e) setNames(character(0), character(0)))
    
    for (i in seq_along(id_vec)) {
      id_i <- id_vec[i]
      if (id_i %in% kept_ids) {
        simplified_desc[i] <- kept_desc[[id_i]] %||% desc_lookup[[id_i]]
      } else {
        # elegir representante con mayor similitud >= 0.25
        if (!is.null(sim_mat) && id_i %in% rownames(sim_mat)) {
          sims <- sim_mat[id_i, kept_ids, drop = FALSE]
          if (!is.null(sims) && length(sims) > 0) {
            best_j <- names(which.max(sims))
            best_val <- as.numeric(sims[1, best_j])
            if (is.finite(best_val) && best_val >= 0.15 && best_j %in% names(kept_desc)) {
              simplified_desc[i] <- kept_desc[[best_j]]
            } else {
              simplified_desc[i] <- desc_lookup[[id_i]]
            }
          }
        }
      }
    }
  }
  
  df <- df %>%
    mutate(Community = community_id, Community_Name = community_name, Ontology = ontology,
           Simplified_Description = simplified_desc) # [ADDED]
  
  return(df) # [MODIFIED] ahora devuelve df con columna extra
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
#write.table(ora_A7, file = "~/CESC_Network/5_Modularity/5_2_ORA/5_2_1_Ora7_louv_max_S15.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(ora_A9, file = "~/CESC_Network/5_Modularity/5_2_ORA/5_2_1_Ora9_louv_max_S15.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

write.table(ora_A7, file = "~/CESC_Network/5_Modularity/5_2_ORA/5_2_1_Ora7_info_max_S15.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ora_A9, file = "~/CESC_Network/5_Modularity/5_2_ORA/5_2_1_Ora9_info_max_S15.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

end_time <- Sys.time()
cat("Tiempo total:", round(end_time - start_time, 2), "minutos\n")
#louvain Tiempo total: 8.98 minutos
#infomap Tiempo total: 9.08 minutos
