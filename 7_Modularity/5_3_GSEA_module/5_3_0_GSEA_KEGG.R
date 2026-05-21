############################################################
## SCRIPT: 5_3_GSEA_KEGG_per_module.R
##
## Goal:
##   Perform KEGG GSEA per module for A7 and A9 Louvain 
##   networks using limma logFC rankings.
##
## Inputs:
##   - Networks:
##       * 5_2_1_Louvain_A7_graph_max.rds
##       * 5_2_1_Louvain_A9_graph_max.rds
##   - limma DE results (all genes, with logFC column):
##       * 2_1_limma_A7_vs_SNT_allGenes.tsv
##       * 2_1_limma_A9_vs_SNT_allGenes.tsv
##
## Outputs (per network):
##   - 5_3_1_GSEA_KEGG_A7_FULL.tsv
##   - 5_3_1_GSEA_KEGG_A7_SIMPLIFIED.tsv
##   - 5_3_1_GSEA_KEGG_A9_FULL.tsv
##   - 5_3_1_GSEA_KEGG_A9_SIMPLIFIED.tsv
##
## Notes:
##   - GSEA is run PER MODULE: for each module, we restrict the
##     ranked vector to genes in that module.
##   - KEGG enrichment is done with gseKEGG() using ENTREZ IDs.
############################################################
start_time <- Sys.time()


suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(clusterProfiler)
  library(vroom)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(future.apply)
  library(readr)
})

plan(multicore, workers = 3)
set.seed(123)

############################################################
## 1. Load networks and limma DE rankings
############################################################

# Louvain community networks
graph_A7 <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A7_graph_max.rds")
graph_A9 <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A9_graph_max.rds")

# limma DE tables (must contain 'Gene' and 'logFC' columns in SYMBOL)
DE_A7 <- readr::read_tsv(
  "~/CESC_Network/4_Network_analisis/4_3_limma_A7_vs_SNT_allGenes_annot_PR_DEG.tsv",
  col_types = cols()
)

DE_A9 <- readr::read_tsv(
  "~/CESC_Network/4_Network_analisis/4_3_limma_A9_vs_SNT_allGenes_annot_PR_DEG.tsv",
  col_types = cols()
)

############################################################
## 2. SYMBOL -> ENTREZID mapping for KEGG
############################################################

# Universo de símbolos que aparecen en DE y en las redes
gene_universe <- unique(c(
  DE_A7$Gene,
  DE_A9$Gene,
  V(graph_A7)$name,
  V(graph_A9)$name
))

# Mapear SYMBOL -> ENTREZID
symbol2entrez <- bitr(
  gene_universe,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# Función para crear ranking logFC pero nombrado por ENTREZID
make_entrez_ranking <- function(DE_tbl, map_tbl) {
  merged <- DE_tbl %>%
    dplyr::select(Gene, logFC) %>%
    inner_join(map_tbl, by = c("Gene" = "SYMBOL")) %>%
    # Si un ENTREZID tiene varios símbolos, nos quedamos con el de |logFC| mayor
    group_by(ENTREZID) %>%
    slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(logFC))
  
  rank_vec <- merged$logFC
  names(rank_vec) <- merged$ENTREZID
  rank_vec <- sort(rank_vec, decreasing = TRUE)
  return(rank_vec)
}

# Ranking KEGG (ENTREZID) para A7 y A9
rank_A7_kegg <- make_entrez_ranking(DE_A7, symbol2entrez)
rank_A9_kegg <- make_entrez_ranking(DE_A9, symbol2entrez)

############################################################
## 3. GSEA function for a single module (KEGG, no simplification)
############################################################
# - Restringe el ranking global (ENTREZID) a los genes del módulo
#   (mapeando los símbolos de la red a ENTREZID).
# - Ejecuta gseKEGG().

gsea_full_module_kegg <- function(ranked_vector,
                                  module_genes_symbol,
                                  map_tbl,
                                  min_module_size = 10) {
  
  # Símbolos del módulo -> ENTREZID
  module_entrez <- map_tbl %>%
    filter(SYMBOL %in% module_genes_symbol) %>%
    pull(ENTREZID) %>%
    unique()
  
  genes_in_module <- intersect(names(ranked_vector), module_entrez)
  
  if (length(genes_in_module) < min_module_size) {
    return(NULL)
  }
  
  module_rank <- ranked_vector[genes_in_module]
  module_rank <- sort(module_rank, decreasing = TRUE)
  
  gseKEGG(
    geneList      = module_rank,
    organism      = "hsa",
    keyType       = "ncbi-geneid",
    minGSSize     = 10,
    maxGSSize     = 500,
    pAdjustMethod = "BH",
    eps           = 1e-10,
    verbose       = FALSE
  )
}

############################################################
## 4. GSEA simplification function (KEGG)
############################################################
# Para KEGG usamos similitud basada en conjunto de genes (JC, Jaccard).

gsea_simplify_kegg <- function(gsea_obj, cutoff = 0.2) {
  
  # Si está vacío, regresar tal cual
  df_tmp <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df_tmp) || nrow(df_tmp) == 0) return(gsea_obj)
  
  gsea_obj <- tryCatch({
    enrichplot::pairwise_termsim(gsea_obj)  # por defecto usa "JC"
  }, error = function(e) gsea_obj)
  
  gsea_obj <- tryCatch({
    simplify(
      gsea_obj,
      cutoff     = cutoff,
      by         = "p.adjust",
      select_fun = min,
      measure    = "JC"
    )
  }, error = function(e) gsea_obj)
  
  return(gsea_obj)
}

############################################################
## 5. Run KEGG GSEA per module for one network
############################################################

run_kegg_gsea_network <- function(graph,
                                  ranked_vector,
                                  map_tbl,
                                  label,
                                  cutoff = 0.2,
                                  min_module_size = 10) {
  
  modules <- sort(unique(V(graph)$Community))
  
  FULL_list  <- list()
  SIMPL_list <- list()
  
  for (m in modules) {
    message("Running KEGG GSEA for module ", m, "...")
    
    module_genes_symbol <- V(graph)[Community == m]$name
    module_name         <- unique(V(graph)[Community == m]$Community_name)
    
    # --- KEGG GSEA para este módulo ---
    gsea_res <- tryCatch({
      gsea_full_module_kegg(
        ranked_vector       = ranked_vector,
        module_genes_symbol = module_genes_symbol,
        map_tbl             = map_tbl,
        min_module_size     = min_module_size
      )
    }, error = function(e) NULL)
    
    if (is.null(gsea_res)) next
    
    df_full <- as.data.frame(gsea_res)
    if (nrow(df_full) == 0) next
    
    df_full$Module_ID   <- m
    df_full$Module_Name <- module_name
    df_full$Ontology    <- "KEGG"   # aquí no hay BP/CC/MF, es una sola colección
    
    FULL_list[[length(FULL_list) + 1]] <- df_full
    
    # --- Versión simplificada ---
    gsea_s <- gsea_simplify_kegg(gsea_res, cutoff = cutoff)
    df_s   <- as.data.frame(gsea_s)
    
    if (!is.null(df_s) && nrow(df_s) > 0) {
      df_s$Module_ID   <- m
      df_s$Module_Name <- module_name
      df_s$Ontology    <- "KEGG"
      SIMPL_list[[length(SIMPL_list) + 1]] <- df_s
    }
  }
  
  FULL  <- dplyr::bind_rows(FULL_list)
  SIMPL <- dplyr::bind_rows(SIMPL_list)
  
  # Output directory
  outdir <- "~/CESC_Network/5_Modularity/5_3_GSEA_module/"
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Save results
  write.table(
    FULL,
    file = file.path(outdir, paste0("5_3_1_GSEA_KEGG_", label, "_FULL_KEGG.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    SIMPL,
    file = file.path(outdir, paste0("5_3_1_GSEA_KEGG_", label, "_SIMPLIFIED_KEGG.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  message("[KEGG ", label, "] FULL pathways: ", nrow(FULL),
          " | SIMPLIFIED pathways: ", nrow(SIMPL))
  
  return(list(FULL = FULL, SIMPLIFIED = SIMPL))
}

############################################################
## 6. Run for A7 and A9
############################################################

kegg_A7 <- run_kegg_gsea_network(
  graph         = graph_A7,
  ranked_vector = rank_A7_kegg,
  map_tbl       = symbol2entrez,
  label         = "A7",
  cutoff        = 0.2
)

kegg_A9 <- run_kegg_gsea_network(
  graph         = graph_A9,
  ranked_vector = rank_A9_kegg,
  map_tbl       = symbol2entrez,
  label         = "A9",
  cutoff        = 0.2
)

end_time <- Sys.time()
cat("Total KEGG GSEA time:", round(end_time - start_time, 2), "minutes\n")
