############################################################
## SCRIPT: 5_3_GSEA_per_module.R
##
## Goal:
##   Perform GO GSEA (BP/CC/MF) per module for A7 and A9 
##   Louvain networks using limma logFC rankings.
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
##   - 5_3_1_GSEA_A7_FULL.tsv
##   - 5_3_1_GSEA_A7_SIMPLIFIED.tsv
##   - 5_3_1_GSEA_A9_FULL.tsv
##   - 5_3_1_GSEA_A9_SIMPLIFIED.tsv
##
## Notes:
##   - GSEA is run PER MODULE: for each module, we restrict the
##     ranked vector to genes in that module.
##   - GO ontology is taken from org.Hs.eg.db via gseGO().
############################################################

start_time <- Sys.time()

BiocManager::install("clusterProfiler")
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

# limma DE tables (must contain 'Gene' and 'logFC' columns)
DE_A7 <- readr::read_tsv(
  "~/CESC_Network/4_Network_analisis/4_3_limma_A7_vs_SNT_allGenes_annot_PR_DEG.tsv",
  col_types = cols()
)

DE_A9 <- readr::read_tsv(
  "~/CESC_Network/4_Network_analisis/4_3_limma_A9_vs_SNT_allGenes_annot_PR_DEG.tsv",
  col_types = cols()
)

# Create named ranking vectors: logFC named by gene symbol
rank_A7 <- setNames(DE_A7$logFC, DE_A7$Gene)
rank_A7 <- sort(rank_A7, decreasing = TRUE)

rank_A9 <- setNames(DE_A9$logFC, DE_A9$Gene)
rank_A9 <- sort(rank_A9, decreasing = TRUE)

############################################################
## 2. GSEA function for a single module (NO simplification)
############################################################
# This function:
#   - Restricts the ranking to genes present in the module.
#   - Runs gseGO() on that module-specific ranking.

gsea_full_module <- function(ranked_vector,
                             module_genes,
                             ont = c("BP", "CC", "MF"),
                             min_module_size = 10) {
  ont <- match.arg(ont)
  
  # Keep only genes that are both in the global ranking and in the module
  genes_in_module <- intersect(names(ranked_vector), module_genes)
  
  # If the module is too small after intersection, skip
  if (length(genes_in_module) < min_module_size) {
    return(NULL)
  }
  
  module_rank <- ranked_vector[genes_in_module]
  module_rank <- sort(module_rank, decreasing = TRUE)
  
  gseGO(
    geneList      = module_rank,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    minGSSize     = 10,
    maxGSSize     = 500,
    pAdjustMethod = "BH",
    eps           = 1e-10,
    verbose       = FALSE
  )
}

############################################################
## 3. GSEA simplification function
############################################################
# Applies pairwise_termsim + simplify() to reduce redundant GO terms.

gsea_simplify <- function(gsea_obj, cutoff = 0.2) {
  
  gsea_obj <- tryCatch({
    enrichplot::pairwise_termsim(gsea_obj)
  }, error = function(e) gsea_obj)
  
  gsea_obj <- tryCatch({
    simplify(
      gsea_obj,
      cutoff     = cutoff,
      by         = "p.adjust",
      select_fun = min,
      measure    = "Wang"
    )
  }, error = function(e) gsea_obj)
  
  return(gsea_obj)
}

############################################################
## 4. Run GSEA per module for one network
############################################################

run_gsea_network <- function(graph,
                             ranked_vector,
                             label,
                             cutoff = 0.2,
                             min_module_size = 10) {
  
  modules    <- unique(V(graph)$Community)
  ontologies <- c("BP", "CC", "MF")
  
  FULL_list  <- list()
  SIMPL_list <- list()
  
  for (m in modules) {
    message("Running module ", m, "...")
    
    module_genes <- V(graph)[Community == m]$name
    module_name  <- unique(V(graph)[Community == m]$Community_name)
    
    for (ont in ontologies) {
      
      gsea_res <- tryCatch({
        gsea_full_module(
          ranked_vector   = ranked_vector,
          module_genes    = module_genes,
          ont             = ont,
          min_module_size = min_module_size
        )
      }, error = function(e) NULL)
      
      if (is.null(gsea_res)) next
      
      df_full <- as.data.frame(gsea_res)
      if (nrow(df_full) == 0) next
      
      df_full$Module_ID   <- m
      df_full$Module_Name <- module_name
      df_full$Ontology    <- ont
      
      FULL_list[[length(FULL_list) + 1]] <- df_full
      
      # ---- Simplified version ----
      gsea_s <- gsea_simplify(gsea_res, cutoff = cutoff)
      df_s   <- as.data.frame(gsea_s)
      
      if (nrow(df_s) > 0) {
        df_s$Module_ID   <- m
        df_s$Module_Name <- module_name
        df_s$Ontology    <- ont
        SIMPL_list[[length(SIMPL_list) + 1]] <- df_s
      }
    }
  }
  
  # Bind all modules together
  FULL  <- dplyr::bind_rows(FULL_list)
  SIMPL <- dplyr::bind_rows(SIMPL_list)
  
  # Output directory
  outdir <- "~/CESC_Network/5_Modularity/5_3_GSEA_module/"
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Save results
  write.table(
    FULL,
    file = file.path(outdir, paste0("5_3_1_GSEA_", label, "_FULL.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  write.table(
    SIMPL,
    file = file.path(outdir, paste0("5_3_1_GSEA_", label, "_SIMPLIFIED.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  message("[", label, "] FULL terms: ", nrow(FULL),
          " | SIMPLIFIED terms: ", nrow(SIMPL))
  
  return(list(FULL = FULL, SIMPLIFIED = SIMPL))
}

############################################################
## 5. Run for A7 and A9
############################################################

gsea_A7 <- run_gsea_network(graph_A7, rank_A7, "A7", cutoff = 0.2)
gsea_A9 <- run_gsea_network(graph_A9, rank_A9, "A9", cutoff = 0.2)

end_time <- Sys.time()
cat("Total GSEA time:", round(end_time - start_time, 2), "minutes\n")
