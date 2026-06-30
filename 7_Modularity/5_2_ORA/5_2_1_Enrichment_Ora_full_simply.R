### ORA with GO enrichment — Louvain modules, HPV-A7 and HPV-A9
### Simplification: "Rel" IC-based similarity with pre-computed semData (fast + correct)
### Outputs: FULL and SIMPLIFIED TSV for A7 and A9
start_time <- Sys.time()

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(clusterProfiler)
  library(enrichplot)
  library(GOSemSim)
  library(vroom)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(purrr)
  library(future.apply)
  library(readr)
})

SIMPLIFY_CUTOFF <- 0.2
plan(multicore, workers = 3)
set.seed(123)

### Load networks ###
graph_A7 <- readRDS("~/CESC_Network/7_Modularity/5_1_Modules/5_1_1_Louvain_A7_graph_max.rds")
graph_A9 <- readRDS("~/CESC_Network/7_Modularity/5_1_Modules/5_1_1_Louvain_A9_graph_max.rds")
cat("[1] A7:", vcount(graph_A7), "nodes |", ecount(graph_A7), "edges\n")
cat("[1] A9:", vcount(graph_A9), "nodes |", ecount(graph_A9), "edges\n")

### Universe ###
universo <- readr::read_rds("~/CESC_Network/2_Prepro_TCGA_GTEx/2_1_Universe_annotation_table.rds") %>%
  pull(hgnc_symbol) %>% unique()
cat("[2] Universe:", length(universo), "genes\n")

### Pre-compute semData once per ontology (avoids rebuilding inside each worker) ###
cat("[3] Pre-computing GO semantic data (Rel, IC-based) for BP, CC, MF...\n")
semData_BP <- GOSemSim::godata(org.Hs.eg.db, ont = "BP", keytype = "SYMBOL", computeIC = TRUE)
semData_CC <- GOSemSim::godata(org.Hs.eg.db, ont = "CC", keytype = "SYMBOL", computeIC = TRUE)
semData_MF <- GOSemSim::godata(org.Hs.eg.db, ont = "MF", keytype = "SYMBOL", computeIC = TRUE)
sem_list <- list(BP = semData_BP, CC = semData_CC, MF = semData_MF)
cat("[3] semData ready.\n")

### ORA for one community x one ontology ###
perform_ora <- function(genes, universe, community_id, community_name, ontology,
                        semData,
                        p_cut = 0.05, q_cut = 0.20) {
  if (length(genes) < 5) return(list(full = data.frame(), simp = data.frame()))

  eg <- tryCatch({
    enrichGO(
      gene          = genes,
      universe      = universe,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cut,
      qvalueCutoff  = q_cut,
      readable      = FALSE
    )
  }, error = function(e) NULL)

  if (is.null(eg) || nrow(as.data.frame(eg)) == 0)
    return(list(full = data.frame(), simp = data.frame()))

  df_full <- as.data.frame(eg) %>%
    arrange(p.adjust, desc(Count)) %>%
    mutate(Community = community_id, Community_Name = community_name, Ontology = ontology)

  # Simplify using pre-computed semData
  eg_sim <- tryCatch({
    eg2 <- enrichplot::pairwise_termsim(eg, method = "Rel", semData = semData)
    clusterProfiler::simplify(eg2, cutoff = SIMPLIFY_CUTOFF, by = "p.adjust",
                               select_fun = min, measure = "Rel", semData = semData)
  }, error = function(e) eg)

  df_simp <- as.data.frame(eg_sim) %>%
    arrange(p.adjust, desc(Count)) %>%
    mutate(Community = community_id, Community_Name = community_name, Ontology = ontology)

  list(full = df_full, simp = df_simp)
}

### Run ORA for all communities of one graph ###
run_network_ora <- function(graph, network_label, universe, sem_list) {
  communities <- sort(unique(V(graph)$Community))
  ontologies  <- c("BP", "CC", "MF")
  cat(sprintf("[ORA] %s — %d communities\n", network_label, length(communities)))

  results <- future_lapply(communities, function(community) {
    genes <- V(graph)[V(graph)$Community == community]$name
    cname <- unique(V(graph)[V(graph)$Community == community]$Community_name)
    csize <- length(genes)
    avg_pr <- mean(V(graph)[V(graph)$Community == community]$PageRank, na.rm = TRUE)
    message(sprintf("  [%s] community %s (%s, %d genes)", network_label, community, cname, csize))

    rows_full <- rows_simp <- list()
    for (ont in ontologies) {
      res <- perform_ora(genes, universe, community, cname, ont, sem_list[[ont]])
      if (nrow(res$full) > 0) rows_full[[ont]] <- res$full %>% mutate(Community_size = csize, Avg_PageRank = avg_pr)
      if (nrow(res$simp) > 0) rows_simp[[ont]] <- res$simp %>% mutate(Community_size = csize, Avg_PageRank = avg_pr)
    }
    list(
      full = dplyr::bind_rows(rows_full),
      simp = dplyr::bind_rows(rows_simp)
    )
  }, future.seed = TRUE)

  full_df <- dplyr::bind_rows(lapply(results, `[[`, "full"))
  simp_df <- dplyr::bind_rows(lapply(results, `[[`, "simp"))
  list(full = full_df, simp = simp_df)
}

### Run for A7 and A9 ###
out_dir    <- "~/CESC_Network/7_Modularity/5_2_ORA/"
cutoff_tag <- gsub("\\.", "_", as.character(SIMPLIFY_CUTOFF))

cat("[4] Running A7...\n")
res_A7 <- run_network_ora(graph_A7, "A7", universo, sem_list)
write.table(res_A7$full, file.path(out_dir, "5_2_1_Ora7_louv_max_FULL.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(res_A7$simp, file.path(out_dir, paste0("5_2_2_Ora7_louv_max_SIMPLIFIED_", cutoff_tag, ".tsv")),
            sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("[4] A7 done. FULL: %d rows | SIMPLIFIED: %d rows\n", nrow(res_A7$full), nrow(res_A7$simp)))
cat(sprintf("    Modules with terms (FULL): %d | (SIMP): %d\n",
            length(unique(res_A7$full$Community_Name)), length(unique(res_A7$simp$Community_Name))))

cat("[5] Running A9...\n")
res_A9 <- run_network_ora(graph_A9, "A9", universo, sem_list)
write.table(res_A9$full, file.path(out_dir, "5_2_1_Ora9_louv_max_FULL.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(res_A9$simp, file.path(out_dir, paste0("5_2_2_Ora9_louv_max_SIMPLIFIED_", cutoff_tag, ".tsv")),
            sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("[5] A9 done. FULL: %d rows | SIMPLIFIED: %d rows\n", nrow(res_A9$full), nrow(res_A9$simp)))
cat(sprintf("    Modules with terms (FULL): %d | (SIMP): %d\n",
            length(unique(res_A9$full$Community_Name)), length(unique(res_A9$simp$Community_Name))))

### Summary ###
cat("\n══════════════ SUMMARY ══════════════\n")
cat(sprintf("[A7] modules: %d  | FULL unique terms: %d  | SIMP unique terms: %d\n",
            length(unique(res_A7$full$Community_Name)),
            length(unique(res_A7$full$Description)),
            length(unique(res_A7$simp$Description))))
cat(sprintf("[A9] modules: %d  | FULL unique terms: %d  | SIMP unique terms: %d\n",
            length(unique(res_A9$full$Community_Name)),
            length(unique(res_A9$full$Description)),
            length(unique(res_A9$simp$Description))))
cat(sprintf("Total time: %.1f min\n", as.numeric(difftime(Sys.time(), start_time, units="mins"))))
