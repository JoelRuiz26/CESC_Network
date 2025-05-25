### Cargar librerías ###
library(igraph)
library(tidyverse)
library(clusterProfiler)
library(vroom)
library(org.Hs.eg.db)

### Cargar redes ###
graph_A7 <- read_graph(file = "~/5_Global_Analisis/6_Modularity/Louvain/3_Filtered_Subnetwork_A7.graphml", format = "graphml")
graph_A9 <- read_graph(file = "~/5_Global_Analisis/6_Modularity/Louvain/3_Filtered_Subnetwork_A9.graphml", format = "graphml")

### Cargar universo de genes ###
universo <- vroom("~/2_Prepro_Data_TCGA/A7_A9/Final_A7_A9_Annot.tsv") %>% pull(hgnc_symbol)

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
    ) %>% as.data.frame()
  }, error = function(e) data.frame())
  
  if (nrow(enrich_results) == 0) return(enrich_results)
  
  enrich_results <- enrich_results %>%
    mutate(Community = community_id, Community_Name = community_name, Ontology = ontology)
  
  return(enrich_results)
}

### Extraer comunidades y realizar ORA ###
extract_communities_and_ora <- function(graph, universe) {
  communities <- unique(V(graph)$Community)
  ontologies <- c("BP", "CC", "MF")
  
  ora_results <- map_dfr(communities, function(community) {
    genes <- V(graph)[V(graph)$Community == community]$name
    community_name <- unique(V(graph)[V(graph)$Community == community]$Community_name)
    
    bind_rows(lapply(ontologies, function(ontology) perform_ora(genes, universe, community, community_name, ontology)))
  })
  
  return(ora_results)
}

### Ejecutar ORA para cada red ###
ora_A7 <- extract_communities_and_ora(graph_A7, universo)
ora_A9 <- extract_communities_and_ora(graph_A9, universo)


write.table(ora_A7, file = "~/7_Enrichment/Enrichment_Louvain/Ora7.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ora_A9, file = "~/7_Enrichment/Enrichment_Louvain/Ora9.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
save.image("~/7_Enrichment/Enrichment_Louvain/Enrichment_community.RData")

#load("~/7_Enrichment/Enrichment_community.RData")




#### Filter data, even if the GO terms are higher that 0.05 in p value-adjusted, its necessary additionally thresholds ###
#Add community size even if are implicit in GenRatio#
Comunidades_A7 <- vroom(file = "~/6_Modularity/Louvain/3_A7Louvain_Filtered.tsv")
Comunidades_A7 <- Comunidades_A7 %>%  
  dplyr::select(Community,size) %>% distinct(Community,.keep_all = TRUE)

Comunidades_A9 <- vroom(file = "~/6_Modularity/Louvain/3_A9Louvain_Filtered.tsv")
Comunidades_A9 <- Comunidades_A9 %>%  
  dplyr::select(Community,size) %>% distinct(Community,.keep_all = TRUE)

ora_filtered_A7 <- ora_A7 %>% left_join(Comunidades_A7, by = "Community") #16 community enriched 
ora_filtered_A9 <- ora_A9 %>% left_join(Comunidades_A9, by = "Community") #23 community enriched
#> colnames(ora_A7)
#"Description" (Go term)     "GeneRatio" (genes in term/genes in community)     
#"BgRatio"  (# genes in background in term / #genes in background)      
#"RichFactor" (GeneRatio/BgRatio) high value means better enrichment
#"FoldEnrichment" (overrepresentation grade, high value means stronger enrichment)
#"zScore"  (higher value means higher random deviation)       
#[8] "pvalue"         "p.adjust"       "qvalue"         "geneID"         "Count"          "Community"      "Community_Name"
#[15] "Ontology" 
#Even most of these are related, they arent redundant 

# Filter by counts #
ora_filtered_A7 <- ora_filtered_A7 %>% dplyr::filter(Count >= 3) #14 community enriched
ora_filtered_A9 <- ora_filtered_A9 %>% dplyr::filter(Count >= 3) # 21 community enriched

# Filter by GeneRatio >10% to get more informative enrichment #
# Although GR consider counts, but in a community of 10 genes where count = 1 , GR = 1/10,
# it would be acceptable in this threshold, but count =1 is not too much informative
#make number the gene ratio
ora_filtered_A7$GeneRatio <- sapply(strsplit(ora_filtered_A7$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
ora_filtered_A9$GeneRatio <- sapply(strsplit(ora_filtered_A9$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
#Filter  GR> 0.05
ora_filtered_A7 <- ora_filtered_A7 %>% filter(GeneRatio > 0.05) #12 community enriched
ora_filtered_A9 <- ora_filtered_A9 %>% filter(GeneRatio > 0.05) # 17 community enriched

#Rich factor
#ora_filtered_A7 <- ora_filtered_A7 %>% filter(RichFactor > 0.2) #5 community
#ora_filtered_A9 <- ora_filtered_A9 %>% filter(RichFactor > 0.2) #5 community

# Fold enrichment >2 (p value only takes significance stadistics but not the enrichment magnitude)
ora_filtered_A7 <- ora_filtered_A7 %>% dplyr::filter(FoldEnrichment >= 2) #12 community enriched
ora_filtered_A9 <- ora_filtered_A9 %>% dplyr::filter(FoldEnrichment >= 2) # 17 community enriched

write.table(ora_filtered_A7, file = "~/7_Enrichment/Enrichment_Louvain/Ora7_filtered.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ora_filtered_A9, file = "~/7_Enrichment/Enrichment_Louvain/Ora9_filtered.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

