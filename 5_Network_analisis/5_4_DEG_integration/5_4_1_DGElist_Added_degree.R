############################################################
## Annotate all genes with:
## - Degree
## - PageRank
## - Degree ranking
## - PageRank ranking
##
## Rank 1 = highest value
############################################################

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(vroom)
})

############################################################
## Load graphs
############################################################
Subgraph_A7 <- readRDS(
  "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds"
)

Subgraph_A9 <- readRDS(
  "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds"
)

############################################################
## Load tables
############################################################
DE_A7 <- vroom(
  "~/CESC_Network/3_DGE/3_1_limma_A7_vs_SNT_allGenes.tsv",
  show_col_types = FALSE
)

DE_A9 <- vroom(
  "~/CESC_Network/3_DGE/3_1_limma_A9_vs_SNT_allGenes.tsv",
  show_col_types = FALSE
)

############################################################
## Function
############################################################
annotate_network <- function(graph, table) {
  
  degree_vec   <- degree(graph)
  pagerank_vec <- page_rank(graph)$vector
  
  network_stats <- tibble(
    Gene = names(degree_vec),
    
    degree = as.numeric(degree_vec),
    
    pagerank = as.numeric(pagerank_vec),
    
    degree_rank =
      rank(-degree_vec, ties.method = "min"),
    
    pagerank_rank =
      rank(-pagerank_vec, ties.method = "min")
  )
  
  table %>%
    left_join(network_stats, by = "Gene")
}

############################################################
## Annotate tables
############################################################
DE_A7_annot <- annotate_network(Subgraph_A7, DE_A7)

DE_A9_annot <- annotate_network(Subgraph_A9, DE_A9)

############################################################
## Save results
############################################################
vroom_write(
  DE_A7_annot,
  "~/CESC_Network/5_Network_analisis/5_4_DEG_integration/5_4_1_DEG_A7_Degree.tsv",
  delim = "\t"
)

vroom_write(
  DE_A9_annot,
  "~/CESC_Network/5_Network_analisis/5_4_DEG_integration/5_4_1_DEG_A9_Degree.tsv",
  delim = "\t"
)


############################################################
## Spearman correlations
############################################################
# Degree importance vs t-statistic
cor_degree <- cor.test(
  -DE_A9_annot$degree_rank,
  DE_A9_annot$t,
  method = "spearman",
  exact = FALSE
)

# PageRank importance vs t-statistic
cor_pagerank <- cor.test(
  -DE_A9_annot$pagerank_rank,
  DE_A9_annot$t,
  method = "spearman",
  exact = FALSE
)

############################################################
## Results
############################################################
cor_degree
# 
# Spearman's rank correlation rho
# 
# data:  -DE_A9_annot$degree_rank and DE_A9_annot$t
# S = 6.9169e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2231052 
# 
cor_pagerank
# 
# 	Spearman's rank correlation rho
# 
# data:  -DE_A9_annot$pagerank_rank and DE_A9_annot$t
# S = 7.0553e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.2075581 