############################################################
## Anotar DGE con PageRank normalizado y Degree (z-score)
## y clasificar en High/Low
############################################################

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(vroom)
})

set.seed(123)

############################################################
## 1) Load graphs
############################################################
Subgraph_A7 <- readRDS('~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds')
Subgraph_A9 <- readRDS('~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds')

############################################################
## 2) Load DGE tables
############################################################
DE_A7_vs_SNT <- vroom::vroom(
  "~/CESC_Network/2_DGE/2_1_limma_A7_vs_SNT_allGenes.tsv",
  show_col_types = FALSE
)

DE_A9_vs_SNT <- vroom::vroom(
  "~/CESC_Network/2_DGE/2_1_limma_A9_vs_SNT_allGenes.tsv",
  show_col_types = FALSE
)

############################################################
## 3) Función: PageRank normalizado + Degree + z-scores
############################################################
get_PR_deg_stats <- function(g) {
  # PageRank (no dirigido; ajusta directed si tu red lo requiere)
  pr_obj <- igraph::page_rank(g, directed = FALSE)
  pr_vec <- pr_obj$vector  # named vector
  
  # Normalización min-max a [0, 1]
  pr_min <- min(pr_vec)
  pr_max <- max(pr_vec)
  pr_norm <- (pr_vec - pr_min) / (pr_max - pr_min)
  
  # Degree
  deg_vec <- igraph::degree(g, mode = "all")
  
  # z-scores
  pr_z  <- as.numeric(scale(pr_norm))
  deg_z <- as.numeric(scale(deg_vec))
  
  tibble(
    gene          = names(pr_vec),
    pagerank      = as.numeric(pr_vec),
    pagerank_norm = as.numeric(pr_norm),
    pagerank_z    = pr_z,
    degree        = as.numeric(deg_vec),
    degree_z      = deg_z
  )
}

############################################################
## 4) Obtener métricas para A7 y A9
############################################################
stats_A7 <- get_PR_deg_stats(Subgraph_A7)
stats_A9 <- get_PR_deg_stats(Subgraph_A9)

############################################################
## 5) Clasificar en High / Low según z-score (> 3)
############################################################
stats_A7 <- stats_A7 %>%
  mutate(
    Degree_group = case_when(
      degree_z > 2           ~ "High_degree",
      TRUE                   ~ "Low_degree"
    ),
    Pagerank_group = case_when(
      pagerank_z > 2         ~ "High_pagerank",
      TRUE                   ~ "Low_pagerank"
    )
  )

stats_A9 <- stats_A9 %>%
  mutate(
    Degree_group = case_when(
      degree_z > 2           ~ "High_degree",
      TRUE                   ~ "Low_degree"
    ),
    Pagerank_group = case_when(
      pagerank_z > 2         ~ "High_pagerank",
      TRUE                   ~ "Low_pagerank"
    )
  )

############################################################
## 6) Anotar tablas DGE con las categorías High/Low
############################################################

DE_A7_vs_SNT_annot <- DE_A7_vs_SNT %>%
  left_join(
    stats_A7 %>%
      select(
        gene,
        pagerank_norm,
        pagerank_z,
        degree,
        degree_z,
        Degree_group,
        Pagerank_group
      ),
    by = c("Gene" = "gene")
  ) %>%
  mutate(
    Degree_group = ifelse(
      is.na(Degree_group), "Not_in_network", Degree_group
    ),
    Pagerank_group = ifelse(
      is.na(Pagerank_group), "Not_in_network", Pagerank_group
    )
  )

DE_A9_vs_SNT_annot <- DE_A9_vs_SNT %>%
  left_join(
    stats_A9 %>%
      select(
        gene,
        pagerank_norm,
        pagerank_z,
        degree,
        degree_z,
        Degree_group,
        Pagerank_group
      ),
    by = c("Gene" = "gene")
  ) %>%
  mutate(
    Degree_group = ifelse(
      is.na(Degree_group), "Not_in_network", Degree_group
    ),
    Pagerank_group = ifelse(
      is.na(Pagerank_group), "Not_in_network", Pagerank_group
    )
  )

############################################################
## 7) Guardar resultados
############################################################
vroom_write(
  DE_A7_vs_SNT_annot,
  "~/CESC_Network/4_Network_analisis/4_3_limma_A7_vs_SNT_allGenes_annot_PR_DEG.tsv"
)

vroom_write(
  DE_A9_vs_SNT_annot,
  "~/CESC_Network/4_Network_analisis/4_3_limma_A9_vs_SNT_allGenes_annot_PR_DEG.tsv"
)
