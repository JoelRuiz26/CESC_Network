#Script for Infomap cluster 
#Input: Networks (.graphml) 
#Output: Networks (graphml) with community named by PageRank algorithm and filterred for enrichment.
#By Joel Ruiz

### Load packages ###
library(igraph)
library(tidyverse)
library(ggplot2)
set.seed(123)

### Load networks ###
graph_A7 <- readRDS('~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds')
graph_A9 <- readRDS('~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds')

###### For A7 network ######
### Infomap ###
Infomap_A7 <- cluster_infomap(graph_A7)        # comunidades
V(graph_A7)$Community <- Infomap_A7$membership 

### PageRank ###
page_ranks_A7 <- page_rank(graph_A7)$vector
V(graph_A7)$PageRank <- page_ranks_A7
#Highest PageRank in community
community_PR_names_A7 <- tapply(V(graph_A7)$name, V(graph_A7)$Community, function(names) {
  ranks <- page_ranks_A7[V(graph_A7)$name %in% names]
  names[which.max(ranks)]
})
V(graph_A7)$Community_name <- community_PR_names_A7[as.character(V(graph_A7)$Community)] #Add to graph

### Filter nodes for enrichment ###
V(graph_A7)$Community_size <- table(V(graph_A7)$Community)[as.character(V(graph_A7)$Community)] #community size
nodes_to_keep_A7 <- V(graph_A7)[V(graph_A7)$Community_size >= 10]
filtered_graph_A7 <- induced_subgraph(graph_A7, vids = nodes_to_keep_A7)


###### For A9 network ######
### Infomap ###
Infomap_A9 <- cluster_infomap(graph_A9)
V(graph_A9)$Community <- Infomap_A9$membership

### PageRank ###
page_ranks_A9 <- page_rank(graph_A9)$vector
V(graph_A9)$PageRank <- page_ranks_A9
#Highest PageRank in community
community_PR_names_A9 <- tapply(V(graph_A9)$name, V(graph_A9)$Community, function(names) {
  ranks <- page_ranks_A9[V(graph_A9)$name %in% names]
  names[which.max(ranks)]
})
V(graph_A9)$Community_name <- community_PR_names_A9[as.character(V(graph_A9)$Community)]

### Filter nodes for enrichment ###
V(graph_A9)$Community_size <- table(V(graph_A9)$Community)[as.character(V(graph_A9)$Community)]
nodes_to_keep_A9 <- V(graph_A9)[V(graph_A9)$Community_size >= 10]
filtered_graph_A9 <- induced_subgraph(graph_A9, vids = nodes_to_keep_A9)


### Save filtered graphs ###
saveRDS(filtered_graph_A9, file = "~/CESC_Network/5_Modularity/5_1_Modules/5_1_1_Infomap_A9_graph_max.rds")
saveRDS(filtered_graph_A7, file = "~/CESC_Network/5_Modularity/5_1_Modules/5_1_1_Infomap_A7_graph_max.rds")


### Verify ###
# Dataframe para A7
A7_Infomap_genes <- tibble(
  Community = as.integer(V(filtered_graph_A7)$Community),
  Community_name = as.character(V(filtered_graph_A7)$Community_name),
  gene = as.character(V(filtered_graph_A7)$name)
) 

A7_sizes <- A7_Infomap_genes %>%
  group_by(Community, Community_name) %>%
  summarise(size = n(), .groups = "drop")

# Dataframe para A9
A9_Infomap_genes <- tibble(
  Community = as.integer(V(filtered_graph_A9)$Community),
  Community_name = as.character(V(filtered_graph_A9)$Community_name),
  gene = as.character(V(filtered_graph_A9)$name)
) 

A9_sizes <- A9_Infomap_genes %>%
  group_by(Community, Community_name) %>%
  summarise(size = n(), .groups = "drop")

# === Plots como objetos ===
pA7 <- ggplot(A7_sizes, aes(x = reorder(Community_name, -size), y = size)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Tamaño de Nodos por Comunidad (A7) - Infomap",
       x = "Comunidad", y = "Tamaño de Nodos") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pA9 <- ggplot(A9_sizes, aes(x = reorder(Community_name, -size), y = size)) +
  geom_bar(stat = "identity", fill = "coral") +
  labs(title = "Tamaño de Nodos por Comunidad (A9) - Infomap",
       x = "Comunidad", y = "Tamaño de Nodos") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# === Combinar en una sola grilla y guardar en PNG ===
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

combined <- pA7 + pA9 + patchwork::plot_layout(ncol = 2)

outdir <- "~/CESC_Network/5_Modularity/5_1_Modules/"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

outfile <- file.path(outdir, "5_1_2_Infomap_community_sizes_A7_A9_max.png")
ggplot2::ggsave(filename = outfile, plot = combined, width = 18, height = 9, dpi = 300, bg = "white")

message("[OK] PNG guardado en: ", outfile)

#A7_SIZES en promedio son de 300 genes, pero hay uno de 3mil, en A9 no estan


