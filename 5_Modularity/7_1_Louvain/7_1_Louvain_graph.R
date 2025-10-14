#Script for Louvain cluster 
#Input: Networks (.graphml) 
#Output: Networks (graphml) with community named by PageRank algorithm and filterred for enrichment.
#By Joel Ruiz

### Load packages ###
library(igraph)
library(tidyverse)
library(ggplot2)
set.seed(123)

### Load networks ###
graph_A7 <- readRDS('~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_elbow.rds')
graph_A9 <- readRDS('~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_elbow.rds')

###### For A7 network ######
### Louvain ###
Louvain_A7 <- cluster_louvain(graph_A7)        #125 communities
V(graph_A7)$Community <- Louvain_A7$membership 

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
### Louvain ###
Louvain_A9 <- cluster_louvain(graph_A9) #235
V(graph_A9)$Community <- Louvain_A9$membership

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
saveRDS(filtered_graph_A9, file = "~/CESC_Network/5_Modularity/7_1_Louvain/7_1_1_Louvain_A9_graph.rds")
saveRDS(filtered_graph_A7, file = "~/CESC_Network/5_Modularity/7_1_Louvain/7_1_1_Louvain_A7_graph.rds")


### Verify ###
# Dataframe para A7
A7_Louvain_genes <- tibble(
  Community = as.integer(V(filtered_graph_A7)$Community),
  Community_name = as.character(V(filtered_graph_A7)$Community_name),
  gene = as.character(V(filtered_graph_A7)$name)
) 

A7_sizes <- A7_Louvain_genes %>%
  group_by(Community, Community_name) %>%
  summarise(size = n(), .groups = "drop")

ggplot(A7_sizes, aes(x = reorder(Community_name, -size), y = size)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Tamaño de Nodos por Comunidad (A7)",
       x = "Comunidad",
       y = "Tamaño de Nodos") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Dataframe para A9
A9_Louvain_genes <- tibble(
  Community = as.integer(V(filtered_graph_A9)$Community),
  Community_name = as.character(V(filtered_graph_A9)$Community_name),
  gene = as.character(V(filtered_graph_A9)$name)
) 

A9_sizes <- A9_Louvain_genes %>%
  group_by(Community, Community_name) %>%
  summarise(size = n(), .groups = "drop")

ggplot(A9_sizes, aes(x = reorder(Community_name, -size), y = size)) +
  geom_bar(stat = "identity", fill = "coral") +
  labs(title = "Tamaño de Nodos por Comunidad (A9)",
       x = "Comunidad",
       y = "Tamaño de Nodos") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
