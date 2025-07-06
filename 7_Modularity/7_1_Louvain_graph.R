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
graph_A7 <- read_graph("~/5_Global_Analisis/2_Parameters_Global/2_Subnetwork_A7_added.graphml", format = "graphml")
graph_A9 <- read_graph("~/5_Global_Analisis/2_Parameters_Global/2_Subnetwork_A9_added.graphml", format = "graphml")
#7968 nodes each network

###### For A7 network ######
### Louvain ###
Louvain_A7 <- cluster_louvain(graph_A7)        #1864 communities
V(graph_A7)$Community <- Louvain_A7$membership #Add to graph

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
nodes_to_keep_A7 <- V(graph_A7)[V(graph_A7)$Community_size >= 6]
filtered_graph_A7 <- induced_subgraph(graph_A7, vids = nodes_to_keep_A7)


###### For A9 network ######
### Louvain ###
Louvain_A9 <- cluster_louvain(graph_A9)
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
nodes_to_keep_A9 <- V(graph_A9)[V(graph_A9)$Community_size >= 6]
filtered_graph_A9 <- induced_subgraph(graph_A9, vids = nodes_to_keep_A9)


### Save filtered graphs ###
write_graph(filtered_graph_A9, file = "~/5_Global_Analisis/3_Modularity/Louvain/3_Filtered_Subnetwork_A9.graphml", format = "graphml")
write_graph(filtered_graph_A7, file = "~/5_Global_Analisis/3_Modularity/Louvain/3_Filtered_Subnetwork_A7.graphml", format = "graphml")


### Verify ###
# Dataframe para A7
#A7_Louvain_genes <- tibble(
#  Community = as.integer(V(filtered_graph_A7)$Community),           # Community number
#  Community_name = as.character(V(filtered_graph_A7)$Community_name), # Highest pagerank Community name
#  size = as.integer(V(filtered_graph_A7)$Community_size),          # Number of nodes by community
#  gene = as.character(V(filtered_graph_A7)$name))                   # Gene name

# Dataframe para A9
#A9_Louvain_genes <- tibble(
#  Community = as.integer(V(filtered_graph_A9)$Community),         
#  Community_name = as.character(V(filtered_graph_A9)$Community_name), 
#  size = as.integer(V(filtered_graph_A9)$Community_size),          
#  gene = as.character(V(filtered_graph_A9)$name))                   

#Save dataframes
#write_tsv(A7_Louvain_genes, "~/5_Global_Analisis/3_Modularity/Louvain/3_A7Louvain_Filtered.tsv")
#write_tsv(A9_Louvain_genes, "~/5_Global_Analisis/3_Modularity/Louvain/3_A9Louvain_Filtered.tsv")

### Module size plot  ###
# For A7
ggplot(A7_Louvain_genes, aes(x = reorder(Community_name, -size), y = size)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Tamaño de Nodos por Comunidad (A7)",
       x = "Comunidad",
       y = "Tamaño de Nodos") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# For A9
ggplot(A9_Louvain_genes, aes(x = reorder(Community_name, -size), y = size)) +
  geom_bar(stat = "identity", fill = "coral") +
  labs(title = "Tamaño de Nodos por Comunidad (A9)",
       x = "Comunidad",
       y = "Tamaño de Nodos") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

