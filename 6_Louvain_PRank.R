#Get Louvain community and pagerank 
#Input: Networks added from 2_Network_topology.R
#Output: Dataframe with Louvain comunities named by highest PageRank and their size

### Load packages ###
library(igraph)
library(tidyverse)

### Read graphs ###
graph_A7 <- read_graph(file = "~/5_Global_Analisis/2_Parameters_Global/2_Subnetwork_A7_added.graphml", 
                       format = "graphml")
graph_A9 <- read_graph(file = "~/5_Global_Analisis/2_Parameters_Global/2_Subnetwork_A9_added.graphml", 
                       format = "graphml")

### Get Louvain community and PageRank  ###
#For A7 network
A7_Louvain_PR<- as_tibble(V(graph_A7))  %>%                #ID number nodes
  mutate(name = V(graph_A7)$name,                          #Node names  
         Community = cluster_louvain(graph_A7)$membership, #Louvain community
         PageRank = page_rank(graph_A7)$vector) %>%        #Pagerank each node un hole network
  group_by(Community) %>%
  summarise(Community_name = name[which.max(PageRank)],    #Highest PageRank
            size = n())                          #Number of nodes by community
        #Comunidad Comunidad_Nombre    Tamano
        #<dbl>      <chr>               <int>
        #1         1 FERMT3              209
        #2         2 RAD1                 70
        #3         3 THEM5                 3
        #4         4 VTA1                300

#For A9 network
A9_Louvain_PR<- as_tibble(V(graph_A9))  %>%     
  mutate(name = V(graph_A9)$name,                      
         Community = cluster_louvain(graph_A9)$membership, 
         PageRank = page_rank(graph_A9)$vector) %>%        
  group_by(Community) %>%
  summarise(Community_name = name[which.max(PageRank)],   
            size = n()) 

### Save df ###
write_tsv(A7_Louvain_PR, "~/5_Global_Analisis/3_Modularity/3_A7Comm_size_PageName.tsv")
write_tsv(A9_Louvain_PR, "~/5_Global_Analisis/3_Modularity/3_A9Comm_size_PageName.tsv")





