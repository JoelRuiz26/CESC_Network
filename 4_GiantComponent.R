#Find minimum giant component with max mutual information in CC networks
#Joel Ruiz Hernandez  
#Input: two coexpression networks sorted by MI
#Output: -Edges with giant component is formed and nodes that contain the subgraph
#        -Graph of growing network and threshold of subraph that contains giant component

#Load libraries
library(igraph)
library(dplyr)
library(vroom)
library(ggplot2)

#load(file = "~/5_Global_Analisis/1_GiantComponent.RData")

### Load networks ###
A7_Full <- vroom(file = "~/5_Global_Analisis/0_1_Final_A7_Annot.sort",
                 col_names = c("GenA", "GenB", "MI")) %>% na.omit()
A9_Full <- vroom(file = "~/5_Global_Analisis/0_1_Final_A9_Annot.sort",
                 col_names = c("GenA", "GenB", "MI")) %>% na.omit()
#Posibles enlaces : 128,845,201
#Full edges in both networks:  64,416,925

### Set the minimum Giant Component in your Networks  ###
total_nodes <- length(unique(c(A7_Full$GenA, A7_Full$GenB))) #11351 Both networks
giant_threshold <- floor(0.5 * total_nodes) + 2  #5677


### Function to calculate total nodes and GC each edge  ###
#(Dataframe output is used for graph)#
calculate_metrics <- function(edgelist, max_links = 100000) { #Explore thresholds as you need
                  #The threshold must have in Largest_Component the giant_threshold
  results <- lapply(1:max_links, function(i) {
    graph <- graph_from_data_frame(edgelist[1:i, ], directed = FALSE)
    c(
      Links = i,
      Total_Genes = vcount(graph),
      Largest_Component = max(components(graph)$csize)
    )
  })
  return(as.data.frame(do.call(rbind, results)))
}
# Calcular métricas para ambas redes
metrics_A7 <- calculate_metrics(A7_Full, max_links = 100000)
metrics_A9 <- calculate_metrics(A9_Full, max_links = 100000)


### Subnetwork that contains the minimum Giant component in the full network  ###

#See cut off
giant_link_A7 <- metrics_A7 %>% 
  filter(Largest_Component >= giant_threshold) %>% 
  slice(1)
# Links  Total_Genes Largest_Component
# 26485  6469        5679
giant_link_A9 <- metrics_A9 %>% 
  filter(Largest_Component >= giant_threshold) %>% 
  slice(1)
# Links  Total_Genes Largest_Component
# 39532  6428        5677

#Slice subnetwork
Subnetwork_A7 <- A7_Full %>% dplyr::slice(1:giant_link_A7$Links)
Subnetwork_A9 <- A9_Full %>% dplyr::slice(1:giant_link_A9$Links)

# Save graphs
graph_A7 <- graph_from_data_frame(d = Subnetwork_A7, directed = FALSE)
E(graph_A7)$weight <- Subnetwork_A7$MI
write_graph(graph_A7, file = "~/5_Global_Analisis/1_Subnetwork_A7.graphml", format = "graphml")

graph_A9 <- graph_from_data_frame(d = Subnetwork_A9, directed = FALSE)
E(graph_A9)$weight <- Subnetwork_A9$MI
write_graph(graph_A9, file = "~/5_Global_Analisis/1_Subnetwork_A9.graphml", format = "graphml")

#save image
save.image(file = "~/5_Global_Analisis/1_GiantComponent.RData")

### Graph growing of conected component in network  ###
GC_graph <- ggplot() +
  geom_line(data = metrics_A7, aes(x = Links, y = Largest_Component, color = "VPH-A7"), size = 1) +
  geom_line(data = metrics_A9, aes(x = Links, y = Largest_Component, color = "VPH-A9"), size = 1) +
  geom_hline(aes(yintercept = giant_threshold, 
                 linetype = "Componente gigante mínimo"), 
             color = "black", size = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c("VPH-A7" = "blue", "VPH-A9" = "red")) +
  scale_linetype_manual(values = c("Componente gigante mínimo" = "dashed")) +
  labs(
    title = "Componente más grande conectado en redes de coexpresión de cáncer cervicouterino",
    x = "Enlaces",
    y = "Genes del componente más grande",
    color = "Red",
    linetype = "") +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "Times New Roman"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.line = element_line(size = 0.5, color = "black"))
#Save
ggsave(filename = "~/5_Global_Analisis/1_GC_A7A9.png", plot = GC_graph)


