# 2_Network_Topology.R
# Input: Two graphs from 1_GiantComponent script ( HPV A7 and A9 CC coexp-networks)
# Output: 1_Two graphs (.graphml) built to be comparable between them 
#         2_TSV with main topological metrics both networks


### Libraries ###
pacman::p_load("igraph", "tidyverse", "ggraph", "tidygraph", "svglite", "ggvenn")

### Load graphs ###
Subgraph_A7 <- readRDS('~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_2_Subnetwork_A7_600k.rds')
Subgraph_A9 <- readRDS('~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_2_Subnetwork_A9_600k.rds')


### Build graphs to be comparable ###
#Get nodes
nodes_A7 <- V(Subgraph_A7)$name
nodes_A9 <- V(Subgraph_A9)$name
#Missing nodes each one
missing_in_A7 <- setdiff(nodes_A9, nodes_A7) #219
missing_in_A9 <- setdiff(nodes_A7, nodes_A9) #526
#Add missing nodes
#graph_A7_add <- add_vertices(Subgraph_A7, nv = length(missing_in_A7), name = missing_in_A7) #7968
#graph_A9_add <- add_vertices(Subgraph_A9, nv = length(missing_in_A9), name = missing_in_A9) #7968
#Verify
#print(all(sort(V(graph_A7_add)$name) == sort(V(graph_A9_add)$name))) #[1] TRUE
#Save graphs (Output 1)
#write_graph(graph_A7_add, '~/5_Global_Analisis/2_Parameters_Global/2_Subnetwork_A7_added.graphml', format = 'graphml')
#write_graph(graph_A9_add, '~/5_Global_Analisis/2_Parameters_Global/2_Subnetwork_A9_added.graphml', format = 'graphml')


### Prepare data and functions  ###
#Make list
graphLists <- list("HPV-A7" = Subgraph_A7, "HPV-A9" = Subgraph_A9)
# Seed
set.seed(123)
#Functions
jaccard_edges <- function(g1, g2) {
  length(E(igraph::intersection(g1, g2))) / length(E(igraph::union(g1, g2)))
}
jaccard_nodes <- function(g1, g2) {
  a <- sort(V(g1)$name)
  b <- sort(V(g2)$name)
  length(intersect(a, b)) / length(union(a, b))
}
topological_metrics <- lapply(graphLists, function(graph) { #Topological metrics
  list(
    Numero_de_nodos = vcount(graph),
    Componentes_conectados = components(graph)$no, # Número de componentes conectados
    Grado_promedio = mean(degree(graph)),
    Coeficiente_de_agrupamiento = transitivity(graph, type = "average"),
    Asortatividad = assortativity_degree(graph),
    Tamano_del_componente_principal = max(components(graph)$csize),
    Porcentaje_componente_principal = max(components(graph)$csize) / vcount(graph) * 100,
    Diametro = diameter(graph),
    Longitud_promedio_caminos_cortos = mean_distance(graph))
})


### Results ###
#Jaccard index
NodesJaccard <- jaccard_nodes(Subgraph_A7, Subgraph_A9) # 0.93454
EdgesJaccard <- jaccard_edges(Subgraph_A7, Subgraph_A9) # 0.19949
#Topological table
topological_metrics_df <- bind_rows(topological_metrics, .id = "network") %>%
  mutate(
    jaccard_nodes = ifelse(network == "HPV-A7", NodesJaccard, NodesJaccard),
    jaccard_edges = ifelse(network == "HPV-A7", EdgesJaccard, EdgesJaccard)
  )
t_topological_metrics <- topological_metrics_df %>% dplyr::select(1:6,jaccard_nodes,jaccard_edges,7:10) %>% t() %>% as.data.frame()
colnames(t_topological_metrics) <- c("VPH-A7","VPH-A9")
t_topological_metrics <- t_topological_metrics %>% dplyr::slice(-1)
          #print(t_topological_metrics)
#                                    VPH-A7    VPH-A9
#Numero_de_nodos                      11162     10855
#Componentes_conectados                   8         9
#Grado_promedio                    107.5076  110.5481
#Coeficiente_de_agrupamiento      0.3087146 0.4170420
#Asortatividad                    0.2869697 0.2666389
#jaccard_nodes                      0.93454   0.93454
#jaccard_edges                    0.1994986 0.1994986
#Tamano_del_componente_principal      11148     10839
#Porcentaje_componente_principal   99.87457  99.85260
#Diametro                                10         9
#Longitud_promedio_caminos_cortos  3.038151  3.019805


# Save tsv (Output 2)
vroom::vroom_write(t_topological_metrics, "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_2_1_Toplogical_features_A7A9.tsv")


### Venn Diagram  ###
#Edges each graph
edges_A7 <- apply(as_edgelist(Subgraph_A7), 1, function(e) paste(sort(e), collapse = "-"))
edges_A9 <- apply(as_edgelist(Subgraph_A9), 1, function(e) paste(sort(e), collapse = "-"))
list_edges <- list("HPV-A7" = edges_A7, "HPV-A9" = edges_A9)
# Plot
VennPlot <- ggvenn(list_edges,
                  fill_color = c("#003366", "#FF0000"), 
                  fill_alpha = 0.5,                    
                  stroke_size = 0.5,                   
                  set_name_size = 8,                   
                  text_size = 10,                       
                  show_percentage = FALSE)
print(VennPlot)
# Save graph (Output 3)
ggsave(filename = "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_2_2_VennDiagram_edges.png",
        plot = VennPlot,
        width = 8, height = 8, dpi = 300)
