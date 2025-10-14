# 2_Network_Topology.R
# Input: Two graphs from 1_GiantComponent script ( HPV A7 and A9 CC coexp-networks)
# Output: 1_Two graphs (.graphml) built to be comparable between them 
#         2_TSV with main topological metrics both networks


### Libraries ###
pacman::p_load("igraph", "tidyverse", "ggraph", "tidygraph", "svglite", "ggvenn")

### Load graphs ###
Subgraph_A7 <- readRDS('~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_elbow.rds')
Subgraph_A9 <- readRDS('~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_elbow.rds')


### Build graphs to be comparable ###
#Get nodes
nodes_A7 <- V(Subgraph_A7)$name
nodes_A9 <- V(Subgraph_A9)$name
#Missing nodes each one
missing_in_A7 <- setdiff(nodes_A9, nodes_A7) #449
missing_in_A9 <- setdiff(nodes_A7, nodes_A9) #835
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
NodesJaccard <- jaccard_nodes(Subgraph_A7, Subgraph_A9) # 0.8851
EdgesJaccard <- jaccard_edges(Subgraph_A7, Subgraph_A9) # 0.1958
#Topological table
topological_metrics_df <- bind_rows(topological_metrics, .id = "network") %>%
  mutate(
    jaccard_nodes = ifelse(network == "HPV-A7", NodesJaccard, NodesJaccard),
    jaccard_edges = ifelse(network == "HPV-A7", EdgesJaccard, EdgesJaccard)
  )
t_topological_metrics <- topological_metrics_df %>% dplyr::select(1:6,jaccard_nodes,jaccard_edges,7:10) %>% t() %>% as.data.frame()
colnames(t_topological_metrics) <- c("HPV-A7","HPV-A9")
t_topological_metrics <- t_topological_metrics %>% dplyr::slice(-1)
          print(t_topological_metrics)
#          HPV-A7    HPV-A9
#          Componentes_conectados                  76       191
#          Grado_promedio                    23.69526  26.60030
#          Coeficiente_de_agrupamiento      0.2594443 0.5693897
#          Asortatividad                    0.3595462 0.5479928
#          jaccard_nodes                    0.7905064 0.7905064
#          jaccard_edges                    0.1271801 0.1271801
#          Tamano_del_componente_principal      10896      9407
#          Porcentaje_componente_principal   98.58849  95.55104
#          Diametro                                12        17
#          Longitud_promedio_caminos_cortos  4.294790  5.508257
          

# Save tsv (Output 2)
saveRDS(t_topological_metrics, "~/CESC_Network/4_Network_analisis/5_2_Network_topology/5_2_1_Toplogical_features_A7A9.rds")

### Venn Diagram  ###
#Edges each graph
# --- Datos para el Venn (edges) ---
edges_A7 <- apply(as_edgelist(Subgraph_A7), 1, function(e) paste(sort(e), collapse = "-"))
edges_A9 <- apply(as_edgelist(Subgraph_A9), 1, function(e) paste(sort(e), collapse = "-"))
list_edges <- list("HPV-A7" = edges_A7, "HPV-A9" = edges_A9)

# --- Venn con número + porcentaje y TÍTULO ---
VennPlot <- ggvenn(
  list_edges,
  fill_color     = c("#003366", "#FF0000"),
  fill_alpha     = 0.5,
  stroke_size    = 0.5,
  set_name_size  = 8,
  text_size      = 10,
  show_percentage= TRUE,   # muestra % (respecto a la unión)
  label_sep      = "\n",   # cuenta en 1a línea, % en 2a
  digits         = 1
) +
  labs(
    title    = "Jaccard edges") +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 35),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    plot.caption  = element_text(hjust = 0.5, size = 9)
  )

print(VennPlot)


# Save graph (Output 3)
ggsave(filename = "~/CESC_Network/4_Network_analisis/5_2_Network_topology/5_2_2_VennDiagram_edges.png",
        plot = VennPlot,
        width = 8, height = 8, dpi = 300)
