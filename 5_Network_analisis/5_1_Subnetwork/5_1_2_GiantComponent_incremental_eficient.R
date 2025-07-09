# --------------------------------------------------------
# Find minimum giant component with max mutual information 
# in Cervical Cancer Co-expression Networks
# Author: Joel Ruiz Hernandez
#
# Input:
#   - Two co-expression networks sorted by Mutual Information (MI)
#
# Output:
#   - Edges required to form the giant component
#   - Nodes included in the resulting subgraph
#   - Growth plot of the largest connected component
#
# Method:
#   We selected the network threshold using the elbow method,
#   identifying the point of maximum curvature in the giant
#   component growth curve (Thorndike, 1953; Barabási, 2016). 
#   This approach identifies the minimal link density needed
#   to achieve global connectivity while limiting noise.
# --------------------------------------------------------

# ===========================
# Load libraries
# ===========================
library(igraph)
library(dplyr)
library(vroom)
library(ggplot2)

# ===========================
# Load precomputed max cutoff
# ===========================
max_cutoff_links <- readRDS(
  "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_1_Max_cutoff_links.rds"
)
cat("\nLoaded max_cutoff_links from RDS:", max_cutoff_links, "\n")

# ===========================
# Load networks
# ===========================
A7_Full <- vroom(
  file = "~/CESC_Network/4_Network/ARACNE-multicore/launch/2_7_Full_counts_A7_annot.sort",
  col_names = c("GenA", "GenB", "MI")
) %>% na.omit()

A9_Full <- vroom(
  file = "~/CESC_Network/4_Network/ARACNE-multicore/launch/2_8_Full_counts_A9_annot.sort",
  col_names = c("GenA", "GenB", "MI")
) %>% na.omit()

# ===========================
# Incremental growth metric calculation
# ===========================
calculate_metrics_incremental <- function(edgelist, max_links) {
  cat("\nStarting incremental calculation up to", max_links, "edges...\n")
  g <- make_empty_graph(directed = FALSE)
  results <- vector("list", max_links)
  
  for (i in 1:max_links) {
    edge_nodes <- c(edgelist$GenA[i], edgelist$GenB[i])
    # Verificar si hay nodos nuevos
    new_nodes <- edge_nodes[!edge_nodes %in% V(g)$name]
    if (length(new_nodes) > 0) {
      g <- add_vertices(g, nv = length(new_nodes), name = new_nodes)
    }
    
    # Agregar la arista
    g <- add_edges(g, edge_nodes)
    
    # Calcular componentes
    comps <- components(g)
    results[[i]] <- c(
      Links = i,
      Total_Genes = vcount(g),
      Largest_Component = max(comps$csize)
    )
    
    if (i %% 100000 == 0) cat(i, "edges processed\n")
  }
  
  as.data.frame(do.call(rbind, results))
}


# ===========================
# Compute growth metrics
# ===========================
metrics_A7 <- calculate_metrics_incremental(A7_Full, max_cutoff_links)

metrics_A9 <- calculate_metrics_incremental(A9_Full, max_cutoff_links)

# ===========================
# Elbow detection function
# ===========================
find_elbow <- function(x, y) {
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  
  line_vec <- c(1, 1) / sqrt(2)
  points_matrix <- cbind(x_norm, y_norm)
  projection_lengths <- points_matrix %*% line_vec
  projection_points <- projection_lengths %*% t(line_vec)
  distances <- sqrt(rowSums((points_matrix - projection_points)^2))
  
  which.max(distances)
}

# ===========================
# Compute elbow thresholds
# ===========================
elbow_idx_A7 <- find_elbow(metrics_A7$Links, metrics_A7$Largest_Component)
elbow_links_A7 <- metrics_A7$Links[elbow_idx_A7]
elbow_genes_A7 <- metrics_A7$Largest_Component[elbow_idx_A7]

elbow_idx_A9 <- find_elbow(metrics_A9$Links, metrics_A9$Largest_Component)
elbow_links_A9 <- metrics_A9$Links[elbow_idx_A9]
elbow_genes_A9 <- metrics_A9$Largest_Component[elbow_idx_A9]

cat("\nDetected thresholds (elbows):\n")
cat("  VPH-A7:", elbow_links_A7, "edges;", elbow_genes_A7, "genes\n")
cat("  VPH-A9:", elbow_links_A9, "edges;", elbow_genes_A9, "genes\n")

# ===========================
# Define common cutoff
# ===========================
common_cutoff_links <- max(elbow_links_A7, elbow_links_A9)
cat("\nCommon threshold (maximum elbow):", common_cutoff_links, "edges\n")

# ===========================
# Plot growth curves
# ===========================
GC_graph_elbow_common <- ggplot() +
  geom_line(data = metrics_A7, aes(x = Links, y = Largest_Component, color = "VPH-A7"), size = 1.2) +
  geom_line(data = metrics_A9, aes(x = Links, y = Largest_Component, color = "VPH-A9"), size = 1.2) +
  geom_vline(xintercept = common_cutoff_links, linetype = "dashed", color = "black", size = 0.8) +
  labs(
    title = "Growth of Largest Connected Component in Cervical Cancer Co-expression Networks",
    subtitle = paste0("Common threshold (maximum elbow): ", format(common_cutoff_links, big.mark = ","), " edges"),
    x = "Number of Edges",
    y = "Number of Genes in Largest Connected Component",
    color = "Network"
  ) +
  scale_color_manual(values = c("VPH-A7" = "blue", "VPH-A9" = "red")) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "top"
  )

ggsave(
  filename = "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_1_GC_graph_elbow_common_incremental_eficient.png",
  plot = GC_graph_elbow_common,
  width = 12,
  height = 8,
  dpi = 300
)

# ===========================
# Export subnetworks at common cutoff
# ===========================
cat("\nExporting subnetworks at common cutoff...\n")

Subnetwork_A7_common <- A7_Full %>% slice(1:common_cutoff_links)
Subnetwork_A9_common <- A9_Full %>% slice(1:common_cutoff_links)

graph_A7_common <- graph_from_data_frame(d = Subnetwork_A7_common, directed = FALSE)
graph_A9_common <- graph_from_data_frame(d = Subnetwork_A9_common, directed = FALSE)

saveRDS(
  graph_A7_common, 
  file = "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_2_Subnetwork_A7_common_incremental_eficient.rds"
)

saveRDS(
  graph_A9_common, 
  file = "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_2_Subnetwork_A9_common_incremental_eficient.rds"
)

save.image("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_3_Image_GiantComponent.RData")
cat("Graphs exported successfully!\n")
