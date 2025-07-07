# --------------------------------------------------------
# Parallel growth-curve calculation with 25 cores using chunks
# Author: Joel Ruiz Hernandez + ChatGPT refinement
# --------------------------------------------------------

# ===========================
# Load libraries
# ===========================
library(igraph)
library(dplyr)
library(vroom)
library(ggplot2)
library(future.apply)

# ===========================
# Parallel setup
# ===========================
plan(multisession, workers = 20)

# ===========================
# Load precomputed cutoff
# ===========================
max_cutoff_links <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_1_Max_cutoff_links.rds")
cat("\nLoaded max_cutoff_links from RDS:", max_cutoff_links, "\n")

# ===========================
# Load networks
# ===========================
A7_Full <- vroom(
  "~/CESC_Network/4_Network/ARACNE-multicore/launch/2_7_Full_counts_A7_annot.sort",
  col_names = c("GenA", "GenB", "MI"),
  show_col_types = FALSE
) %>% na.omit()

A9_Full <- vroom(
  "~/CESC_Network/4_Network/ARACNE-multicore/launch/2_8_Full_counts_A9_annot.sort",
  col_names = c("GenA", "GenB", "MI"),
  show_col_types = FALSE
) %>% na.omit()

# ===========================
# Parallel chunked computation
# ===========================
calculate_metrics_chunked <- function(edgelist, max_links, chunk_size = 10000) {
  chunks <- split(1:max_links, ceiling(seq_along(1:max_links)/chunk_size))
  
  results <- future_lapply(chunks, function(indices) {
    local_results <- vector("list", length(indices))
    for (j in seq_along(indices)) {
      i <- indices[j]
      g <- graph_from_data_frame(edgelist[1:i, ], directed = FALSE)
      local_results[[j]] <- c(
        Links = i,
        Total_Genes = vcount(g),
        Largest_Component = max(components(g)$csize)
      )
    }
    do.call(rbind, local_results)
  })
  
  as.data.frame(do.call(rbind, results))
}

# ===========================
# Compute growth curves
# ===========================
cat("\nComputing metrics for A7 in parallel...\n")
metrics_A7 <- calculate_metrics_chunked(A7_Full, max_cutoff_links)

cat("\nComputing metrics for A9 in parallel...\n")
metrics_A9 <- calculate_metrics_chunked(A9_Full, max_cutoff_links)

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
# Find elbows
# ===========================
elbow_idx_A7 <- find_elbow(metrics_A7$Links, metrics_A7$Largest_Component)
elbow_links_A7 <- metrics_A7$Links[elbow_idx_A7]
elbow_genes_A7 <- metrics_A7$Largest_Component[elbow_idx_A7]

elbow_idx_A9 <- find_elbow(metrics_A9$Links, metrics_A9$Largest_Component)
elbow_links_A9 <- metrics_A9$Links[elbow_idx_A9]
elbow_genes_A9 <- metrics_A9$Largest_Component[elbow_idx_A9]

cat("\nDetected thresholds (elbows):\n")
cat("  VPH-A7: ", elbow_links_A7, " edges; ", elbow_genes_A7, " genes\n")
cat("  VPH-A9: ", elbow_links_A9, " edges; ", elbow_genes_A9, " genes\n")

# ===========================
# Define common cutoff
# ===========================
common_cutoff_links <- max(elbow_links_A7, elbow_links_A9)
cat("\nCommon threshold (maximum elbow): ", common_cutoff_links, " edges\n")

# ===========================
# Plot
# ===========================
GC_graph_elbow_common <- ggplot() +
  geom_line(data = metrics_A7, aes(x = Links, y = Largest_Component, color = "VPH-A7"), size = 1.5) +
  geom_line(data = metrics_A9, aes(x = Links, y = Largest_Component, color = "VPH-A9"), size = 1.5) +
  geom_vline(xintercept = common_cutoff_links, linetype = "dashed", color = "black", size = 0.6) +
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

# ===========================
# Save figure
# ===========================
ggsave(
  "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_1_GC_graph_elbow_paralel.png",
  plot = GC_graph_elbow_common,
  width = 12,
  height = 8,
  dpi = 300
)

# ===========================
# Export subnetworks
# ===========================
cat("\nExporting subnetworks at common cutoff...\n")
Subnetwork_A7_common <- A7_Full %>% slice(1:common_cutoff_links)
Subnetwork_A9_common <- A9_Full %>% slice(1:common_cutoff_links)

graph_A7_common <- graph_from_data_frame(d = Subnetwork_A7_common, directed = FALSE)
graph_A9_common <- graph_from_data_frame(d = Subnetwork_A9_common, directed = FALSE)

saveRDS(
  graph_A7_common,
  "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_2_Subnetwork_A7_paralel.rds"
)
saveRDS(
  graph_A9_common,
  "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_2_Subnetwork_A9_paralel.rds"
)

cat("Graphs exported successfully!\n")
