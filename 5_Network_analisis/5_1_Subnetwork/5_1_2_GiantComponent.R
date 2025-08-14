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

library(igraph)
library(dplyr)
library(vroom)
library(ggplot2)

# --------------------------------------------------------
# Load precomputed maximum cutoff value (if available)
# --------------------------------------------------------
# NOTE: If you don't have an RDS, we will cap by the file length later.
max_cutoff_links <- 1000000
#max_cutoff_links <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_1_Max_cutoff_links.rds")
cat("\n[INFO] Requested max_cutoff_links:", max_cutoff_links, "\n")

# --------------------------------------------------------
# Load networks (already sorted by MI)
# --------------------------------------------------------
A7_Full <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_A7_edgelist_allnodes.rds")
A9_Full <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_A9_edgelist_allnodes.rds")

# Effective limits to avoid out-of-range indexing
imax_A7 <- min(max_cutoff_links, nrow(A7_Full))
imax_A9 <- min(max_cutoff_links, nrow(A9_Full))
cat("[INFO] Usable edges -> A7:", imax_A7, "| A9:", imax_A9, "\n")

# --------------------------------------------------------
# Fast incremental growth metrics (sampling)
# --------------------------------------------------------
# Instead of rebuilding the graph at each i, we keep a single graph that grows.
# We record metrics every 'step' edges (coarse curve) to be fast,
# and later we refine at step=1 only around the elbow window.

calculate_metrics_incremental <- function(edgelist, imax, step = 1000L) {
  # Build an empty graph and add edges one by one
  g <- make_empty_graph(directed = FALSE)
  out <- vector("list", length = ceiling(imax/step) + 1L)
  ptr <- 1L
  
  # Helper to add vertices if new
  add_vertices_if_new <- function(g, nodes) {
    new_nodes <- nodes[!(nodes %in% V(g)$name)]
    if (length(new_nodes) > 0L) {
      g <- add_vertices(g, nv = length(new_nodes), name = new_nodes)
    }
    g
  }
  
  for (i in seq_len(imax)) {
    e <- c(edgelist$GenA[i], edgelist$GenB[i])
    g <- add_vertices_if_new(g, e)
    g <- add_edges(g, e)
    
    if (i %% step == 0L || i == imax) {
      comps <- components(g)
      out[[ptr]] <- c(
        Links = i,
        Total_Genes = vcount(g),
        Largest_Component = max(comps$csize)
      )
      ptr <- ptr + 1L
    }
  }
  as.data.frame(do.call(rbind, out))
}

# --------------------------------------------------------
# Elbow detection (normalized distance to the diagonal)
# --------------------------------------------------------
find_elbow <- function(x, y) {
  # Normalize to [0,1] to make the scale comparable
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  # Project to the (1,1) direction and measure the orthogonal distance
  line_vec <- c(1, 1) / sqrt(2)
  points_matrix <- cbind(x_norm, y_norm)
  proj_len <- points_matrix %*% line_vec
  proj_pts <- proj_len %*% t(line_vec)
  dists <- sqrt(rowSums((points_matrix - proj_pts)^2))
  which.max(dists)
}

# --------------------------------------------------------
# Compute growth metrics (coarse + refine around elbow)
# --------------------------------------------------------
cat("\n[INFO] Computing coarse curves (step = 1000)...\n")
metrics_A7_coarse <- calculate_metrics_incremental(A7_Full, imax = imax_A7, step = 1000L)

metrics_A9_coarse <- calculate_metrics_incremental(A9_Full, imax = imax_A9, step = 1000L)


# Elbow on coarse curves
elbow_idx_A7_c <- find_elbow(metrics_A7_coarse$Links, metrics_A7_coarse$Largest_Component)
elbow_idx_A9_c <- find_elbow(metrics_A9_coarse$Links, metrics_A9_coarse$Largest_Component)

elbow_links_A7_c <- metrics_A7_coarse$Links[elbow_idx_A7_c]
elbow_links_A9_c <- metrics_A9_coarse$Links[elbow_idx_A9_c]
cat("[INFO] Coarse elbows at ~", elbow_links_A7_c, "(A7) and ~", elbow_links_A9_c, "(A9)\n")

# Refine +/- window around the coarse elbow with step = 1
refine_window <- 2000L
ref_range_A7 <- c(max(1L, elbow_links_A7_c - refine_window), min(imax_A7, elbow_links_A7_c + refine_window))
ref_range_A9 <- c(max(1L, elbow_links_A9_c - refine_window), min(imax_A9, elbow_links_A9_c + refine_window))

cat("[INFO] Refining A7 in [", ref_range_A7[1], ", ", ref_range_A7[2], "] (step = 1)\n", sep = "")
metrics_A7_refine <- calculate_metrics_incremental(A7_Full, imax = ref_range_A7[2], step = 1L) %>%
  dplyr::filter(Links >= ref_range_A7[1])

cat("[INFO] Refining A9 in [", ref_range_A9[1], ", ", ref_range_A9[2], "] (step = 1)\n", sep = "")
metrics_A9_refine <- calculate_metrics_incremental(A9_Full, imax = ref_range_A9[2], step = 1L) %>%
  dplyr::filter(Links >= ref_range_A9[1])


# Merge: keep coarse outside the refine windows, and dense inside
metrics_A7 <- dplyr::bind_rows(
  dplyr::filter(metrics_A7_coarse, Links < ref_range_A7[1] | Links > ref_range_A7[2]),
  metrics_A7_refine
) %>% dplyr::arrange(Links)

metrics_A9 <- dplyr::bind_rows(
  dplyr::filter(metrics_A9_coarse, Links < ref_range_A9[1] | Links > ref_range_A9[2]),
  metrics_A9_refine
) %>% dplyr::arrange(Links)

# --------------------------------------------------------
# Compute elbow thresholds (final, on refined curves)
# --------------------------------------------------------
elbow_idx_A7 <- find_elbow(metrics_A7$Links, metrics_A7$Largest_Component)
elbow_links_A7 <- metrics_A7$Links[elbow_idx_A7]
elbow_genes_A7 <- metrics_A7$Largest_Component[elbow_idx_A7]

elbow_idx_A9 <- find_elbow(metrics_A9$Links, metrics_A9$Largest_Component)
elbow_links_A9 <- metrics_A9$Links[elbow_idx_A9]
elbow_genes_A9 <- metrics_A9$Largest_Component[elbow_idx_A9]

cat("\nDetected thresholds (elbows):\n")
cat("  HPV-A7: ", elbow_links_A7, " edges; ", elbow_genes_A7, " genes\n", sep = "")
cat("  HPV-A9: ", elbow_links_A9, " edges; ", elbow_genes_A9, " genes\n", sep = "")

# --------------------------------------------------------
# Define common cutoff (maximum of the two elbows)
# --------------------------------------------------------
common_cutoff_links <- max(elbow_links_A7, elbow_links_A9)
cat("\nCommon threshold (maximum elbow): ", common_cutoff_links, " edges\n", sep = "")

# --------------------------------------------------------
# Real-scale geometry to draw chord/perpendiculars for each curve
# --------------------------------------------------------
elbow_geom <- function(x, y, idx) {
  # Returns: chord endpoints, elbow point, and its orthogonal projection
  A <- c(x[1], y[1]); B <- c(x[length(x)], y[length(y)])
  AB <- B - A; AB_len2 <- sum(AB^2)
  P  <- c(x[idx], y[idx])
  t  <- ((P[1]-A[1])*AB[1] + (P[2]-A[2])*AB[2]) / AB_len2
  Pproj <- c(A[1] + t*AB[1], A[2] + t*AB[2])
  list(x_start=A[1], y_start=A[2], x_end=B[1], y_end=B[2],
       x_elbow=P[1], y_elbow=P[2], x_proj=Pproj[1], y_proj=Pproj[2])
}

geom_A7 <- elbow_geom(metrics_A7$Links, metrics_A7$Largest_Component, elbow_idx_A7)
geom_A9 <- elbow_geom(metrics_A9$Links, metrics_A9$Largest_Component, elbow_idx_A9)


# --------------------------------------------------------
# Plot growth curves with full elbow annotations
# --------------------------------------------------------
max_x <- max(metrics_A7$Links, metrics_A9$Links)
max_y <- max(metrics_A7$Largest_Component, metrics_A9$Largest_Component)

GC_graph_elbow_common <- ggplot() +
  geom_line(data = metrics_A7, aes(x = Links, y = Largest_Component, color = "HPV-A7"), linewidth = 1.5) +
  geom_line(data = metrics_A9, aes(x = Links, y = Largest_Component, color = "HPV-A9"), linewidth = 1.5) +
  
  # --- Chords (reference straight line) for each network ---
  annotate("segment",
           x = geom_A7$x_start, y = geom_A7$y_start,
           xend = geom_A7$x_end, yend = geom_A7$y_end,
           linetype = "dotdash", linewidth = 0.6, alpha = 0.7, color = "blue") +
  annotate("segment",
           x = geom_A9$x_start, y = geom_A9$y_start,
           xend = geom_A9$x_end, yend = geom_A9$y_end,
           linetype = "dotdash", linewidth = 0.6, alpha = 0.7, color = "red") +
  
  # --- Perpendiculars from elbow to chord ---
  annotate("segment",
           x = geom_A7$x_elbow, y = geom_A7$y_elbow,
           xend = geom_A7$x_proj,  yend = geom_A7$y_proj,
           linetype = "dotted", linewidth = 0.7, color = "grey30") +
  annotate("segment",
           x = geom_A9$x_elbow, y = geom_A9$y_elbow,
           xend = geom_A9$x_proj,  yend = geom_A9$y_proj,
           linetype = "dotted", linewidth = 0.7, color = "grey30") +
  
  # --- Elbow points ---
  annotate("point", x = geom_A7$x_elbow, y = geom_A7$y_elbow,
           size = 2.6, shape = 21, stroke = 1, fill = "white", color = "blue") +
  annotate("point", x = geom_A9$x_elbow, y = geom_A9$y_elbow,
           size = 2.6, shape = 21, stroke = 1, fill = "white", color = "red") +
  
  # --- Vertical lines at each elbow ---
  geom_vline(xintercept = elbow_links_A7, color = "blue", linetype = "longdash", linewidth = 0.3) +
  geom_vline(xintercept = elbow_links_A9, color = "red",  linetype = "longdash", linewidth = 0.3) +
  
  # --- Common cutoff (black) ---
  geom_vline(xintercept = common_cutoff_links, linetype = "dashed", color = "black", linewidth = 0.7) +
  
  # --- Labels ---
  annotate("text", x = geom_A7$x_elbow, y = geom_A7$y_elbow,
           label = paste0("A7 elbow\nx=", format(elbow_links_A7, big.mark=",")),
           hjust = -0.05, vjust = -0.8, size = 3.6, color = "blue") +
  annotate("text", x = geom_A9$x_elbow, y = geom_A9$y_elbow,
           label = paste0("A9 elbow\nx=", format(elbow_links_A9, big.mark=",")),
           hjust = -0.05, vjust = -0.8, size = 3.6, color = "red") +
  annotate("text", x = common_cutoff_links, y = max_y*0.06,
           label = paste0("Common cutoff = ", format(common_cutoff_links, big.mark=",")),
           angle = 90, vjust = -0.5, size = 4.0, color = "black") +
  labs(
    title = "Growth of Largest Connected Component in Cervical Cancer Co-expression Networks",
#    subtitle = "Elbow method: chord (dot-dash), perpendicular (dotted), elbows (circles), common cut (black line)",
    x = "Number of Edges",
    y = "Number of Genes in Largest Connected Component",
    color = "Network"
  ) +
  scale_color_manual(values = c("HPV-A7" = "blue", "HPV-A9" = "red")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max_x)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max_y)) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "top",
    legend.spacing.x = grid::unit(0.5, 'cm'),
    panel.border = element_rect(linewidth = 1),
    panel.grid.major = element_line(linewidth = 0.4, linetype = 'dotted', color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 1, color = "black")
  )

# Print + Save
#print(GC_graph_elbow_common)


ggsave(filename = "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_1_GC_graph_elbow.png",
        plot = GC_graph_elbow_common, width = 12, height = 8, dpi = 600)

# --------------------------------------------------------
# Export subnetworks at common cutoff
# --------------------------------------------------------
cat("\nExporting subnetworks at common cutoff...\n")
 Subnetwork_A7_common <- A7_Full %>% slice(1:common_cutoff_links)
 Subnetwork_A9_common <- A9_Full %>% slice(1:common_cutoff_links)

graph_A7_common <- graph_from_data_frame(d = Subnetwork_A7_common, directed = FALSE)
graph_A9_common <- graph_from_data_frame(d = Subnetwork_A9_common, directed = FALSE)

saveRDS(graph_A7_common, file = "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_2_Subnetwork_A7_elbow.rds")
saveRDS(graph_A9_common, file = "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_2_Subnetwork_A9_elbow.rds")

save.image("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_3_GiantComponent.RData")
cat("Graphs exported successfully!\n")


#load("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_GiantComponent_FAST.RData")



