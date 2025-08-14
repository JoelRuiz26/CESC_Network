# =========================================================
# Giant Component + Elbow (rápido y correcto: igraph coarse + DSU refine con fallback)
# Autor: Joel Ruiz Hernandez
# =========================================================

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(ggplot2)
  library(grid)
})

# --------------------------------------------------------
# Parámetros
# --------------------------------------------------------
max_cutoff_links <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_1_Max_cutoff_links.rds")
cat("\n[INFO] Requested max_cutoff_links:", max_cutoff_links, "\n")

A7_Full <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_A7_edgelist_allnodes.rds")
A9_Full <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_A9_edgelist_allnodes.rds")

imax_A7 <- min(max_cutoff_links, nrow(A7_Full))
imax_A9 <- min(max_cutoff_links, nrow(A9_Full))
cat("[INFO] Usable edges -> A7:", imax_A7, "| A9:", imax_A9, "\n")

# ========================================================
# Utilidades
# ========================================================
.as_mat_edges <- function(df, upto) {
  as.matrix(df[seq_len(upto), c("GenA","GenB")])
}

# Curva coarse correcta con igraph (checkpoints espaciados)
calculate_metrics_coarse_igraph <- function(edgelist, imax, step = 5000L) {
  imax <- min(imax, nrow(edgelist))
  take <- unique(sort(c(seq(from = step, to = imax, by = step), imax)))
  Links <- integer(length(take))
  Total_Genes <- integer(length(take))
  Largest_Component <- integer(length(take))
  
  for (j in seq_along(take)) {
    i <- take[j]
    E <- .as_mat_edges(edgelist, i)
    g <- graph_from_edgelist(E, directed = FALSE)
    cs <- components(g)$csize
    Links[j] <- i
    Total_Genes[j] <- vcount(g)
    Largest_Component[j] <- max(cs)
  }
  data.frame(Links, Total_Genes, Largest_Component)
}

# DSU refine warm-start SOLO en la ventana [start, end]
refine_window_dsu <- function(edgelist, start, end) {
  stopifnot(start >= 1L, end >= start)
  a <- trimws(as.character(edgelist$GenA[1:end]))
  b <- trimws(as.character(edgelist$GenB[1:end]))
  nodes <- unique(c(a, b))
  idA <- match(a, nodes); idB <- match(b, nodes)
  n_nodes <- length(nodes)
  
  parent <- seq_len(n_nodes)
  comp_sz <- rep.int(1L, n_nodes)
  seen <- rep.int(FALSE, n_nodes)
  seen_count <- 0L
  largest <- 0L
  
  dsu_find <- function(x) {
    while (parent[x] != x) {
      parent[x] <- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  dsu_union <- function(x, y) {
    rx <- dsu_find(x); ry <- dsu_find(y)
    if (rx == ry) return(comp_sz[rx])
    # union simple (sin heurística; suficiente y robusto)
    parent[ry] <- rx
    comp_sz[rx] <- comp_sz[rx] + comp_sz[ry]
    comp_sz[ry] <- 0L
    comp_sz[rx]
  }
  
  # Warm-start hasta start-1
  if (start > 1L) {
    for (i in 1:(start-1L)) {
      u <- idA[i]; v <- idB[i]
      if (!seen[u]) { seen[u] <- TRUE; seen_count <- seen_count + 1L }
      if (!seen[v]) { seen[v] <- TRUE; seen_count <- seen_count + 1L }
      if (u != v) {
        s <- dsu_union(u, v)
        if (s > largest) largest <- s
      }
    }
  }
  
  # Pasada fina SOLO en [start, end]
  n <- end - start + 1L
  Links <- integer(n)
  Total_Genes <- integer(n)
  Largest_Component <- integer(n)
  
  ptr <- 0L
  for (i in start:end) {
    u <- idA[i]; v <- idB[i]
    if (!seen[u]) { seen[u] <- TRUE; seen_count <- seen_count + 1L }
    if (!seen[v]) { seen[v] <- TRUE; seen_count <- seen_count + 1L }
    if (u != v) {
      s <- dsu_union(u, v)
      if (s > largest) largest <- s
    }
    ptr <- ptr + 1L
    Links[ptr] <- i
    Total_Genes[ptr] <- seen_count
    Largest_Component[ptr] <- largest
  }
  
  data.frame(Links, Total_Genes, Largest_Component)
}

# Fallback: refine con igraph incremental SOLO en la ventana (correcto; algo más pesado)
refine_window_igraph <- function(edgelist, start, end) {
  E <- .as_mat_edges(edgelist, end)
  # grafo hasta start-1
  if (start > 1L) {
    g <- graph_from_edgelist(E[1:(start-1L), , drop = FALSE], directed = FALSE)
  } else {
    g <- make_empty_graph(directed = FALSE)
  }
  n <- end - start + 1L
  Links <- integer(n); Total_Genes <- integer(n); Largest_Component <- integer(n)
  
  ptr <- 0L
  for (i in start:end) {
    # asegurar vértices
    e <- E[i, ]
    for (v in e) {
      if (!(v %in% V(g)$name)) g <- add_vertices(g, 1, name = v)
    }
    g <- add_edges(g, as.vector(e))
    cs <- components(g)$csize
    ptr <- ptr + 1L
    Links[ptr] <- i
    Total_Genes[ptr] <- vcount(g)
    Largest_Component[ptr] <- max(cs)
  }
  data.frame(Links, Total_Genes, Largest_Component)
}

# Detección del elbow (distancia ortogonal a la diagonal normalizada)
find_elbow <- function(x, y) {
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  line_vec <- c(1, 1) / sqrt(2)
  P <- cbind(x_norm, y_norm)
  proj <- (P %*% line_vec) %*% t(line_vec)
  dists <- sqrt(rowSums((P - proj)^2))
  which.max(dists)
}

# Geometría para dibujar cuerda y perpendicular
elbow_geom <- function(x, y, idx) {
  A <- c(x[1], y[1]); B <- c(x[length(x)], y[length(y)])
  AB <- B - A; AB2 <- sum(AB^2)
  P  <- c(x[idx], y[idx])
  t  <- ((P[1]-A[1])*AB[1] + (P[2]-A[2])*AB[2]) / AB2
  Pp <- c(A[1] + t*AB[1], A[2] + t*AB[2])
  list(x_start=A[1], y_start=A[2], x_end=B[1], y_end=B[2],
       x_elbow=P[1], y_elbow=P[2], x_proj=Pp[1], y_proj=Pp[2])
}

# ========================================================
# 1) Curvas coarse (correctas y rápidas)
# ========================================================
cat("\n[INFO] Computing coarse curves (step = 5000)...\n")
metrics_A7_coarse <- calculate_metrics_coarse_igraph(A7_Full, imax = imax_A7, step = 5000L)
metrics_A9_coarse <- calculate_metrics_coarse_igraph(A9_Full, imax = imax_A9, step = 5000L)

elbow_idx_A7_c <- find_elbow(metrics_A7_coarse$Links, metrics_A7_coarse$Largest_Component)
elbow_idx_A9_c <- find_elbow(metrics_A9_coarse$Links, metrics_A9_coarse$Largest_Component)
elbow_links_A7_c <- metrics_A7_coarse$Links[elbow_idx_A7_c]
elbow_links_A9_c <- metrics_A9_coarse$Links[elbow_idx_A9_c]
cat("[INFO] Coarse elbows ~", elbow_links_A7_c, "(A7), ~", elbow_links_A9_c, "(A9)\n")

# ========================================================
# 2) Refinamiento (DSU warm-start) + verificación; si falla -> igraph
# ========================================================
refine_window <- 2500L
ref_range_A7 <- c(max(1L, elbow_links_A7_c - refine_window), min(imax_A7, elbow_links_A7_c + refine_window))
ref_range_A9 <- c(max(1L, elbow_links_A9_c - refine_window), min(imax_A9, elbow_links_A9_c + refine_window))

cat("[INFO] Refining A7 in [", ref_range_A7[1], ", ", ref_range_A7[2], "] ...\n", sep = "")
metrics_A7_refine <- refine_window_dsu(A7_Full, start = ref_range_A7[1], end = ref_range_A7[2])

# Chequeo de consistencia en el punto medio de la ventana A7
mid_A7 <- floor(mean(ref_range_A7))
E_mid <- .as_mat_edges(A7_Full, mid_A7)
g_mid <- graph_from_edgelist(E_mid, directed = FALSE)
lc_mid_igraph <- max(components(g_mid)$csize)
lc_mid_dsu <- metrics_A7_refine$Largest_Component[which(metrics_A7_refine$Links == mid_A7)]

if (length(lc_mid_dsu)==0L || abs(lc_mid_dsu - lc_mid_igraph) > 1L) {
  cat("[WARN] DSU mismatch on A7 (midpoint). Falling back to igraph refine.\n")
  metrics_A7_refine <- refine_window_igraph(A7_Full, start = ref_range_A7[1], end = ref_range_A7[2])
}

cat("[INFO] Refining A9 in [", ref_range_A9[1], ", ", ref_range_A9[2], "] ...\n", sep = "")
metrics_A9_refine <- refine_window_dsu(A9_Full, start = ref_range_A9[1], end = ref_range_A9[2])

mid_A9 <- floor(mean(ref_range_A9))
E_mid9 <- .as_mat_edges(A9_Full, mid_A9)
g_mid9 <- graph_from_edgelist(E_mid9, directed = FALSE)
lc_mid9_igraph <- max(components(g_mid9)$csize)
lc_mid9_dsu <- metrics_A9_refine$Largest_Component[which(metrics_A9_refine$Links == mid_A9)]

if (length(lc_mid9_dsu)==0L || abs(lc_mid9_dsu - lc_mid9_igraph) > 1L) {
  cat("[WARN] DSU mismatch on A9 (midpoint). Falling back to igraph refine.\n")
  metrics_A9_refine <- refine_window_igraph(A9_Full, start = ref_range_A9[1], end = ref_range_A9[2])
}

# Fusionar: coarse fuera de la ventana + refine dentro
metrics_A7 <- bind_rows(
  filter(metrics_A7_coarse, Links < ref_range_A7[1] | Links > ref_range_A7[2]),
  metrics_A7_refine
) %>% arrange(Links)

metrics_A9 <- bind_rows(
  filter(metrics_A9_coarse, Links < ref_range_A9[1] | Links > ref_range_A9[2]),
  metrics_A9_refine
) %>% arrange(Links)

# ========================================================
# 3) Elbows finales y corte común
# ========================================================
elbow_idx_A7 <- find_elbow(metrics_A7$Links, metrics_A7$Largest_Component)
elbow_links_A7 <- metrics_A7$Links[elbow_idx_A7]
elbow_genes_A7 <- metrics_A7$Largest_Component[elbow_idx_A7]

elbow_idx_A9 <- find_elbow(metrics_A9$Links, metrics_A9$Largest_Component)
elbow_links_A9 <- metrics_A9$Links[elbow_idx_A9]
elbow_genes_A9 <- metrics_A9$Largest_Component[elbow_idx_A9]

cat("\nDetected thresholds (elbows):\n")
cat("  HPV-A7: ", elbow_links_A7, " edges; ", elbow_genes_A7, " genes\n", sep = "")
cat("  HPV-A9: ", elbow_links_A9, " edges; ", elbow_genes_A9, " genes\n", sep = "")

common_cutoff_links <- max(elbow_links_A7, elbow_links_A9)
cat("\nCommon threshold (maximum elbow): ", common_cutoff_links, " edges\n", sep = "")

# ========================================================
# 4) Geometría y figura con todos los señalamientos
# ========================================================
geom_A7 <- elbow_geom(metrics_A7$Links, metrics_A7$Largest_Component, elbow_idx_A7)
geom_A9 <- elbow_geom(metrics_A9$Links, metrics_A9$Largest_Component, elbow_idx_A9)

max_x <- max(metrics_A7$Links, metrics_A9$Links)
max_y <- max(metrics_A7$Largest_Component, metrics_A9$Largest_Component)

GC_graph_elbow_common <- ggplot() +
  geom_line(data = metrics_A7, aes(x = Links, y = Largest_Component, color = "HPV-A7"), linewidth = 1.4) +
  geom_line(data = metrics_A9, aes(x = Links, y = Largest_Component, color = "HPV-A9"), linewidth = 1.4) +
  annotate("segment", x = geom_A7$x_start, y = geom_A7$y_start, xend = geom_A7$x_end, yend = geom_A7$y_end,
           linetype = "dotdash", linewidth = 0.6, alpha = 0.7, color = "blue") +
  annotate("segment", x = geom_A9$x_start, y = geom_A9$y_start, xend = geom_A9$x_end, yend = geom_A9$y_end,
           linetype = "dotdash", linewidth = 0.6, alpha = 0.7, color = "red") +
  annotate("segment", x = geom_A7$x_elbow, y = geom_A7$y_elbow, xend = geom_A7$x_proj, yend = geom_A7$y_proj,
           linetype = "dotted", linewidth = 0.7, color = "grey30") +
  annotate("segment", x = geom_A9$x_elbow, y = geom_A9$y_elbow, xend = geom_A9$x_proj, yend = geom_A9$y_proj,
           linetype = "dotted", linewidth = 0.7, color = "grey30") +
  annotate("point", x = geom_A7$x_elbow, y = geom_A7$y_elbow, size = 2.6, shape = 21, stroke = 1, fill = "white", color = "blue") +
  annotate("point", x = geom_A9$x_elbow, y = geom_A9$y_elbow, size = 2.6, shape = 21, stroke = 1, fill = "white", color = "red") +
  geom_vline(xintercept = elbow_links_A7, color = "blue", linetype = "longdash", linewidth = 0.3) +
  geom_vline(xintercept = elbow_links_A9, color = "red",  linetype = "longdash", linewidth = 0.3) +
  geom_vline(xintercept = common_cutoff_links, linetype = "dashed", color = "black", linewidth = 0.7) +
  annotate("text", x = geom_A7$x_elbow, y = geom_A7$y_elbow, label = paste0("A7 elbow\nx=", format(elbow_links_A7, big.mark=",")),
           hjust = -0.05, vjust = -0.8, size = 3.6, color = "blue") +
  annotate("text", x = geom_A9$x_elbow, y = geom_A9$y_elbow, label = paste0("A9 elbow\nx=", format(elbow_links_A9, big.mark=",")),
           hjust = -0.05, vjust = -0.8, size = 3.6, color = "red") +
  annotate("text", x = common_cutoff_links, y = max_y*0.06, label = paste0("Common cutoff = ", format(common_cutoff_links, big.mark=",")),
           angle = 90, vjust = -0.5, size = 4.0, color = "black") +
  labs(title = "Growth of Largest Connected Component in Cervical Cancer Co-expression Networks",
       x = "Number of Edges", y = "Number of Genes in Largest Connected Component", color = "Network") +
  scale_color_manual(values = c("HPV-A7" = "blue", "HPV-A9" = "red")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max_x)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max_y)) +
  theme_bw(base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = "top",
        legend.spacing.x = unit(0.5, 'cm'),
        panel.border = element_rect(linewidth = 1),
        panel.grid.major = element_line(linewidth = 0.4, linetype = 'dotted', color = "grey80"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 1, color = "black"))

print(GC_graph_elbow_common)







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


