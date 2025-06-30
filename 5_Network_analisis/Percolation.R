library(igraph)
library(vroom)

# ---------------------- CONFIG ----------------------- #
path_A7 <- "/home/jruiz/CESC_Network/4_Network/ARACNE-multicore/launch/2_7_Full_counts_A7_annot.sort"
path_A9 <- "/home/jruiz/CESC_Network/4_Network/ARACNE-multicore/launch/2_8_Full_counts_A9_annot.sort"

max_edges <- 5000000
out_folder <- "/home/jruiz/CESC_Network/5_Network_analisis/percolation_plots/"
dir.create(out_folder, showWarnings = FALSE)

# ---------------------- FUNCTION ------------------------- #
analyze_percolation <- function(path, net_name) {
  cat("\nReading:", net_name, "\n")
  edges <- vroom(path, col_names = c("Gene1", "Gene2", "MI"), n_max = max_edges)
  
  all_genes <- unique(c(edges$Gene1, edges$Gene2))
  total_nodes <- length(all_genes)
  cat("Total unique genes:", total_nodes, "\n")
  
  g <- make_empty_graph(n = total_nodes, directed = FALSE)
  V(g)$name <- all_genes
  
  giant_sizes <- numeric(nrow(edges))
  
  for (i in seq_len(nrow(edges))) {
    g <- add_edges(g, c(edges$Gene1[i], edges$Gene2[i]))
    comps <- components(g)
    giant_sizes[i] <- max(comps$csize) / total_nodes
  }
  
  edge_fraction <- seq_along(giant_sizes) / nrow(edges)
  deriv <- c(0, diff(giant_sizes) / diff(edge_fraction))
  
  # Percolation point = max derivative
  percol_idx <- which.max(deriv)
  percol_MI <- edges$MI[percol_idx]
  
  cat("\n===== Percolation Results =====\n")
  cat("Network:", net_name, "\n")
  cat("Percolation edge index:", percol_idx, "\n")
  cat("Fraction of edges added:", round(edge_fraction[percol_idx], 4), "\n")
  cat("Giant component fraction:", round(giant_sizes[percol_idx], 4), "\n")
  cat("Mutual Information threshold:", percol_MI, "\n")
  
  # -------------------- PLOT -------------------- #
  pdf(file = paste0(out_folder, net_name, "_Percolation.pdf"), width=8, height=6)
  plot(edge_fraction, giant_sizes, type='l', lwd=2,
       main=paste("Percolation Analysis:", net_name),
       xlab="Fraction of Edges Added",
       ylab="Fraction of Nodes in Giant Component")
  abline(v=edge_fraction[percol_idx], col="red", lty=2)
  points(edge_fraction[percol_idx], giant_sizes[percol_idx], col="red", pch=19)
  legend("bottomright", legend=paste("MI threshold:", round(percol_MI, 3)), bty="n")
  dev.off()
  
  png(file = paste0(out_folder, net_name, "_Percolation.png"), width=800, height=600)
  plot(edge_fraction, giant_sizes, type='l', lwd=2,
       main=paste("Percolation Analysis:", net_name),
       xlab="Fraction of Edges Added",
       ylab="Fraction of Nodes in Giant Component")
  abline(v=edge_fraction[percol_idx], col="red", lty=2)
  points(edge_fraction[percol_idx], giant_sizes[percol_idx], col="red", pch=19)
  legend("bottomright", legend=paste("MI threshold:", round(percol_MI, 3)), bty="n")
  dev.off()
  
  return(list(
    network_name = net_name,
    total_nodes = total_nodes,
    percolation_index = percol_idx,
    edge_fraction = edge_fraction[percol_idx],
    giant_fraction = giant_sizes[percol_idx],
    MI_threshold = percol_MI
  ))
}

# --------------------- RUN BOTH NETWORKS ---------------------- #
result_A7 <- analyze_percolation(path_A7, "A7")
result_A9 <- analyze_percolation(path_A9, "A9")

# --------------------- SAVE RESULTS --------------------------- #
save(result_A7, result_A9, file = paste0(out_folder, "PercolationResults.RData"))
print(result_A7)
print(result_A9)
