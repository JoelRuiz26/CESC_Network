# SCRIPT: DeltaCon Similarity Analysis between A7 and A9 Graphs (Robust + Parallel)
# Author: Jose Joel Ruiz Hernandez | Date: 2025-06-05
# REFERENCE:
# Koutra D, Vogelstein JT, Faloutsos C. DeltaCon: A Principled Massive-Graph Similarity Function.
# In: Proceedings of the 2013 SIAM International Conference on Data Mining. 2013.

start_time <- Sys.time()

# -------------------------- LIBRARIES ----------------------------------------
library(igraph)
library(Matrix)       
library(furrr)        # for parallel processing with purrr
library(future)       # for configuring parallel plans

# ---------------------- PARALLEL PLAN ----------------------------------------
# Configure parallel plan: (adjust as needed).
plan(multisession, workers = 2)

# --------------------------- LOAD GRAPHS -------------------------------------
g_A7 <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds")
g_A9 <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds")

# Ensure both graphs have the same node set
nodes_order <- sort(union(V(g_A7)$name, V(g_A9)$name))
g_A7 <- igraph::permute(g_A7, match(nodes_order, V(g_A7)$name))
g_A9 <- igraph::permute(g_A9, match(nodes_order, V(g_A9)$name))

# -------------------- DELTACON MATRIX COMPUTATION ----------------------------
# DeltaCon requires computing the **affinity matrix** S for each graph:
# Let A be the adjacency matrix and D the degree matrix.
# The affinity matrix S is derived from:
#
#     S = (I + ε² D - ε A)⁻¹
#
# where:
#     - I is the identity matrix
#     - ε (epsilon) is a small regularization constant (1 / 1 + kmax)

compute_similarity_matrix_parallel <- function(graph) {
  message("🧠 Calculando matriz de afinidad con ΔCon...")
  
  epsilon <- 1 / (1 + max(igraph::degree(graph)))
  
  # Build adjacency matrix A (sparse)
  A <- igraph::as_adjacency_matrix(graph, sparse = TRUE)
  
  # Build degree matrix D (diagonal)
  D <- Diagonal(x = Matrix::rowSums(A))
  
  # Identity matrix
  I <- Diagonal(n = nrow(A))
  
  # Compute the matrix to invert: M = I + ε² D - ε A
  M <- I + (epsilon^2) * D - epsilon * A
  
  # Solve the full system: S = M⁻¹
  # solve() internally chooses the best method (LU or Cholesky)
  # based on the properties of M.
  # It measures how strongly connected two nodes are in a network,
  # considering direct and indirect paths.
  S <- solve(M, I)
  
  rownames(S) <- colnames(S) <- V(graph)$name
  return(S)
}

# -------------------- COMPUTE AFFINITY MATRICES ------------------------------
# Compute S_A7 and S_A9 for the A7 and A9 graphs
S_A7 <- compute_similarity_matrix_parallel(g_A7)
S_A9 <- compute_similarity_matrix_parallel(g_A9)

# -------------------- ALIGN MATRICES TO SAME ORDER ---------------------------
# Get node names and common order
genes <- sort(intersect(rownames(S_A7), rownames(S_A9)))
S_A7 <- S_A7[genes, genes]
S_A9 <- S_A9[genes, genes]

# -------------------- DELTACON DISTANCE --------------------------------------
# Compute DeltaCon distance based on the two similarity matrices:
#
#     δ(S1, S2) = sqrt(Σᵢⱼ (√S1ᵢⱼ - √S2ᵢⱼ)²)
#
# Which is a squared Euclidean distance between the **square-rooted** affinities.

stopifnot(all(rownames(S_A7) == rownames(S_A9)))
stopifnot(identical(rownames(S_A7), rownames(S_A9)))
stopifnot(all(colnames(S_A7) == colnames(S_A9)))

delta_con_distance <- function(S1, S2) {
  sqrt(sum((sqrt(S1) - sqrt(S2))^2, na.rm = TRUE))
}

# Compute distance and transform into similarity
delta_con <- delta_con_distance(S_A7, S_A9)
delta_con_similarity <- 1 / (1 + delta_con)

# ------------------------- OUTPUT --------------------------------------------
cat("✅ DeltaCon Distance (A7 vs A9): ", round(delta_con, 4), "\n")
cat("✅ DeltaCon Similarity (A7 vs A9): ", round(delta_con_similarity, 6), "\n")

# ------------------------- SAVE OUTPUT ---------------------------------------
out_dir <- "~/CESC_Network/5_Network_analisis/5_3_DeltaCon/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(S_A7, file.path(out_dir, "S_A7_matrix.rds"))
saveRDS(S_A9, file.path(out_dir, "S_A9_matrix.rds"))
saveRDS(delta_con, file.path(out_dir, "DeltaCon_distance_A7_A9.rds"))
saveRDS(delta_con_similarity,
        file.path(out_dir, "DeltaCon_similarity_A7_A9.rds"))

save.image(file.path(out_dir, "DeltaCon_A7_A9_image.RData"))
#load(file.path(out_dir, "DeltaCon_A7_A9_image.RData"))

end_time <- Sys.time()
cat("🕒 Total execution time:", round(end_time - start_time, 2), "seconds\n")

#readRDS(file.path(out_dir, "DeltaCon_distance_A7_A9.rds"))
# [1] 65.09439
# > readRDS(file.path(out_dir, "DeltaCon_similarity_A7_A9.rds"))
# [1] 0.01512988