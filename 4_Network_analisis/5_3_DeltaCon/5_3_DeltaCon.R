# SCRIPT: DeltaCon Similarity Analysis between A7 and A9 Graphs (Robust + Parallel)
# Author: Jose Joel Ruiz Hernandez | Date: 2025-06-05
#REFERENCE:
# Koutra D, Vogelstein JT, Faloutsos C. DeltaCon: A Principled Massive-Graph Similarity Function.
# In: Proceedings of the 2013 SIAM International Conference on Data Mining. 2013.

start_time <- Sys.time()

# -------------------------- LIBRARIES ----------------------------------------
library(igraph)
library(Matrix)
library(furrr)
library(future)

# ---------------------- PARALLEL PLAN (safe) ---------------------------------
plan(multisession, workers = 2)

# --------------------------- LOAD GRAPHS -------------------------------------
g_A7 <- readRDS("~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds")
g_A9 <- readRDS("~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds")


# Align both graphs to same node order
nodes_order <- sort(union(V(g_A7)$name, V(g_A9)$name))

g_A7 <- igraph::permute(g_A7, match(nodes_order, V(g_A7)$name))
g_A9 <- igraph::permute(g_A9, match(nodes_order, V(g_A9)$name))

# ----------------- FUNCTION: Compute Affinity Matrix (DeltaCon) --------------
compute_similarity_matrix_parallel <- function(graph) {
  message("🧠 Calculando matriz de afinidad con ΔCon...")
  
  epsilon <- 1 / (1 + max(igraph::degree(graph)))
  
  A <- igraph::as_adjacency_matrix(graph, sparse = TRUE)
  D <- Diagonal(x = Matrix::rowSums(A))
  I <- Diagonal(n = nrow(A))
  M <- I + (epsilon^2) * D - epsilon * A
  
  n <- ncol(M)
  workers <- future::nbrOfWorkers()
  block_size <- ceiling(n / workers)
  blocks <- split(seq_len(n), ceiling(seq_len(n) / block_size))
  
  solve_block <- function(cols) {
    tryCatch({
      solve(M, I[, cols, drop = FALSE])
    }, error = function(e) {
      message("❌ Error resolviendo columnas: ", paste(head(cols, 5), "..."))
      message("→ ", e$message)
      return(matrix(NA, nrow = nrow(M), ncol = length(cols)))
    })
  }
  
  result_blocks <- future_map(blocks, solve_block,
                              .options = furrr_options(seed = TRUE))
  
  S <- do.call(cbind, result_blocks)
  rownames(S) <- colnames(S) <- V(graph)$name
  return(S)
}

# -------------------- COMPUTE AFFINITY MATRICES ------------------------------
S_A7 <- compute_similarity_matrix_parallel(g_A7)
S_A9 <- compute_similarity_matrix_parallel(g_A9)

# -------------------- ALIGN MATRICES TO SAME ORDER ---------------------------
genes <- sort(intersect(rownames(S_A7), rownames(S_A9)))
S_A7 <- S_A7[genes, genes]
S_A9 <- S_A9[genes, genes]

# -------------------- DELTACON DISTANCE --------------------------------------
delta_con_distance <- function(S1, S2) {
  sqrt(sum((sqrt(S1) - sqrt(S2))^2, na.rm = TRUE))
}

delta_con <- delta_con_distance(S_A7, S_A9)
delta_con_similarity <- 1 / (1 + delta_con)

# ------------------------- OUTPUT --------------------------------------------
cat("✅ DeltaCon Distance (A7 vs A9): ", round(delta_con, 4), "\n")
cat("✅ DeltaCon Similarity (A7 vs A9): ", round(delta_con_similarity, 6), "\n")

# ------------------------- SAVE OUTPUT ---------------------------------------
out_dir <- "~/CESC_Network/4_Network_analisis/5_3_DeltaCon/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(S_A7, file.path(out_dir, "S_A7_matrix.rds"))
saveRDS(S_A9, file.path(out_dir, "S_A9_matrix.rds"))
saveRDS(delta_con, file.path(out_dir, "DeltaCon_distance_A7_A9.rds"))
saveRDS(delta_con_similarity,
        file.path(out_dir, "DeltaCon_similarity_A7_A9.rds"))

save.image(file.path(out_dir, "DeltaCon_A7_A9_image.RData"))

end_time <- Sys.time()
cat("🕒 Total execution time:", round(end_time - start_time, 2), "seconds\n")
