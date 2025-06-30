# SCRIPT: DeltaCon Similarity Analysis between Two Graphs (Parallelized)
# Author: Jose Joel Ruiz Hernandez Contact: joelruizhernandez406@gmail.com
# Date: [2025-04-24]
# 
# DESCRIPTION:
# This script computes the DeltaCon similarity between two co-expression 
# networks using parallel matrix operations for scalability.
#
# INPUTS:
#   - Edgelist per graph (both most have same nodes)
#
# OUTPUTS:
#   - DeltaCon Distance and Similarity values printed to console
#   - RData saved environment
#
# REFERENCE:
# Koutra D, Vogelstein JT, Faloutsos C. DeltaCon: A Principled Massive-Graph Similarity Function.
# In: Proceedings of the 2013 SIAM International Conference on Data Mining. 2013.
################################################################################
start_time <- Sys.time()
# Load required libraries
library(tidyverse)   
library(igraph)      
library(Matrix)     
library(furrr)       # for parallel processing with purrr
library(future)      # for configuring parallel plans

#-------------------------- PARALLELIZATION PLAN -------------------------------
# Configure parallel plan: (adjust as needed).
plan(multisession, workers = 4)
#---------------------------- DATA LOADING -------------------------------------
# Load top 600,000 edges from A7 and A9 networks. 
# These files should be sorted by mutual information (MI).
A7 <- vroom::vroom("~/4_MI/ARACNE-multicore/launch/Final_A7_Annot.sort",
                   col_names = c("GenA", "GenB", "MI"),
                   show_col_types = FALSE) %>% 
  slice(1:600000)

A9 <- vroom::vroom("~/4_MI/ARACNE-multicore/launch/Final_A9_Annot.sort",
                   col_names = c("GenA", "GenB", "MI"),
                   show_col_types = FALSE) %>%
  slice(1:600000)

#--------------------------- GRAPH CREATION ------------------------------------
# Create igraph objects from edgelists (ignoring weights).
g_A7 <- igraph::graph_from_data_frame(A7[, 1:2], directed = FALSE)
g_A9 <- igraph::graph_from_data_frame(A9[, 1:2], directed = FALSE)

# Ensure both graphs have the same node set
nodes_order <- sort(union(V(g_A7)$name, V(g_A9)$name))
g_A7 <- igraph::permute(g_A7, match(nodes_order, V(g_A7)$name))
g_A9 <- igraph::permute(g_A9, match(nodes_order, V(g_A9)$name))

#-------------------- DELTACON MATRIX COMPUTATION ------------------------------
# DeltaCon requires computing the **affinity matrix** S for each graph:
# Let A be the adjacency matrix and D the degree matrix.
# The affinity matrix S is derived from:
#
#     S = (I + ε² D - ε A)⁻¹
#
# where:
#     - I is the identity matrix
#     - ε (epsilon) is a small regularization constant (e.g., 0.001, or 1 / 1 + kmax)
#
# This function solves the linear system in **blocks** (for memory efficiency),
# using **parallel computation** for speed.
#epsilon <- 1 / (1 + max(degree(graph)))
#block_size <- ceiling(ncol(A) / workers)  # más dinámico

compute_similarity_matrix_parallel <- function(graph, epsilon = 0.001) {
  # Build adjacency matrix A (sparse)
  A <- as_adjacency_matrix(graph, sparse = TRUE)
  
  # Build degree matrix D (diagonal)
  D <- Diagonal(x = Matrix::rowSums(A))
  
  # Identity matrix
  I <- Diagonal(n = nrow(A))
  
  # Compute the matrix to invert: M = I + ε² D - ε A
  M <- I + (epsilon^2) * D - epsilon * A
  
  # Number of columns (nodes)
  n <- ncol(M)

  workers <- future::nbrOfWorkers()
  block_size <- ceiling(n / workers)
  
  # Define column blocks to invert M in pieces (for memory efficiency)
  blocks <- split(seq_len(n), ceiling(seq_len(n) / block_size))
  
  # Solve M⁻¹ I[:, cols] for each block of columns in parallel
  solve_block <- function(cols) {
    solve(M, I[, cols, drop = FALSE])
  }   #Use solve(), which internally chooses the best method (LU or Cholesky) 
      #based on the properties of MM.
  
  # Parallel solve using future_map (from furrr)
  result_blocks <- future_map(blocks, solve_block, .options = furrr_options(seed = TRUE))
  
  # Combine block solutions into full affinity matrix S
  S <- do.call(cbind, result_blocks)
  #It measures how strongly connected two nodes are in a network,
  #considering direct and indirect paths.
  rownames(S) <- colnames(S) <- V(graph)$name
  return(S)
}

#-------------------- COMPUTE AFFINITY MATRICES ------------------------------
# Compute S_A7 and S_A9 for the A7 and A9 graphs
S_A7 <- compute_similarity_matrix_parallel(g_A7)
S_A9 <- compute_similarity_matrix_parallel(g_A9)

# Get node names and common order
genes <- intersect(rownames(S_A7), rownames(S_A9))
genes <- sort(genes)
# Reorder matrices
S_A7 <- S_A7[genes, genes]
S_A9 <- S_A9[genes, genes]

colnames(S_A7) <- rownames(S_A7)
colnames(S_A9) <- rownames(S_A9)
#------------------------ DELTACON DISTANCE ------------------------------------
# Compute DeltaCon distance based on the two similarity matrices:
#
#     δ(S1, S2) = sqrt(Σᵢⱼ (√S1ᵢⱼ - √S2ᵢⱼ)²)
#
# Which is a squared Euclidean distance between the **square-rooted** affinities.
stopifnot(all(rownames(S_A7) == rownames(S_A9)))
stopifnot(identical(rownames(S_A7), rownames(S_A9)))  
stopifnot(all(colnames(S_A7) == colnames(S_A9)))

delta_con_distance <- function(S1, S2) {
  sqrt(sum((sqrt(S1) - sqrt(S2))^2))
}

# Compute distance and transform into similarity
delta_con <- delta_con_distance(S_A7, S_A9)
delta_con_similarity <- 1 / (1 + delta_con)

#------------------------- OUTPUT ----------------------------------------------
# Print distance and similarity
cat("DeltaCon Distance:", delta_con, "\n")
cat("DeltaCon Similarity:", delta_con_similarity, "\n")

#DeltaCon Distance: 41.96834 
#DeltaCon Similarity: 0.02327295 (near 1: same networks, near 0: very different networks)
#Total execution time: 1.036897 

#-------------------------- SAVE SESSION ---------------------------------------
save.image("~/5_Global_Analisis/3_DeltaCon/DeltaCon.RData")
end_time <- Sys.time()
cat("Total execution time:", end_time - start_time, "\n")

#load("~/5_Global_Analisis/3_DeltaCon/DeltaCon_parallel.RData")




