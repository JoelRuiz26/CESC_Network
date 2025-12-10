############################################################
## SCRIPT: 5_3_1_DeltaCon_full_A7_vs_A9.R
##
## Goal:
##   Compute FULL DeltaCon affinity matrices and distance
##   between A7 and A9 subnetworks (~17k genes), evitando
##   sesgos numéricos ligados al orden alfabético de genes.
##
## Inputs:
##   - 5_2_4_Subnetwork_A7_max.rds
##   - 5_2_4_Subnetwork_A9_max.rds
##
## Outputs:
##   - 5_3_1_S_A7_full.rds
##   - 5_3_1_S_A9_full.rds
##   - 5_3_1_DeltaCon_distance_A7_A9_full.rds
##   - 5_3_1_DeltaCon_similarity_A7_A9_full.rds
############################################################

rm(list = ls())
gc()

library(igraph)
library(Matrix)

start_time <- Sys.time()

#--------------------------- LOAD GRAPHS ---------------------------------------

g_A7 <- readRDS(
  "~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds"
)
g_A9 <- readRDS(
  "~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds"
)

#---------------------- ALINEAR MISMO UNIVERSO DE GENES -----------------------

genes_A7 <- V(g_A7)$name
genes_A9 <- V(g_A9)$name

cat("Genes A7:", length(genes_A7), "\n")
cat("Genes A9:", length(genes_A9), "\n")

# Usamos SOLO la intersección (genes presentes en ambos)
common_genes <- intersect(genes_A7, genes_A9)
cat("Genes comunes (intersección):", length(common_genes), "\n")

# IMPORTANTE:
#   intersect() devuelve common_genes ordenados, típicamente alfabéticamente.
#   Para evitar sesgos numéricos ligados al orden lexicográfico, imponemos
#   un ORDEN ALEATORIO FIJO (reproducible con seed).
set.seed(12345)
perm <- sample(length(common_genes))
common_genes <- common_genes[perm]

# Índices en cada grafo que corresponden a common_genes (en este NUEVO orden)
idx_A7 <- match(common_genes, genes_A7)
idx_A9 <- match(common_genes, genes_A9)

stopifnot(!any(is.na(idx_A7)))
stopifnot(!any(is.na(idx_A9)))

# Subgrafos inducidos EN EL MISMO ORDEN de common_genes
g_A7 <- igraph::induced_subgraph(g_A7, vids = idx_A7)
g_A9 <- igraph::induced_subgraph(g_A9, vids = idx_A9)

# Renombramos explícitamente por seguridad
V(g_A7)$name <- common_genes
V(g_A9)$name <- common_genes

stopifnot(identical(V(g_A7)$name, V(g_A9)$name))

cat("✅ A7 y A9 alineados con el mismo universo y orden de genes (aleatorio).\n")
cat("Número de genes en DeltaCon full:", length(common_genes), "\n\n")

#------------------ FUNCIÓN: MATRIZ DE AFINIDAD COMPLETA -----------------------
# S = (I + eps^2 D - eps A)^(-1)
# - Un solo factorizado (Cholesky) de M
# - Resolución por bloques para no explotar RAM

compute_affinity_full <- function(graph,
                                  epsilon = NULL,
                                  block_size = 1500) {
  cat(">> Construyendo matrices A, D, I...\n")
  
  A <- igraph::as_adjacency_matrix(graph, sparse = TRUE)
  D <- Diagonal(x = Matrix::rowSums(A))
  n <- nrow(A)
  I <- Diagonal(n)
  
  if (is.null(epsilon)) {
    epsilon <- 1 / (1 + max(igraph::degree(graph)))
  }
  cat("   epsilon elegido:", epsilon, "\n")
  
  cat(">> Construyendo M = I + eps^2 D - eps A ...\n")
  M <- I + (epsilon^2) * D - epsilon * A
  
  cat(">> Factorizando M (Cholesky)...\n")
  M_chol <- try(Cholesky(M, LDL = FALSE, perm = TRUE), silent = TRUE)
  
  if (inherits(M_chol, "try-error")) {
    cat("!! Cholesky falló, usando solve(M, RHS) directo (más lento)...\n")
    use_chol <- FALSE
  } else {
    use_chol <- TRUE
  }
  
  cols <- seq_len(n)
  blocks <- split(cols, ceiling(seq_along(cols) / block_size))
  S_blocks <- vector("list", length(blocks))
  
  cat(">> Resolviendo M^{-1} por bloques (", length(blocks), " bloques)...\n", sep = "")
  
  for (i in seq_along(blocks)) {
    idx <- blocks[[i]]
    cat("   - Bloque", i, "de", length(blocks),
        " (cols:", min(idx), "-", max(idx), ")\n")
    
    RHS <- I[, idx, drop = FALSE]  # columnas de la identidad
    
    if (use_chol) {
      S_blocks[[i]] <- Matrix::solve(M_chol, RHS)
    } else {
      S_blocks[[i]] <- solve(M, RHS)
    }
    
    gc()
  }
  
  cat(">> Combinando bloques en matriz S completa...\n")
  S <- do.call(cbind, S_blocks)
  
  rownames(S) <- V(graph)$name
  colnames(S) <- V(graph)$name
  
  S
}

#-------------------- COMPUTE AFFINITY MATRICES -------------------------------

cat("=== Calculando S_A7 (full) ===\n")
S_A7 <- compute_affinity_full(g_A7, epsilon = NULL, block_size = 1500)
gc()

cat("\n=== Calculando S_A9 (full) ===\n")
S_A9 <- compute_affinity_full(g_A9, epsilon = NULL, block_size = 1500)
gc()

#-------------------- ASEGURAR ÓRDENES (SIN REORDENAR) ------------------------

# Aquí SOLO comprobamos que el orden de genes sea idéntico.
# NO volvemos a usar sort(), para NO reintroducir orden alfabético.

stopifnot(identical(rownames(S_A7), rownames(S_A9)))
stopifnot(identical(colnames(S_A7), colnames(S_A9)))

genes <- rownames(S_A7)
cat("Primeros 10 genes en el orden interno de S_A7/S_A9:\n")
print(genes[1:10])

#------------------------ DELTACON DISTANCE (POR BLOQUES) ---------------------
# δ(S1,S2) = sqrt( Σ_ij ( sqrt(S1_ij) - sqrt(S2_ij) )^2 )
# Bloques de filas para no duplicar matrices enormes.

delta_con_distance_blocked <- function(S1, S2, block_size = 800) {
  stopifnot(all(dim(S1) == dim(S2)))
  
  n <- nrow(S1)
  blocks <- split(seq_len(n), ceiling(seq_along(seq_len(n)) / block_size))
  
  total <- 0
  
  cat(">> Calculando DeltaCon por bloques (", length(blocks), " bloques)...\n", sep = "")
  
  for (i in seq_along(blocks)) {
    idx <- blocks[[i]]
    cat("   - Bloque", i, "de", length(blocks),
        " (filas:", min(idx), "-", max(idx), ")\n")
    
    B1 <- as.matrix(S1[idx, , drop = FALSE])
    B2 <- as.matrix(S2[idx, , drop = FALSE])
    
    B1 <- pmax(B1, 0)
    B2 <- pmax(B2, 0)
    
    diff_sqrt <- sqrt(B1) - sqrt(B2)
    total <- total + sum(diff_sqrt^2)
    
    rm(B1, B2, diff_sqrt)
    gc()
  }
  
  sqrt(total)
}

cat("\n=== Calculando DeltaCon Distance & Similarity (full) ===\n")
delta_con <- delta_con_distance_blocked(S_A7, S_A9, block_size = 800)
delta_con_similarity <- 1 / (1 + delta_con)

cat("✅ DeltaCon (FULL) Distance (A7 vs A9): ", round(delta_con, 4), "\n")
cat("✅ DeltaCon (FULL) Similarity (A7 vs A9): ", round(delta_con_similarity, 6), "\n")

#------------------------- SAVE OUTPUTS ----------------------------------------

out_dir <- "~/CESC_Network/4_Network_analisis/5_3_DeltaCon/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(S_A7, file.path(out_dir, "5_3_1_S_A7_full.rds"))
saveRDS(S_A9, file.path(out_dir, "5_3_1_S_A9_full.rds"))

saveRDS(delta_con,
        file.path(out_dir, "5_3_1_DeltaCon_distance_A7_A9_full.rds"))
saveRDS(delta_con_similarity,
        file.path(out_dir, "5_3_1_DeltaCon_similarity_A7_A9_full.rds"))

end_time <- Sys.time()
cat("\nTiempo total de ejecución:", end_time - start_time, "\n")
