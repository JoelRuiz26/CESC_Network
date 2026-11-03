### ============================================================
### 01_ase_procrustes_distance.R
###
### Represent the full HPV-A7 and HPV-A9 co-expression networks (CESC) as
### points in a shared latent space, and measure per-gene + global distance
### between clades. Method: Adjacency Spectral Embedding (ASE) estimated
### separately per network, then aligned with an orthogonal Procrustes
### rotation -- the two-sample Random Dot Product Graph (RDPG) comparison
### of:
###
###   Sussman DL, Tang M, Fishkind DE, Priebe CE (2012). "A Consistent
###   Adjacency Spectral Embedding for Stochastic Blockmodel Graphs."
###   Journal of the American Statistical Association, 107(499):1119-1128.
###   doi:10.1080/01621459.2012.699795
###
###   Tang M, Athreya A, Sussman DL, Lyzinski V, Park Y, Priebe CE (2017).
###   "A Semiparametric Two-Sample Hypothesis Testing Problem for Random
###   Graphs." Journal of Computational and Graphical Statistics,
###   26(2):344-354. doi:10.1080/10618600.2016.1193505
###
### Restricted to the GENES SHARED by both networks (98.1% of nodes), with
### identical row/column order in both adjacency matrices, so A7 and A9
### embeddings land in a space where "row i" means the same gene in both --
### a precondition for Procrustes (and for the per-gene distance to mean
### anything). Networks are unweighted (no edge attribute survives the
### A7/A9 matching step in 05_topology), matching the plain RDPG/Bernoulli-
### edge model these methods assume.
###
### Kept deliberately to ONE free parameter: the embedding dimension d,
### chosen automatically via the Zhu & Ghodsi (2006) profile-likelihood
### scree-plot elbow (Computational Statistics & Data Analysis 51:918-930)
### -- no manually-picked threshold.
###
### Matrices stay sparse throughout (igraph/Matrix dgCMatrix) and only a
### truncated SVD (irlba) is ever computed -- the dense n x n matrix is
### never built, so this comfortably fits in RAM for n ~ 15,000 nodes.
###
### -> results/01_shared_gene_order.tsv        (canonical gene order used)
### -> results/01_ase_X7.rds, 01_ase_X9_aligned.rds  (n x d embeddings)
### -> results/01_gene_latent_distance.tsv     (per-gene distance table)
### -> results/01_global_distance_summary.tsv  (mean/median/Frobenius)
### ============================================================
suppressPackageStartupMessages({ library(igraph); library(Matrix); library(irlba); library(dplyr); library(readr) })

set.seed(123)
out_dir  <- "~/CESC_Network/11_latent_embedding/results"
topo_dir <- "~/CESC_Network/05_topology/results"
K_SCREE  <- 50   # how many singular values to inspect for the elbow

## -- 1) Load full matched networks, restrict to shared genes, same order --
g7 <- readRDS(file.path(topo_dir, "5_2_4_Subnetwork_A7_max.rds"))
g9 <- readRDS(file.path(topo_dir, "5_2_4_Subnetwork_A9_max.rds"))
shared <- sort(intersect(V(g7)$name, V(g9)$name))
cat(sprintf("Shared genes: %d (A7 total %d, A9 total %d)\n", length(shared), vcount(g7), vcount(g9)))
write_tsv(tibble(gene = shared), file.path(out_dir, "01_shared_gene_order.tsv"))

adj_in_order <- function(g, genes) {
  gs <- induced_subgraph(g, V(g)[name %in% genes])
  A  <- as_adjacency_matrix(gs, sparse = TRUE)
  A[genes, genes]
}
A7 <- adj_in_order(g7, shared)
A9 <- adj_in_order(g9, shared)

## -- 2) Pick ONE embedding dimension d for both, via the shared/average
## structure -- Zhu & Ghodsi (2006) profile-likelihood elbow on the scree
## of singular values of (A7+A9)/2 --
get_elbow <- function(sv) {
  d <- sort(sv, decreasing = TRUE); q <- length(d)
  ll <- function(i) {
    mu1 <- mean(d[1:i]); mu2 <- mean(d[(i + 1):q])
    s2  <- (sum((d[1:i] - mu1)^2) + sum((d[(i + 1):q] - mu2)^2)) / (q - 2)
    sum(dnorm(d[1:i], mu1, sqrt(s2), log = TRUE)) + sum(dnorm(d[(i + 1):q], mu2, sqrt(s2), log = TRUE))
  }
  which.max(vapply(seq_len(q - 1), ll, numeric(1)))
}
sv_bar <- irlba(( A7 + A9 ) / 2, nv = K_SCREE, nu = K_SCREE)$d
d <- get_elbow(sv_bar)
cat(sprintf("Embedding dimension (profile-likelihood elbow of the %d-singular-value scree): d = %d\n", K_SCREE, d))
write_tsv(tibble(k = seq_along(sv_bar), singular_value = sv_bar, elbow_d = d), file.path(out_dir, "01_scree_singular_values.tsv"))

## -- 3) ASE per network: X = U_d %*% diag(sqrt(S_d)) --
ase <- function(A, d) { s <- irlba(A, nv = d, nu = d); s$u %*% diag(sqrt(s$d[1:d])) }
X7 <- ase(A7, d)
X9 <- ase(A9, d)

## -- 4) Orthogonal Procrustes: rotate X9 onto X7's frame --
## argmin_R || X7 - X9 R ||_F  s.t. R orthogonal  ->  R = UV' from svd(X9'X7)
sv <- svd(crossprod(X9, X7))
R  <- sv$u %*% t(sv$v)
X9_aligned <- X9 %*% R

saveRDS(X7,         file.path(out_dir, "01_ase_X7.rds"))
saveRDS(X9,         file.path(out_dir, "01_ase_X9_raw.rds"))        # pre-alignment, needed by 02's permutation null
saveRDS(X9_aligned, file.path(out_dir, "01_ase_X9_aligned.rds"))

## -- 5) Distances --
gene_dist <- tibble(gene = shared, distance = sqrt(rowSums((X7 - X9_aligned)^2))) |> arrange(desc(distance))
write_tsv(gene_dist, file.path(out_dir, "01_gene_latent_distance.tsv"))

global_dist <- tibble(
  d = d, n_genes = length(shared),
  mean_distance = mean(gene_dist$distance), median_distance = median(gene_dist$distance),
  frobenius_distance = sqrt(sum((X7 - X9_aligned)^2)),
  procrustes_fit_ss = sum(sv$d)   # sum of singular values = alignment goodness-of-fit (higher = better)
)
write_tsv(global_dist, file.path(out_dir, "01_global_distance_summary.tsv"))

cat("\n== Global A7 vs A9 latent distance (CESC, shared genes) ==\n")
print(as.data.frame(global_dist))
cat("\nTop 10 genes with the largest shift in latent position between A7 and A9:\n")
print(as.data.frame(head(gene_dist, 10)))
message("\n[OK] Saved embeddings + distance tables to ", out_dir)
