################################################################################
# Find the depth N at which both networks yield the same number of shared nodes
# AND the same number of edges in the induced shared-node subgraph.
#
# Strategy: precompute the "activation depth" of every edge in each network
# (the minimum N at which the edge is counted in the shared subgraph), then
# scan the difference function d(N) = a9_edges(N) - a7_edges(N) for zeros
# or sign changes. No output files are written.
################################################################################

suppressPackageStartupMessages({ library(vroom); library(dplyr) })

DIR <- "/STORAGE/csbig/jruiz/cesc_network/cervix_ARACNE"
read_net <- function(f)
  vroom(f, col_names = c("GenA","GenB","MI","pvalue"),
        skip = 1, show_col_types = FALSE) %>% select(GenA, GenB, MI)

cat("[1] Loading networks...\n")
a7 <- read_net(file.path(DIR, "A7_network_sort.txt"))
a9 <- read_net(file.path(DIR, "A9_network_sort.txt"))
cat("  A7:", nrow(a7), "edges |  A9:", nrow(a9), "edges\n\n")

# ── First-appearance rank per node ──────────────────────────────────────────
cat("[2] Computing first-appearance ranks...\n")
first_rank <- function(net) {
  genes <- c(rbind(net$GenA, net$GenB))
  ranks <- rep(seq_len(nrow(net)), each = 2)
  tapply(ranks, genes, min)
}
fa7 <- first_rank(a7)
fa9 <- first_rank(a9)

all_nodes <- union(names(fa7), names(fa9))
fa7_v <- setNames(rep(Inf, length(all_nodes)), all_nodes); fa7_v[names(fa7)] <- fa7
fa9_v <- setNames(rep(Inf, length(all_nodes)), all_nodes); fa9_v[names(fa9)] <- fa9

# ── Activation depth per edge ────────────────────────────────────────────────
# A7 edge at rank r contributes when N >= max(r, fa9[GenA], fa9[GenB])
cat("[3] Computing activation depths...\n")
safe <- function(tbl, genes) { v <- tbl[genes]; v[is.na(v)] <- Inf; unname(v) }

a7_act <- pmax(seq_len(nrow(a7)), safe(fa9, a7$GenA), safe(fa9, a7$GenB))
a9_act <- pmax(seq_len(nrow(a9)), safe(fa7, a9$GenA), safe(fa7, a9$GenB))

# Remove edges that never activate (Inf — endpoint never appears in other net)
a7_act <- a7_act[is.finite(a7_act)]
a9_act <- a9_act[is.finite(a9_act)]

# ── Build sorted activation lists for O(log n) cumulative counts ─────────────
a7_sorted <- sort(a7_act)
a9_sorted <- sort(a9_act)

edge_count <- function(sorted_acts, N) findInterval(N, sorted_acts)

shared_nodes_at <- function(N) sum(fa7_v <= N & fa9_v <= N)

# ── Scan: find where d(N) = a9_edges(N) - a7_edges(N) = 0 or changes sign ──
cat("[4] Scanning for convergence...\n\n")

max_N <- max(nrow(a7), nrow(a9))

# Coarse pass: every 500
coarse <- seq(0, max_N, by = 500)
e7c <- sapply(coarse, function(N) edge_count(a7_sorted, N))
e9c <- sapply(coarse, function(N) edge_count(a9_sorted, N))
dc  <- e9c - e7c

# Find sign-change brackets and exact zeros in coarse pass
exact_coarse  <- coarse[dc == 0]
bracket_idx   <- which(diff(sign(dc[dc != 0])) != 0)  # brackets ignoring zeros

cat("=== Coarse scan (step 500) ===\n")
if (length(exact_coarse) > 0) {
  cat("Exact zeros at N =", paste(exact_coarse, collapse=", "), "\n")
} else {
  cat("No exact zeros found at step=500\n")
}

# Fine-grained search in each bracket
fine_results <- list()

if (length(bracket_idx) > 0) {
  for (bi in bracket_idx) {
    n_lo <- coarse[bi]
    n_hi <- coarse[bi + 1]
    cat(sprintf("  Bracket: N=[%d, %d]  d=[%d, %d] — refining...\n",
                n_lo, n_hi, dc[bi], dc[bi+1]))
    fine <- seq(n_lo, n_hi, by = 1)
    e7f  <- sapply(fine, function(N) edge_count(a7_sorted, N))
    e9f  <- sapply(fine, function(N) edge_count(a9_sorted, N))
    df   <- e9f - e7f
    exact_fine <- fine[df == 0]
    if (length(exact_fine) > 0) {
      for (N in exact_fine) {
        sn <- shared_nodes_at(N)
        cat(sprintf("  ✓ CONVERGENCE at N=%d: shared_nodes=%d  A7_edges=%d  A9_edges=%d\n",
                    N, sn, edge_count(a7_sorted, N), edge_count(a9_sorted, N)))
        fine_results[[length(fine_results)+1]] <- data.frame(
          N=N, shared_nodes=sn,
          a7_edges=edge_count(a7_sorted,N),
          a9_edges=edge_count(a9_sorted,N))
      }
    } else {
      # Report the two boundary points of the sign change
      sc <- which(diff(sign(df)) != 0)[1]
      cat(sprintf("  No exact zero — sign changes between N=%d (d=%d) and N=%d (d=%d)\n",
                  fine[sc], df[sc], fine[sc+1], df[sc+1]))
    }
  }
}

# Also check very small N (1:500) with step=1 in case crossing happens early
cat("\n=== Fine scan N=1 to 500 (step 1) ===\n")
tiny <- 1:500
e7t <- sapply(tiny, function(N) edge_count(a7_sorted, N))
e9t <- sapply(tiny, function(N) edge_count(a9_sorted, N))
dt  <- e9t - e7t
exact_tiny <- tiny[dt == 0 & e7t > 0]   # non-trivial zeros
if (length(exact_tiny) > 0) {
  cat("Convergence points (non-trivial) at small N:\n")
  for (N in exact_tiny) {
    sn <- shared_nodes_at(N)
    cat(sprintf("  ✓ N=%d: shared_nodes=%d  A7_edges=%d  A9_edges=%d\n",
                N, sn, e7t[N], e9t[N]))
  }
} else {
  cat("No non-trivial convergence in N=1..500\n")
  # Show the minimum absolute difference
  best_tiny <- tiny[which.min(abs(dt))]
  cat(sprintf("  Closest: N=%d  A7=%d  A9=%d  diff=%d\n",
              best_tiny, e7t[best_tiny], e9t[best_tiny], dt[best_tiny]))
}

cat("\n=== Summary ===\n")
if (length(fine_results) > 0 || length(exact_coarse) > 0) {
  cat("Convergence EXISTS — see points above.\n")
} else {
  # Report the global minimum |diff| point
  all_N  <- seq(0, max_N, by = 100)
  e7a    <- sapply(all_N, function(N) edge_count(a7_sorted, N))
  e9a    <- sapply(all_N, function(N) edge_count(a9_sorted, N))
  da     <- e9a - e7a
  best_i <- which.min(abs(da))
  sn     <- shared_nodes_at(all_N[best_i])
  cat(sprintf("NO convergence found. Global minimum |diff| = %d at N=%d\n",
              abs(da[best_i]), all_N[best_i]))
  cat(sprintf("  shared_nodes=%d  A7_edges=%d  A9_edges=%d\n",
              sn, e7a[best_i], e9a[best_i]))
}

cat("\nDone.\n")
