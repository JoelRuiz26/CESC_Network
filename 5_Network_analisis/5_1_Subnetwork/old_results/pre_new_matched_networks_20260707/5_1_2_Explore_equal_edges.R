################################################################################
# Exploratory script: find the MI-rank cutoff at which both networks share
# the same node identity AND the same number of edges.
#
# Logic:
#   Both networks are pre-sorted by MI (highest first) and already deduplicated
#   and restricted to a shared 17,480-gene universe.
#
#   For a given depth N (top-N edges considered in each network):
#     - shared_nodes(N) = genes that have appeared in BOTH A7[1:N] AND A9[1:N]
#     - a7_edges(N)     = edges in A7[1:N] where BOTH endpoints are in shared_nodes(N)
#     - a9_edges(N)     = edges in A9[1:N] where BOTH endpoints are in shared_nodes(N)
#
#   Key insight: an A7 edge at rank r contributes to a7_edges(N) only when
#     N >= max(r, first_appearance_in_A9[GenA], first_appearance_in_A9[GenB])
#   This lets us precompute each edge's "activation depth" and build cumulative
#   counts efficiently without looping over all N.
################################################################################

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

DIR_NET <- "/STORAGE/csbig/jruiz/cesc_network/cervix_ARACNE"
DIR_OUT <- "~/CESC_Network/5_Network_analisis/5_1_Subnetwork"

# ── 1) Load unique sorted networks ────────────────────────────────────────────
cat("[1] Loading networks...\n")
read_net <- function(path) {
  vroom(path,
        col_names = c("GenA", "GenB", "MI", "pvalue"),
        skip = 1, show_col_types = FALSE) %>%
    select(GenA, GenB, MI)
}

a7 <- read_net(file.path(DIR_NET, "A7_network_sort_unique.txt"))
a9 <- read_net(file.path(DIR_NET, "A9_network_sort_unique.txt"))

cat("  A7 unique edges:", nrow(a7), "\n")
cat("  A9 unique edges:", nrow(a9), "\n\n")

# ── 2) First-appearance rank of each node in each network ─────────────────────
cat("[2] Computing first-appearance ranks...\n")

first_rank <- function(net) {
  genes  <- c(rbind(net$GenA, net$GenB))         # interleaved: A, B, A, B, ...
  ranks  <- rep(seq_len(nrow(net)), each = 2)
  tapply(ranks, genes, min)
}

fa7 <- first_rank(a7)   # named int vector: gene -> first rank in A7
fa9 <- first_rank(a9)   # named int vector: gene -> first rank in A9

all_nodes <- union(names(fa7), names(fa9))
cat("  Total unique genes across both networks:", length(all_nodes), "\n\n")

# Safe lookup: return Inf if gene absent (never appears in that network)
lookup <- function(table, genes) {
  vals <- table[genes]
  vals[is.na(vals)] <- Inf
  unname(vals)
}

# ── 3) Precompute "activation depth" for every edge in each network ───────────
# An A7 edge at rank r is counted in a7_edges(N) when:
#   N >= max(r, fa9[GenA], fa9[GenB])   <- both endpoints must be in A9[1:N] too
cat("[3] Computing activation depths...\n")

a7_activation <- pmax(
  seq_len(nrow(a7)),
  lookup(fa9, a7$GenA),
  lookup(fa9, a7$GenB)
)

a9_activation <- pmax(
  seq_len(nrow(a9)),
  lookup(fa7, a9$GenA),
  lookup(fa7, a9$GenB)
)

# ── 4) Cumulative edge counts at checkpoints ──────────────────────────────────
cat("[4] Building cumulative profile...\n")

max_N <- max(nrow(a7), nrow(a9))

# Checkpoints: coarse (every 2000) + key landmarks
checkpoints <- unique(sort(c(
  seq(2000, max_N, by = 2000),
  115784, 117407, 149929,          # known cutoffs
  nrow(a7), nrow(a9)
)))
checkpoints <- checkpoints[checkpoints <= max_N]

# Precompute sorted activation depths for fast cumulative counts
a7_act_sorted <- sort(a7_activation[is.finite(a7_activation)])
a9_act_sorted <- sort(a9_activation[is.finite(a9_activation)])

# shared_nodes(N): genes where fa7 <= N AND fa9 <= N
shared_at <- function(N) {
  sum(lookup(fa7, all_nodes) <= N & lookup(fa9, all_nodes) <= N)
}

results <- tibble(
  N            = checkpoints,
  shared_nodes = vapply(checkpoints, shared_at, integer(1)),
  a7_edges     = findInterval(checkpoints, a7_act_sorted),
  a9_edges     = findInterval(checkpoints, a9_act_sorted),
) %>%
  mutate(diff = a9_edges - a7_edges)

cat("\n════════════════════════════════════════════════════════\n")
cat(" EXPLORATORY PROFILE (shared-node induced subgraph)\n")
cat("════════════════════════════════════════════════════════\n")
print(results, n = Inf)

# ── 5) Find crossover / minimum-difference point ──────────────────────────────
cat("\n[5] Locating crossover and minimum-difference points...\n")

# Point of minimum absolute difference
min_diff_idx <- which.min(abs(results$diff))
cat(sprintf("  Min |diff| at N = %d  →  A7: %d  A9: %d  diff: %d  shared nodes: %d\n",
            results$N[min_diff_idx],
            results$a7_edges[min_diff_idx],
            results$a9_edges[min_diff_idx],
            results$diff[min_diff_idx],
            results$shared_nodes[min_diff_idx]))

# Sign changes (crossover points)
sign_changes <- which(diff(sign(results$diff)) != 0)
if (length(sign_changes) > 0) {
  cat("  Crossover points (diff changes sign) between checkpoints:\n")
  for (i in sign_changes) {
    cat(sprintf("    N = %d (diff=%d) → N = %d (diff=%d)\n",
                results$N[i],   results$diff[i],
                results$N[i+1], results$diff[i+1]))
  }
} else {
  cat("  No exact crossover found in the checkpoint grid.\n")
}

# ── 6) Fine-grained search around the minimum-diff region ────────────────────
cat("\n[6] Fine-grained search (step = 100) around minimum-diff region...\n")

n_window <- results$N[max(1, min_diff_idx - 1)]
n_end    <- results$N[min(nrow(results), min_diff_idx + 1)]
fine_N   <- seq(n_window, n_end, by = 100)

fine_results <- tibble(
  N            = fine_N,
  shared_nodes = vapply(fine_N, shared_at, integer(1)),
  a7_edges     = findInterval(fine_N, a7_act_sorted),
  a9_edges     = findInterval(fine_N, a9_act_sorted)
) %>%
  mutate(diff = a9_edges - a7_edges)

# Exact equality
exact <- fine_results %>% filter(diff == 0)
if (nrow(exact) > 0) {
  cat("  ✓ EXACT EQUAL EDGE COUNTS FOUND:\n")
  print(exact)
} else {
  cat("  No exact match at step=100. Closest points:\n")
  best <- fine_results %>% slice(which.min(abs(diff)))
  print(best)
  sign2 <- which(diff(sign(fine_results$diff)) != 0)
  if (length(sign2) > 0) {
    cat("  Crossover between:\n")
    cat(sprintf("    N=%d (diff=%d)  and  N=%d (diff=%d)\n",
                fine_results$N[sign2[1]],   fine_results$diff[sign2[1]],
                fine_results$N[sign2[1]+1], fine_results$diff[sign2[1]+1]))
  }
}

# ── 7) Plot ───────────────────────────────────────────────────────────────────
cat("\n[7] Generating plot...\n")

plot_df <- results %>%
  select(N, a7_edges, a9_edges) %>%
  tidyr::pivot_longer(cols = c(a7_edges, a9_edges),
                      names_to = "Network", values_to = "Edges") %>%
  mutate(Network = recode(Network, a7_edges = "HPV-A7", a9_edges = "HPV-A9"))

p <- ggplot(plot_df, aes(x = N, y = Edges, color = Network, linetype = Network)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = results$N[min_diff_idx], linetype = "dotted",
             color = "grey40", linewidth = 0.8) +
  annotate("text", x = results$N[min_diff_idx], y = max(plot_df$Edges) * 0.5,
           label = paste0("N=", results$N[min_diff_idx]),
           hjust = -0.1, size = 5, color = "grey30") +
  scale_color_manual(values = c("HPV-A7" = "#003366", "HPV-A9" = "#CC0000")) +
  scale_linetype_manual(values = c("HPV-A7" = "solid", "HPV-A9" = "dashed")) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Shared-node induced edges vs. MI-rank depth",
    subtitle = "Edges counted only for genes present in both networks at depth N",
    x        = "Depth N (top-N edges per network)",
    y        = "Edges in shared-node subgraph",
    color    = "Network", linetype = "Network"
  ) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    axis.title    = element_text(size = 13, face = "bold"),
    axis.text     = element_text(size = 11),
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = alpha("white", 0.9), color = "grey60"),
    legend.title  = element_text(size = 12, face = "bold"),
    legend.text   = element_text(size = 11)
  )

ggsave(file.path(DIR_OUT, "5_1_2_Explore_equal_edges.png"),
       p, width = 10, height = 6, dpi = 300, bg = "white")

cat("Saved: 5_1_2_Explore_equal_edges.png\n")
cat("\nDone.\n")
