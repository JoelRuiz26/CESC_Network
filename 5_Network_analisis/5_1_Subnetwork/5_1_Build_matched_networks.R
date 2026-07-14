################################################################################
# 5_1_Build_matched_networks.R
#
# Builds the two comparable network graphs used downstream, using ONLY the new
# MI-sorted, deduplicated edgelists in 4_Network:
#   4_Network/A7/network_A7_AP_sort_unique.txt
#   4_Network/A9/network_A9_AP_sort_unique.txt
#
# Matching strategy (edge- AND connectivity-matched):
#   1) Identify the SMALLER network (fewer edges). Take its giant connected
#      component; its edge count defines the target N. Save it as-is.
#   2) For the LARGER network, add edges in decreasing-MI order until its giant
#      connected component reaches exactly N edges (merely cutting at N leaves it
#      disconnected, so we keep adding top-MI edges until the giant component
#      itself holds N edges). Save that giant component.
#
# Result: both graphs are single connected components with the same edge count N.
#
# Outputs (RDS names preserved so all downstream scripts pick them up unchanged):
#   5_2_4_Subnetwork_A7_max.rds , 5_2_4_Subnetwork_A9_max.rds
#   network_AX_AP_matched.graphml  (in 4_Network)
################################################################################

suppressPackageStartupMessages(library(igraph))

net_dir <- "~/CESC_Network/4_Network"
out_dir <- "~/CESC_Network/5_Network_analisis/5_1_Subnetwork"

## ── Load the new deduplicated, MI-sorted edgelists ────────────────────────────
read_el <- function(net)
  read.delim(file.path(net_dir, net, sprintf("network_%s_AP_sort_unique.txt", net)),
             header = TRUE, stringsAsFactors = FALSE)          # Regulator Target MI pvalue (MI desc)

el <- list(A7 = read_el("A7"), A9 = read_el("A9"))
n_edges <- sapply(el, nrow)

small <- names(which.min(n_edges))
large <- names(which.max(n_edges))
cat(sprintf("Edges  ->  A7: %d | A9: %d\n", n_edges["A7"], n_edges["A9"]))
cat(sprintf("Smaller network = %s | Larger network = %s\n\n", small, large))

## ── Helpers ───────────────────────────────────────────────────────────────────
graph_topk    <- function(d, k) simplify(graph_from_data_frame(d[seq_len(k), 1:2], directed = FALSE))
giant_of      <- function(g) { cm <- components(g); induced_subgraph(g, which(cm$membership == which.max(cm$csize))) }
giant_nedges  <- function(g) ecount(giant_of(g))

## ── 1) Smaller network: giant component defines target N ──────────────────────
g_small <- giant_of(graph_topk(el[[small]], nrow(el[[small]])))
N <- ecount(g_small)
cat(sprintf("[1] %s giant component: %d nodes / %d edges  ->  target N = %d\n",
            small, vcount(g_small), N, N))

## ── 2) Larger network: grow top-MI edges until giant component has N edges ─────
k <- N
repeat {
  e <- giant_nedges(graph_topk(el[[large]], k))
  if (e >= N) break
  k <- k + (N - e)                        # jump forward by the current deficit
}
while (k > N && giant_nedges(graph_topk(el[[large]], k - 1)) >= N) k <- k - 1   # smallest such k

g_large <- giant_of(graph_topk(el[[large]], k))
cat(sprintf("[2] %s needs top-%d edges; giant component: %d nodes / %d edges\n",
            large, k, vcount(g_large), ecount(g_large)))
if (ecount(g_large) != N)
  warning(sprintf("Giant component of %s has %d edges, not exactly N=%d (edge merge overshoot).",
                  large, ecount(g_large), N))

## ── Save (map back to A7 / A9 names) ──────────────────────────────────────────
graphs <- setNames(list(g_small, g_large), c(small, large))
for (net in c("A7", "A9")) {
  g <- graphs[[net]]
  write_graph(g, file.path(net_dir, sprintf("network_%s_AP_matched.graphml", net)), format = "graphml")
  saveRDS(g,   file.path(out_dir, sprintf("5_2_4_Subnetwork_%s_max.rds", net)))
}

cat("\n[OK] Saved. Both graphs are connected and edge-matched:\n")
for (net in c("A7", "A9")) {
  g <- graphs[[net]]
  cat(sprintf("  %s: %d nodes / %d edges / %d component(s)\n",
              net, vcount(g), ecount(g), components(g)$no))
}
