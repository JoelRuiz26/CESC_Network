### ============================================================
### 03_key_gene_process_network.R  (EXPLORATORY / TEMPORARY)
###
### For the cross-cohort "consistent" selective hub genes
### (01_consistent_selective_hubs.tsv), builds two REAL network plots (one
### per clade). Since these key genes are rarely directly linked to each
### other (they were picked for degree+dependency, not mutual adjacency),
### each pair of key genes is connected via a minimum-spanning tree over
### shortest-path distances in the full ARACNe network, and every gene that
### lies ON one of those shortest paths is added to the plot as a small
### "connector" node -- so the plot shows an actual connected subnetwork,
### not just isolated dots. Organic force-directed (fr) layout, deep jewel-
### tone palette. No biological-process nodes (gene-gene network only).
###
### -> figures/03_key_gene_process_network.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(igraph)
  library(tidygraph); library(ggraph); library(patchwork)
})

out_dir <- "~/CESC_Network/08_cross_cohort_validation"
mod_dir <- "~/CESC_Network/06_modules/results"
topo_dir <- "~/CESC_Network/05_topology/results"
FF <- "Times New Roman"
set.seed(42)

genes_tbl <- vroom(file.path(out_dir, "results/01_consistent_selective_hubs.tsv"), show_col_types = FALSE)

## ── Steiner-tree-like connector subgraph: MST over shortest-path distances
## between key genes, then union of the actual shortest paths + any
## incidental direct edges among all involved nodes
build_connector_graph <- function(key_genes, full_g) {
  key_genes <- intersect(key_genes, V(full_g)$name)
  dmat <- distances(full_g, v = key_genes, to = key_genes)
  cg <- graph_from_adjacency_matrix(dmat, mode = "undirected", weighted = TRUE, diag = FALSE)
  V(cg)$name <- key_genes
  tree <- mst(cg, weights = E(cg)$weight)
  ed <- igraph::as_data_frame(tree, what = "edges")
  allnodes <- key_genes
  for (i in seq_len(nrow(ed))) {
    sp <- shortest_paths(full_g, from = ed$from[i], to = ed$to[i], output = "vpath")$vpath[[1]]
    allnodes <- union(allnodes, names(sp))
  }
  induced_subgraph(full_g, allnodes)
}

build_clade_data <- function(clade, own_label, own_col, both_col) {
  key <- genes_tbl %>% filter(consistent_group %in% c(own_label, "Both (A7 & A9)"))
  full_g <- readRDS(file.path(topo_dir, sprintf("5_2_4_Subnetwork_%s_max.rds", clade)))
  net_g <- build_connector_graph(key$gene, full_g)

  louv <- readRDS(file.path(mod_dir, sprintf("5_1_1_Louvain_%s_graph_max.rds", clade)))
  louv_tbl <- tibble(gene = V(louv)$name, pagerank = V(louv)$PageRank)

  node_tbl <- tibble(gene = V(net_g)$name) %>%
    left_join(key %>% select(gene, consistent_group), by = "gene") %>%
    left_join(louv_tbl, by = "gene") %>%
    mutate(category = case_when(consistent_group == "Both (A7 & A9)" ~ "Both (A7 & A9)",
                                consistent_group == own_label ~ own_label,
                                TRUE ~ "Connector gene"),
           pagerank = ifelse(is.na(pagerank), min(pagerank, na.rm = TRUE), pagerank))

  list(net_g = net_g, node_tbl = node_tbl, own_col = own_col, both_col = both_col, own_label = own_label)
}

build_plot <- function(d, clade_title) {
  edges <- igraph::as_data_frame(d$net_g, what = "edges") %>% transmute(from, to)
  nodes <- d$node_tbl %>% transmute(name = gene, category, size = pagerank, key = category != "Connector gene")

  g <- tbl_graph(nodes = nodes, edges = edges, node_key = "name", directed = FALSE)
  lay <- create_layout(g, layout = "kk")

  cat_cols <- c(setNames(d$own_col, d$own_label), "Both (A7 & A9)" = d$both_col, "Connector gene" = "#D9D3C8")

  ggraph(lay) +
    geom_edge_arc(strength = 0.12, colour = "#A8A398", linewidth = 0.75, alpha = 0.65) +
    geom_node_point(data = ~ filter(.x, !key), fill = "#D9D3C8", colour = "#8F887A",
                    shape = 21, size = 5.2, stroke = 0.6) +
    geom_node_point(data = ~ filter(.x, key), aes(fill = category, size = size),
                    shape = 21, colour = "white", stroke = 1.3) +
    scale_fill_manual(values = cat_cols, name = "Key gene category",
                      breaks = c(d$own_label, "Both (A7 & A9)")) +
    scale_size_continuous(range = c(9, 16), guide = "none") +
    geom_node_text(data = ~ filter(.x, !key), aes(label = name), family = FF,
                   size = 3.2, colour = "#7A7264", repel = TRUE, max.overlaps = Inf, segment.size = 0.2) +
    geom_node_text(data = ~ filter(.x, key), aes(label = name), family = FF, fontface = "bold",
                   size = 5.2, colour = "grey12", repel = TRUE, max.overlaps = Inf, segment.size = 0.3,
                   bg.color = "white", bg.r = 0.13) +
    labs(title = clade_title) +
    theme_void(base_family = FF) +
    theme(plot.title = element_text(size = 21, face = "bold", hjust = 0.5, margin = margin(b = 12)),
          legend.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 12),
          legend.position = "bottom",
          plot.margin = margin(20, 20, 10, 20))
}

dA7 <- build_clade_data("A7", "HPV-A7", "#1F4E8C", "#5B3A8E")
dA9 <- build_clade_data("A9", "HPV-A9", "#8C2F39", "#5B3A8E")

pA <- build_plot(dA7, "HPV-A7 network")
pB <- build_plot(dA9, "HPV-A9 network")

fig <- (pA | pB) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FF, size = 22, face = "bold"))

ggsave(file.path(out_dir, "figures/03_key_gene_process_network.png"), fig, width = 20, height = 11, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "figures/03_key_gene_process_network.pdf"), fig, width = 20, height = 11, device = cairo_pdf, bg = "white")
message("Saved: figures/03_key_gene_process_network.png / .pdf")
