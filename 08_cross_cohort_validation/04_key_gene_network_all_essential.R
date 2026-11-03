### ============================================================
### 04_key_gene_network_all_essential.R  (EXPLORATORY / TEMPORARY)
###
### FULL ("complete") essential hub sets -- i.e. including pan-essential
### genes, not just the selective/non-pan-essential subset -- restricted to
### genes whose category (A7-selective hub / A9-selective hub / hub in
### both) is IDENTICAL in both CESC and HTMCP-CC ("consistent" across
### cohorts). NO synthetic shortest-path connector nodes this time -- only
### REAL direct ARACNe edges among these genes are shown, so every node on
### the plot is one of the actual key genes and is labelled; genes with no
### direct link to another key gene appear as small separate components.
### Non-pan-essential ("selective") genes get a gold ring.
###
### -> figures/04_key_gene_network_all_essential.png / .pdf
### ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(vroom); library(igraph)
  library(tidygraph); library(ggraph); library(patchwork)
})

out_dir <- "~/CESC_Network/08_cross_cohort_validation"
dep_dir <- "~/CESC_Network/07_dependency/results"
dep_htmcp_dir <- "~/CESC_Network/07_dependency/validation_uganda/results"
topo_dir <- "~/CESC_Network/05_topology/results"
FF <- "Times New Roman"
set.seed(42)

pan <- vroom(file.path(dep_dir, "6_1_common_essential_genes.tsv"), delim = "\t", show_col_types = FALSE)$gene_u %>% toupper()

elbow <- function(v) {
  n <- length(v); x <- (seq_len(n) - 1) / (n - 1); y <- (v - v[1]) / (v[n] - v[1])
  list(idx = which.max(abs(y - x)), val = v[which.max(abs(y - x))])
}

## complete essential-hub classification: NO pan-essential filtering here
classify_cohort <- function(res_dir) {
  tab <- vroom(file.path(res_dir, "6_2_gene_dep_degree_table.tsv"), show_col_types = FALSE) %>%
    filter(!is.na(degree), !is.na(mean_dep)) %>% distinct(network, gene, .keep_all = TRUE)
  nets <- sort(unique(tab$network))
  info <- map(nets, function(nk) {
    d <- tab %>% filter(network == nk)
    ek <- elbow(sort(d$degree, decreasing = TRUE)); ed <- elbow(sort(d$mean_dep))
    list(net = nk, ek = ek, ed = ed)
  }) %>% set_names(nets)
  flags <- tab %>% left_join(map_dfr(info, ~ tibble(network = .x$net, k_cut = .x$ek$val, d_cut = .x$ed$val)), by = "network") %>%
    mutate(hub = degree >= k_cut, ess = hub & mean_dep <= d_cut, net = str_remove(network, "HPV-")) %>%
    select(gene, net, hub, ess) %>% pivot_wider(names_from = net, values_from = c(hub, ess)) %>%
    mutate(across(starts_with("hub_") | starts_with("ess_"), ~ replace_na(.x, FALSE)))
  dep <- vroom(file.path(res_dir, "6_7_hub_dependency_medians.tsv"), show_col_types = FALSE) %>% select(gene, gene_u)
  dep %>% inner_join(flags, by = "gene") %>% filter(ess_A7 | ess_A9) %>%
    mutate(grp = case_when(ess_A7 & ess_A9 ~ "Both", ess_A7 ~ "A7", ess_A9 ~ "A9", TRUE ~ "None"),
           is_pan = toupper(gene_u) %in% pan)
}

cesc  <- classify_cohort(dep_dir)
htmcp <- classify_cohort(dep_htmcp_dir)
consistent <- inner_join(cesc %>% select(gene, grp_cesc = grp, is_pan), htmcp %>% select(gene, grp_htmcp = grp), by = "gene") %>%
  mutate(consistent_group = case_when(grp_cesc == "Both" & grp_htmcp == "Both" ~ "Both (A7 & A9)",
                                      grp_cesc == "A7" & grp_htmcp == "A7" ~ "HPV-A7",
                                      grp_cesc == "A9" & grp_htmcp == "A9" ~ "HPV-A9",
                                      TRUE ~ NA_character_)) %>%
  filter(!is.na(consistent_group)) %>% select(gene, consistent_group, is_pan)

cat("Consistent COMPLETE essential hubs (cross-cohort):\n")
print(table(consistent$consistent_group, ifelse(consistent$is_pan, "pan-essential", "selective")))
write_tsv(consistent, file.path(out_dir, "results/04_consistent_complete_essential_hubs.tsv"))

build_clade_data <- function(clade, own_label, own_col, both_col) {
  key <- consistent %>% filter(consistent_group %in% c(own_label, "Both (A7 & A9)"))
  full_g <- readRDS(file.path(topo_dir, sprintf("5_2_4_Subnetwork_%s_max.rds", clade)))
  key_genes <- intersect(key$gene, V(full_g)$name)
  net_g0 <- induced_subgraph(full_g, key_genes)

  ## Integrate EVERY disconnected component (isolated single genes AND
  ## small already-connected clusters, e.g. a 2-3 gene chain) into the main
  ## network via a shortest path (in the full ARACNe network) from that
  ## component's closest member to the giant component's closest member,
  ## adding whichever intermediate "connector" genes are needed -- so
  ## nothing is left stranded on its own.
  comps0 <- igraph::components(net_g0)
  giant_id <- which.max(comps0$csize)
  giant_names <- names(comps0$membership)[comps0$membership == giant_id]
  other_ids <- setdiff(seq_len(comps0$no), giant_id)

  edge_extra <- list(); connector_genes <- character(0)
  for (cid in other_ids) {
    comp_names <- names(comps0$membership)[comps0$membership == cid]
    dv <- distances(full_g, v = comp_names, to = giant_names)
    best <- which(dv == min(dv), arr.ind = TRUE)[1, ]
    src <- comp_names[best[1]]; target <- giant_names[best[2]]
    sp <- names(shortest_paths(full_g, from = src, to = target, output = "vpath")$vpath[[1]])
    if (length(sp) > 1) {
      edge_extra[[length(edge_extra) + 1]] <- tibble(from = sp[-length(sp)], to = sp[-1])
      connector_genes <- union(connector_genes, setdiff(sp, key_genes))
    }
  }
  edges0 <- igraph::as_data_frame(net_g0, what = "edges")
  edges_all <- bind_rows(edges0, bind_rows(edge_extra)) %>%
    mutate(a = pmin(from, to), b = pmax(from, to)) %>% distinct(a, b) %>% transmute(from = a, to = b)
  all_names <- union(key_genes, connector_genes)
  net_g <- graph_from_data_frame(edges_all, directed = FALSE, vertices = tibble(name = all_names))

  louv <- readRDS(file.path("~/CESC_Network/06_modules/results", sprintf("5_1_1_Louvain_%s_graph_max.rds", clade)))
  louv_tbl <- tibble(gene = V(louv)$name, pagerank = V(louv)$PageRank)

  node_tbl <- tibble(gene = V(net_g)$name) %>%
    left_join(key %>% select(gene, consistent_group, is_pan), by = "gene") %>%
    left_join(louv_tbl, by = "gene") %>%
    mutate(category = case_when(consistent_group == "Both (A7 & A9)" ~ "Both (A7 & A9)",
                                consistent_group == own_label ~ own_label,
                                TRUE ~ "Gen conector"),
           selective = !is.na(is_pan) & !is_pan,
           is_connector = category == "Gen conector",
           pagerank = ifelse(is.na(pagerank), min(pagerank, na.rm = TRUE), pagerank))

  list(net_g = net_g, node_tbl = node_tbl, own_col = own_col, both_col = both_col, own_label = own_label)
}

build_plot <- function(d, clade_title) {
  edges <- igraph::as_data_frame(d$net_g, what = "edges") %>% transmute(from, to)
  nodes <- d$node_tbl %>%
    transmute(name = gene, category, is_connector,
              size = ifelse(is_connector, min(pagerank) * 0.4, pagerank),
              hub_type = if_else(selective, "Selectivo (no pan-esencial)", "Pan-esencial"))

  g <- tbl_graph(nodes = nodes, edges = edges, node_key = "name", directed = FALSE)

  ## Piecewise layout: a single global "stress" call places disconnected
  ## components with no shared reference frame, which can fling a small
  ## component far away and stretch the whole coordinate range into mostly
  ## empty space. Instead: lay out the giant component on its own (the
  ## anchor), lay out every other multi-node component on its own and slot
  ## it into a tidy row underneath, and spread single, unconnected genes
  ## into their own row below that -- then combine into one x/y vector for
  ## a manual layout, so edges still resolve correctly by node name.
  comps <- igraph::components(d$net_g)
  giant_id <- which.max(comps$csize)
  giant_names <- names(comps$membership)[comps$membership == giant_id]
  other_ids <- setdiff(seq_len(comps$no), giant_id)

  xy_giant <- graphlayouts::layout_with_stress(induced_subgraph(d$net_g, giant_names))
  pos <- tibble(name = giant_names, x = xy_giant[, 1], y = xy_giant[, 2])
  y_bottom <- min(pos$y); x_cursor <- min(pos$x); row_y <- y_bottom - 2.0

  small_ids <- other_ids[comps$csize[other_ids] > 1]
  for (cid in small_ids) {
    nm <- names(comps$membership)[comps$membership == cid]
    sub_g <- induced_subgraph(d$net_g, nm)
    xy <- if (vcount(sub_g) > 2) graphlayouts::layout_with_stress(sub_g) else layout_as_star(sub_g)
    sx <- xy[, 1] - min(xy[, 1]) + x_cursor
    sy <- xy[, 2] - mean(range(xy[, 2])) + row_y
    pos <- bind_rows(pos, tibble(name = nm, x = sx, y = sy))
    x_cursor <- max(sx) + 1.4
  }

  iso_names <- names(comps$membership)[comps$csize[comps$membership] == 1]
  if (length(iso_names) > 0) {
    pos <- bind_rows(pos, tibble(name = iso_names,
                                 x = seq(min(pos$x), max(x_cursor, max(pos$x)), length.out = length(iso_names)),
                                 y = row_y - 2.2))
  }
  pos <- pos[match(nodes$name, pos$name), ]
  ## scale up inter-node spacing (in DATA coordinates) so much bigger text/
  ## nodes still have room -- this is what actually makes labels look bigger
  ## relative to the page, since font size is fixed in points regardless of
  ## how large the saved canvas is
  pos$x <- pos$x * 2.5; pos$y <- pos$y * 2.5
  lay <- create_layout(g, layout = "manual", x = pos$x, y = pos$y)

  cat_cols <- c(setNames(d$own_col, d$own_label), "Both (A7 & A9)" = d$both_col, "Gen conector" = "#D9D3C8")
  ring_cols <- c("Selectivo (no pan-esencial)" = "#E8A33D", "Pan-esencial" = "white")

  ggraph(lay) +
    geom_edge_diagonal(colour = "#B7B2A6", linewidth = 0.6, alpha = 0.8) +
    geom_node_point(data = ~ filter(.x, is_connector), fill = "#D9D3C8", colour = "#8F887A",
                    shape = 21, size = 9, stroke = 0.7) +
    geom_node_point(data = ~ filter(.x, !is_connector), aes(fill = category, colour = hub_type, size = size),
                    shape = 21, stroke = 2.6) +
    scale_fill_manual(values = cat_cols, name = "Categoría del gen", breaks = c(d$own_label, "Both (A7 & A9)")) +
    scale_colour_manual(values = ring_cols, name = "Tipo de hub",
                        guide = guide_legend(override.aes = list(fill = "grey70", size = 11))) +
    scale_size_continuous(range = c(16, 32), guide = "none") +
    geom_node_text(data = ~ filter(.x, is_connector), aes(label = name), family = FF,
                   size = 9.5, colour = "grey35", repel = TRUE, max.overlaps = 12, segment.size = 0.3,
                   box.padding = 0.35, force = 1, max.time = 2, max.iter = 20000) +
    geom_node_text(data = ~ filter(.x, !is_connector), aes(label = name), family = FF, fontface = "bold",
                   size = 13, colour = "grey5", repel = TRUE, max.overlaps = 18, segment.size = 0.4,
                   segment.colour = "grey55", bg.color = alpha("white", 0.85), bg.r = 0.16,
                   box.padding = 0.45, force = 1, max.time = 2, max.iter = 20000) +
    labs(title = clade_title) +
    theme_void(base_family = FF) +
    theme(plot.title = element_text(size = 42, face = "bold", hjust = 0.5, margin = margin(b = 18)),
          legend.title = element_text(size = 26, face = "bold"), legend.text = element_text(size = 25),
          legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0.1, "cm"),
          plot.margin = margin(24, 24, 14, 24))
}

dA7 <- build_clade_data("A7", "HPV-A7", "#1F4E8C", "#5B3A8E")
dA9 <- build_clade_data("A9", "HPV-A9", "#8C2F39", "#5B3A8E")

pA <- build_plot(dA7, "Red HPV-A7")
pB <- build_plot(dA9, "Red HPV-A9")

fig <- (pA | pB) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = FF, size = 22, face = "bold"))

ggsave(file.path(out_dir, "figures/04_key_gene_network_all_essential.png"), fig, width = 48, height = 28, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "figures/04_key_gene_network_all_essential.pdf"), fig, width = 48, height = 28, device = cairo_pdf, bg = "white", limitsize = FALSE)
message("Saved: figures/04_key_gene_network_all_essential.png / .pdf")
