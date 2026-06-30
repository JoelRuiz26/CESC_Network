################################################################################
# Topological metrics evolution as edges increase from 20k to 117,407 (step 500)
################################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(vroom)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

DIR_NET <- "/STORAGE/csbig/jruiz/cesc_network/cervix_ARACNE"
DIR_OUT <- "~/CESC_Network/5_Network_analisis/5_2_Network_topology"

# ── Load full networks ─────────────────────────────────────────────────────────
cat("[1] Loading networks...\n")
read_net <- function(f)
  vroom(f, col_names = c("GenA","GenB","MI","pvalue"),
        skip = 1, show_col_types = FALSE) %>% select(GenA, GenB, MI)

a7_full <- read_net(file.path(DIR_NET, "A7_network_sort.txt"))
a9_full <- read_net(file.path(DIR_NET, "A9_network_sort.txt"))

# ── Steps ─────────────────────────────────────────────────────────────────────
steps <- seq(20000, 115784, by = 500)
cat("[2] Computing metrics at", length(steps), "steps (step=500)...\n")

# ── Metric computation ─────────────────────────────────────────────────────────
compute_metrics <- function(el, name, N) {
  g    <- graph_from_data_frame(el, directed = FALSE)
  g    <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  comp <- components(g)
  deg  <- igraph::degree(g)

  list(
    network       = name,
    N             = N,
    Nodos         = vcount(g),
    Densidad      = edge_density(g),
    Componentes   = comp$no,
    Grado_medio   = mean(deg),
    Clustering    = transitivity(g, type = "average"),
    Transitividad = transitivity(g, type = "global"),
    Asortatividad = assortativity_degree(g),
    GC_size       = max(comp$csize),
    GC_pct        = max(comp$csize) / vcount(g) * 100,
    Diametro      = diameter(g, unconnected = TRUE),
    Mean_distance = mean_distance(g, unconnected = FALSE)
  )
}

jaccard_nodes_fn <- function(g1, g2) {
  a <- V(g1)$name; b <- V(g2)$name
  length(intersect(a, b)) / length(union(a, b))
}
jaccard_edges_fn <- function(g1, g2) {
  e1 <- apply(as_edgelist(g1), 1, function(e) paste(sort(e), collapse="-"))
  e2 <- apply(as_edgelist(g2), 1, function(e) paste(sort(e), collapse="-"))
  length(intersect(e1, e2)) / length(union(e1, e2))
}

# ── Main loop ──────────────────────────────────────────────────────────────────
results <- vector("list", length(steps))

for (i in seq_along(steps)) {
  N <- steps[i]
  if (i %% 20 == 0 || i == 1 || i == length(steps))
    cat(sprintf("  Step %d/%d  (N=%d)\n", i, length(steps), N))

  el7 <- a7_full %>% slice(1:min(N, nrow(a7_full)))
  el9 <- a9_full %>% slice(1:min(N, nrow(a9_full)))

  m7 <- compute_metrics(el7, "HPV-A7", N)
  m9 <- compute_metrics(el9, "HPV-A9", N)

  g7 <- graph_from_data_frame(el7, directed=FALSE) %>%
        igraph::simplify(remove.multiple=TRUE, remove.loops=TRUE)
  g9 <- graph_from_data_frame(el9, directed=FALSE) %>%
        igraph::simplify(remove.multiple=TRUE, remove.loops=TRUE)

  jn <- jaccard_nodes_fn(g7, g9)
  je <- jaccard_edges_fn(g7, g9)

  m7$Jaccard_nodos <- jn; m7$Jaccard_enlaces <- je
  m9$Jaccard_nodos <- jn; m9$Jaccard_enlaces <- je

  results[[i]] <- bind_rows(as.data.frame(m7), as.data.frame(m9))
}

df <- bind_rows(results)
saveRDS(df, file.path(DIR_OUT, "5_2_6_topology_evolution_data.rds"))
cat("[3] Data saved. Building plots...\n")

# ── Plot ───────────────────────────────────────────────────────────────────────
pal <- c("HPV-A7" = "#003366", "HPV-A9" = "#CC0000")

metrics_meta <- list(
  list(col="Nodos",         label="N° nodos"),
  list(col="Densidad",      label="Densidad"),
  list(col="Componentes",   label="Componentes"),
  list(col="Grado_medio",   label="Grado medio"),
  list(col="Clustering",    label="Clustering local"),
  list(col="Transitividad", label="Transitividad global"),
  list(col="Asortatividad", label="Asortatividad"),
  list(col="GC_pct",        label="% Componente gigante"),
  list(col="Diametro",      label="Diámetro"),
  list(col="Mean_distance", label="Long. media caminos"),
  list(col="Jaccard_nodos", label="Jaccard nodos"),
  list(col="Jaccard_enlaces",label="Jaccard enlaces")
)

make_panel <- function(col, label) {
  ggplot(df, aes(x = N, y = .data[[col]], color = network, linetype = network)) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = pal) +
    scale_linetype_manual(values = c("HPV-A7"="solid","HPV-A9"="dashed")) +
    scale_x_continuous(labels = comma, breaks = seq(20000, 117407, by = 20000)) +
    labs(title = label, x = NULL, y = NULL, color = NULL, linetype = NULL) +
    theme_bw(base_family = "Times New Roman") +
    theme(
      plot.title       = element_text(size = 28, face = "bold", hjust = 0.5),
      axis.text        = element_text(size = 20),
      axis.text.x      = element_text(angle = 30, hjust = 1),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position  = "none"
    )
}

panels <- lapply(metrics_meta, function(m) make_panel(m$col, m$label))

# Shared legend
legend_plot <- ggplot(df %>% filter(N == 20000),
                      aes(x=N, y=Grado_medio, color=network, linetype=network)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = pal) +
  scale_linetype_manual(values = c("HPV-A7"="solid","HPV-A9"="dashed")) +
  labs(color=NULL, linetype=NULL) +
  theme_bw() +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 30),
        legend.key.width = unit(2.5, "cm"),
        legend.key.size  = unit(1.5, "cm"))

get_legend <- function(p) {
  g <- ggplotGrob(p)
  g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
}
leg <- get_legend(legend_plot)

combined <- wrap_plots(panels, ncol = 3) +
  plot_annotation(
    title = "Evolución de métricas topológicas — CC Co-expression Networks",
    subtitle = "Top-N enlaces ordenados por MI (paso 500), N = 20,000 → 117,407",
    theme = theme(
      plot.title    = element_text(size = 40, face = "bold", hjust = 0.5,
                                   family = "Times New Roman"),
      plot.subtitle = element_text(size = 28, hjust = 0.5, color = "grey30",
                                   family = "Times New Roman")
    )
  )

final <- combined / patchwork::wrap_elements(leg) +
  plot_layout(heights = c(20, 1))

ggsave(file.path(DIR_OUT, "5_2_6_Topology_evolution.png"),
       final, width = 24, height = 28, dpi = 300, bg = "white")

cat("[4] Done. Saved: 5_2_6_Topology_evolution.png\n")
