################################################################################
# Network Stability Analysis
#
# Objective: Demonstrate that the topological differences between HPV-A7 and
# HPV-A9 co-expression networks persist across all edge cutoffs (N = 20,000
# to 117,407, step 500), and justify the choice of N = 117,184 as the cutoff
# that maximises the number of shared nodes between both networks.
#
# Output (all saved to this directory):
#   5_2_1_A_cutoff_justification.png  — shared nodes vs N; highlights N=117,184
#   5_2_1_B_topology_evolution.png    — all topology metrics vs N (both nets)
#   5_2_1_C_metric_delta.png          — per-metric A9−A7 difference vs N
#   5_2_1_stability_data.rds          — raw per-step data frame
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

DIR_NET  <- "/STORAGE/csbig/jruiz/cesc_network/cervix_ARACNE"
DIR_OUT  <- "~/CESC_Network/5_Network_analisis/5_2_1_Net_stability"
CUTOFF   <- 117184L    # chosen cutoff (maximises shared nodes within A7 range)

PALETTE  <- c("HPV-A7" = "#003366", "HPV-A9" = "#CC0000")
STEPS    <- seq(20000L, 117407L, by = 500L)
BASE_SZ  <- 18          # base font size for panels

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Load complete sorted networks
# ─────────────────────────────────────────────────────────────────────────────
cat("[1] Loading networks...\n")
read_net <- function(f)
  vroom(f, col_names = c("GenA","GenB","MI","pvalue"),
        skip = 1, show_col_types = FALSE) %>% select(GenA, GenB, MI)

a7 <- read_net(file.path(DIR_NET, "A7_network_sort.txt"))
a9 <- read_net(file.path(DIR_NET, "A9_network_sort.txt"))
cat(sprintf("  A7: %d edges | A9: %d edges\n", nrow(a7), nrow(a9)))

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Compute per-step metrics
# ─────────────────────────────────────────────────────────────────────────────
cat(sprintf("[2] Computing metrics at %d steps...\n", length(STEPS)))

build_graph <- function(el, n)
  graph_from_data_frame(slice(el, 1:min(n, nrow(el))), directed = FALSE) %>%
  igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE)

fast_metrics <- function(g, name, N) {
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
    GC_pct        = max(comp$csize) / vcount(g) * 100
  )
}

jaccard_nodes_fn <- function(g1, g2) {
  a <- V(g1)$name; b <- V(g2)$name
  length(intersect(a, b)) / length(union(a, b))
}

# Jaccard on edges: use sorted integer-pair keys (fast, avoids R-level string apply)
jaccard_edges_fn <- function(g1, g2) {
  edge_key <- function(g) {
    el  <- as_edgelist(g, names = FALSE)
    paste(pmin(el[,1], el[,2]), pmax(el[,1], el[,2]), sep = "-")
  }
  # node names may differ between graphs; work with canonical name strings
  edge_key_name <- function(g) {
    el <- as_edgelist(g, names = TRUE)
    ifelse(el[,1] < el[,2],
           paste(el[,1], el[,2], sep = "\t"),
           paste(el[,2], el[,1], sep = "\t"))
  }
  e1 <- edge_key_name(g1)
  e2 <- edge_key_name(g2)
  length(intersect(e1, e2)) / length(union(e1, e2))
}

results <- vector("list", length(STEPS))

for (i in seq_along(STEPS)) {
  N <- STEPS[i]
  if (i %% 20 == 0 || i == 1 || i == length(STEPS))
    cat(sprintf("  %d/%d  N=%d\n", i, length(STEPS), N))

  g7 <- build_graph(a7, N)
  g9 <- build_graph(a9, N)

  m7 <- fast_metrics(g7, "HPV-A7", N)
  m9 <- fast_metrics(g9, "HPV-A9", N)

  jn <- jaccard_nodes_fn(g7, g9)
  je <- jaccard_edges_fn(g7, g9)

  m7$Jaccard_nodos   <- jn;  m7$Jaccard_enlaces <- je
  m9$Jaccard_nodos   <- jn;  m9$Jaccard_enlaces <- je

  results[[i]] <- bind_rows(as.data.frame(m7), as.data.frame(m9))
}

df <- bind_rows(results)
saveRDS(df, file.path(DIR_OUT, "5_2_1_stability_data.rds"))
cat("[3] Data saved.\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Plot A — Cutoff justification: shared nodes vs N
# ─────────────────────────────────────────────────────────────────────────────
cat("[4] Plot A: cutoff justification...\n")

shared_df <- df %>%
  filter(network == "HPV-A7") %>%
  select(N, Jaccard_nodos, Nodos) %>%
  left_join(
    df %>% filter(network == "HPV-A9") %>% select(N, Nodos_A9 = Nodos),
    by = "N"
  ) %>%
  mutate(Shared_nodes = round(Jaccard_nodos * (Nodos + Nodos_A9) / (1 + Jaccard_nodos)))

# approximate shared nodes directly from Jaccard definition:
# J = |I| / |U| = |I| / (|A7| + |A9| - |I|)  => |I| = J*(|A7|+|A9|)/(1+J)
# Use node-level counts from df
node_df <- df %>%
  select(N, network, Nodos) %>%
  pivot_wider(names_from = network, values_from = Nodos) %>%
  rename(A7 = `HPV-A7`, A9 = `HPV-A9`) %>%
  left_join(df %>% filter(network == "HPV-A7") %>% select(N, Jaccard_nodos), by = "N") %>%
  mutate(
    Shared = round(Jaccard_nodos * (A7 + A9) / (1 + Jaccard_nodos)),
    pct_A7 = Shared / A7 * 100
  )

# N=117,184 is not in the 500-step sequence; use the known value from 5_1_1 analysis
SHARED_AT_CUTOFF <- 17436L

p_just <- ggplot(node_df, aes(x = N)) +
  geom_line(aes(y = Shared), color = "#444444", linewidth = 1.2) +
  geom_vline(xintercept = CUTOFF, linetype = "dashed", color = "#CC0000", linewidth = 1.1) +
  annotate("point", x = CUTOFF, y = SHARED_AT_CUTOFF,
           color = "#CC0000", size = 5, shape = 21, fill = "white", stroke = 2) +
  annotate("text",
           x = CUTOFF - 2000, y = 7500,
           label = sprintf("N = %s\n%s nodos compartidos",
                           format(CUTOFF, big.mark=","),
                           format(SHARED_AT_CUTOFF, big.mark=",")),
           hjust = 1, color = "#CC0000", size = BASE_SZ * 0.6,
           fontface = "bold") +
  scale_x_continuous(labels = comma, breaks = seq(20000, 117407, by = 20000)) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Nodos compartidos entre HPV-A7 y HPV-A9 en función del corte",
    subtitle = "Se elige N = 117,184 porque maximiza los nodos compartidos dentro del rango válido de A7",
    x        = "N (top-N enlaces por red, ordenados por MI)",
    y        = "Nodos en la intersección"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(size = BASE_SZ * 2,   face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = BASE_SZ * 1.3, hjust = 0.5, color = "grey30"),
    axis.title    = element_text(size = BASE_SZ * 1.4, face = "bold"),
    axis.text     = element_text(size = BASE_SZ * 1.2),
    axis.text.x   = element_text(angle = 30, hjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

ggsave(file.path(DIR_OUT, "5_2_1_A_cutoff_justification.png"),
       p_just, width = 16, height = 8, dpi = 300, bg = "white")

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Plot B — Topology metric evolution (multi-panel)
# ─────────────────────────────────────────────────────────────────────────────
cat("[5] Plot B: topology evolution...\n")

metrics_meta <- list(
  list(col = "Nodos",         label = "N° de nodos"),
  list(col = "Densidad",      label = "Densidad"),
  list(col = "Componentes",   label = "N° de componentes"),
  list(col = "Grado_medio",   label = "Grado medio"),
  list(col = "Clustering",    label = "Clustering local\n(transitividad promedio)"),
  list(col = "Transitividad", label = "Transitividad global"),
  list(col = "Asortatividad", label = "Asortatividad"),
  list(col = "GC_pct",        label = "% Componente gigante"),
  list(col = "Jaccard_nodos", label = "Jaccard — nodos"),
  list(col = "Jaccard_enlaces",label = "Jaccard — enlaces")
)

make_panel <- function(col, label) {
  ggplot(df, aes(x = N, y = .data[[col]],
                 color = network, linetype = network)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = CUTOFF, linetype = "dotted",
               color = "grey40", linewidth = 0.7) +
    scale_color_manual(values = PALETTE) +
    scale_linetype_manual(values = c("HPV-A7" = "solid", "HPV-A9" = "dashed")) +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       breaks = seq(20000, 117407, by = 30000)) +
    labs(title  = label, x = NULL, y = NULL,
         color  = NULL, linetype = NULL) +
    theme_bw() +
    theme(
      plot.title       = element_text(size = BASE_SZ * 1.5, face = "bold", hjust = 0.5),
      axis.text        = element_text(size = BASE_SZ),
      axis.text.x      = element_text(angle = 30, hjust = 1, size = BASE_SZ * 0.85),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
      legend.position  = "none"
    )
}

panels <- lapply(metrics_meta, function(m) make_panel(m$col, m$label))

# shared legend
leg_data <- df %>% filter(N == min(STEPS)) %>%
  select(network, Nodos) %>%
  mutate(x = 1)

leg_p <- ggplot(leg_data, aes(x=x, y=Nodos, color=network, linetype=network)) +
  geom_line(linewidth = 1.8) +
  scale_color_manual(values = PALETTE) +
  scale_linetype_manual(values = c("HPV-A7"="solid","HPV-A9"="dashed")) +
  labs(color=NULL, linetype=NULL) +
  theme_void() +
  theme(
    legend.position   = "bottom",
    legend.text       = element_text(size = BASE_SZ * 2),
    legend.key.width  = unit(3, "cm"),
    legend.key.height = unit(1.2, "cm"),
    legend.spacing.x  = unit(1, "cm")
  )

get_legend <- function(p) {
  gt <- ggplotGrob(p)
  gt$grobs[[which(sapply(gt$grobs, function(x) x$name) == "guide-box")]]
}
leg <- tryCatch(get_legend(leg_p),
                error = function(e) patchwork::guide_area())

combined <- wrap_plots(panels, ncol = 2) +
  plot_annotation(
    title    = "Estabilidad topológica — HPV-A7 vs HPV-A9",
    subtitle = sprintf("Top-N enlaces por red (paso 500) · línea punteada = corte elegido (N = %s)",
                       format(CUTOFF, big.mark=",")),
    theme = theme(
      plot.title    = element_text(size = BASE_SZ * 2.5, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = BASE_SZ * 1.5, hjust = 0.5, color = "grey30")
    )
  )

final_B <- combined /
  patchwork::wrap_elements(leg) +
  plot_layout(heights = c(25, 1))

ggsave(file.path(DIR_OUT, "5_2_1_B_topology_evolution.png"),
       final_B, width = 22, height = 40, dpi = 300, bg = "white")

# ─────────────────────────────────────────────────────────────────────────────
# 5.  Plot C — Metric delta (A9 − A7) vs N
# ─────────────────────────────────────────────────────────────────────────────
cat("[6] Plot C: metric deltas...\n")

delta_df <- df %>%
  select(N, network, Nodos, Densidad, Componentes, Grado_medio,
         Clustering, Transitividad, Asortatividad, GC_pct) %>%
  pivot_wider(names_from = network,
              values_from = c(Nodos, Densidad, Componentes, Grado_medio,
                              Clustering, Transitividad, Asortatividad, GC_pct)) %>%
  mutate(
    Delta_Nodos         = `Nodos_HPV-A9`         - `Nodos_HPV-A7`,
    Delta_Densidad      = `Densidad_HPV-A9`      - `Densidad_HPV-A7`,
    Delta_Componentes   = `Componentes_HPV-A9`   - `Componentes_HPV-A7`,
    Delta_Grado         = `Grado_medio_HPV-A9`   - `Grado_medio_HPV-A7`,
    Delta_Clustering    = `Clustering_HPV-A9`    - `Clustering_HPV-A7`,
    Delta_Transitividad = `Transitividad_HPV-A9` - `Transitividad_HPV-A7`,
    Delta_Asortatividad = `Asortatividad_HPV-A9` - `Asortatividad_HPV-A7`,
    Delta_GC_pct        = `GC_pct_HPV-A9`        - `GC_pct_HPV-A7`
  ) %>%
  select(N, starts_with("Delta_"))

delta_long <- delta_df %>%
  pivot_longer(-N, names_to = "Metrica", values_to = "Delta") %>%
  mutate(Metrica = recode(Metrica,
    Delta_Nodos         = "Nodos",
    Delta_Densidad      = "Densidad",
    Delta_Componentes   = "Componentes",
    Delta_Grado         = "Grado medio",
    Delta_Clustering    = "Clustering local",
    Delta_Transitividad = "Transitividad global",
    Delta_Asortatividad = "Asortatividad",
    Delta_GC_pct        = "% Comp. gigante"
  ))

make_delta_panel <- function(metric) {
  sub_df  <- delta_long %>% filter(Metrica == metric)
  yrange  <- range(sub_df$Delta, na.rm = TRUE)
  pos_neg <- if (all(sub_df$Delta >= 0)) "#CC0000" else if (all(sub_df$Delta <= 0)) "#003366" else "black"

  ggplot(sub_df, aes(x = N, y = Delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_line(color = "#555555", linewidth = 1) +
    geom_vline(xintercept = CUTOFF, linetype = "dotted",
               color = "#CC0000", linewidth = 0.8) +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       breaks = seq(20000, 117407, by = 30000)) +
    labs(title = metric,
         x = NULL,
         y = "A9 − A7") +
    theme_bw() +
    theme(
      plot.title       = element_text(size = BASE_SZ * 1.5, face = "bold", hjust = 0.5),
      axis.title.y     = element_text(size = BASE_SZ * 1.1),
      axis.text        = element_text(size = BASE_SZ),
      axis.text.x      = element_text(angle = 30, hjust = 1, size = BASE_SZ * 0.85),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7)
    )
}

delta_panels <- lapply(unique(delta_long$Metrica), make_delta_panel)

combined_C <- wrap_plots(delta_panels, ncol = 2) +
  plot_annotation(
    title    = "Diferencia de métricas topológicas  (HPV-A9 − HPV-A7)",
    subtitle = sprintf("Diferencias positivas = A9 > A7 · línea punteada roja = corte elegido (N = %s)",
                       format(CUTOFF, big.mark = ",")),
    theme = theme(
      plot.title    = element_text(size = BASE_SZ * 2.5, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = BASE_SZ * 1.5, hjust = 0.5, color = "grey30")
    )
  )

ggsave(file.path(DIR_OUT, "5_2_1_C_metric_delta.png"),
       combined_C, width = 22, height = 28, dpi = 300, bg = "white")

cat("[7] All done. Outputs in", DIR_OUT, "\n")
cat("    5_2_1_A_cutoff_justification.png\n")
cat("    5_2_1_B_topology_evolution.png\n")
cat("    5_2_1_C_metric_delta.png\n")
cat("    5_2_1_stability_data.rds\n")
