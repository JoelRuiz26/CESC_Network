# =========================================================
# LCC + Independent Elbows (fast & faithful: igraph coarse + window refine)
# Autor: Joel Ruiz Hernandez
# =========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(scales)
  library(grid)
  library(igraph)
})

# --------------------------------------------------------
# Cutoffs of each network which contains the maximum links (MI-sorted)  that contains all conected nodes
# --------------------------------------------------------
xmax_A7 <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_1_A7_cutoff_links.rds")
xmax_A9 <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_1_A9_cutoff_links.rds")

cat("\n[INFO] xmax per network -> A7:", xmax_A7, "| A9:", xmax_A9, "\n")

# --------------------------------------------------------
# Cargar edgelists (ordenados por MI)
# --------------------------------------------------------
A7_Full <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_A7_edgelist_allnodes.rds")
A9_Full <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_A9_edgelist_allnodes.rds")

imax_A7 <- min(xmax_A7, nrow(A7_Full))
imax_A9 <- min(xmax_A9, nrow(A9_Full))
cat("[INFO] Usable edges -> A7:", imax_A7, "| A9:", imax_A9, "\n")

# ========================================================
# Utilidades
# ========================================================
.as_mat_edges <- function(df, upto) {
  as.matrix(df[seq_len(upto), c("GenA","GenB")])
}

calculate_metrics_coarse_igraph <- function(edgelist, imax, step = 10000L) {
  imax <- min(imax, nrow(edgelist))
  take <- unique(sort(c(seq(from = step, to = imax, by = step), imax)))
  Links <- integer(length(take))
  Total_Genes <- integer(length(take))
  Largest_Component <- integer(length(take))
  for (j in seq_along(take)) {
    i <- take[j]
    E <- .as_mat_edges(edgelist, i)
    g <- graph_from_edgelist(E, directed = FALSE)
    cs <- components(g)$csize
    Links[j] <- i
    Total_Genes[j] <- vcount(g)
    Largest_Component[j] <- max(cs)
  }
  data.frame(Links, Total_Genes, Largest_Component)
}

refine_window_igraph <- function(edgelist, start, end) {
  E <- .as_mat_edges(edgelist, end)
  if (start > 1L) g <- graph_from_edgelist(E[1:(start-1L), , drop = FALSE], directed = FALSE)
  else g <- make_empty_graph(directed = FALSE)
  n <- end - start + 1L
  Links <- integer(n); Total_Genes <- integer(n); Largest_Component <- integer(n)
  ptr <- 0L
  for (i in start:end) {
    e <- E[i, ]
    for (v in e) if (!(v %in% V(g)$name)) g <- add_vertices(g, 1, name = v)
    g <- add_edges(g, as.vector(e))
    cs <- components(g)$csize
    ptr <- ptr + 1L
    Links[ptr] <- i
    Total_Genes[ptr] <- vcount(g)
    Largest_Component[ptr] <- max(cs)
  }
  data.frame(Links, Total_Genes, Largest_Component)
}

# --- NUEVO: elbow con ventana + desempate configurable ---
find_elbow <- function(x, y, window = NULL, tie = c("nearest","middle","first"), guess = NULL, eps = 1e-12) {
  tie <- match.arg(tie)
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  P <- cbind(x_norm, y_norm)
  v <- c(1,1) / sqrt(2)
  dists <- sqrt(rowSums((P - (P %*% v) %*% t(v))^2))
  
  keep <- rep(TRUE, length(x))
  if (!is.null(window)) keep <- (x >= window[1] & x <= window[2])
  
  # si la ventana no retiene nada, cae al dominio completo
  if (!any(keep)) keep <- rep(TRUE, length(x))
  
  m <- max(dists[keep])
  idxs <- which(keep & (abs(dists - m) <= eps))
  
  if (length(idxs) == 1L) return(idxs)
  
  if (tie == "nearest" && !is.null(guess)) {
    return(idxs[ which.min(abs(x[idxs] - guess)) ])
  } else if (tie == "middle") {
    return(idxs[ ceiling(length(idxs)/2) ])
  } else { # "first"
    return(idxs[1])
  }
}

elbow_geom <- function(x, y, idx) {
  A <- c(x[1], y[1]); B <- c(x[length(x)], y[length(y)])
  AB <- B - A; AB2 <- sum(AB^2)
  P  <- c(x[idx], y[idx])
  t  <- ((P[1]-A[1])*AB[1] + (P[2]-A[2])*AB[2]) / AB2
  Pp <- c(A[1] + t*AB[1], A[2] + t*AB[2])
  list(x_start=A[1], y_start=A[2], x_end=B[1], y_end=B[2],
       x_elbow=P[1], y_elbow=P[2], x_proj=Pp[1], y_proj=Pp[2])
}

# ========================================================
# 1) Coarse por red (hasta su propio xmax)
# ========================================================
coarse_step <- 10000L
cat("\n[INFO] Computing coarse curves (step = ", coarse_step, ")...\n", sep = "")
metrics_A7_coarse <- calculate_metrics_coarse_igraph(A7_Full, imax = imax_A7, step = coarse_step)
metrics_A9_coarse <- calculate_metrics_coarse_igraph(A9_Full, imax = imax_A9, step = coarse_step)

# Elbow coarse (aquí puede salir múltiplo del step; es sólo guía)
elbow_idx_A7_c <- find_elbow(metrics_A7_coarse$Links, metrics_A7_coarse$Largest_Component)
elbow_idx_A9_c <- find_elbow(metrics_A9_coarse$Links, metrics_A9_coarse$Largest_Component)
elbow_links_A7_c <- metrics_A7_coarse$Links[elbow_idx_A7_c]
elbow_links_A9_c <- metrics_A9_coarse$Links[elbow_idx_A9_c]
cat("[INFO] Coarse elbows ~", elbow_links_A7_c, "(A7), ~", elbow_links_A9_c, "(A9)\n")

# ========================================================
# 2) Refine independiente por red (±1.5*step alrededor del elbow coarse)
#    y fusión con coarse fuera de la ventana
# ========================================================
refine_window <- as.integer(coarse_step * 1.5)

ref_range_A7 <- c(max(1L, elbow_links_A7_c - refine_window),
                  min(imax_A7, elbow_links_A7_c + refine_window))
ref_range_A9 <- c(max(1L, elbow_links_A9_c - refine_window),
                  min(imax_A9, elbow_links_A9_c + refine_window))

cat("[INFO] Refining A7 in [", ref_range_A7[1], ", ", ref_range_A7[2], "] ...\n", sep = "")
metrics_A7_refine <- refine_window_igraph(A7_Full, start = ref_range_A7[1], end = ref_range_A7[2])

cat("[INFO] Refining A9 in [", ref_range_A9[1], ", ", ref_range_A9[2], "] ...\n", sep = "")
metrics_A9_refine <- refine_window_igraph(A9_Full, start = ref_range_A9[1], end = ref_range_A9[2])

metrics_A7 <- bind_rows(
  filter(metrics_A7_coarse, Links < ref_range_A7[1] | Links > ref_range_A7[2]),
  metrics_A7_refine
) %>% distinct(Links, .keep_all=TRUE) %>% arrange(Links)

metrics_A9 <- bind_rows(
  filter(metrics_A9_coarse, Links < ref_range_A9[1] | Links > ref_range_A9[2]),
  metrics_A9_refine
) %>% distinct(Links, .keep_all=TRUE) %>% arrange(Links)

# --- ANCLAR AL ORIGEN ---
metrics_A7 <- bind_rows(tibble(Links=0L, Total_Genes=0L, Largest_Component=0L), metrics_A7) %>%
  distinct(Links, .keep_all = TRUE) %>% arrange(Links)

metrics_A9 <- bind_rows(tibble(Links=0L, Total_Genes=0L, Largest_Component=0L), metrics_A9) %>%
  distinct(Links, .keep_all = TRUE) %>% arrange(Links)

# ========================================================
# 3) Elbows finales independientes (RESTRINGIDOS a la ventana refinada)
#    + desempate: el más cercano al elbow coarse
# ========================================================
elbow_idx_A7 <- find_elbow(metrics_A7$Links, metrics_A7$Largest_Component,
                           window = ref_range_A7, tie = "nearest", guess = elbow_links_A7_c)
elbow_idx_A9 <- find_elbow(metrics_A9$Links, metrics_A9$Largest_Component,
                           window = ref_range_A9, tie = "nearest", guess = elbow_links_A9_c)

elbow_links_A7 <- metrics_A7$Links[elbow_idx_A7]
elbow_genes_A7 <- metrics_A7$Largest_Component[elbow_idx_A7]

elbow_links_A9 <- metrics_A9$Links[elbow_idx_A9]
elbow_genes_A9 <- metrics_A9$Largest_Component[elbow_idx_A9]

cat("\nIndependent elbows:\n")
cat("  HPV-A7 elbow:", format(elbow_links_A7, big.mark=","), "edges; LCC =", elbow_genes_A7, "\n")
cat("  HPV-A9 elbow:", format(elbow_links_A9, big.mark=","), "edges; LCC =", elbow_genes_A9, "\n")

# (Opcional) Si por alguna razón el elbow cayó fuera (no debería), re-centrar y refinar una vez más:
recenter_refine <- function(edgelist, imax, metrics_coarse, elbow_links, old_range, coarse_step) {
  if (elbow_links >= old_range[1] && elbow_links <= old_range[2]) return(NULL)
  new_range <- c(max(1L, elbow_links - as.integer(coarse_step*1.5)),
                 min(imax, elbow_links + as.integer(coarse_step*1.5)))
  cat("[INFO] Re-centering refine to [", new_range[1], ", ", new_range[2], "] ...\n", sep = "")
  met_ref <- refine_window_igraph(edgelist, start = new_range[1], end = new_range[2])
  met <- bind_rows(filter(metrics_coarse, Links < new_range[1] | Links > new_range[2]),
                   met_ref) %>%
    distinct(Links, .keep_all=TRUE) %>% arrange(Links)
  list(metrics = met, window = new_range)
}

recA7 <- recenter_refine(A7_Full, imax_A7, metrics_A7_coarse, elbow_links_A7, ref_range_A7, coarse_step)
if (!is.null(recA7)) {
  metrics_A7 <- recA7$metrics
  ref_range_A7 <- recA7$window
  elbow_idx_A7 <- find_elbow(metrics_A7$Links, metrics_A7$Largest_Component,
                             window = ref_range_A7, tie = "nearest", guess = elbow_links_A7)
  elbow_links_A7 <- metrics_A7$Links[elbow_idx_A7]
  elbow_genes_A7 <- metrics_A7$Largest_Component[elbow_idx_A7]
}

recA9 <- recenter_refine(A9_Full, imax_A9, metrics_A9_coarse, elbow_links_A9, ref_range_A9, coarse_step)
if (!is.null(recA9)) {
  metrics_A9 <- recA9$metrics
  ref_range_A9 <- recA9$window
  elbow_idx_A9 <- find_elbow(metrics_A9$Links, metrics_A9$Largest_Component,
                             window = ref_range_A9, tie = "nearest", guess = elbow_links_A9)
  elbow_links_A9 <- metrics_A9$Links[elbow_idx_A9]
  elbow_genes_A9 <- metrics_A9$Largest_Component[elbow_idx_A9]
}



# ========================================================
# 4) Figura (revista, desde el origen, sin diagonal)
# ========================================================

# 1) Asegura que las curvas arrancan en (0,0)
ensure_origin <- function(df) {
  if (!any(df$Links == 0L)) {
    df <- bind_rows(tibble(Links=0L, Total_Genes=0L, Largest_Component=0L), df)
  } else {
    df$Largest_Component[df$Links == 0L] <- 0L
    df$Total_Genes[df$Links == 0L]       <- 0L
  }
  arrange(df, Links)
}
metrics_A7p <- ensure_origin(metrics_A7)
metrics_A9p <- ensure_origin(metrics_A9)

# 2) Datos de codos
fmt_int <- function(z) formatC(as.integer(z), format = "d", big.mark = ",")
elbows <- tibble(
  Network = c("HPV-A7","HPV-A9"),
  x = c(elbow_links_A7, elbow_links_A9),
  y = c(elbow_genes_A7,  elbow_genes_A9),
  label = c(
    paste0("HPV-A7 elbow\nx=", fmt_int(elbow_links_A7)),
    paste0("HPV-A9 elbow\nx=", fmt_int(elbow_links_A9))
  )
)

# 3) Límites del panel (curvas intactas)
x_end <- 2e6                  # cámbialo a 3e6 si lo prefieres
x_lim <- c(0, x_end)
y_max <- max(metrics_A7p$Largest_Component, metrics_A9p$Largest_Component, na.rm = TRUE)
y_lim <- c(0, y_max * 1.03)   # holgura mínima arriba

# 4) Posiciones de ETIQUETAS: franja inferior derecha, APILADAS sin solape
dy <- diff(y_lim)
y_band_top <- min(elbows$y) - 0.02 * dy       # justo debajo del codo más bajo
step       <- 0.07 * dy                        # separación vertical entre etiquetas

elbows_pos <- elbows %>%
  arrange(desc(y)) %>%                         # orden: primero la más alta
  mutate(
    x_lab = pmin(x + 0.12 * diff(x_lim), x_end * 0.94),  # a la derecha del codo
    y_lab = y_band_top - (row_number() - 1) * step        # apilado vertical
  )

# 5) Paleta
pal <- c("HPV-A7" = "blue", "HPV-A9" = "red")

# 6) Plot
p <- ggplot() +
  geom_line(data = metrics_A7p, aes(Links, Largest_Component, color = "HPV-A7"), linewidth = 1.3) +
  geom_line(data = metrics_A9p, aes(Links, Largest_Component, color = "HPV-A9"), linewidth = 1.3) +
  
  geom_vline(data = elbows, aes(xintercept = x, color = Network),
             linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  geom_point(data = elbows, aes(x, y, fill = Network),
             shape = 21, size = 3, stroke = 0.9, color = "white", show.legend = FALSE) +
  
  # Conectores y etiquetas (borde negro, fondo transparente)
  geom_segment(data = elbows_pos,
               aes(x = x, y = y, xend = x_lab, yend = y_lab),
               linetype = "dotted", linewidth = 0.45, color = "grey40") +
  geom_label(data = elbows_pos,
             aes(x = x_lab, y = y_lab, label = label),
             fill = NA,            # transparente: no tapa líneas
             label.colour = "black",
             color = "black",
             label.size = 0.3,
             size = 3.8,
             label.padding = unit(3, "pt"),
             show.legend = FALSE) +
  
  scale_color_manual(values = pal, name = "Network") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(limits = x_lim, expand = c(0, 0), labels = label_comma()) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0), labels = label_comma()) +
  
  labs(
    title = "Growth of the Largest Connected Component in Cervical Cancer Co-expression Networks",
    x     = "Number of Edges (MI-sorted)",
    y     = "Number of Genes in Largest Connected Component"
  ) +
  
  theme_bw() +
  theme(
    plot.title.position = "plot",
    plot.title  = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(t = 12, b = 10)),
    axis.title  = element_text(size = 16, face = "bold"),
    axis.text   = element_text(size = 13),
    legend.position = "top",
    legend.title    = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 13),
    legend.box.margin = margin(t = 6, b = 8),
    panel.grid.major = element_line(linewidth = 0.35, linetype = "dotted", color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 22, r = 18, b = 14, l = 16)
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 2))) +
  coord_cartesian(clip = "off")

#print(p)
ggsave("GC_elbows_labels_stacked.png", p, width = 14, height = 6.5, dpi = 600)

# --------------------------------------------------------
# 5) (Opcional) Exportar subredes al elbow independiente
# --------------------------------------------------------
save.image("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_2_GiantComponent.RData")

ggsave("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_3_GC_graph_elbow.png",
        plot = p, width = 12, height = 8, dpi = 600)

Subnetwork_A7_elbow <- A7_Full %>% slice(1:elbow_links_A7)
Subnetwork_A9_elbow <- A9_Full %>% slice(1:elbow_links_A7)

graph_A7_elbow <- graph_from_data_frame(d = Subnetwork_A7_elbow, directed = FALSE)
graph_A9_elbow <- graph_from_data_frame(d = Subnetwork_A9_elbow, directed = FALSE)

saveRDS(graph_A7_elbow, "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_elbow.rds") 
saveRDS(graph_A9_elbow, "~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_elbow.rds")

cat("Graphs exported successfully!\n")

#load("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_2_GiantComponent.RData")

