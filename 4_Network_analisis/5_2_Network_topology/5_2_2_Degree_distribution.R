################################################################################
# Degree Distribution Analysis for Cervical Cancer Co-expression Networks
# Author: Joel Ruiz Hernandez
# Description:
#   This script generates two  plots:
#     (1) Degree distribution
#     (2) Complementary cumulative degree distribution (CCDF)
#   Both graphs highlight the heavy-tailed, heterogeneous nature typical of
#   biological networks, showing hubs with high degree.
################################################################################

library(igraph)
library(ggplot2)
library(scales)

################################################################################
# Load Graphs
################################################################################

A7_HPV <- readRDS('~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_elbow.rds')
A9_HPV <- readRDS('~/CESC_Network/4_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_elbow.rds')

################################################################################
# Degree Distribution (PDF)
################################################################################
#   The degree distribution shows the probability p(k) that a gene node
#   has exactly k connections. Biological co-expression networks typically
#   display heavy-tailed distributions where most nodes have low degree
#   but a few "hub" genes are highly connected.

get_degree_dist_df <- function(graph, group_name) {
  stopifnot(inherits(graph, "igraph"))
  deg <- igraph::degree(graph, mode = "all")     # << clave
  deg <- as.numeric(deg)
  deg_dist <- table(deg) / length(deg)
  data.frame(
    Degree = as.numeric(names(deg_dist)),
    Probability = as.numeric(deg_dist),
    Network = group_name,
    row.names = NULL, check.names = FALSE
  )
}

df_A7 <- get_degree_dist_df(A7_HPV, "HPV_A7")
df_A9 <- get_degree_dist_df(A9_HPV, "HPV_A9")


df_A7 <- get_degree_dist_df(A7_HPV, "HPV_A7")
df_A9 <- get_degree_dist_df(A9_HPV, "HPV_A9")
df_combined <- rbind(df_A7, df_A9)

# Professional-style plot
p_pdf <- ggplot(df_combined, aes(x = Degree, y = Probability, 
                                 color = Network, shape = Network, linetype = Network)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3, size = 1) +
  scale_color_manual(values = c(HPV_A7 = "#003366", HPV_A9 = "#FF0000")) +
  scale_shape_manual(values = c(HPV_A7 = 16, HPV_A9 = 17)) +
  scale_linetype_manual(values = c(HPV_A7 = "solid", HPV_A9 = "dashed")) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bl") +
  labs(
    title = "Degree Distribution in CC Co-expression Networks",
    x = "Degree (k)",
    y = "Probability p(k)",
    color = "Network",
    shape = "Network",
    linetype = "Network"
  ) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.title         = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title         = element_text(size = 12, face = "bold"),
    axis.text          = element_text(size = 10),
    panel.border       = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major   = element_line(color = "grey90"),
    panel.grid.minor   = element_line(color = "grey95"),
    axis.line          = element_blank(),
    legend.position    = c(0.95, 0.95),         # top-right inside
    legend.justification = c("right", "top"),
    legend.direction   = "vertical",
    legend.background  = element_rect(fill = alpha("white", 0.8), color = "black"),
    legend.key         = element_blank(),
    legend.title       = element_text(size = 11, face = "bold"),
    legend.text        = element_text(size = 10)
  )

ggsave("~/CESC_Network/4_Network_analisis/5_2_Network_topology/5_2_2_Degree_distribution.png",
       p_pdf, width = 7, height = 6, dpi = 300)
print(p_pdf)

################################################################################
#  Complementary Cumulative Distribution (CCDF)
################################################################################
#   P(K ≥ k) is the complementary cumulative distribution function (CCDF).
#   It indicates the probability that a gene has at least k connections.
#   This smooths the noise in the tail and emphasizes the presence of hubs.
#   Heavy-tailed CCDF curves in log-log suggest heterogeneous network
#   organization, typical in biological systems.
# ========= CCDF desde el eje Y (k=0), SOLO PUNTOS =========
library(ggplot2)
library(dplyr)
library(scales)
library(igraph)

# Escala log10(1+x), permite k = 0
log10_1p_trans <- function() {
  scales::trans_new(
    name = "log10_1p",
    transform = function(x) log10(1 + x),
    inverse   = function(x) 10^x - 1
  )
}
library(igraph)
library(dplyr)

build_ccdf_df <- function(graph, label) {
  stopifnot(inherits(graph, "igraph"))
  # Usa explícitamente igraph::degree y garantiza vector numérico
  deg <- as.integer(igraph::degree(graph, mode = "all"))
  deg <- deg[!is.na(deg)]
  if (length(deg) == 0L) {
    return(tibble::tibble(Degree = integer(0), CCDF = numeric(0), Network = label))
  }
  kmax <- max(deg)
  h    <- tabulate(deg + 1L, nbins = kmax + 1L)  # frecuencias 0..kmax
  ccdf <- rev(cumsum(rev(h))) / length(deg)      # P(K >= k)
  tibble::tibble(Degree = 0:kmax, CCDF = ccdf, Network = label)
}

# (opcional pero recomendable) Evita conflictos definitivamente:
# if (exists("degree", inherits = FALSE)) rm(degree)
# library(conflicted); conflict_prefer("degree", "igraph")

df_ccdf <- dplyr::bind_rows(
  build_ccdf_df(A7_HPV, "HPV_A7"),
  build_ccdf_df(A9_HPV, "HPV_A9")
)


x_max <- max(df_ccdf$Degree)
y_min <- min(df_ccdf$CCDF[df_ccdf$CCDF > 0])

p_ccdf <- ggplot(df_ccdf, aes(Degree, CCDF, color = Network)) +
  # SOLO puntos (sin líneas que unan)
  geom_point(size = 1.7, alpha = 0.7, na.rm = TRUE) +
  
  scale_color_manual(values = c(HPV_A7 = "#003366", HPV_A9 = "#FF0000")) +
  
  # X en log10(1+x) para incluir 0; Y en log10 clásico
  scale_x_continuous(
    trans  = log10_1p_trans(),
    limits = c(0, x_max),
    expand = c(0, 0),
    breaks = { br <- c(0, 1, 10^(0:floor(log10(max(1, x_max))))); br[br <= x_max] },
    labels = label_comma()
  ) +
  scale_y_log10(
    limits = c(y_min, 1),
    expand = c(0, 0),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  annotation_logticks(sides = "l") +
  
  labs(
    title = "CCDF of Degree Distribution in Co-expression Networks",
    x = "Degree (k)",
    y = expression(P(K >= k)),
    color = "Network"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10),
    legend.position     = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background   = element_rect(fill = alpha("white", 0.8), color = "black")
  )

print(p_ccdf)


ggsave("~/CESC_Network/4_Network_analisis/5_2_Network_topology/5_1_2_Degree_CCDF_distribution.png",
       p_ccdf, width = 7, height = 6, dpi = 300)
