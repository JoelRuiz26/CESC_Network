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
A7_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_2_Subnetwork_A7_600k.rds")
A9_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_1_2_2_Subnetwork_A9_600k.rds")

################################################################################
# Degree Distribution (PDF)
################################################################################
#   The degree distribution shows the probability p(k) that a gene node
#   has exactly k connections. Biological co-expression networks typically
#   display heavy-tailed distributions where most nodes have low degree
#   but a few "hub" genes are highly connected.

get_degree_dist_df <- function(graph, group_name) {
  deg <- degree(graph, mode = "all")
  deg_dist <- table(deg) / length(deg)
  data.frame(Degree = as.numeric(names(deg_dist)), 
             Probability = as.numeric(deg_dist), 
             Network = group_name)
}

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

ggsave("~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_2_2_Degree_distribution.png",
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

build_ccdf_df <- function(graph, label) {
  deg <- degree(graph)
  xs <- sort(unique(deg))
  ys <- sapply(xs, function(k) sum(deg >= k) / length(deg))
  data.frame(Degree = xs, CCDF = ys, Network = label)
}

df_ccdf <- rbind(
  build_ccdf_df(A7_HPV, "HPV_A7"),
  build_ccdf_df(A9_HPV, "HPV_A9")
)

p_ccdf <- ggplot(df_ccdf, aes(x = Degree, y = CCDF, 
                              color = Network, shape = Network, linetype = Network)) +
  geom_point(size = 2.5, alpha = 0.7) +
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
    title = "CCDF of Degree Distribution in Co‑expression Networks",
    x = "Degree (k)", 
    y = expression(P(K >= k)),
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
    legend.position    = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.direction   = "vertical",
    legend.background  = element_rect(fill = alpha("white", 0.8), color = "black"),
    legend.key         = element_blank(),
    legend.title       = element_text(size = 11, face = "bold"),
    legend.text        = element_text(size = 10)
  )

ggsave("~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_1_2_Degree_CCDF_distribution.png",
       p_ccdf, width = 7, height = 6, dpi = 300)
print(p_ccdf)
