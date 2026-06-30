################################################################################
# Delta-degree z-score — Differential Degree Connectivity
# HPV-A7 vs HPV-A9 Co-expression Networks
# Author: Joel Ruiz Hernandez
################################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(ggrepel)
})

################################################################################
# 1) Load graphs
################################################################################

A7_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A7_max.rds")
A9_HPV <- readRDS("~/CESC_Network/5_Network_analisis/5_1_Subnetwork/5_2_4_Subnetwork_A9_max.rds")

################################################################################
# 2) Shared nodes
################################################################################

common_nodes <- intersect(V(A7_HPV)$name, V(A9_HPV)$name)

deg_A7_common <- as.numeric(igraph::degree(A7_HPV, v = common_nodes))
deg_A9_common <- as.numeric(igraph::degree(A9_HPV, v = common_nodes))

degree_df <- data.frame(
  Gene   = common_nodes,
  deg_A7 = deg_A7_common,
  deg_A9 = deg_A9_common,
  stringsAsFactors = FALSE
)

degree_df <- degree_df[degree_df$deg_A7 > 0 & degree_df$deg_A9 > 0, ]

################################################################################
# 3) Delta-degree and z-score
################################################################################

degree_df$delta <- degree_df$deg_A7 - degree_df$deg_A9
degree_df$z     <- as.numeric(scale(degree_df$delta))

################################################################################
# 4) Color category
################################################################################

z_thr <- 1.96

degree_df$color_group <- ifelse(
  degree_df$z >  z_thr, "A7 higher",
  ifelse(degree_df$z < -z_thr, "A9 higher", "Similar")
)

degree_df$color_group <- factor(
  degree_df$color_group,
  levels = c("A7 higher", "Similar", "A9 higher")
)

################################################################################
# 5) Count per group
################################################################################

group_counts <- as.data.frame(table(degree_df$color_group))
colnames(group_counts) <- c("color_group", "n")

legend_labels <- setNames(
  sprintf("%s  (n = %d)", group_counts$color_group, group_counts$n),
  group_counts$color_group
)

################################################################################
# 6) Rank-order x-axis
################################################################################

degree_df <- degree_df[order(degree_df$z), ]
degree_df$rank <- seq_len(nrow(degree_df))

################################################################################
# 7) Top hub genes in each tail
################################################################################

top_genes <- rbind(
  head(degree_df[order(-degree_df$z), ], 8),
  head(degree_df[order( degree_df$z), ], 8)
)

################################################################################
# 8) Plot
################################################################################

pal <- c("A7 higher" = "#003366",
         "Similar"   = "grey70",
         "A9 higher" = "#CC0000")

p_delta <- ggplot(degree_df, aes(x = rank, y = z, color = color_group)) +
  
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin =  z_thr, ymax = Inf,
           fill = "#003366", alpha = 0.06) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = -z_thr,
           fill = "#CC0000", alpha = 0.06) +
  
  geom_hline(yintercept = 0,      linetype = "solid",  color = "grey40",
             linewidth = 0.5) +
  geom_hline(yintercept =  z_thr, linetype = "dashed", color = "#003366",
             linewidth = 0.45) +
  geom_hline(yintercept = -z_thr, linetype = "dashed", color = "#CC0000",
             linewidth = 0.45) +
  
  geom_point(
    aes(size  = ifelse(color_group == "Similar", 0.6, 1.1),
        alpha = ifelse(color_group == "Similar", 0.25, 0.80))
  ) +
  scale_size_identity() +
  scale_alpha_identity() +
  
  geom_text_repel(
    data               = top_genes,
    aes(label          = Gene),
    size               = 2.8,
    family             = "Times New Roman",
    fontface           = "italic",
    max.overlaps       = 20,
    min.segment.length = 0.2,
    segment.size       = 0.3,
    segment.color      = "grey50",
    box.padding        = 0.4,
    show.legend        = FALSE
  ) +
  
  scale_color_manual(values = pal, labels = legend_labels) +
  
  labs(
    title    = "Differential Degree Connectivity",
    subtitle = sprintf(
      "Z-score of \u0394degree (A7\u2212A9) per shared gene  |  threshold |z| = %.2f",
      z_thr),
    x     = "Gene rank (ordered by z-score)",
    y     = "Z-score of \u0394degree",
    color = NULL
  ) +
  
  theme_classic(base_family = "Times New Roman") +
  theme(
    plot.title         = element_text(size = 16, face = "bold",    hjust = 0.5),
    plot.subtitle      = element_text(size = 10, color = "grey35", hjust = 0.5),
    axis.title         = element_text(size = 12, face = "bold"),
    axis.text          = element_text(size = 10, color = "black"),
    axis.line          = element_line(linewidth = 0.5),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = c(0.5, 0.97),
    legend.justification = c("center", "top"),
    legend.direction   = "horizontal",
    legend.background  = element_rect(fill  = alpha("white", 0.85),
                                      color = "grey70", linewidth = 0.3),
    legend.key         = element_blank(),
    legend.text        = element_text(size = 9.5)
  )

p_delta

ggsave(
  "~/CESC_Network/5_Network_analisis/5_2_Network_topology/5_2_5_Delta_Degree_Zscore.png",
  p_delta,
  width = 10, height = 6, dpi = 600, bg = "white"
)

message("Saved: Delta-degree z-score")