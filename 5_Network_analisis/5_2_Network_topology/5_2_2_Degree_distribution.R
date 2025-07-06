library(igraph)
library(ggplot2)
library(scales)

# Load graphs
A7_HPV <- read_graph("~/5_Global_Analisis/1_Subnetworks/A7_subgraph_600mil.graphml", format = "graphml")
A9_HPV <- read_graph("~/5_Global_Analisis/1_Subnetworks/A9_subgraph_600mil.graphml", format = "graphml")

# Degree distribution function
get_degree_dist_df <- function(graph, group_name) {
  deg <- degree(graph, mode = "all")
  deg_dist <- table(deg)/length(deg)
  data.frame(k = as.numeric(names(deg_dist)), pk = as.numeric(deg_dist), Red = group_name)
}

# Get data
df_A7 <- get_degree_dist_df(A7_HPV, "A7")
df_A9 <- get_degree_dist_df(A9_HPV, "A9")
df_combined <- rbind(df_A7, df_A9)

# Calculate x-axis limit
x_upper_limit <- quantile(df_combined$k, 0.95)

# Plot with adjusted axes
ggplot(df_combined, aes(x = k, y = pk, color = Red)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(NA, x_upper_limit)
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  annotation_logticks() +
  labs(
    title = "Distribución de grado en redes de coexpresión de CC",
    x = "Grado (k)",
    y = "Probabilidad p(k)",
    color = "Red"  # Changed from "Grupo" to "Red"
  ) +
  scale_color_manual(values = c("A7" = "#003366", "A9" = "#FF0000")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 20, margin = margin(r = 15)),
    axis.title.x = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    aspect.ratio = 0.8
  ) +
  coord_cartesian(ylim = c(min(df_combined$pk), max(df_combined$pk)*1.1))

save.image("~/5_Global_Analisis/2_Parameters_Global/3_Image_Degree_distribution.RData")
load("~/5_Global_Analisis/2_Parameters_Global/3_Image_Degree_distribution.RData")
