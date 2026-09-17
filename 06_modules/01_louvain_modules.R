#Script for Louvain cluster
#Input: Networks (.graphml)
#Output: Networks (graphml) with communities (size >= 10) named M1, M2, M3, ...
#        by descending size (M1 = largest module), filtered for enrichment.
#By Joel Ruiz

### Load packages ###
library(igraph)
library(tidyverse)
library(ggplot2)
set.seed(123)

### Load networks ###
graph_A7 <- readRDS('~/CESC_Network/05_topology/results/5_2_4_Subnetwork_A7_max.rds')
graph_A9 <- readRDS('~/CESC_Network/05_topology/results/5_2_4_Subnetwork_A9_max.rds')

###### For A7 network ######
### Louvain ###
Louvain_A7 <- cluster_louvain(graph_A7)        #125 communities
V(graph_A7)$Community <- Louvain_A7$membership 

### Normalized degree (degree / max_degree_in_network), kept as node attribute ###
deg_A7     <- degree(graph_A7)
norm_deg_A7 <- deg_A7 / max(deg_A7)
V(graph_A7)$NormDegree <- norm_deg_A7

### PageRank (kept as node attribute) ###
V(graph_A7)$PageRank <- page_rank(graph_A7)$vector

### Filter nodes for enrichment ###
V(graph_A7)$Community_size <- table(V(graph_A7)$Community)[as.character(V(graph_A7)$Community)] #community size
nodes_to_keep_A7 <- V(graph_A7)[V(graph_A7)$Community_size >= 10]
filtered_graph_A7 <- induced_subgraph(graph_A7, vids = nodes_to_keep_A7)

### Name modules M1, M2, M3, ... by descending size (M1 = largest) ###
mod_rank_A7 <- names(sort(table(V(filtered_graph_A7)$Community), decreasing = TRUE))
mod_name_A7 <- setNames(paste0("M", seq_along(mod_rank_A7)), mod_rank_A7)
V(filtered_graph_A7)$Community_name <- unname(mod_name_A7[as.character(V(filtered_graph_A7)$Community)])


###### For A9 network ######
### Louvain ###
Louvain_A9 <- cluster_louvain(graph_A9) #235
V(graph_A9)$Community <- Louvain_A9$membership

### Normalized degree, kept as node attribute ###
deg_A9      <- degree(graph_A9)
norm_deg_A9 <- deg_A9 / max(deg_A9)
V(graph_A9)$NormDegree <- norm_deg_A9

### PageRank (kept as node attribute) ###
V(graph_A9)$PageRank <- page_rank(graph_A9)$vector

### Filter nodes for enrichment ###
V(graph_A9)$Community_size <- table(V(graph_A9)$Community)[as.character(V(graph_A9)$Community)]
nodes_to_keep_A9 <- V(graph_A9)[V(graph_A9)$Community_size >= 10]
filtered_graph_A9 <- induced_subgraph(graph_A9, vids = nodes_to_keep_A9)

### Name modules M1, M2, M3, ... by descending size (M1 = largest) ###
mod_rank_A9 <- names(sort(table(V(filtered_graph_A9)$Community), decreasing = TRUE))
mod_name_A9 <- setNames(paste0("M", seq_along(mod_rank_A9)), mod_rank_A9)
V(filtered_graph_A9)$Community_name <- unname(mod_name_A9[as.character(V(filtered_graph_A9)$Community)])


### Save filtered graphs ###
saveRDS(filtered_graph_A9, file = "~/CESC_Network/06_modules/results/5_1_1_Louvain_A9_graph_max.rds")
saveRDS(filtered_graph_A7, file = "~/CESC_Network/06_modules/results/5_1_1_Louvain_A7_graph_max.rds")


### Verify ###
# Dataframe para A7
A7_Louvain_genes <- tibble(
  Community = as.integer(V(filtered_graph_A7)$Community),
  Community_name = as.character(V(filtered_graph_A7)$Community_name),
  gene = as.character(V(filtered_graph_A7)$name)
) 

A7_sizes <- A7_Louvain_genes %>%
  group_by(Community, Community_name) %>%
  summarise(size = n(), .groups = "drop")

# Dataframe para A9
A9_Louvain_genes <- tibble(
  Community = as.integer(V(filtered_graph_A9)$Community),
  Community_name = as.character(V(filtered_graph_A9)$Community_name),
  gene = as.character(V(filtered_graph_A9)$name)
) 

A9_sizes <- A9_Louvain_genes %>%
  group_by(Community, Community_name) %>%
  summarise(size = n(), .groups = "drop")


library(patchwork)

A7_sizes <- A7_sizes %>% arrange(desc(size)) %>%
  mutate(Community_name = factor(Community_name, levels = Community_name))
A9_sizes <- A9_sizes %>% arrange(desc(size)) %>%
  mutate(Community_name = factor(Community_name, levels = Community_name))

base_theme <- theme_bw(base_family = "Times New Roman") +
  theme(
    axis.text.x        = element_text(size = 17, face = "bold", angle = 45, hjust = 1),
    axis.text.y        = element_text(size = 20, face = "bold"),
    axis.title.x       = element_text(size = 22, face = "bold", margin = margin(t = 8)),
    axis.title.y       = element_text(size = 22, face = "bold", margin = margin(r = 8)),
    plot.title         = element_text(size = 26, face = "bold", hjust = 0.5,
                                      margin = margin(b = 8)),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border       = element_rect(color = "black", linewidth = 0.8)
  )

pA7 <- ggplot(A7_sizes, aes(x = Community_name, y = size)) +
  geom_bar(stat = "identity", fill = "#003366", alpha = 0.88, width = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::comma) +
  labs(title = sprintf("HPV-A7  (%d modules)", nrow(A7_sizes)), x = NULL, y = "Number of genes") +
  base_theme

pA9 <- ggplot(A9_sizes, aes(x = Community_name, y = size)) +
  geom_bar(stat = "identity", fill = "#CC2200", alpha = 0.88, width = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::comma) +
  labs(title = sprintf("HPV-A9  (%d modules)", nrow(A9_sizes)), x = NULL, y = "Number of genes") +
  base_theme

combined <- pA7 + pA9 + plot_layout(ncol = 2) +
  plot_annotation(
    title = "Louvain Module Size Distribution — HPV Co-expression Networks",
    theme = theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5,
                                             family = "Times New Roman",
                                             margin = margin(b = 12)))
  )

outfile <- "~/CESC_Network/06_modules/results/5_1_2_Louvain_community_sizes_A7_A9_max.png"
ggsave(filename = outfile, plot = combined, width = 28, height = 12, dpi = 300, bg = "white")
message("[OK] Saved: ", outfile)

