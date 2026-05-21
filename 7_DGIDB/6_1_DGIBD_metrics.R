# ====================================================
# Network-based Drug Repurposing using DGIdb
# ====================================================
# Author: Joel 
# Fecha: 2025
# Referencia: Guney et al. 2016 - Nature Communications
# ----------------------------
# SETUP
# ----------------------------
setwd("~/CESC_Network/6_DGIDB/")
library(tidyverse)
library(vroom)
library(igraph)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggrepel)

set.seed(123)

# ----------------------------
# LOAD DGIdb DATA
# ----------------------------
url_interactions <- "https://dgidb.org/data/latest/interactions.tsv"
dgidb_interactions <- read_tsv(url_interactions, show_col_types = FALSE)
glimpse(dgidb_interactions)

table(dgidb_interactions$interaction_type,
      dgidb_interactions$interaction_source_db_name)

# ----------------------------
# FILTER APPROVED DRUGS
# ----------------------------
dgidb_interactions_fda <- dgidb_interactions %>%
  filter(approved == "TRUE")

# Score numérico
dgidb_interactions_fda <- dgidb_interactions_fda %>%
  mutate(
    interaction_score = ifelse(interaction_score == "NULL", "0", interaction_score),
    interaction_score = as.numeric(interaction_score),
    drug_name = as.character(drug_name),
    gene_name = as.character(gene_name)
  )

length(unique(dgidb_interactions_fda$drug_name))  # ~ 3146

# ----------------------------
# LOAD CERVIX NETWORK GENES
# ----------------------------
graph_A7 <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A7_graph_max.rds")
A7_genes <- V(graph_A7)$name

graph_A9 <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A9_graph_max.rds")
A9_genes <- V(graph_A9)$name

genes_universe <- unique(c(A7_genes, A9_genes))
length(genes_universe)

# ----------------------------
# KEEP DRUGS WITH ≥1 GENE IN THE UNIVERSE
# ----------------------------
drug_stats <- dgidb_interactions_fda %>%
  group_by(drug_name) %>%
  summarise(
    n_original_interactions = n(),
    n_in_universe = sum(gene_name %in% genes_universe),
    .groups = "drop"
  )

valid_drugs <- drug_stats %>%
  filter(n_in_universe >= 1)

dgidb_filtered <- dgidb_interactions_fda %>%
  inner_join(valid_drugs, by = "drug_name")

table(dgidb_filtered$interaction_type)

# ----------------------------
# GROUP interaction_type INTO 3 CLASSES
# ----------------------------
dgidb_filtered <- dgidb_filtered %>%
  mutate(
    interaction_type = ifelse(
      interaction_type == "NULL" | interaction_type == "" | is.na(interaction_type),
      "other/unknown",
      interaction_type
    )
  )

interaction_to_group <- function(x) {
  dplyr::case_when(
    # INHIBITS
    x %in% c(
      "antagonist", "antibody", "antisense oligonucleotide",
      "blocker", "cleavage", "inhibitor",
      "inhibitory allosteric modulator",
      "inverse agonist", "negative modulator",
      "partial antagonist", "suppressor"
    ) ~ "inhibits",
    # ACTIVATES
    x %in% c(
      "activator", "agonist", "chaperone", "cofactor",
      "inducer", "partial agonist", "positive modulator",
      "stimulator", "vaccine"
    ) ~ "activates",
    # REGULATES / UNKNOWN
    x %in% c(
      "other/unknown", "adduct", "allosteric modulator",
      "binder", "ligand", "modulator",
      "multitarget", "potentiator", "product of", "substrate",
      "NA", "None", "n/a", "immunotherapy"
    ) ~ "regulates",
    TRUE ~ "regulates"
  )
}

dgidb_filtered <- dgidb_filtered %>%
  mutate(
    interaction_group = interaction_to_group(interaction_type),
    interaction_group = factor(
      interaction_group,
      levels = c("inhibits", "activates", "regulates")
    )
  )

table(dgidb_filtered$interaction_group)

# ----------------------------
# 1) Unique edges restricted to gene universe
# ----------------------------
edges_all <- dgidb_filtered %>%
  filter(gene_name %in% genes_universe) %>%
  distinct(drug_name, gene_name, interaction_group)

# ----------------------------
# 2) Top N drugs and genes
# ----------------------------
N_TOP_DRUGS <- 20
N_TOP_GENES <- 20

drug_stats_global <- edges_all %>%
  group_by(drug_name) %>%
  summarise(n_genes = n_distinct(gene_name), .groups = "drop") %>%
  arrange(desc(n_genes), drug_name)

top_drugs <- drug_stats_global %>%
  slice_head(n = N_TOP_DRUGS)

gene_stats <- edges_all %>%
  group_by(gene_name) %>%
  summarise(n_drugs = n_distinct(drug_name), .groups = "drop") %>%
  arrange(desc(n_drugs), gene_name)

top_genes <- gene_stats %>%
  slice_head(n = N_TOP_GENES)

# ----------------------------
# 3) Edges to plot
# ----------------------------
edges_plot <- edges_all %>%
  mutate(
    interaction_group = factor(
      interaction_group,
      levels = c("inhibits", "activates", "regulates")
    )
  )

# ----------------------------
# 4) Node layout (x-y coordinates)
# ----------------------------
LEFT_X  <- 0
RIGHT_X <- 0.22

## 4a) Drug nodes (left side)
drug_nodes <- edges_plot %>%
  distinct(drug_name) %>%
  left_join(drug_stats_global, by = "drug_name") %>%
  mutate(
    is_top = drug_name %in% top_drugs$drug_name
  ) %>%
  arrange(desc(is_top), desc(n_genes), drug_name) %>%
  mutate(
    y = if_else(is_top, row_number(), row_number() + 15),
    x = LEFT_X,
    label = if_else(
      is_top,
      paste0(drug_name, " (", n_genes, ")"),
      NA_character_
    )
  )

## 4b) Gene nodes (right side)
gene_nodes <- edges_plot %>%
  distinct(gene_name) %>%
  left_join(gene_stats, by = "gene_name") %>%
  arrange(desc(n_drugs), gene_name) %>%
  mutate(
    x = RIGHT_X,
    y = row_number(),
    label = if_else(
      gene_name %in% top_genes$gene_name,
      paste0(gene_name, " (", n_drugs, ")"),
      NA_character_
    )
  )

## 4c) All nodes
nodes_all <- bind_rows(
  drug_nodes %>% mutate(node_type = "Drug"),
  gene_nodes %>% mutate(node_type = "Gene")
)

# ----------------------------
# 5) Edges with coordinates
# ----------------------------
edges_coord <- edges_plot %>%
  left_join(
    drug_nodes %>% select(drug_name, x_drug = x, y_drug = y),
    by = "drug_name"
  ) %>%
  left_join(
    gene_nodes %>% select(gene_name, x_gene = x, y_gene = y),
    by = "gene_name"
  )

# ----------------------------
# 6) Stats for subtitle
# ----------------------------
total_drugs        <- length(unique(dgidb_filtered$drug_name))
total_genes        <- length(unique(dgidb_filtered$gene_name))
total_interactions <- nrow(edges_all)

subtitle_text <- paste0(
  "Total dataset: ", format(total_drugs, big.mark = ","), " drugs, ",
  format(total_genes, big.mark = ","), " genes, ",
  format(total_interactions, big.mark = ","), " interactions"
)

# ----------------------------
# 7) Colors
# ----------------------------
group_colors <- c(
  inhibits  = "#E15759", 
  activates = "#59A14F",  
  regulates = "#4E79A7"   
)

# ----------------------------
# 8) FINAL PLOT
# ----------------------------
p <- ggplot() +
  # Títulos arriba
  annotate(
    "text",
    x = LEFT_X,
    y = max(nodes_all$y) + 20,
    label = paste0("FDA approved drugs\n(n = ", format(total_drugs, big.mark = ","), ")"),
    hjust = 0.5, vjust = 0,
    size = 3.5, fontface = "bold",
    color = "#2D3047", lineheight = 0.9
  ) +
  annotate(
    "text",
    x = RIGHT_X,
    y = max(nodes_all$y) + 20,
    label = paste0("Genes in cervix cancer network\n(n = ", format(total_genes, big.mark = ","), ")"),
    hjust = 0.5, vjust = 0,
    size = 3.5, fontface = "bold",
    color = "#2D3047", lineheight = 0.9
  ) +
  
  # Edges
  geom_segment(
    data = edges_coord,
    aes(
      x = x_drug, y = y_drug,
      xend = x_gene, yend = y_gene,
      colour = interaction_group
    ),
    alpha = 0.02,
    linewidth = 0.3
  ) +
  
  # Nodes
  geom_point(
    data = nodes_all,
    aes(x = x, y = y, shape = node_type),
    size = 2.5,
    colour = "grey20",
    fill = ifelse(nodes_all$node_type == "Drug", "#FF9D9A", "#76B7B2"),
    stroke = 0
  ) +
  
  # Etiquetas fármacos top
  geom_text_repel(
    data = filter(drug_nodes, !is.na(label)),
    aes(x = x, y = y, label = label),
    direction = "y",
    hjust = 1,
    nudge_x = -0.03,
    nudge_y = 0.6,          # constante, sin idx
    size = 3.2,
    fontface = "bold",
    color = "#2D3047",
    segment.size = 0.2,
    segment.color = "grey50",
    segment.alpha = 3,
    box.padding = 0.3,
    point.padding = 0.25,
    max.overlaps = 50,
    min.segment.length = 0.1,
    seed = 123
  ) +
  
  # Etiquetas genes top
  geom_text_repel(
    data = filter(gene_nodes, !is.na(label)),
    aes(x = x, y = y, label = label),
    direction = "y",
    hjust = 0,
    nudge_x = 0.03,
    nudge_y = 0.8,          # constante
    size = 3.2,
    fontface = "bold",
    color = "#2D3047",
    segment.size = 0.2,
    segment.color = "grey50",
    segment.alpha = 3,
    box.padding = 0.3,
    point.padding = 0.25,
    max.overlaps = 50,
    min.segment.length = 0.1,
    seed = 456
  ) +
  
  scale_colour_manual(
    name = "Interaction group (edges)",
    values = group_colors,
    guide = guide_legend(
      override.aes = list(
        linewidth = 2,
        alpha = 0.5
      ),
      keywidth = unit(1.2, "cm"),
      keyheight = unit(0.3, "cm")
    )
  ) +
  scale_shape_manual(
    name = "Node type",
    values = c(Drug = 21, Gene = 19),
    guide = guide_legend(
      override.aes = list(
        size = 3,
        fill = c("#FF6B6B", "#4ECDC4")
      )
    )
  ) +
  scale_x_continuous(
    breaks = NULL,
    limits = c(LEFT_X - 0.15, RIGHT_X + 0.15),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.10, 0.25))
  ) +
  labs(
    title = "Bipartite map of drug-gene interactions (DGIdb)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "top",
    legend.box       = "vertical",
    legend.margin    = margin(t = 10, b = 5),
    legend.spacing.x = unit(0.4, "cm"),
    legend.spacing.y = unit(0, "cm"),
    panel.grid       = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks       = element_blank(),
    axis.text.x      = element_text(
      face = "bold",
      size = 11,
      color = "#2D3047",
      margin = margin(t = 25, b = 5)
    ),
    plot.title = element_text(
      face = "bold",
      size = 16,
      color = "#2D3047",
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 10.5,
      color = "#555555",
      hjust = 0.5,
      lineheight = 1.3,
      margin = margin(b = 20, t = 5)
    ),
    legend.title = element_text(
      face = "bold",
      size = 10,
      color = "#2D3047"
    ),
    legend.text = element_text(
      size = 9,
      color = "#555555"
    ),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.margin       = margin(15, 20, 15, 20)
  )

# Guardar con dimensiones óptimas
ggsave(
  "drug_gene_bipartite_top20_final.png",
  p,
  width  = 18,
  height = 10,
  dpi    = 800,
  bg     = "white"
)



# Guardar con dimensiones óptimas
ggsave(
  "drug_gene_bipartite_top20_final.pdf",
  p,
  width  = 18,
  height = 10,
  dpi    = 800,
  bg     = "white"
)


cat("\n✓ Gráfico guardado como: drug_gene_bipartite_top20_final.png\n")
cat("✓ Estadísticas incluidas:\n")
cat("  - Total drugs (drug_name):", format(total_drugs, big.mark = ","), "\n")
cat("  - Total genes (gene_name):", format(total_genes, big.mark = ","), "\n")
cat("  - Total interactions:", format(total_interactions, big.mark = ","), "\n")

