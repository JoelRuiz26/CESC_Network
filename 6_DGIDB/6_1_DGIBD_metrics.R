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

# ----------------------------
# LOAD DGIdb DATA
# ----------------------------
url_interactions <- "https://dgidb.org/data/latest/interactions.tsv"
dgidb_interactions <- read_tsv(url_interactions, show_col_types = FALSE)
glimpse(dgidb_interactions)
#Rows: 98,239
#Columns: 13
#$ gene_claim_name               <chr> "CYP2D6", "PPARG", "ATAD5", "…
#$ gene_concept_id               <chr> "hgnc:2625", "hgnc:9236", "hg…
#$ gene_name                     <chr> "CYP2D6", "PPARG", "ATAD5", "…
#$ interaction_source_db_name    <chr> "DTC", "DTC", "DTC", "DTC", "…
#$ interaction_source_db_version <chr> "9/2/20", "9/2/20", "9/2/20",…
#$ interaction_type              <chr> "NULL", "NULL", "NULL", "NULL…
#$ interaction_score             <chr> "0.017709164", "0.84012274", …
#$ drug_claim_name               <chr> "RACLOPRIDE", "KALOPANAX-SAPO…
#$ drug_concept_id               <chr> "ncit:C152139", "chembl:CHEMB…
#$ drug_name                     <chr> "RACLOPRIDE", "CHEMBL:CHEMBL1…
#$ approved                      <chr> "FALSE", "FALSE", "FALSE", "F…
#$ immunotherapy                 <chr> "FALSE", "FALSE", "FALSE", "F…
#$ anti_neoplastic               <chr> "FALSE", "FALSE", "FALSE", "F…

table(dgidb_interactions$interaction_type, dgidb_interactions$interaction_source_db_name)

#CancerCommons   CGI ChEMBL CIViC ClearityFoundationBiomarkers ClearityFoundationClinicalTrial COSMIC  DoCM   DTC   FDA
#activator                             0     0    134     0                            0                               4      0     0     0     0
#agonist                               0     0   1443     0                            0                               0      0     0     0     0
#antibody                              6     0      0     0                            0                               0      0     0     0     0
#antisense oligonucleotide             0     0      0     0                            0                               0      0     0     0     0
#binder                                0     0    259     0                            0                               0      0     0     0     0
#blocker                               0     0   1149     0                            0                               1      0     0     0     0
#cleavage                              0     0     82     0                            0                               0      0     0     0     0
#immunotherapy                         0     0      0     0                            0                               1      0     0     0     0
#inhibitor                           103     0  11910     0                            0                             232      0     0    94     0
#inverse agonist                       0     0     36     0                            0                               0      0     0     0     0
#modulator                             0     0    609     0                            0                               0      0     0     0     0
#negative modulator                    0     0    142     0                            0                               0      0     0     0     0
#NULL                                  0   368      0  1120                          163                               0     37    76 23785   427
#other/unknown                         0     0     68     0                            0                               2      0     0     0     0
#positive modulator                    0     0   1039     0                            0                               0      0     0     0     0
#potentiator                           0     0     51     0                            0                               0      0     0     0     0
#vaccine                               0     0     19     0                            0                               0      0     0     0     0

#GuideToPharmacology JAX-CKB MyCancerGenome MyCancerGenomeClinicalTrial   NCI OncoKB PharmGKB  TALC TdgClinicalTrial  TEND   TTD
#activator                                 439       0              0                           1     0      0        0     0                0     0     0
#agonist                                  4397       0              0                           0     0      0        0     9                0     0     0
#antibody                                  211       0             57                           0     0      0        0    21                0     0     0
#antisense oligonucleotide                   0       0              0                           0     0      0        0     4                0     0     0
#binder                                      0       0              0                           0     0      0        0     0                0     0     0
#blocker                                   669       0              0                           0     0      0        0     0                0     0     0
#cleavage                                    0       0              0                           0     0      0        0     0                0     0     0
#immunotherapy                               0       0             10                           0     0      0        0     0                0     0     0
#inhibitor                                9447       0            613                         316     0      0        0   522                0     0     0
#inverse agonist                             0       0              0                           0     0      0        0     0                0     0     0
#modulator                                 630       0              0                           0     0      0        0     0                0     0     0
#negative modulator                          0       0              0                           0     0      0        0     0                0     0     0
#NULL                                      277    2148              0                           0  6231    148     5976    22             5071  2242 15245
#other/unknown                               0       0            149                           0     0      0        0     0                0     0     0
#positive modulator                          0       0              0                           0     0      0        0     0                0     0     0
#potentiator                                 0       0              0                           0     0      0        0     0                0     0     0
#vaccine                                     0       0              0                           0     0      0        0    24                0     0     0


# ----------------------------
# FILTER APPROVED DRUGS
# ----------------------------
dgidb_interactions_fda <- dgidb_interactions %>%
  filter(approved == "TRUE")
dim(dgidb_interactions) #[1] 37970    13
table(dgidb_interactions_fda$interaction_type, dgidb_interactions_fda$interaction_source_db_name)

# Make score numeric type
# Filtro mínimo de confianza por interacción (ajusta el umbral)
dgidb_interactions_fda <- dgidb_interactions_fda %>%
  mutate(
    interaction_score = ifelse(interaction_score == "NULL", "0", interaction_score),
    interaction_score = as.numeric(interaction_score)
  )
summary(dgidb_interactions_fda$interaction_score)
#> summary(dgidb_interactions_fda$interaction_score)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.00000   0.02012   0.09299   0.71798   0.38467 157.52301 

dgidb_interactions_fda <- dgidb_interactions_fda %>%
  mutate(
    drug_name = as.character(drug_name),
    gene_name = as.character(gene_name)
  )

length(unique(dgidb_interactions_fda$drug_name)) #[1] 3146


graph_A7 <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A7_graph_max.rds")
A7_genes <- V(graph_A7)$name
length(A7_genes) #[1] 17481


graph_A9 <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A9_graph_max.rds")
A9_genes <- V(graph_A9)$name
length(A9_genes) #[1] 17481

# Gene universe
genes_universe <- unique(c(A7_genes, A9_genes))
length(genes_universe) #[1] 174812
# Gene universe (you already have this)

#Compute, per drug, how many interactions it has overall
#    and how many are with genes in your universe
drug_stats <- dgidb_interactions_fda %>%
  group_by(drug_name) %>%
  summarise(
    n_original_interactions = n(),
    n_in_universe = sum(gene_name %in% genes_universe),
    .groups = "drop"
  )

# Keep only drugs with at least 1 interaction in your gene universe
valid_drugs <- drug_stats %>%
  filter(n_in_universe >= 1)

# Filter the main DGIdb table to those drugs
#    and attach the informative counts
dgidb_filtered <- dgidb_interactions_fda %>%
  inner_join(valid_drugs, by = "drug_name")

table(dgidb_filtered$interaction_type)
#activator            agonist           antibody 
#154               1590                 81 
#binder            blocker           cleavage 
#80               1053                 69 
#immunotherapy          inhibitor    inverse agonist 
#7               8709                 11 
#modulator negative modulator               NULL 
#301                 49              23822 
#other/unknown positive modulator        potentiator 
#157                939                 42 


# ----------------------------
# Group interaction types into broader classes
# ----------------------------

# Recode NULL / empty / NA as "other/unknown" for clarity
dgidb_filtered <- dgidb_filtered %>%
  mutate(
    interaction_type = ifelse(
      interaction_type == "NULL" | interaction_type == "" | is.na(interaction_type),
      "other/unknown",
      interaction_type
    )
  )

# Map fine-grained interaction_type to 3 broad groups:
#   - inhibits
#   - activates
#   - regulates (includes unknown / mixed / non-directional)
interaction_to_group <- function(x) {
  dplyr::case_when(
    # INHIBITS group
    x %in% c(
      "antagonist", "antibody", "antisense oligonucleotide",
      "blocker", "cleavage", "inhibitor",
      "inhibitory allosteric modulator",
      "inverse agonist", "negative modulator",
      "partial antagonist", "suppressor"
    ) ~ "inhibits",
    
    # ACTIVATES group
    x %in% c(
      "activator", "agonist", "chaperone", "cofactor",
      "inducer", "partial agonist", "positive modulator",
      "stimulator", "vaccine"
    ) ~ "activates",
    
    # REGULATES / UNKNOWN group
    x %in% c(
      "other/unknown", "adduct", "allosteric modulator",
      "binder", "ligand", "modulator",
      "multitarget", "potentiator", "product of", "substrate",
      "NA", "None", "n/a", "immunotherapy"
    ) ~ "regulates",
    
    # Any other residual value is treated as regulates
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

# Quick check
table(dgidb_filtered$interaction_type)
table(dgidb_filtered$interaction_group)

#inhibits activates regulates 
#9972      2683     24409 



library(tidyverse)
library(ggrepel)

set.seed(123)

# ----------------------------
# 1) Unique edges restricted to gene universe
# ----------------------------
edges_all <- dgidb_filtered %>%
  filter(gene_name %in% genes_universe) %>%
  distinct(drug_name, gene_name, interaction_group)

# ----------------------------
# 2) Top 10 drugs (global) and top 10 genes
# ----------------------------

N_TOP_DRUGS <- 20
N_TOP_GENES <- 20

# nº de genes únicos por fármaco (todas las interacciones)
drug_stats_global <- edges_all %>%
  group_by(drug_name) %>%
  summarise(n_genes = n_distinct(gene_name), .groups = "drop") %>%
  arrange(desc(n_genes), drug_name)

top_drugs <- drug_stats_global %>%
  slice_head(n = N_TOP_DRUGS)

# nº de fármacos únicos por gen (global)
gene_stats <- edges_all %>%
  group_by(gene_name) %>%
  summarise(n_drugs = n_distinct(drug_name), .groups = "drop") %>%
  arrange(desc(n_drugs), gene_name)

top_genes <- gene_stats %>%
  slice_head(n = N_TOP_GENES)

# ----------------------------
# 3) Subset of edges to plot
#    (keep all edges touching top drugs OR top genes)
# ----------------------------
# ----------------------------
# 3) Edges to plot
#    (ALL edges; top lists are only for labels)
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

# Make plot narrower by reducing X-distance between groups
LEFT_X  <- 0
RIGHT_X <- 0.35   # más delgado el “cuerpo” de líneas

## 4a) Drug nodes (left side), ordered by n_genes desc
## 4a) Drug nodes (left side) — clean ordered block for top 10
## 4a) Drug nodes (left side) — clean ordered block for top 10
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
    # texto base
    label_raw = paste0(drug_name, " (", n_genes, ")"),
    # wrap solo si es top, con cortes en espacios
    label = if_else(
      is_top,
      paste0(drug_name, " (", n_genes, ")"),
      NA_character_
    )
  )



## 4b) Gene nodes (right side)
## 4a) Drug nodes (left side) — clean ordered block for top N
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
    
    # Etiqueta final SIN wraps (una sola línea)
    label = if_else(
      is_top,
      paste0(drug_name, " (", n_genes, ")"),
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
# 6) Colors (only for edges = interaction_group)
#    Use a highly contrasting ColorBrewer-like palette
# ====================================================
# CALCULAR ESTADÍSTICAS PARA EL SUBTÍTULO
# ====================================================
# Calcular estadísticas únicas del dataset filtrado
total_drugs <- length(unique(dgidb_filtered$drug_name))  # Usar drug_name en lugar de drug_claim_name
total_genes <- length(unique(dgidb_filtered$gene_name))  # Usar gene_name en lugar de gene_claim_name
total_interactions <- nrow(edges_all)

# Crear texto del subtítulo con estadísticas
subtitle_text <- paste0(
  "Total dataset: ", format(total_drugs, big.mark = ","), " drugs, ", 
  format(total_genes, big.mark = ","), " genes, ", 
  format(total_interactions, big.mark = ","), " interactions"
)

# ====================================================
# 6) Colors (only for edges = interaction_group)
#    Use high-contrast, colorblind-friendly palette
# ====================================================
group_colors <- c(
  inhibits  = "#D73027",  # Strong red
  activates = "#1A9850",  # Strong green
  regulates = "#4575B4"   # Strong blue
)

# ====================================================
# 7) FINAL PLOT con etiquetas superiores mejoradas
# ====================================================
p <- ggplot() +
  
  # 1. AÑADIR TEXTO EN LA PARTE SUPERIOR DE LAS BARRAS
  # Texto para fármacos (lado izquierdo)
  annotate(
    "text",
    x = LEFT_X,
    y = max(nodes_all$y) + 20,  # Posicionado arriba
    label = paste0("FDA approved drugs\n(n = ", format(total_drugs, big.mark = ","), ")"),
    hjust = 0.5,
    vjust = 0,
    size = 3.5,
    fontface = "bold",
    color = "#2D3047",
    lineheight = 0.9
  ) +
  
  # Texto para genes (lado derecho)
  annotate(
    "text",
    x = RIGHT_X,
    y = max(nodes_all$y) + 20,  # Posicionado arriba
    label = paste0("Genes in cervix cancer network\n(n = ", format(total_genes, big.mark = ","), ")"),
    hjust = 0.5,
    vjust = 0,
    size = 3.5,
    fontface = "bold",
    color = "#2D3047",
    lineheight = 0.9
  ) +
  
  # 2. EDGES (conexiones)
  geom_segment(
    data = edges_coord,
    aes(
      x = x_drug, y = y_drug,
      xend = x_gene, yend = y_gene,
      colour = interaction_group
    ),
    alpha = 0.05,
    linewidth = 0.3
  ) +
  
  # 3. NODES
  geom_point(
    data = nodes_all,
    aes(x = x, y = y, shape = node_type),
    size = 2.5,
    colour = "grey20",
    fill = ifelse(nodes_all$node_type == "Drug", "#FF6B6B", "#4ECDC4"),
    stroke = 0.4
  ) +
  
  # 4. ETIQUETAS PARA FÁRMACOS TOP
  # 4. ETIQUETAS PARA FÁRMACOS TOP
  geom_text_repel(
    data = filter(drug_nodes, !is.na(label)),
    aes(x = x, y = y, label = label),
    direction = "y",
    hjust = 1,
    nudge_x = -0.04,
    nudge_y = ifelse(row_number(filter(drug_nodes, !is.na(label))) %% 2 == 0, 0.8, -0.8),
    size = 3.2,
    fontface = "bold",
    color = "#2D3047",
    segment.size = 0.2,
    segment.color = "grey50",
    segment.alpha = 2.5,
    box.padding = 0.3,
    point.padding = 0.25,
    max.overlaps = 50,
    min.segment.length = 0.1,
    seed = 123
  ) +
  
  # 5. ETIQUETAS PARA GENES TOP
  # 5. ETIQUETAS PARA GENES TOP
  geom_text_repel(
    data = filter(gene_nodes, !is.na(label)),
    aes(x = x, y = y, label = label),
    direction = "y",
    hjust = 0,
    nudge_x = 0.04,
    nudge_y = ifelse(row_number(filter(gene_nodes, !is.na(label))) %% 2 == 0, 1.0, -1.0),
    size = 3.2,
    fontface = "bold",
    color = "#2D3047",
    segment.size = 0.2,
    segment.color = "grey50",
    segment.alpha = 2.5,
    box.padding = 0.3,
    point.padding = 0.25,
    max.overlaps = 50,
    min.segment.length = 0.1,
    seed = 456
  ) +
  
  # 6. ESCALAS DE COLOR Y FORMA
  scale_colour_manual(
    name = "Interaction group (edges)",
    values = group_colors,
    guide = guide_legend(
      override.aes = list(
        linewidth = 2,
        alpha = 0.8
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
  
  # 7. ESCALA X (con títulos de grupo movidos más abajo)
  scale_x_continuous(
    breaks = NULL,   # <- sin marcas numéricas en el eje X
    limits = c(LEFT_X - 0.15, RIGHT_X + 0.15),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  
  # 8. ESCALA Y (con más espacio en la parte superior para las nuevas etiquetas)
  scale_y_continuous(
    expand = expansion(mult = c(0.10, 0.25))   # más espacio arriba SIN hacer más alto el plot
  ) +
  
  # 9. ETIQUETAS DEL GRÁFICO
  labs(
    title = "Bipartite map of drug-gene interactions (DGIdb)",
    x = NULL,
    y = NULL
  ) +
  
  # 10. TEMA MEJORADO
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.margin = margin(t = 10, b = 5),
    legend.spacing.x = unit(0.4, "cm"),
    legend.spacing.y = unit(0, "cm"),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      face = "bold",
      size = 11,
      color = "#2D3047",
      margin = margin(t = 25, b = 5)  # Más margen superior para los títulos
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
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 20, 15, 20)  # Márgenes ajustados
  )

# Mostrar el gráfico
print(p)

# Guardar con dimensiones óptimas
ggsave("drug_gene_bipartite_top10_final.png", p,
       width = 18,   # <-- más ancho
       height = 10,
       dpi = 400,
       bg = "white")


# Mensaje de confirmación
cat("\n✓ Gráfico guardado como: drug_gene_bipartite_top10_final.png\n")
cat("✓ Estadísticas incluidas:\n")
cat("  - Total drugs (drug_name):", format(total_drugs, big.mark = ","), "\n")
cat("  - Total genes (gene_name):", format(total_genes, big.mark = ","), "\n")
cat("  - Total interactions:", format(total_interactions, big.mark = ","), "\n")