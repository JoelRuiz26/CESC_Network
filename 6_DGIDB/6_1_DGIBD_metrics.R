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
length(A7_genes) #[1] 19502

graph_A9 <- readRDS("~/CESC_Network/5_Modularity/5_1_Modules/5_2_1_Louvain_A9_graph_max.rds")
A9_genes <- V(graph_A9)$name
length(A9_genes) #[1] 19502

# Gene universe
genes_universe <- unique(c(A7_genes, A9_genes))
length(genes_universe) #[1] 19502
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



#############plot#############}
# Cargar librerías
library(tidyverse)
library(ggrepel)
library(patchwork)

# CORREGIR: Usar gene_claim_name en lugar de gene_name y mantener todos los registros
dgidb_filtered_corrected <- dgidb_filtered %>%
  mutate(gene_display = ifelse(gene_claim_name == "NULL" | is.na(gene_claim_name), 
                               gene_name, 
                               gene_claim_name))

# Preparar datos para los tres paneles
# Panel 1: Distribución de genes por fármaco (usando gene_claim_name)
farmacos_genes <- dgidb_filtered_corrected %>%
  dplyr::select(drug_name, gene_display) %>%
  distinct() %>%
  group_by(drug_name) %>%
  summarise(n_genes = n()) %>%
  ungroup() %>%
  arrange(desc(n_genes))

# Panel 2: Distribución de fármacos por gen (usando gene_claim_name)
genes_farmacos <- dgidb_filtered_corrected %>%
  dplyr::select(gene_display, drug_name) %>%
  distinct() %>%
  group_by(gene_display) %>%
  summarise(n_drugs = n()) %>%
  ungroup() %>%
  arrange(desc(n_drugs))

# Panel 3: Distribución del interaction score promedio por fármaco
farmacos_score <- dgidb_filtered_corrected %>%
  group_by(drug_name) %>%
  summarise(
    mean_score = mean(interaction_score, na.rm = TRUE),
    median_score = median(interaction_score, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(desc(mean_score))

# Identificar los top para etiquetar (reducido para evitar solapamiento)
top_farmacos <- farmacos_genes %>% head(6)
top_genes <- genes_farmacos %>% head(6)
top_scores <- farmacos_score %>% 
  filter(mean_score > quantile(farmacos_score$mean_score, 0.95, na.rm = TRUE)) %>%
  head(6)

# Estadísticas
stats_farmacos <- farmacos_genes %>%
  summarise(
    total = n(),
    media = mean(n_genes),
    mediana = median(n_genes),
    maximo = max(n_genes)
  )

stats_genes <- genes_farmacos %>%
  summarise(
    total = n(),
    media = mean(n_drugs),
    mediana = median(n_drugs),
    maximo = max(n_drugs)
  )

stats_score <- farmacos_score %>%
  summarise(
    total = n(),
    media_mean = mean(mean_score, na.rm = TRUE),
    mediana_mean = median(mean_score, na.rm = TRUE),
    maximo_mean = max(mean_score, na.rm = TRUE)
  )

# FUNCIÓN MEJORADA para calcular posiciones de etiquetas
calcular_posiciones_inteligentes <- function(data, n_labels, type = "histogram") {
  if (type == "histogram") {
    hist_info <- hist(data, breaks = 50, plot = FALSE)
    # Usar diferentes niveles de densidad para distribuir etiquetas
    niveles <- seq(0.1, 0.9, length.out = n_labels)
    posiciones <- quantile(hist_info$density, niveles, na.rm = TRUE)
  } else {
    # Para densidad simple
    dens <- density(data, na.rm = TRUE)
    niveles <- seq(0.1, 0.9, length.out = n_labels)
    posiciones <- quantile(dens$y, niveles, na.rm = TRUE)
  }
  return(posiciones)
}

# PANEL 1: Distribución de genes por fármaco (MEJORADO)
y_positions1 <- calcular_posiciones_inteligentes(farmacos_genes$n_genes, nrow(top_farmacos))
top_farmacos_labeled <- top_farmacos %>%
  mutate(
    y_position = y_positions1[1:n()],
    label = paste0(drug_name, "\n(", n_genes, " genes)")
  )

panel1 <- ggplot(farmacos_genes, aes(x = n_genes)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    fill = "#1f77b4",
    alpha = 0.7,
    color = "white",
    size = 0.2
  ) +
  geom_density(
    color = "#ff7f0e",
    size = 1.2,
    alpha = 0.8
  ) +
  geom_point(
    data = top_farmacos_labeled,
    aes(y = y_position),
    color = "#d62728",
    size = 3,
    alpha = 0.8
  ) +
  geom_text_repel(
    data = top_farmacos_labeled,
    aes(y = y_position, label = label),
    size = 3.0,
    color = "black",
    fontface = "bold",
    segment.color = "gray40",
    segment.size = 0.4,
    box.padding = 0.8,
    point.padding = 0.5,
    max.overlaps = 15,
    min.segment.length = 0.1,
    force = 1,
    nudge_x = 10,
    nudge_y = 0.01,
    direction = "y"
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Distribution of Unique Target Genes per Drug",
    subtitle = "Top 6 multi-target drugs",
    x = "Number of Unique Target Genes per Drug",
    y = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "gray30", size = 9),
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank()
  )

# PANEL 2: Distribución de fármacos por gen (MEJORADO)
y_positions2 <- calcular_posiciones_inteligentes(genes_farmacos$n_drugs, nrow(top_genes))
top_genes_labeled <- top_genes %>%
  mutate(
    y_position = y_positions2[1:n()],
    label = paste0(gene_display, "\n(", n_drugs, " drugs)")
  )

panel2 <- ggplot(genes_farmacos, aes(x = n_drugs)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    fill = "#2ca02c",
    alpha = 0.7,
    color = "white",
    size = 0.2
  ) +
  geom_density(
    color = "#9467bd",
    size = 1.2,
    alpha = 0.8
  ) +
  geom_point(
    data = top_genes_labeled,
    aes(y = y_position),
    color = "#e377c2",
    size = 3,
    alpha = 0.8
  ) +
  geom_text_repel(
    data = top_genes_labeled,
    aes(y = y_position, label = label),
    size = 3.0,
    color = "black",
    fontface = "bold",
    segment.color = "gray40",
    segment.size = 0.4,
    box.padding = 0.8,
    point.padding = 0.5,
    max.overlaps = 15,
    min.segment.length = 0.1,
    force = 1,
    nudge_x = 15,
    nudge_y = 0.01,
    direction = "y"
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Distribution of Unique Drugs per Target Gene",
    subtitle = "Top 6 druggable genes",
    x = "Number of Unique Drugs per Target Gene",
    y = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "gray30", size = 9),
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank()
  )

# PANEL 3: Distribución del interaction score promedio por fármaco (NUEVO)
y_positions3 <- calcular_posiciones_inteligentes(farmacos_score$mean_score, nrow(top_scores), "density")
top_scores_labeled <- top_scores %>%
  mutate(
    y_position = y_positions3[1:n()],
    label = paste0(drug_name, "\n(", round(mean_score, 2), ")")
  )

panel3 <- ggplot(farmacos_score, aes(x = mean_score)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    fill = "#ff7f0e",
    alpha = 0.7,
    color = "white",
    size = 0.2
  ) +
  geom_density(
    color = "#1f77b4",
    size = 1.2,
    alpha = 0.8
  ) +
  geom_point(
    data = top_scores_labeled,
    aes(y = y_position),
    color = "#d62728",
    size = 3,
    alpha = 0.8
  ) +
  geom_text_repel(
    data = top_scores_labeled,
    aes(y = y_position, label = label),
    size = 3.0,
    color = "black",
    fontface = "bold",
    segment.color = "gray40",
    segment.size = 0.4,
    box.padding = 0.8,
    point.padding = 0.5,
    max.overlaps = 15,
    min.segment.length = 0.1,
    force = 1,
    nudge_x = 2,
    nudge_y = 0.01,
    direction = "y"
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.1))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Distribution of Mean Interaction Score per Drug",
    subtitle = "Top 6 drugs with highest mean interaction scores",
    x = "Mean Interaction Score per Drug",
    y = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "gray30", size = 9),
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank()
  )

# COMBINAR LOS TRES PANELES
plot_combinado <- (panel1 + panel2) / panel3 +
  plot_annotation(
    title = "Comprehensive Drug-Target Interaction Analysis",
    subtitle = paste(
      "Based on", format(stats_farmacos$total, big.mark = ","), 
      "FDA-approved drugs"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 10, color = "gray50", hjust = 0.5)
    )
  )

# Mostrar el gráfico combinado
print(plot_combinado)

# Guardar en alta calidad
ggsave("comprehensive_drug_gene_analysis.png", plot_combinado,
       width = 16, height = 12, dpi = 300, bg = "white")
ggsave("comprehensive_drug_gene_analysis.pdf", plot_combinado,
       width = 16, height = 12)

