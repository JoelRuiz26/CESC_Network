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

# -----




# =========================================================
# 2) Analizar fármacos significativos por z-score
#    y su tipo de interacción (inhibits / activates / regulates)
# =========================================================

# --- cargar resultados de proximidad (zc) ---
prox_A7 <- readr::read_tsv(
  "8_5_ClosestProximity_A7_FULL.tsv",
  show_col_types = FALSE
)

prox_A9 <- readr::read_tsv(
  "8_5_ClosestProximity_A9_FULL.tsv",
  show_col_types = FALSE


# Umbral clásico en network medicine (~p = 0.05 one-sided)
z_thresh <- -1.65

sig_A7 <- prox_A7 %>%
  filter(!is.na(zc), zc <= z_thresh) %>%
  mutate(Cohort = "A7")

sig_A9 <- prox_A9 %>%
  filter(!is.na(zc), zc <= z_thresh) %>%
  mutate(Cohort = "A9")

sig_all <- bind_rows(sig_A7, sig_A9)

cat("Umbral zc <=", z_thresh, "\n")
cat("  A7: fármacos significativos =", n_distinct(sig_A7$drug), "\n")
cat("  A9: fármacos significativos =", n_distinct(sig_A9$drug), "\n\n")

# =========================================================
# 3) Cruzar con DGIdb (edges_all) y contar por tipo de interacción
# =========================================================

# edges_all viene de arriba:
#   distinct(drug_name, gene_name, interaction_group)

edges_sig <- edges_all %>%
  inner_join(sig_all,
             by = c("drug_name" = "drug"))

# Resumen: cuántos FÁRMACOS y cuántas INTERACCIONES por cohorte y tipo
summary_drugs <- edges_sig %>%
  group_by(Cohort, interaction_group) %>%
  summarise(
    n_drugs = n_distinct(drug_name),
    n_edges = n(),
    .groups = "drop"
  ) %>%
  arrange(Cohort, interaction_group)

print(summary_drugs)

# Guardar tabla para revisar después
readr::write_tsv(summary_drugs,
                 "8_8_Significant_drugs_by_interaction_group.tsv")
library(vroom)

sig_drug_prox <- vroom("~/CESC_Network/6_DGIDB/8_8_Significant_drugs_by_interaction_group.tsv")
# =========================================================
# 4) Gráficas fancy
# =========================================================

# Colores consistentes para interaction_group
group_cols <- c(
  inhibits  = "#d73027",  # rojo
  activates = "#1a9850",  # verde
  regulates = "#4575b4"   # azul
)

# --- (A) Barras: número de fármacos por tipo de interacción y cohorte ---
p_drugs <- summary_drugs %>%
  ggplot(aes(x = interaction_group,
             y = n_drugs,
             fill = interaction_group)) +
  geom_col(width = 0.7, color = "black") +
  facet_wrap(~ Cohort, nrow = 1) +
  scale_fill_manual(values = group_cols, name = "Interaction group") +
  labs(
    x = "Tipo de interacción DGIdb",
    y = "Número de fármacos (zc ≤ -1.65)",
    title = "Fármacos significativos por tipo de interacción",
    subtitle = "Sólo fármacos FDA con al menos una diana en la red CESC"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1, face = "bold"),
    strip.text   = element_text(face = "bold"),
    plot.title   = element_text(face = "bold", size = 16),
    legend.position = "right"
  )

ggsave("8_8_Bar_drugs_by_interaction_group.png",
       p_drugs, width = 10, height = 5, dpi = 300)

# --- (B) Barras apiladas: número de interacciones por cohorte ---
p_edges <- summary_drugs %>%
  ggplot(aes(x = Cohort,
             y = n_edges,
             fill = interaction_group)) +
  geom_col(width = 0.7, color = "black") +
  scale_fill_manual(values = group_cols,
                    name = "Interaction group") +
  labs(
    x = "Cohorte HPV",
    y = "Número de interacciones fármaco–gen",
    title = "Carga de interacciones DGIdb entre fármacos significativos y la red",
    subtitle = paste0("Umbral zc ≤ ", z_thresh,
                      " (proximidad menor de lo esperado por azar)")
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x  = element_text(face = "bold"),
    plot.title   = element_text(face = "bold", size = 16),
    legend.position = "right"
  )

ggsave("8_8_Stacked_edges_by_interaction_group.png",
       p_edges, width = 7, height = 5, dpi = 300)

# =========================================================
# 5) (Opcional) tabla de fármacos con su mezcla de tipos de interacción
# =========================================================

drug_profile <- edges_sig %>%
  group_by(Cohort, drug_name) %>%
  summarise(
    has_inhibits   = any(interaction_group == "inhibits"),
    has_activates  = any(interaction_group == "activates"),
    has_regulates  = any(interaction_group == "regulates"),
    n_targets      = n_distinct(gene_name),
    .groups = "drop"
  )

readr::write_tsv(drug_profile,
                 "8_8_Significant_drug_profiles.tsv")
