# ====================================================
# Network-based Drug Repurposing using DGIdb
# ====================================================
# Author: Joel 
# Fecha: 2025
# Referencia: Guney et al. 2016 - Nature Communications
# ----------------------------------------------------
# Este script:
# 1. Usa TODOS los fármacos aprobados de DGIdb
# 2. Asocia targets via DGIdb
# 3. Filtra targets en red
# 4. Calcula distancias a DGE
# 5. Rankea fármacos por proximidad
# 6. Visualiza resultados

# ====================================================

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
head(dgidb_interactions)
# A tibble: 6 × 13
#gene_claim_name                gene_concept_id gene_name interaction_source_d…¹ interaction_source_d…² interaction_type interaction_score drug_claim_name
#<chr>                          <chr>           <chr>     <chr>                  <chr>                  <chr>            <chr>             <chr>          
#1 ANGIOTENSIN II RECEPTOR TYPE-1 hgnc:336        AGTR1     TTD                    2020.06.01             NULL             0.66747039706391… ANGIOTENSIN II 
#2 VASCULAR ENDOTHELIAL GROWTH F… hgnc:6307       KDR       TTD                    2020.06.01             NULL             0.026412309       HKI-272        
#3 KIT                            hgnc:6342       KIT       DoCM                   10-Apr-24              NULL             1.029066006546944 IMATINIB       
#4 CCL4                           hgnc:10630      CCL4      NCI                    14-Sep-17              NULL             3.182283105193606 CLODRONATE     
#5 BAZ2B                          hgnc:963        BAZ2B     DTC                    9/2/20                 NULL             0.00767656        OXYTETRACYCLINE
#6 CSNK2A3                        hgnc:2458       CSNK2A3   PharmGKB               4/5/24                 NULL             0.72424374118199… cisplatin      
# ℹ abbreviated names: ¹​interaction_source_db_name, ²​interaction_source_db_version
# ℹ 5 more variables: drug_concept_id <chr>, drug_name <chr>, approved <chr>, immunotherapy <chr>, anti_neoplastic <chr>
colnames(dgidb_interactions)
#[1] "gene_claim_name"               "gene_concept_id"               "gene_name"                     "interaction_source_db_name"   
#[5] "interaction_source_db_version" "interaction_type"              "interaction_score"             "drug_claim_name"              
#[9] "drug_concept_id"               "drug_name"                     "approved"                      "immunotherapy"                
#[13] "anti_neoplastic"   

# ----------------------------
# FILTER APPROVED DRUGS
# ----------------------------
dgidb_interactions <- dgidb_interactions %>%
  filter(approved == "TRUE")
dim(dgidb_interactions) #[1] 37970    13

# Make score numeric type
# Filtro mínimo de confianza por interacción (ajusta el umbral)
dgidb_interactions <- dgidb_interactions %>%
  mutate(
    interaction_score = ifelse(interaction_score == "NULL", "0", interaction_score),
    interaction_score = as.numeric(interaction_score)
  )

dgidb_interactions <- dgidb_interactions %>% filter(interaction_score >= 20)
dim(dgidb_interactions) #[1] 136    13 with 20 filter

summary(dgidb_interactions$interaction_score)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.00000   0.02012   0.09299   0.71798   0.38467 157.52301 

# ----------------------------
# DEFINE STANDARDIZATION FUNCTION
# ----------------------------
normalize_drug_name <- function(name) {
  name %>%
    tolower() %>%
    str_replace_all("[[:punct:]]", " ") %>%
    str_squish()
}

# ----------------------------
# APPLY STANDARDIZATION (DGIdb)
# ----------------------------
dgidb_interactions <- dgidb_interactions %>%
  mutate(norm_drug_name = normalize_drug_name(drug_name)) %>%
  mutate(norm_drug_name_claim = normalize_drug_name(drug_claim_name))

# ----------------------------
# GET UNIQUE DGIdb DRUG NAMES
# ----------------------------
dgidb_drug_names <- union(
  unique(dgidb_interactions$norm_drug_name),
  unique(dgidb_interactions$norm_drug_name_claim)
)
length(dgidb_drug_names) #[1] ~5k

# ----------------------------
# (REPLACING OCTAD) USE ALL DGIdb DRUGS
# ----------------------------
all_matches_A7 <- dgidb_drug_names
all_matches_A9 <- dgidb_drug_names

cat("✅ Total final A7 matches (DGIdb-only):", length(all_matches_A7), "\n")
cat("✅ Total final A9 matches (DGIdb-only):", length(all_matches_A9), "\n")

# ----------------------------
# SAVE RESULTS (mismos nombres de salida)
# ----------------------------

write_tsv(tibble(A7_matches = all_matches_A7), "~/CESC_Network/6_DGIDB/8_1_Drugs_A7_matches_OCTAD_DGIDB.tsv")
write_tsv(tibble(A9_matches = all_matches_A9), "~/CESC_Network/6_DGIDB/8_1_Drugs_A9_matches_OCTAD_DGIDB.tsv")

# ----------------------------
# GET DRUG-TARGET INTERACTIONS
# ----------------------------

interactions_A7 <- dgidb_interactions %>%
  filter(norm_drug_name %in% all_matches_A7 | norm_drug_name_claim %in% all_matches_A7) %>%
  dplyr::select(drug_name, drug_claim_name, gene_name, interaction_score) %>%
  distinct()

interactions_A9 <- dgidb_interactions %>%
  filter(norm_drug_name %in% all_matches_A9 | norm_drug_name_claim %in% all_matches_A9) %>%
  dplyr::select(drug_name, drug_claim_name, gene_name, interaction_score) %>%
  distinct()

# ----------------------------
# LOAD NETWORK GENES
# ----------------------------
graph_A7 <- readRDS("~/CESC_Network/5_Modularity/7_1_Louvain/7_1_1_Louvain_A7_graph.rds")
graph_A9 <- readRDS("~/CESC_Network/5_Modularity/7_1_Louvain/7_1_1_Louvain_A9_graph.rds")

network_genes_A7 <- V(graph_A7)$name
network_genes_A9 <- V(graph_A9)$name

summary(graph_A7)
#+ attr: name (v/c), Community (v/n), PageRank (v/n), Community_name (v/c), Community_size (v/n), MI (e/n)

# ----------------------------
# FILTER DRUG TARGETS IN NETWORK
# ----------------------------
interactions_A7 <- interactions_A7 %>%
  filter(gene_name %in% network_genes_A7)

interactions_A9 <- interactions_A9 %>%
  filter(gene_name %in% network_genes_A9)

# ----------------------------
# LOAD DGE (key genes)
# ----------------------------
DGE_A7 <- readRDS("~/OCTAD_Cervical_Cancer/1_DGE_signature/1_0_Output_rds/1_5_CoreGenes_min/core_genes_A7.rds")
DGE_A9 <- readRDS("~/OCTAD_Cervical_Cancer/1_DGE_signature/1_0_Output_rds/1_5_CoreGenes_min/core_genes_A9.rds")

dge_genes_A7 <- unique(DGE_A7)
dge_genes_A9 <- unique(DGE_A9)

# ----------------------------
# NETWORK PROXIMITY CALCULATION
# ----------------------------
#Proximity(T,D)= mean over t∈T of min distance(t, D)
#Selectivity score S = mean_dist_to_nonDGE / mean_dist_to_DGE

# Load libraries
library(tidyverse); library(igraph); library(ComplexHeatmap); library(circlize); library(ggplot2); library(viridis)

#----------------------------------------------------
# 1. Normalization and selectivity + proximity
#----------------------------------------------------
interactions_A7 <- interactions_A7 %>% mutate(interaction_score_scaled = (interaction_score - min(interaction_score, na.rm=TRUE)) / (max(interaction_score, na.rm=TRUE) - min(interaction_score, na.rm=TRUE)))
interactions_A9 <- interactions_A9 %>% mutate(interaction_score_scaled = (interaction_score - min(interaction_score, na.rm=TRUE)) / (max(interaction_score, na.rm=TRUE) - min(interaction_score, na.rm=TRUE)))

mean_min_distance <- function(graph, from_nodes, to_nodes) {
  from_nodes <- intersect(from_nodes, V(graph)$name)
  to_nodes <- intersect(to_nodes, V(graph)$name)
  if(length(from_nodes) == 0 || length(to_nodes) == 0) return(NA)
  d <- distances(graph, v=from_nodes, to=to_nodes)
  mean(apply(d, 1, min), na.rm=TRUE)
}

get_selectivity_scores <- function(graph, drug_targets, DGE, nonDGE) {
  prox_DGE <- mean_min_distance(graph, drug_targets, DGE)
  prox_nonDGE <- mean_min_distance(graph, drug_targets, nonDGE)
  if (is.na(prox_DGE) || prox_DGE == 0) return(tibble(prox_to_DGE=prox_DGE, prox_to_nonDGE=prox_nonDGE, selectivity_score=NA))
  tibble(prox_to_DGE=prox_DGE, prox_to_nonDGE=prox_nonDGE, selectivity_score=prox_nonDGE/prox_DGE)
}

compute_zscore <- function(graph, targets, DGE, nperm=500) {
  obs <- mean_min_distance(graph, targets, DGE)
  deg <- igraph::degree(graph)
  bins <- cut(deg, breaks = quantile(deg, probs = seq(0,1,0.1), na.rm=TRUE), include.lowest=TRUE)
  bins_t <- bins[targets]; bins_d <- bins[DGE]
  dist_rand <- replicate(nperm, {
    if (any(is.na(bins_t)) || any(is.na(bins_d))) {
      samp_t <- sample(V(graph)$name, length(targets))
      samp_d <- sample(V(graph)$name, length(DGE))
    } else {
      samp_t <- unlist(lapply(unique(bins_t), function(b) sample(names(bins[bins==b]), sum(bins_t==b))))
      samp_d <- unlist(lapply(unique(bins_d), function(b) sample(names(bins[bins==b]), sum(bins_d==b))))
    }
    mean_min_distance(graph, samp_t, samp_d)
  })
  (obs - mean(dist_rand, na.rm=TRUE)) / sd(dist_rand, na.rm=TRUE)
}

compute_drug_scores <- function(graph, interactions, dge_genes) {
  dge_genes <- intersect(dge_genes, V(graph)$name)
  nonDGE <- setdiff(V(graph)$name, dge_genes)
  all_drugs <- unique(interactions$drug_name)
  
  results <- lapply(all_drugs, function(drug) {
    targets <- interactions %>% filter(drug_name==drug) %>% pull(gene_name) %>% intersect(V(graph)$name)
    if (length(targets)==0) return(NULL)
    sel <- get_selectivity_scores(graph, targets, dge_genes, nonDGE)
    mi <- mean(interactions %>% filter(drug_name==drug) %>% pull(interaction_score_scaled), na.rm=TRUE)
    z <- compute_zscore(graph, targets, dge_genes)
    tibble(drug=drug, mean_interaction_score=mi, prox_to_DGE=sel$prox_to_DGE,
           selectivity_score=sel$selectivity_score, combined_score=sel$selectivity_score*mi,
           z_score=z)
  }) %>% bind_rows() %>% arrange(desc(combined_score))
  results
}

# Run for A7 and A9
results_A7 <- compute_drug_scores(graph_A7, interactions_A7, dge_genes_A7) #slow

results_A9 <- compute_drug_scores(graph_A9, interactions_A9, dge_genes_A9) #slow

save.image("~/CESC_Network/6_DGIDB/8_4_Image_DGIDB.RData")
# load("~/CESC_Network/8_DGIDB/8_4_Image_DGIDB.RData")

#----------------------------------------------------

# =============================================================
# 4️⃣ PREPARE COMMUNITY × DRUG MATRIX
# =============================================================
prepare_comm_matrix <- function(graph, interactions, dge_genes, results_df) {
  gene2comm <- setNames(V(graph)$Community_name, V(graph)$name)
  dge <- intersect(dge_genes, V(graph)$name)
  comms <- unique(gene2comm[dge])
  
  genes_by_comm <- lapply(comms, function(c) intersect(names(gene2comm)[gene2comm == c], dge))
  names(genes_by_comm) <- comms
  
  genes_all <- intersect(c(interactions$gene_name, dge), V(graph)$name)
  dist_mat_all <- distances(graph, v = genes_all, to = genes_all)
  
  drugs <- results_df$drug
  M <- sapply(drugs, function(drug) {
    targets <- intersect(interactions %>% filter(drug_name == drug) %>% pull(gene_name), genes_all)
    sapply(comms, function(c) {
      g <- genes_by_comm[[c]]
      if (length(targets) == 0 || length(g) == 0) return(NA)
      mean(sapply(g, function(x) min(dist_mat_all[x, targets], na.rm = TRUE)), na.rm = TRUE)
    })
  })
  rownames(M) <- comms
  colnames(M) <- drugs
  M
}

# =============================================================
# 5️⃣ BUILD MATRICES
# =============================================================
mat_A7 <- prepare_comm_matrix(graph_A7, interactions_A7, dge_genes_A7, results_A7)
mat_A9 <- prepare_comm_matrix(graph_A9, interactions_A9, dge_genes_A9, results_A9)

# Remove all-NA rows
mat_A7 <- mat_A7[rowSums(is.na(mat_A7)) < ncol(mat_A7), , drop=FALSE]
mat_A9 <- mat_A9[rowSums(is.na(mat_A9)) < ncol(mat_A9), , drop=FALSE]

# =============================================================
# 5️⃣ Annotate Heatmap Function (Clean and Professional)
# =============================================================
annotate_heatmap <- function(mat, results, cohort_label) {
  require(ComplexHeatmap)
  require(circlize)
  
  # Map of z-scores
  z_map <- setNames(results$z_score, results$drug)
  
  # Column labels with asterisks
  new_colnames <- colnames(mat)
  new_colnames <- ifelse(!is.na(z_map[new_colnames]) & z_map[new_colnames] < -1.5,
                         paste0(new_colnames, " *"), new_colnames)
  colnames(mat) <- new_colnames
  
  # Single-hue color scale (e.g., navy blue)
  proximity_palette <- colorRamp2(
    c(min(mat, na.rm=TRUE), max(mat, na.rm=TRUE)),
    c("#08306B", "#FFFAEC")
  )
  
  Heatmap(
    mat,
    name = "Mean Proximity",
    col = proximity_palette,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 14, fontface = "bold"), 
    column_names_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Mean Proximity Distance",
      legend_height = unit(4, "cm")
    ),
    column_title = paste0(
      "HPV Clade ", cohort_label, 
      ": Drug–Gene Community Proximity\n"
    ),
    column_title_gp = gpar(fontsize = 16, fontface = "bold")
  )
}

# =============================================================
# 6️⃣ DRAW HEATMAPS
# =============================================================
ht_A7 <- annotate_heatmap(
  mat_A7,
  results_A7 %>% filter(drug %in% colnames(mat_A7)),
  "A7"
)
draw(ht_A7)

ht_A9 <- annotate_heatmap(
  mat_A9,
  results_A9 %>% filter(drug %in% colnames(mat_A9)),
  "A9"
)
draw(ht_A9)

# 7️⃣ SAVE AS PDF and PNG
# =============================================================
pdf("~/CESC_Network/6_DGIDB/8_3_A7_Heatmap_with_Zscore.pdf", width=12, height=8)
draw(ht_A7)
dev.off()

pdf("~/CESC_Network/6_DGIDB/8_3_A9_Heatmap_with_Zscore.pdf", width=12, height=8)
draw(ht_A9)
dev.off()

png("~/CESC_Network/6_DGIDB/8_3_A7_Heatmap_with_Zscore.png", width=1200, height=800, res=150)
draw(ht_A7)
dev.off()

png("~/CESC_Network/6_DGIDB/8_3_A9_Heatmap_with_Zscore.png", width=1200, height=800, res=150)
draw(ht_A9)
dev.off()
