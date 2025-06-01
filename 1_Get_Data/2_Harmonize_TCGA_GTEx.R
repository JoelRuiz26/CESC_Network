################################################################################
# Script: Integrate and Harmonize TCGA and GTEx Expression Data (Counts + Metadata)
# Author: Joel Ruiz Hernandez contact: josejoelruizhernandez@gmail.com
# Date: 20/12/2024
# Description:
#   This script performs a secure and traceable integration of RNA-seq gene 
#   expression count matrices from TCGA (CESC) and GTEx (Cervix Uteri), using
#   GENCODE v26 as reference for gene ID validation. The output includes a unified
#   count matrix and harmonized sample metadata.

# INPUT:
#   - TCGA metadata:          "~/1_Get_Data_TCGA/1_2_Metadata.tsv"
#   - TCGA counts (unstranded): "~/1_Get_Data_TCGA/1_1_unstranded_counts.tsv"
#   - GTEx raw counts (RDS): "~/1_Get_Data_TCGA/GTEx_Cervix_raw_counts.rds"
#
# OUTPUT:
#   - Combined count matrix:   "~/1_Get_Data_TCGA/2_combined_counts_TCGA_GTEx.rds"
#   - Combined sample metadata:"~/1_Get_Data_TCGA/2_combined_metadata_TCGA_GTEx.rds"
################################################################################

# ========================= #
# Load Required Libraries   #
# ========================= #

library(tidyverse)
library(biomaRt)
library(vroom)

# ============================ #
# Load Input Data             #
# ============================ #

metadata_TCGA <- vroom(file = "~/1_Get_Data_TCGA/1_2_Metadata.tsv")                 # [278 x 8]
unstranded_counts_TCGA <- vroom(file = "~/1_Get_Data_TCGA/1_1_unstranded_counts.tsv") # [60660 x 279]
raw_counts_GTEx <- readRDS("~/1_Get_Data_TCGA/GTEx_Cervix_raw_counts.rds")           # [63856 x 19]

# NOTE:
# - GTEx uses: GRCh38 + GENCODE v26
# - TCGA uses: GRCh38 + GENCODE v22
# - Although both use GRCh38, using different GENCODE versions may cause 
#   mismatches in gene IDs, pseudogenes, and annotation discrepancies.

# ========================================= #
# Step 1: Connect to Ensembl v100 (GENCODE v26)
# ========================================= #

ensembl_v26 <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl", version = 100)

# ===================================================== #
# Step 2: Extract all valid GENCODE v26 gene IDs (with version)
# ===================================================== #

genes_validos_v26 <- getBM(
  attributes = c("ensembl_gene_id_version"),
  mart = ensembl_v26
) %>%
  pull(ensembl_gene_id_version) %>%
  unique()

# ================================================================= #
# Step 3: Filter TCGA genes to keep only those present in GENCODE v26
# ================================================================= #

tcga_genes_full <- unstranded_counts_TCGA$gene
tcga_valid <- tcga_genes_full %in% genes_validos_v26
tcga_counts_filtered <- unstranded_counts_TCGA[tcga_valid, ]
# Result: From 60660 to 58185 rows (~4% excluded)

# ================================================================= #
# Step 4: Filter GTEx genes to keep only those present in GENCODE v26
# ================================================================= #

gtex_genes_full <- rownames(raw_counts_GTEx)
gtex_valid <- gtex_genes_full %in% genes_validos_v26
gtex_counts_filtered <- raw_counts_GTEx[gtex_valid, ]
# Result: 43596 genes retained

# ======================================================== #
# Step 5: Identify the exact intersecting genes (with version)
# ======================================================== #

common_genes_exact <- intersect(
  tcga_counts_filtered$gene,
  rownames(gtex_counts_filtered)
)
# Result: 36919 common genes

# ========================================================== #
# Step 6: Subset both matrices by common genes (same order)
# ========================================================== #

tcga_final <- tcga_counts_filtered %>%
  filter(gene %in% common_genes_exact) %>%
  column_to_rownames("gene") %>%
  as.matrix()

gtex_final <- gtex_counts_filtered[common_genes_exact, ]

# ============================================== #
# Step 7: Combine the filtered matrices by column
# ============================================== #

stopifnot(identical(rownames(tcga_final), rownames(gtex_final)))  # Ensure gene order matches
combined_counts <- cbind(tcga_final, gtex_final)
dim(combined_counts)  # [36919 x 297]

# ====================================== #
# Step 8: Save the combined count matrix
# ====================================== #

#saveRDS(combined_counts, "~/1_Get_Data_TCGA/2_combined_counts_TCGA_GTEx.rds")

####Get universe for enrichment
#raw_counts_GTEx <- raw_counts_GTEx %>% as.data.frame()
#genes <- rownames(raw_counts_GTEx) %>% as.character()
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 100)

#annot_universe <- getBM(
#  attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "hgnc_symbol"),
#  filters = "ensembl_gene_id_version",
#  values = genes,
#  mart = ensembl)%>% filter(hgnc_symbol != "") #VOLVER A CORRER

#dim(annot_universe) #[1] 24546     3
#length(unique(annot_universe$hgnc_symbol))
#[1] 21435
#saveRDS(annot_universe,file = "~/3_DGE/DESeq2/annot_universe_24546.rds")

