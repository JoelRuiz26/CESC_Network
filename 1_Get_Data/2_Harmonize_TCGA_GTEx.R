# ============================================================ #
# Title: Harmonize Raw RNA-Seq Counts from TCGA-CESC and GTEx
# Description:
#   This script merges raw gene expression data from TCGA-CESC and GTEx 
#   (cervix uteri), standardizes gene identifiers, filters, and 
#   constructs a unified expression matrix and metadata table.
# Author: Joel Ruiz Hernández
# Contact: josejoelruizhernandez@gmail.com
# Date: 2025-04-05
# ============================================================ #
setwd("~/CESC_Network_local/1_Get_Data/")
# ===================== LOAD REQUIRED LIBRARIES ===================== #
library(recount3)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)
library(vroom)
library(stringr)

# ===================== LOAD AND PREPARE COUNT MATRICES ===================== #

# Load TCGA and GTEx raw counts
tcga <- vroom("~/CESC_Network_local/1_Get_Data/1_1_unstranded_counts.tsv")
gtex <- readRDS("~/CESC_Network_local/1_Get_Data/1_3_GTEx_Cervix_raw_counts.rds") %>%
  as.data.frame() %>%
  rownames_to_column("gene")

# Standardize gene IDs (remove Ensembl version suffix)
# Remove Ensembl version suffix
tcga <- tcga %>% mutate(gene = str_remove(gene, "\\..*$"))
gtex <- gtex %>% mutate(gene = str_remove(gene, "\\..*$"))

# Collapse duplicated genes by summing counts
tcga <- tcga %>%
  group_by(gene) %>%
  summarise(across(everything(), sum), .groups = "drop")

gtex <- gtex %>%
  group_by(gene) %>%
  summarise(across(everything(), sum), .groups = "drop")

# Keep only intersecting genes
genes_comunes <- intersect(tcga$gene, gtex$gene) 
#length(genes_comunes) #[1] 57562
tcga <- tcga %>% filter(gene %in% genes_comunes) %>% arrange(gene)
dim(tcga) #[1] 57562   279
gtex <- gtex %>% filter(gene %in% genes_comunes) %>% arrange(gene)
dim(gtex) #[1] 57562    20

# Ensure same gene order
stopifnot(all(tcga$gene == gtex$gene))

# Combine count matrices
tcga_mat <- tcga %>% column_to_rownames("gene")
gtex_mat <- gtex %>% column_to_rownames("gene")
combined_counts <- cbind(tcga_mat, gtex_mat) %>%
  as.data.frame() %>%
  rownames_to_column("gene")
dim(combined_counts) #[1] 57562   298

# Remove genes with zero expression across all samples
combined_counts <- combined_counts[rowSums(combined_counts[,-1]) > 0, ]
dim(combined_counts) #[1] 56033   298

# ===================== BUILD COMBINED METADATA ===================== #

# Load TCGA metadata and remove unused columns
metadata_tcga <- vroom("~/CESC_Network_local/1_Get_Data/1_2_Metadata.tsv") %>%
  dplyr::select(-race, -frec, -percent)

#Identify GTEx samples from count matrix
combined_samples <- colnames(combined_counts)[-1]  # remove gene column
gtex_samples <- setdiff(combined_samples, metadata_tcga$specimenID)

#Append GTEx metadata rows
metadata_tcga <- metadata_tcga %>%
  bind_rows(
    tibble(specimenID = gtex_samples, cases.submitter_id = gtex_samples) %>%
      mutate(across(-c(specimenID, cases.submitter_id), ~ "Solid Tissue Normal"))
  )

# Filter HPV and add 'source' column
metadata_tcga <- metadata_tcga %>%
  filter(!HPV_clade %in% "otro", !HPV_type %in% "HPV70") %>%
  mutate(
    HPV_clade = recode(HPV_clade, "A7" = "A7_clade", "A9" = "A9_clade"),
    source = if_else(specimenID %in% gtex_samples, "GTEX", "TCGA")
  ) %>%
  dplyr::select(specimenID, source, everything())

metadata_full <- metadata_tcga %>%
  mutate(across(everything(), ~replace_na(.x, "Solid Tissue Normal")))

# ===================== FINALIZE AND SAVE ===================== #

# Filter count matrix by valid metadata samples
combined_counts <- combined_counts %>%
  dplyr::select(gene, all_of(metadata_full$specimenID))

# Save results
saveRDS(metadata_full, file = "~/CESC_Network_local/1_Get_Data/2_Harmoni_metadata_TCGA_GTEx.rds")
saveRDS(combined_counts, file = "~/CESC_Network_local/1_Get_Data/2_Harmoni_counts_TCGA_GTEx.rds")
