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

# ===================== LOAD REQUIRED LIBRARIES ===================== #
library(recount3)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)
library(vroom)
library(stringr)

# ===================== PART 1: LOAD AND PREPARE COUNT MATRICES ===================== #

# Load TCGA and GTEx raw counts
tcga <- vroom("~/1_Get_Data_TCGA/1_1_unstranded_counts.tsv")
gtex <- readRDS("~/1_Get_Data_TCGA/1_GTEx_Cervix_raw_counts.rds") %>%
  as.data.frame() %>%
  rownames_to_column("gene")

# Standardize gene IDs (remove Ensembl version suffix)
tcga <- tcga %>% mutate(gene = str_remove(gene, "\\..*$"))
gtex <- gtex %>% mutate(gene = str_remove(gene, "\\..*$"))

# Remove genes with zero expression across all samples
tcga <- tcga[rowSums(tcga[,-1]) > 0, ]
gtex <- gtex[rowSums(gtex[,-1]) > 0, ]

# Keep only intersecting genes
genes_comunes <- intersect(tcga$gene, gtex$gene)
tcga <- tcga %>% filter(gene %in% genes_comunes) %>% arrange(gene)
gtex <- gtex %>% filter(gene %in% genes_comunes) %>% arrange(gene)

# Ensure same gene order
stopifnot(all(tcga$gene == gtex$gene))

# Combine count matrices
tcga_mat <- tcga %>% column_to_rownames("gene")
gtex_mat <- gtex %>% column_to_rownames("gene")
combined_counts <- cbind(tcga_mat, gtex_mat) %>%
  as.data.frame() %>%
  rownames_to_column("gene")

# ===================== PART 2: BUILD COMBINED METADATA ===================== #

# Step 1: Load TCGA metadata and remove unused columns
metadata_tcga <- vroom("~/1_Get_Data_TCGA/1_2_Metadata.tsv") %>%
  select(-race, -frec, -percent)

# Step 2: Identify GTEx samples from count matrix
combined_samples <- colnames(combined_counts)[-1]  # remove gene column
gtex_samples <- setdiff(combined_samples, metadata_tcga$specimenID)

# Step 3: Append GTEx metadata rows
metadata_tcga <- metadata_tcga %>%
  bind_rows(
    tibble(specimenID = gtex_samples, cases.submitter_id = gtex_samples) %>%
      mutate(across(-c(specimenID, cases.submitter_id), ~ "Solid Tissue Normal"))
  )

# Step 4: Filter HPV and add 'source' column
metadata_tcga <- metadata_tcga %>%
  filter(!HPV_clade %in% "otro", !HPV_type %in% "HPV70") %>%
  mutate(
    HPV_clade = recode(HPV_clade, "A7" = "A7_clade", "A9" = "A9_clade"),
    source = if_else(specimenID %in% gtex_samples, "GTEX", "TCGA")
  ) %>%
  select(specimenID, source, everything())

metadata_tcga <- metadata_tcga %>%
  mutate(across(everything(), ~replace_na(.x, "Solid Tissue Normal")))

# ===================== PART 3: FINALIZE AND SAVE ===================== #

# Filter count matrix by valid metadata samples
combined_counts <- combined_counts %>%
  select(gene, all_of(metadata_tcga$specimenID))

# Save results
saveRDS(metadata_tcga, file = "~/1_Get_Data_TCGA/2_Harmoni_metadata_TCGA_GTEx.rds")
saveRDS(combined_counts, file = "~/1_Get_Data_TCGA/2_Harmoni_counts_TCGA_GTEx.rds")
