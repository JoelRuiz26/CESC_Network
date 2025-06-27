# ===================== SCRIPT OVERVIEW ===================== #
# Title: Get Raw Gene Expression Data from GTEx and TCGA-CESC (Cervix)
# Description:
#   This script retrieves raw RNA-Seq counts and metadata for cervix uteri samples
#   from both GTEx and TCGA-CESC. It includes downloading, processing, and saving 
#   of both expression matrices and sample-level metadata for downstream analyses.
# Author: Joel Ruiz Hernandez
# Contact: josejoelruizhernandez@gmail.com
# Date: 2025-04-05
# =========================================================== #
setwd("~/CESC_Network_local/1_Get_Data/")
# ===================== LOAD REQUIRED LIBRARIES ===================== #
library(recount3)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)
library(vroom)

# ===================== PART 1: GTEx (CERVIX UTERI) ===================== #

# Manually download RSE object from GTEx
rse_cervix <- create_rse_manual(
  project = "CERVIX_UTERI",
  project_home = "data_sources/gtex"
)

# View available assays (important to know the structure)
assayNames(rse_cervix)
# [1] "raw_counts"

# Extract raw count matrix (similar to TCGAbiolinks)
counts_matrix <- assays(rse_cervix)[["raw_counts"]] 
# Check dimensions
dim(counts_matrix)
# [1] 63856    19

# Check annotation version
metadata(rse_cervix)$annotation
# [1] "gencode_v26"

# Extract sample metadata
gtex_meta <- as.data.frame(colData(rse_cervix))
colnames(gtex_meta)

# ===================== PART 2: TCGA-CESC ===================== #

# Define sample types
sample_types <- c("Primary Tumor", "Solid Tissue Normal")

# Query RNA-Seq data from TCGA-CESC
query_TCGA_CESC <- GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  access = "open",
  sample.type = sample_types
)

# Extract basic info from query
Output_query_TCGA <- getResults(query_TCGA_CESC) %>%
  dplyr::select(cases, cases.submitter_id, sample_type)
# head(Output_query_TCGA)
#    cases                            cases.submitter_id   sample_type
# 1  TCGA-C5-A1M5-01A-11R-A13Y-07     TCGA-C5-A1M5         Primary Tumor
# 2  TCGA-EK-A2R9-01A-11R-A18M-07     TCGA-EK-A2R9         Primary Tumor

# Download all RNA-seq data files (~307 files)
GDCdownload(query_TCGA_CESC)

# Prepare the SummarizedExperiment object
TCGA_prepare <- GDCprepare(query_TCGA_CESC, summarizedExperiment = TRUE)

# View available assays
# => unstranded
# => stranded_first
# => stranded_second
# => tpm_unstrand
# => fpkm_unstrand
# => fpkm_uq_unstrand

# Extract "unstranded" assay and convert to tibble
unstranded <- assay(TCGA_prepare, "unstranded")
unstranded <- as_tibble(unstranded, rownames = "gene")

# Check dimensions
dim(unstranded)
# [1] 60660   308

# Preview values
head(unstranded[1:10, 1:3])
# Example:
# gene                 TCGA-C5-A1M5...  TCGA-EK-A2R9...
# ENSG00000000003.15         6914             10140
# ENSG00000000005.6             5                 1

# ===================== BUILD TCGA METADATA ===================== #

# Load HPV distribution data
HPV_TCGA <- vroom(file = "~/CESC_Network_local/0_HPV_Distribution/1_1_HPV_Clade_Distrib.tsv")
# dim(HPV_TCGA)  # [1] 275   5

# Preview HPV metadata
# head(HPV_TCGA)
# `TCGA Case ID` tipos_HPV Clado_filogenetico frec percent

# Load clinical metadata
Clinical_info <- GDCquery_clinic("TCGA-CESC", "clinical") %>%
  dplyr::select("submitter_id", "figo_stage", "race", "primary_diagnosis")

# Merge all metadata components
Metadata <- Output_query_TCGA %>%
  left_join(Clinical_info, by = c("cases.submitter_id" = "submitter_id")) %>%
  right_join(HPV_TCGA, by = c("cases.submitter_id" = "TCGA Case ID")) %>%
  mutate_at(vars(figo_stage:Clado_filogenetico, -race),
            ~ ifelse(sample_type == "Solid Tissue Normal", "Solid Tissue Normal", .)) %>%
  dplyr::rename(
    specimenID = cases,
    HPV_type = tipos_HPV,
    HPV_clade = Clado_filogenetico
  )

# Check dimensions
#dim(Metadata)  # [1] 278   10

# Check column names
# colnames(Metadata)
# "specimenID" "cases.submitter_id" "sample_type" "figo_stage"
# "race" "primary_diagnosis" "HPV_type" "HPV_clade" "frec" "percent"

# Filter the count matrix by valid sample IDs in metadata
Cases_metadata <- Metadata %>% pull(specimenID)
unstranded_counts <- unstranded %>% dplyr::select(gene, all_of(Cases_metadata))

# Check filtered dimensions
# dim(unstranded)         # [1] 60660 308
#dim(unstranded_counts)  # [1] 60660 279
# Note: 29 samples likely HPV-negative and excluded

# Optional check
# all(Cases_metadata %in% colnames(unstranded))  # TRUE

# Count samples by HPV clade
CladeHPV_specimenID <- Metadata %>%
  group_by(HPV_clade) %>%
  summarise(num_samples = n()) %>%
  pivot_wider(names_from = HPV_clade, values_from = num_samples)
# A7    A9  negative  otro
# 68   202     3        5

# ===================== SAVE FINAL OUTPUTS ===================== #
#Save GTEX counts
saveRDS(counts_matrix, file = "~/CESC_Network_local/1_Get_Data/1_3_GTEx_Cervix_raw_counts.rds")

# Save TCGA metadata
Metadata <- as_tibble(Metadata)
vroom_write(Metadata, "~/CESC_Network_local/1_Get_Data/1_2_Metadata.tsv", delim = "\t")

# Save count matrix (filtered)
unstranded_counts <- as_tibble(unstranded_counts)
vroom_write(unstranded_counts, "~/CESC_Network_local/1_Get_Data/1_1_unstranded_counts.tsv", delim = "\t")

# Optionally save session
#save.image(file = "~/CESC_Network_local/1_Get_Data/0_Image_data_RNA.RData")
#load(file = "~/CESC_Network_local/1_Get_Data_TCGA/0_Image_data_RNA.RData")

