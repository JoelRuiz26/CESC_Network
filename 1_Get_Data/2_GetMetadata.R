############# 1) GET METADATA (TCGA-CESC clinicalMatrix from Xena) #############

suppressPackageStartupMessages({
  library(UCSCXenaTools)
  library(dplyr)
  library(vroom)
})

## Output directory
out_dir <- "~/CESC_Network/1_Get_Data/Octad_matrix/1_Output_rds/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## 1) Load XenaData and locate TCGA-CESC clinicalMatrix in tcgaHub
data("XenaData")

cesc_clin_info <- XenaData %>%
  filter(
    XenaHostNames == "tcgaHub",
    grepl("Cervical", XenaCohorts),
    grepl("CESC_clinicalMatrix", XenaDatasets)
  )

## 2) Build XenaHub ONLY with that dataset
df_cesc <- XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = cesc_clin_info$XenaDatasets[1])

## 3) Download and prepare clinical matrix
xe_cesc <- XenaQuery(df_cesc) %>%
  XenaDownload(destdir = out_dir)

clin_tcga_cesc <- XenaPrepare(xe_cesc)

## 4) Save full clinical matrix
saveRDS(
  clin_tcga_cesc,
  file = file.path(out_dir, "2_0_TCGA_CESC_clinicalMatrix_full.rds")
)


############# 2) Add HPV, summarized stage / histology + GTEx normals ##########
## HPV_genotypes.tsv from supplementary PDF of doi: 10.1038/s41598-020-71300-7
## PMID: 32868873; PMCID: PMC7459114.

HPV_TCGA <- vroom::vroom(
  "~/CESC_Network/0_HPV_Distribution/1_1_HPV_Clade_Distrib.tsv",
  show_col_types = FALSE
)

## OCTAD cervix pheno (TCGA + GTEx)
pheno <- readRDS(
  "~/CESC_Network/1_Get_Data/Octad_matrix/1_Output_rds/1_0_Cervix_pheno_TCGA_GTEx.rds"
)

## 2A) Start from TCGA clinical and tag source
meta_all <- clin_tcga_cesc %>%
  mutate(source = "TCGA")

## 2B) Append GTEx rows with only key info; rest will be NA and filled later
gtex_rows <- pheno %>%
  filter(data.source == "GTEX") %>%
  transmute(
    sampleID          = sample.id,
    sample_type       = "Solid Tissue Normal",
    tumor_tissue_site = biopsy.site,
    source            = "GTEX"
  )

meta_all <- bind_rows(
  meta_all,
  # Bind with GTEx rows; union of columns, TCGA cols NA for GTEx
  gtex_rows
)

############# 3) Add HPV, Stage_group and histological_group ###################

meta_all <- meta_all %>%
  # HPV from literature (only TCGA rows will match bcr_patient_barcode)
  left_join(HPV_TCGA, by = "bcr_patient_barcode") %>%
  mutate(
    # Mark normals based on sample_type
    HPV_type  = ifelse(sample_type == "Solid Tissue Normal", "Solid Tissue Normal", HPV_type),
    HPV_clade = ifelse(sample_type == "Solid Tissue Normal", "Solid Tissue Normal", HPV_clade),
    # Recode A7/A9 to A7_clade / A9_clade for consistency
    HPV_clade = dplyr::recode(HPV_clade,
                              "A7" = "A7_clade",
                              "A9" = "A9_clade",
                              .default = HPV_clade),
    # Summarize clinical_stage into Stage_group (Stage I–IV)
    Stage_group = gsub("Stage ([IV]+).*", "Stage \\1", clinical_stage),
    Stage_group = ifelse(
      grepl("^Stage [IV]+$", Stage_group),
      Stage_group,
      "(NOT REPORTED)"
    ),
    # Summarize histological_type into histological_group
    histological_group = dplyr::case_when(
      sample_type == "Solid Tissue Normal"                          ~ "Solid Tissue Normal",
      histological_type == "Cervical Squamous Cell Carcinoma"       ~ "Squamous",
      histological_type == "Adenosquamous"                          ~ "Adenosquamous",
      histological_type %in% c(
        "Endocervical Adenocarcinoma of the Usual Type",
        "Endocervical Type of Adenocarcinoma",
        "Mucinous Adenocarcinoma of Endocervical Type",
        "Endometrioid Adenocarcinoma of Endocervix"
      )                                                             ~ "Adenocarcinoma",
      TRUE                                                          ~ "(NOT REPORTED)"
    )
  )

############# 4) Compact table + filters (A7, A9, normals) #####################

metadata_cesc_all <- meta_all %>%
  dplyr::select(
    sampleID,
    source,
    bcr_patient_barcode,
    sample_type,
    HPV_type,
    HPV_clade,
    clinical_stage,
    Stage_group,
    histological_type,
    histological_group,
    tumor_tissue_site,
    vital_status,
    days_to_last_followup,
    days_to_death
  ) %>%
  filter(
    !grepl("metast", sample_type %||% "", ignore.case = TRUE),              # remove metastasis
    !(HPV_clade %in% "otro"),                                              # drop "otro"
    !(HPV_type  %in% "HPV70"),                                             # drop HPV70
    HPV_clade %in% c("A7_clade", "A9_clade", "Solid Tissue Normal") |      # keep A7/A9/normals
      is.na(HPV_clade)                                                     # (por si quieres inspeccionar NA)
  )

dim(metadata_cesc_all) #[1] 314  15
table(metadata_cesc_all$source)
#GTEX TCGA 
#10  304 
table(metadata_cesc_all$histological_group, metadata_cesc_all$HPV_clade, useNA = "ifany")

#                      A7_clade A9_clade Solid Tissue Normal <NA>
#(NOT REPORTED)             0        0                   0    1
#Adenocarcinoma            11       26                   0   10
#Adenosquamous              2        1                   0    3
#Solid Tissue Normal        0        0                  13    0
#Squamous                  53      175                   0   19

## Save final summarized metadata (TCGA + GTEx, filtered)
saveRDS(
  metadata_cesc_all,
  file = file.path(out_dir, "2_1_CESC_TCGA_GTEx_metadata_HPV_stage_histology.rds")
)

# ============================
#remove rows with NA HPV_clade

metadata_cesc_all_noNA <- metadata_cesc_all %>%
  filter(!is.na(HPV_clade))

saveRDS(
  metadata_cesc_all_noNA,
  file = file.path(out_dir, "2_2_CESC_TCGA_GTEx_metadata_clade.rds")
)

#=================== CLade samples for OCTAD

out_dir_OCTAD <- "/home/jruiz/OCTAD_Cervical_Cancer/0_Get_Data/"

# A7 samples
A7_samples <- metadata_cesc_all_noNA %>%
  filter(HPV_clade == "A7_clade") %>%
  pull(sampleID)

saveRDS(
  A7_samples,
  file = file.path(out_dir_OCTAD, "0_A7_Samples.rds")
)

# A9 samples
A9_samples <- metadata_cesc_all_noNA %>%
  filter(HPV_clade == "A9_clade") %>%
  pull(sampleID)

saveRDS(
  A9_samples,
  file = file.path(out_dir_OCTAD, "0_A9_Samples.rds")
)
